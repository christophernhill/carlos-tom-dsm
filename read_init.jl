using Pkg
Pkg.add("Plots")
Pkg.add("AbstractPlotting")
Pkg.add("GLMakie")
Pkg.add("Printf")
Pkg.add("LaTeXStrings")
using Plots
using AbstractPlotting
using GLMakie
using Printf
using LaTeXStrings

# Grid
dxf="delX.bin"; nx=35;
dyf="delY.bin"; ny=35;
dzf="delR.bin"; nz=47;
dx=bswap.( reinterpret(Float64, read(dxf,sizeof(Float64)*nx ) ) )
dy=bswap.( reinterpret(Float64, read(dyf,sizeof(Float64)*ny ) ) )
dz=bswap.( reinterpret(Float64, read(dzf,sizeof(Float64)*nz ) ) )
zf=[0;cumsum(dz)];zc=0.5*(zf[2:end]+zf[1:end-1])
xf=[0;cumsum(dx)];xc=0.5*(xf[2:end]+xf[1:end-1])
yf=[0;cumsum(dy)];yc=0.5*(yf[2:end]+yf[1:end-1])

# Bathymetry
bfile="topog.mercator"
bathy=bswap.( reinterpret(Float64, read(bfile,sizeof(Float64)*nx*ny ) ) )
bathy=reshape(bathy,(nx,ny))

# Initial fields
uinitf="u0.bin"
vinitf="v0.bin"
tinitf="t0.bin"
sinitf="s0.bin"
u0=bswap.( reinterpret(Float64, read(uinitf,sizeof(Float64)*nx*ny*nz ) ) )
u0=reshape(u0,(nx,ny,nz))
v0=bswap.( reinterpret(Float64, read(vinitf,sizeof(Float64)*nx*ny*nz ) ) )
v0=reshape(v0,(nx,ny,nz))
t0=bswap.( reinterpret(Float64, read(tinitf,sizeof(Float64)*nx*ny*nz ) ) )
t0=reshape(t0,(nx,ny,nz))
s0=bswap.( reinterpret(Float64, read(sinitf,sizeof(Float64)*nx*ny*nz ) ) )
s0=reshape(s0,(nx,ny,nz))

# Read in output tiles
function get_tiled_file(fpref,nsx,snx,nsy,sny,nr)
    fld=ones(snx*nsx,sny*nsy,nr);
    for j in 1:nsy, i in 1:nsx
        fn=fpref * @sprintf(".%3.3d",i) * @sprintf(".%3.3d",j) * ".data"
        phi=bswap.( reinterpret(Float32, read(fn,sizeof(Float32)*snx*sny*nr ) ) )
        phi=reshape(phi,(snx,sny,nr))
        i0=(i-1)*snx+1;i1=i0-1+snx
        j0=(j-1)*sny+1;j1=j0-1+sny
        fld[i0:i1,j0:j1,:].=phi
    end
    return fld 
end
# Read in output Depth map
nsx=5; nsy=5;
snx=7; sny=7;
nr=1;
fpref="../run08/Depth"
dmap=get_tiled_file(fpref,nsx,snx,nsy,sny,nr)
nr=47;
fpref="../run08/hFacC"
hfc=get_tiled_file(fpref,nsx,snx,nsy,sny,nr)
fpref="../run08/hFacW"
hfw=get_tiled_file(fpref,nsx,snx,nsy,sny,nr)
fpref="../run08/hFacS"
hfs=get_tiled_file(fpref,nsx,snx,nsy,sny,nr)

function w_z(u,v,dx,dy,dz,fw,fs,i,j,k)
        # Compute dw/dz for each cell in column
        vol=dx[i]*dy[j]*dz[k]
        uflxp=u[i+1,j,k]*dy[j]*fw[i+1,j,k]*dz[k]
        uflxm=u[i  ,j,k]*dy[j]*fw[i  ,j,k]*dz[k]
        vflxp=v[i,j+1,k]*dx[i]*fs[i,j+1,k]*dz[k]
        vflxm=v[i,j  ,k]*dx[i]*fs[i,j  ,k]*dz[k]
        rv=1. / vol
        wz=rv*(
                - uflxp + uflxm 
                - vflxp + vflxm 
               )
end
nx=nsx*snx;
ny=nsy*sny;
wzdz=[w_z(u0,v0,dx,dy,dz,hfw,hfs,i,j,k) for i in 1:nx-1 , j in 1:ny-1, k in 1:nr]

# Test function for checking formula
## wz=[sin(zc[k]/sum(dz)*2π) for i in 1:nx-1 , j in 1:ny-1, k in 1:nr]
## wzdz=[wz[i,j,k]*dz[k]*2π/sum(dz) for i in 1:nx-1 , j in 1:ny-1, k in 1:nr]

# Compute integral from bottom, then put k=1 back at the top!
# This needs a dz...
# w0=cat(cumsum(wzdz[:,:,end:-1:1],dims=3)[:,:,end:-1:1], wzdz[:,:,1]*0, dims=3)
# w0z0=w0[:,:,1];
w0=cumsum([wzdz[:,:,k].*dz[k] for k=nr:-1:1])[end:-1:1]
# w0z0=sum([wzdz[:,:,k].*dz[k] for k=nr:-1:1])
w0z0=w0[1]
Plots.heatmap(xc/1000,yc/1000,w0z0'*24*3600,aspect_ratio=:equal,xlabel="X (km)",ylabel="Y (km)",title=L"Surface w0 ($\frac{\mathrm{m}}{\mathrm{day}})$")
savefig("w0_surface_initial.png")

# Make some open boundary stats and plots
nbch=35
nbcz=47
nbct=1461
function read_obcs(fpref,nbch,nbcz,nbct)
        fn=fpref * @sprintf(".bin")
        phi=bswap.( reinterpret(Float64, read(fn,sizeof(Float64)*nbch*nbcz*nbct ) ) )
        phi=reshape(phi,(nbch,nbcz,nbct) )
        return phi
end
uobcs=read_obcs("OBSu",nbch,nbcz,nbct);
uobcn=read_obcs("OBNu",nbch,nbcz,nbct);
uobce=read_obcs("OBEu",nbch,nbcz,nbct);
uobcw=read_obcs("OBWu",nbch,nbcz,nbct);

vobcs=read_obcs("OBSv",nbch,nbcz,nbct);
vobcn=read_obcs("OBNv",nbch,nbcz,nbct);
vobce=read_obcs("OBEv",nbch,nbcz,nbct);
vobcw=read_obcs("OBWv",nbch,nbcz,nbct);

# make some plots
p1=Plots.plot(xc/1000,u0[1,:,1],title="u0 on Western Edge",label="")
p2=Plots.plot(xc/1000,uobcw[:,1,1],title="uw obcs @ t=0",label="") # west, level 1, time 1
p3=Plots.plot(xc/1000,u0[end,:,1],title="u0 on Eastern Edge",label="")
p4=Plots.plot(xc/1000,uobce[:,1,1],title="ue obcs @ t=0",label="") # east, level 1, time 1

p5=Plots.plot(xc/1000,u0[:,1,1],title="u0 on Southern Edge",label="")
p6=Plots.plot(xc/1000,uobcs[:,1,1],title="us obcs @ t=0",label="") # south, level 1, time 1
p7=Plots.plot(xc/1000,u0[:,end,1],title="u0 on Northern Edge",label="")
p8=Plots.plot(xc/1000,uobcn[:,1,1],title="un obcs @ t=0",label="") # north, level 1, time 1
pl=[p1,p2,p3,p4,p5,p6,p7,p8];
pgrid=Plots.plot(pl...,layout=(4,2))

pt=Plots.scatter(ones(3),annotation=(2,2,"U initial field and OBCS"),size=(200,100),
           marker=0,markeralpha=0,axis=false, grid=false, leg=false, label=false, border=:none)
Plots.plot(pt,pgrid,layout=grid(2,1,heights=[0.1,0.9]))
Plots.plot!(size=(800,1200))

savefig("u_initial_bcs.png")

p1=Plots.plot(xc/1000,v0[1,:,1],title="v0 (tangential) on Western Edge",label="")
p2=Plots.plot(xc/1000,vobcw[:,1,1],title="vw (tangential) obcs @ t=0",label="") # west, level 1, time 1
p3=Plots.plot(xc/1000,v0[end,:,1],title="v0 (tangential) on Eastern Edge",label="")
p4=Plots.plot(xc/1000,vobce[:,1,1],title="ve (tangential) obcs @ t=0",label="") # east, level 1, time 1

p5=Plots.plot(xc/1000,v0[:,1,1],title="v0 (normal) on Southern Edge",label="")
p6=Plots.plot(xc/1000,vobcs[:,1,1],title="vs (normal) obcs @ t=0",label="") # south, level 1, time 1
p7=Plots.plot(xc/1000,v0[:,end,1],title="v0 (normal) on Northern Edge",label="")
p8=Plots.plot(xc/1000,vobcn[:,1,1],title="vn (normal) obcs @ t=0",label="") # north, level 1, time 1
pl=[p1,p2,p3,p4,p5,p6,p7,p8];
pgrid=Plots.plot(pl...,layout=(4,2))

pt=Plots.scatter(ones(3),annotation=(2,2,"V initial field and OBCS"),size=(200,100),
           marker=0,markeralpha=0,axis=false, grid=false, leg=false, label=false, border=:none)
Plots.plot(pt,pgrid,layout=grid(2,1,heights=[0.1,0.9]))
Plots.plot!(size=(800,1200))

savefig("v_initial_bcs.png")

# Maps of vertically integrated and top layer surface horizontal flow
fubx(f)=0.5*(f[1:end-1,1:end-1].+f[2:end,1:end-1])
fvby(f)=0.5*(f[1:end-1,1:end-1].+f[1:end-1,2:end])
# Barotropic for U, V and cell center
dw=sum([hfw[:,:,k].*dz[k] for k=1:nr])
ds=sum([hfs[:,:,k].*dz[k] for k=1:nr])
dc=sum([hfc[:,:,k].*dz[k] for k=1:nr])
UB(f)=sum([f[:,:,k].*hfw[:,:,k].*dz[k] for k=1:nr])./dw
VB(f)=sum([f[:,:,k].*hfs[:,:,k].*dz[k] for k=1:nr])./ds
fbaro(f)=sum([f[:,:,k].hfc[:,:,k].*dz[k] for k=1:nr])./dc

U0=fubx(UB(u0))
V0=fvby(VB(v0))
# U0=fubx(u0[:,:,1])
# V0=fvby(v0[:,:,1])
f = GLMakie.Figure(resolution = (1600, 1600));GLMakie.Axis(f[1, 1],title="Vertically averaged flow, m/s");f
asize=7
lscale=100
GLMakie.arrows!(xc[1:end-1]/1000,yc[1:end-1]/1000,U0,V0,arrowsize=asize,lengthscale=lscale)
xref=minimum(xc[1:end-1]/1000)-10
yref=minimum(yc[1:end-1]/1000)-10
refarr=round(abs((maximum(U0[:])^2+maximum(V0[:])^2)^0.5),sigdigits=1)
GLMakie.arrows!([xref],[yref],[refarr],[0],arrowsize=asize,lengthscale=lscale)
lstr=latexstring(@sprintf("%2.1f",refarr))
GLMakie.text!(lstr, position = (xref, yref), align = (:left, :center), offset=(0,-20))
f
save("barotropic_flow_t_0.png",f)

# Now do top layer currents with and without barotropic component
usurf=fubx(u0[:,:,1])
vsurf=fvby(v0[:,:,1])
f = GLMakie.Figure(resolution = (1600, 1600));GLMakie.Axis(f[1, 1],title="Total surface layer flow, m/s");f
asize=7
lscale=50
GLMakie.arrows!(xc[1:end-1]/1000,yc[1:end-1]/1000,usurf,vsurf,arrowsize=asize,lengthscale=lscale)
xref=minimum(xc[1:end-1]/1000)-10
yref=minimum(yc[1:end-1]/1000)-10
refarr=round(abs((maximum(usurf[:])^2+maximum(vsurf[:])^2)^0.5),sigdigits=1)
GLMakie.arrows!([xref],[yref],[refarr],[0],arrowsize=asize,lengthscale=lscale)
lstr=latexstring(@sprintf("%2.1f",refarr))
GLMakie.text!(lstr, position = (xref, yref), align = (:left, :center), offset=(0,-20))
f
save("total_surface_layer_flow_t_0.png",f)

usurf=usurf-U0
vsurf=vsurf-V0
f = GLMakie.Figure(resolution = (1600, 1600));GLMakie.Axis(f[1, 1],title="Surface layer flow minus vertically averaged flow, m/s");f
asize=7
lscale=50
GLMakie.arrows!(xc[1:end-1]/1000,yc[1:end-1]/1000,usurf,vsurf,arrowsize=asize,lengthscale=lscale)
xref=minimum(xc[1:end-1]/1000)-10
yref=minimum(yc[1:end-1]/1000)-10
refarr=round(abs((maximum(usurf[:])^2+maximum(vsurf[:])^2)^0.5),sigdigits=1)
GLMakie.arrows!([xref],[yref],[refarr],[0],arrowsize=asize,lengthscale=lscale)
lstr=latexstring(@sprintf("%2.1f",refarr))
GLMakie.text!(lstr, position = (xref, yref), align = (:left, :center), offset=(0,-20))
f
save("surface_layer_flow_nobaro_t_0.png",f)
# GLMakie.destroy!(GLMakie.global_gl_screen())


# Let try correcting the barotropic divergence
#  First calculate the barotropic transport
txbaro=UB(u0).*dw
tybaro=VB(v0).*ds
psi=zeros(1:nx+1,1:ny+1)
psi[2:end,1].=-cumsum(tybaro[:,1].*dx)
psi[1:end-1,2:end].=cumsum([txbaro[i,j].*dy[j] for i=1:nx, j=1:ny],dims=2)
for j=2:size(psi,2)
    psi[:,j]=psi[:,j]+psi[:,1]
end
psi[end,2:end-1]=psi[end-1,2:end-1]-tybaro[end,2:end].*dx[2:end]
Plots.contourf(xf/1000,yf/1000, psi'/1e6,title="Barotropic stream function, Sv")
savefig("barotropic_stream_function.png")
# to print in xy orientation with x across screen, y up and 1,1 in bottome left
# psi[:,end:-1:1]'

# Compute boundary flow imbalance (net convergence, divergence)
# volume in at south, out at east, out at north, in at west
sv_add_tot=sum(tybaro[1:end-1,1].*dx[1:end-1])
ev_sub_tot=sum(txbaro[end,1:end-1].*dy[1:end-1])
nv_sub_tot=sum(tybaro[1:end-1,end].*dx[1:end-1])
wv_add_tot=sum(txbaro[1,1:end-1].*dy[1:end-1])
net_imbalance_Sv=(wv_add_tot+sv_add_tot-ev_sub_tot-nv_sub_tot)/1.e6
net_h_per_day=net_imbalance_Sv*1e6/(sum(dx[1:end-1])*sum(dy[1:end-1]))*24*3600

# Now lets try and compute a numerically non-divergent barotropic transport and flow 
# field 
# approx. barotropic divergence
# bd=txbaro[2:end,1:end-1]-txbaro[1:end-1,1:end-1]+tybaro[1:end-1,2:end]-tybaro[1:end-1,1:end-1]
#
# -psi_x is meridional transport
#  psi_y is zonal transport
nd_tx= (psi[1:end-1,2:end]-psi[1:end-1,1:end-1])./dy
nd_ty=-(psi[2:end,1:end-1]-psi[1:end-1,1:end-1])./dx
# see what vectors look like
txb_a=fubx(txbaro)
tyb_a=fvby(tybaro)
f = GLMakie.Figure(resolution = (1600, 1600));GLMakie.Axis(f[1, 1],title="Surface transport uncorrected, m^2/s");f
asize=7
lscale=5e-2
GLMakie.arrows!(xc[1:end-1]/1000,yc[1:end-1]/1000,txb_a,tyb_a,arrowsize=asize,lengthscale=lscale)
xref=minimum(xc[1:end-1]/1000)-10
yref=minimum(yc[1:end-1]/1000)-10
refarr=round(abs((maximum(txb_a[:])^2+maximum(tyb_a[:])^2)^0.5),sigdigits=1)
GLMakie.arrows!([xref],[yref],[refarr],[0],arrowsize=asize,lengthscale=lscale)
lstr=latexstring(@sprintf("%2.1f",refarr))
GLMakie.text!(lstr, position = (xref, yref), align = (:left, :center), offset=(0,-20))
f
save("surface_transport_uncorrected.png",f)

txb_a=fubx(nd_tx)
tyb_a=fvby(nd_ty)
f = GLMakie.Figure(resolution = (1600, 1600));GLMakie.Axis(f[1, 1],title="Surface transport corrected, m^2/s");f
asize=7
lscale=5e-2
GLMakie.arrows!(xc[1:end-1]/1000,yc[1:end-1]/1000,txb_a,tyb_a,arrowsize=asize,lengthscale=lscale)
xref=minimum(xc[1:end-1]/1000)-10
yref=minimum(yc[1:end-1]/1000)-10
refarr=round(abs((maximum(txb_a[:])^2+maximum(tyb_a[:])^2)^0.5),sigdigits=1)
GLMakie.arrows!([xref],[yref],[refarr],[0],arrowsize=asize,lengthscale=lscale)
lstr=latexstring(@sprintf("%2.1f",refarr))
GLMakie.text!(lstr, position = (xref, yref), align = (:left, :center), offset=(0,-20))
f
save("surface_transport_corrected.png",f)

txb_a=fubx(nd_tx-txbaro)
tyb_a=fvby(nd_ty-tybaro)
f = GLMakie.Figure(resolution = (1600, 1600));GLMakie.Axis(f[1, 1],title="Surface transport corrections, m^2/s");f
asize=7
lscale=5e-2
GLMakie.arrows!(xc[1:end-1]/1000,yc[1:end-1]/1000,txb_a,tyb_a,arrowsize=asize,lengthscale=lscale)
xref=minimum(xc[1:end-1]/1000)-10
yref=minimum(yc[1:end-1]/1000)-10
refarr=round(abs((maximum(txb_a[:])^2+maximum(tyb_a[:])^2)^0.5),sigdigits=1)
GLMakie.arrows!([xref],[yref],[refarr],[0],arrowsize=asize,lengthscale=lscale)
lstr=latexstring(@sprintf("%2.1f",refarr))
GLMakie.text!(lstr, position = (xref, yref), align = (:left, :center), offset=(0,-20))
f
save("surface_transport_correction.png",f)

GLMakie.destroy!(GLMakie.global_gl_screen())

bd=txbaro[2:end,1:end-1]-txbaro[1:end-1,1:end-1]+tybaro[1:end-1,2:end]-tybaro[1:end-1,1:end-1];
p1=Plots.heatmap(bd[:,end:-1:1]'*24*3600/4000,aspect_ratio=:equal)
txb=nd_tx;
tyb=nd_ty;
bd=(txb[2:end,1:end-1]-txb[1:end-1,1:end-1])./dx[1:end-1]+(tyb[1:end-1,2:end]-tyb[1:end-1,1:end-1])./dy[1:end-1];
p2=Plots.heatmap(bd[:,end:-1:1]'*24*3600/4000,aspect_ratio=:equal)
Plots.plot(p1,p2)
savefig("divergences.png")
