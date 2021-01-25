;....basic version of MIRO for simple simulations of Farday Rotation Measure 
;...derived from Bonafede, Vazza et al. 2013 MNRAS

pro miro,input,dir_out=dir_out,res=res,n0=n0,rc=rc,bmean=bmean,kin=kin,kmax=kmax,alphak=alphak,seed=seed,name=name,ncentre=ncentre,alphaB=alphaB,nclump=nclump,write_disk=write_disk,cc_ncc=cc_ncc,t500=t500,temp=temp 
 if not IS_DEF(cc_ncc) then cc_ncc='ncc'
  t0=systime(1)
  print,'>>STARTING MIRO<<'
  CPU,TPOOL_NTHREADS=8

   if not IS_DEF(dir_out)   then dir_out='./'   ;...folder for outputs
   if not IS_DEF(res)       then res=2 ;...cell resolution in kpc
   if not IS_DEF(n0)        then n0=128 ;..grid size
   if not IS_DEF(kin)       then  kin=10     ;...minimum k for the power-law [0:n0*0.5-1]
   if not IS_DEF(kmax)      then  kmax=(n0*0.5)-1   ;...maximum k  (n0*0.5-1 is the Nyquist frequency)
   if not IS_DEF(seed)      then seed=systime(1) ;..random seed number
   if not IS_DEF(alphak)    then  alphak=1.6666 ;...spectral index for 1D power specrtra (1.666 is for kolmogorov)
   if not IS_DEF(bmean)     then  bmean=1 ;...[muG] wanted normalisation of the rms B-field within the volume
   if not IS_DEF(ncentre)   then ncentre=3.3e-3 ;...[part/cm^3] particle density in the cluster core
   if not IS_DEF(ncentre)   then rc=290 ;...[kpc] core radius of cluster
   if not IS_DEF(name)      then name="test" ;...identifier for this run
   if not IS_DEF(alphaB)    then alphaB=0.5  ;...assumed exponent of the B(n) \propto n^alfaB scaling
   if not IS_DEF(nclump)    then nclump=0  ;..number of clumps (randomly distributed
   if not IS_DEF(write_disk) then write_disk=0   ;...if write_disk=1 write all generated 3D grids on disk 
   if not IS_DEF(cc_ncc) then cc_ncc='ncc' ;...
     
   ;...3D density field  -> initialised empty density field. 
  d0=fltarr(n0,n0,n0)
  d0(*,*,*)=1 ;...values assumed to be in 1/cm^3
  temp=fltarr(n0,n0,n0)
  temp(*,*,*)=1.
                   ;a=build_coma(d0,rc,beta,ncentre,n0,res,centre)  ;....this creates a beta-model profile, otherwise uniform density
 
  a=build_cluster(d0,rc,beta,n0,ncentre,res,centre,cc_ncc,r500,t500,temp)

  print,"galaxy cluster model created"
 


   if nclump gt 0 then begin  ;....this adds nclump clumps of gas in the cluster atmosphere
   a=add_clumps(nclump,n0,dclump,seed,res)
   idd=where(dclump gt 0.1)
   d0(idd)*=(1+dclump(idd))
   dclump=0
   print,"gas clumps created"
   endif

 
   t0=systime(1)
   a=badd(n0,d0,bx,by,bz,alphak,kin,kmax,seed,alphaB)  ;...generation of 3-D magnetic fields
   print,'3D magnetic field generated'  
   ;print,systime(1)-t0
   a=normb(bx,by,bz,bmean,n0)   ;..magnetic fields are renormalized to physical units
   
print,'3D magnetic field renormalized'  
print,'min/max Bx',min(bx),max(bx)
print,'min/max By',min(by),max(by)
print,'min/max Bz',min(bz),max(bz)

;....computing the 3D power spectrum (just one component for brevity)

pk3d=powerspectrum3d(bx)
nk=size(pk3d)
nyq=nk(1)
kk=indgen(nk(1))+1
set_plot,'ps'
device,filename=dir_out+'statistics_'+name+'.eps',xsize=16,ysize=14,/color,font_size=7,/encapsul
loadct,13
!p.font=1
!p.multi=[0,2,2]

plot,kk,pk3d*kk^2.,/xlog,/ylog,yrange=[1e-10,1],title='Magnetic power spectrum',ytitle='E!dB!n(k)',xtitle='k'
k1=kin
k2=kmax-1
pk1=alog10(pk3d(k1)*kk(k1)^2.)
plots,kk(k1),10.^(pk1)
plots,kk(k2),10.^(pk1-alphak*alog10(kk(k2)/float(kk(k1)))),/continue,linestyle=2,thick=5,col=80


;.........



print,'generating outputs'

;....various outputs

   ima=fltarr(n0,n0,3) ;...map making

 for c=0,n0-1 do begin ;...projection along z-axis
 ima(*,*,0)+=d0(*,*,c)  ;...projected density
 ima(*,*,1)+=sqrt(bx(*,*,c)^2.+by(*,*,c)^2.+bz(*,*,c)^2.)/float(n0)   ;...projected volume-weighted mean B
 ima(*,*,2)+=bz(*,*,c)*d0(*,*,c)*812.*res     ;....Rotation Measure for each pixel
 endfor
  d0=0
 
contour,(ima(*,*,2)),nlevels=128,title='Rotation Measure',xrange=[0,n0-1],yrange=[0,n0-1],/xstyle,/ystyle,/fill

minb=-2*bmean
maxb=2*bmean
nb=100
db=(maxb-minb)/float(nb)
h1=histogram((bx),min=minb,max=maxb,nbins=nb)
h2=histogram((by),min=minb,max=maxb,nbins=nb)
h3=histogram((bz),min=minb,max=maxb,nbins=nb)
xb=minb+db*indgen(nb)
plot,xb,h1,linestyle=0,xtitle='B-field[!9m!1!G]',ytitle='dN/dB',title='distribution functions for 3 components',yrange=[1,max(h1)],/ystyle,psym=10
oplot,xb,h2,linestyle=1,psym=10
oplot,xb,h3,linestyle=5,psym=10

plot,bx(*,n0*0.5-1,n0*0.5-1),linestyle=0,xtitle='x-axis',ytitle='B-field[!9m!1!G]',title='1-D trend of 3 components',xrange=[0,n0-1],/xstyle
 oplot,by(*,n0*0.5-1,n0*0.5-1),linestyle=1
 oplot,bz(*,n0*0.5-1,n0*0.5-1),linestyle=5
device,/close
device,/encapsul

  bx=0
  by=0
  bz=0
  
;...writing maps on disk

tag=string(n0,'(i3)')+'_'+string(kin,'(i3)')+'_'+string(kmax,'(i3)') ;...name label that records spectral information, for multiple outputs

print,'generating .fits file'
 writefits,dir_out+'maps_'+tag+name+'.fits',ima   ;...fits datacbue

print,'computing RM profile with percentile distribution'
nn=size(im)
dist=fltarr(n0,n0)
  for  y=0,n0-1 do begin
  for  x=0,n0-1 do begin
    dist(x,y)= sqrt( (x-(n0*0.5-0.5))^2.+(y-(n0*0.5-0.5))^2.)
  endfor
  endfor
  irm=reform(ima(*,*,2))
  dist=uint(dist)
  mad=max(dist)
  rm_prof=fltarr(mad,9)

  for r=0L,mad-1 do begin
   ri=where(dist gt r and dist lt r+2 ,nr)
   per=percentiles(abs(irm(ri)),value=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
   rm_prof(r,0:8)=per(0:8)
  endfor

 set_plot,'ps'
 !p.multi=0
 !p.font=1
 device,filename=dir_out+"RM_sim_vs_obs"+name+".eps",/encapsul,/color,xsize=15,ysize=13,font_size=10
 loadct,13
 ro=res*(0.5+indgen(max(dist)))
 plot,ro,rm_prof(*,0),xrange=[0,mad*res],yrange=[0.1,1e4],/ylog,/ystyle,/xstyle,/nodata,xtitle='r[kpc]',ytitle='<|RM|>[rad/m!e2!n]',title='Observed(points) vs simulated (colors)'

   x1=1
   x2=mad-1

   loadct,1
   l1=[ro,ro(x2-1-indgen(x2-x1))]
   l2=[rm_prof(x1:x2,0),rm_prof(x2-1-indgen(x2-x1),8)]
   l3=[rm_prof(x1:x2,1),rm_prof(x2-1-indgen(x2-x1),7)]
   l4=[rm_prof(x1:x2,2),rm_prof(x2-1-indgen(x2-x1),6)]
   l5=[rm_prof(x1:x2,3),rm_prof(x2-1-indgen(x2-x1),5)]
   l6=rm_prof(x1:x2,4)

     polyfill,l1,l2,col=200,/fill,noclip=0
     polyfill,l1,l3,col=150,/fill,noclip=0
     polyfill,l1,l4,col=100,/fill,noclip=0
     polyfill,l1,l5,col=50,/fill,noclip=0

     oplot,l1,l6,col=250,thick=6,linestyle=2
   

;...overplotting real RM data of COMA
 r_coma=[51,124,372,532,919,1250,1489]
 rm_coma=[-256,-120,372,51,21,6,12]
 erm_coma=[50,22,51,4,30,12,27]
 ;  erm_coma=[50,22,51,4,20.5,5.5,11.5]
 disp_coma=[303,166,154,16,65,56,37]
 beam_coma=[35,56,10,16,7,33,4]
 res_beam=1 ;...kpc
 beam_equiv=sqrt(beam_coma)*res_beam ;..

 r_coma2=[2451,2543,2075,1982,1824,1667,1113]
 rm_coma2=[4.8,2.8,13,15,18,-22,69]
 erm_coma2=[0.7,2,4,4,5,4,9]

 loadct,13

 errplot,r_coma,abs(rm_coma)-erm_coma,abs(rm_coma)+erm_coma,col=250,thick=5;,psym=4;,col=70,psym=6,symsize=1.1,thick=4
 errplot,r_coma2,abs(rm_coma2)-erm_coma2,abs(rm_coma2)+erm_coma2,col=230,thick=5;,psym=1;,col=70,psym=6,symsize=1.1,thick=4

 oplot,r_coma,abs(rm_coma),col=250,thick=5,psym=4,symsize=1.5;,col=70,psym=6,symsize=1.1,thick=4
 oplot,r_coma2,abs(rm_coma2),col=230,thick=5,psym=1,symsize=1.5;,col=70,psym=6,symsize=1.1,thick=4


 device,/close
 device,/encapsul
 
 ; print,systime(1)-t0
  
if write_disk eq 1 then begin
print,'writing datacubes on disk'
    openw,5,dir_out+'3D_fields_'+tag+name+'.bin'   ;....binary output of 3D B-components and gas density
     writeu,5,bx
     writeu,5,by
     writeu,5,bz
     writeu,5,d0
    close,5
    endif
      
  end



function Pk,k,index=index
    p=1.*(k+1)^(index)
return, p
end

function badd,n0,d0,bx,by,bz,alphak,kin,kmax,seed,alphaB
;...routine that creates the 3D model of magnetic field

nng=size(d0)
ng1=nng(1)
ng2=nng(2)
ng3=nng(3)

 ng=ng3      
  ik2=-(alphak+2.)       ;...spectral index in A space, n in real space n=ik-2
  pok=fltarr(3*ng)      ;....1D power spectrum
  pok(*)=1e-10          ;...normalization of the spectrum
  k=1+indgen(3*ng)      ;...array of wavenumbers


  ;this creates the PS
  for kk=kin,kmax-1 do begin
     pok(kk)=Pk(kk,index=ik2)
  endfor
  
  

;....array of k-distances from the centre of the box is built here
;....this procedure is done instead of a more trivial 3D loop, and gives a speedup of a factor 25%

  krx=indgen(ng1)-ng1*0.5+1
  kry=indgen(ng2)-ng2*0.5+1
  krz=indgen(ng3)-ng3*0.5+1
  krx=reform(krx,ng1,1,1)
  kry=reform(kry,1,ng2,1)
  krz=reform(krz,1,1,ng3)
  kra=rebin(krx,ng1,ng2,ng3)^2
  krx=0
  kra=temporary(kra)+rebin(kry,ng1,ng2,ng3)^2                 ; conserve memory
  kry=0
  kra=temporary(kra)+rebin(krz,ng1,ng2,ng3)^2
  krz=0
  kra=fix(temporary(sqrt(kra)))
 
 
 ;extract random numbers for amplitude and phase
  ;amplitude from raylegh, phase uniform between 0 and 2 pi
  ; we extract amplitude and 3 angles (polar coordinates) 
  ; the components of A are the projection on the 3 axis
 
 a1_seed=83243.*seed
 a2_seed=22342.*seed
 a3_seed=22973.*seed
 p_seed=34213.*seed


 x_a=randomu(a1_seed,ng,ng,ng,/uniform)                  ;.....array of random numbers to extract the amplitudes
 phi2=randomu(a2_seed,ng,ng,ng,/uniform)*!pi    ;.....array of phase from 0 and pi
 phi1=randomu(p_seed,ng,ng,ng,/uniform)*2.*!pi ;.....array of phase from 0 and 2pi
 phi3=randomu(a3_seed,ng,ng,ng,/uniform)*2.*!pi ;.....array of phase from 0 and 2pi
 
 
 ;....arrays of complex B-field components
 
 Bkz=complexarr(ng1,ng2,ng3)
 Bky=complexarr(ng1,ng2,ng3)
 Bkx=complexarr(ng1,ng2,ng3)
 
  for k=kin,kmax do begin
  krk=where(kra eq k,nr)
  kk=array_indices(kra,krk)
  for i=0L,nr-1 do begin
  kx=kk(0,i)
  ky=kk(1,i)
  kz=kk(2,i)
 
   am=rayleigh(sqrt(pok(k)),float(x_a(kx,ky,kz))) ;extract amplitude from rayleigh ;distr with sigma=pow(kA)
 
    
      ph2=phi2(kx,ky,kz)
      ph3=phi3(kx,ky,kz)  
      ph1=phi1(kx,ky,kz)  
 
        c2=cos(ph2)
        s2=sin(ph2)
        c3=cos(ph3)
        s3=sin(ph3)
        c1=cos(ph1)
        s1=sin(ph1)
                ampx=abs(am*s2*c3) ; x component of Ak (modulus)
                ampy=abs(am*s2*s3) ; y component of Ak (modulus)
                ampz=abs(am*c2)                     ; z component of Ak (modulus)

                Akx=ampx*complex(c1,s1)
                Aky=ampy*complex(c1,s1)
                Akz=ampz*complex(c1,s1)

                k_vec=[kx,ky,kz]
                A_vec= [Akx,Aky,Akz]
   
                rot1=cross(real_part(k_vec),real_part(A_vec))
                rot2=cross(imaginary(k_vec),imaginary(A_vec))

                Bkx(kx,ky,kz)=complex(rot1(0),rot2(0))  ; B component z direction in fourier space
                Bky(kx,ky,kz)=complex(rot1(1),rot2(1))
                Bkz(kx,ky,kz)=complex(rot1(2),rot2(2))
      endfor
      endfor
    

 phi=0
 phi2=0
 x_a=0
 y_a=0
 z_a=0
 ka=0


  Bkz=temporary(shift(Bkz,-ng/2.+1,-ng/2.+1,-ng/2.+1))
  Bz=float(fft(real_part(Bkz),1))   ; Bz in the real space
  bkz=0
  Bky=temporary(shift(Bky,-ng/2.+1,-ng/2.+1,-ng/2.+1))
  By=float(fft(real_part(Bky),1))   ; By in the real space
  bky=0
  Bkx=temporary(shift(Bkx,-ng/2.+1,-ng/2.+1,-ng/2.+1))
  Bx=float(fft(real_part(Bkx),1))   ; Bx in the real space
  bkx=0
 
 ;...we impose a scaling with the gas density field 
   bx=temporary(bx)*d0^(alphaB)
   by=temporary(by)*d0^(alphaB)
   bz=temporary(bz)*d0^(alphaB)
   
 
  end
  
  
;returns a random variable extracted from the rayleigh distribution
;with parameter sigma
function rayleigh,sigma,x
    y=sigma*sqrt(-2.*alog(float(x)))
return,y
end


function normB,bx,by,bz,bmean,n0

iw1=where(finite(bx) ne 1)
iw2=where(finite(bx) eq 1)
bx(iw1)=min(bx(iw2))
by(iw1)=min(by(iw2))
bz(iw1)=min(bz(iw2))

totmean=(n0^3.*bmean^2.)/float(3.)

tbx=total(bx^2.,/double)
tby=total(bx^2.,/double)
tbz=total(bx^2.,/double)

bx=temporary(bx)*sqrt(totmean/float(tbx))
by=temporary(by)*sqrt(totmean/float(tby))
bz=temporary(bz)*sqrt(totmean/float(tbz))

end


function cross,vector1,vector2,result
result=[0.0,0.0,0.0]
     result (0) = vector1 (1) * vector2 (2) - vector1 (2) * vector2 (1)
     result (1) = - (vector1 (0) * vector2 (2) - vector1 (2) * vector2 (0))
     result (2) = vector1 (0) * vector2 (1) - vector1 (1) * vector2 (0)
return,result
end

function build_coma,d0,rc,beta,ncentre,n0,res,centre
  nn=size(d0)
  ng=nn(1)
  IF not IS_DEF(rc) then rc=290.  ;...core radius in kpc
  if not IS_DEF(beta) then beta=0.75   ;...beta parameter of the beta profile
  if not IS_DEF(n0) then ncentre=3.44e-3  ;...central density in 1/cm^3
  if not IS_DEF(ng) then ng=256 
  if not IS_DEF(centre) then centre=[0.5,0.5,0.5]
  
  ic=centre(0)*ng-1
  jc=centre(1)*ng-1
  lc=centre(2)*ng-1

           ; resolution in kpc (same of the sim)
  for l=0,ng-1 do begin
     for j=0,ng-1 do begin
        for i=0,ng-1 do begin
           r=sqrt((l-lc)^2.+(j-jc)^2.+(i-ic)^2.)*res
           d0(i,j,l)=ncentre*(1.+r^2./float(rc^2.))^(-1.5*beta)
        endfor 
     endfor
  endfor
 
 end
 
 
;...function to create a specific galaxy cluster atmosphere

function build_cluster,d0,rc,beta,n0,ncentre,res,centre,cc_ncc,r500,t500,temp  
  nn=size(d0)
  ng=nn(1)
  IF not IS_DEF(rc) then rc=290.  ;...core radius in kpc
  if not IS_DEF(beta) then beta=0.75   ;...beta parameter of the beta profile
  if not IS_DEF(ncentre) then ncentre=3.44e-3  ;...central density in 1/cm^3
  if not IS_DEF(ng) then ng=256
  if not IS_DEF(centre) then centre=[0.5,0.5,0.5]
  if not IS_DEF(gamm) then gamm=3.
  if not IS_DEF(r500) then r500=rc*5.
  if not IS_DEF(t500) then t500=5. ;...keV

  ic=centre(0)*ng-1
  jc=centre(1)*ng-1
  lc=centre(2)*ng-1

  ;.....best-fit formulas from Ghirardini+2018, for average CC or NCC clusters

  if cc_ncc eq 'ncc' then begin
    ncentre=exp(-4.9)
    rc=exp(-2.7)*r500
    rs=exp(-0.51)*r500
    alf=0.7
    bet=0.39
    eps=2.6
    t0=1.09
    logr=exp(-4.4)*r500
    rt=0.45*r500
    tm0=0.66
    acool=1.33
    c2=0.3
  endif

  if cc_ncc eq 'cc' then begin
    ncentre=exp(-3.9)
    rc=exp(-3.2)*r500
    rs=exp(0.17)*r500
    alf=0.8
    bet=0.49
    eps=4.67
    t0=1.32
    logr=exp(-2.8)*r500
    rt=0.4*r500
    tm0=0.22
    acool=0.74
    c2=0.33
  endif


  ; resolution in kpc (same of the sim)
  for l=0,ng-1 do begin
    for j=0,ng-1 do begin
      for i=0,ng-1 do begin
        r=sqrt((l-lc-0.5)^2.+(j-jc-0.5)^2.+(i-ic-0.5)^2.)*res
        ;        d0(i,j,l)=ncentre*(1.+r^2./float(rc^2.))^(-1.5*beta)
        d0(i,j,l)=sqrt(ncentre^2.*((r/float(rc))^(-alf))/float((1.+r^2./float(rc^2.))^(3.*bet-alf*0.5))*1./float(1.+(r^gamm)/float(rs^gamm))^(eps/float(gamm)))
        temp(i,j,l)=t0*t500*((tm0+(r/float(logr))^acool)/float(1.+(r/float(logr))^acool))*1/float(1.+r^2./float(rt^2.))^(0.5*c2)
      endfor
    endfor
  endfor

 ;...screen plotting of density profile
  prof=prof_radial(d0)
  proft=prof_radial(temp)
  nr=size(prof)
  rr=0.5+indgen(nr(1))
  set_plot,'x';
  !p.multi=[0,2,0]
  plot,rr,prof,/ylog,/xlog,xrange=[1.,1e2],charsize=2.2,title=cc_ncc
  plot,rr,proft,/ylog,/xlog,yrange=[1e7,1e8],xrange=[1.,1e2],charsize=2.2,title=cc_ncc



end

 
function build_fila,d0,rc,beta,ncentre,n0,res,centre
  nn=size(d0)
  ng=nn(1)
  IF not IS_DEF(rc) then rc=290.  ;...core radius in kpc
  if not IS_DEF(beta) then beta=0.75   ;...beta parameter of the beta profile
  if not IS_DEF(n0) then ncentre=3.44e-3  ;...central density in 1/cm^3
  if not IS_DEF(ng) then ng=256
  if not IS_DEF(centre) then centre=[0.5,0.5,0.5]

  ic=centre(0)*ng-1
  jc=centre(1)*ng-1
  lc=centre(2)*ng-1

  ; resolution in kpc (same of the sim)
  for l=0,ng-1 do begin
    for j=0,ng-1 do begin
      for i=0,ng-1 do begin
        r=sqrt((l-lc)^2.+(i-ic)^2.)*res
        d0(i,j,l)=ncentre*(1.+r^2./float(rc^2.))^(-1.5*beta)
      endfor
    endfor
  endfor

end



function add_clumps,nclump,ng,dclump,seed,res


  xc=10+(2*ng*0.5-20)*randomu(121*seed,nclump) ;...random extraction of points and densities
  yc=10+(2*ng*0.5-20)*randomu(566*seed,nclump)
  zc=10+(2*ng*0.5-20)*randomu(67*seed,nclump)


  ddc=100*randomu(141*seed,nclump,/poisson) ;...overdensity with respect to the cluster gas model
  rrc=1+10*randomu(31*seed,nclump)
  dclump=fltarr(ng,ng,ng)
  dclump(*,*,*)=0.
  for c=0,nclump-1 do begin
  
    irr=uint(rrc(c))
    ;....each clump is assumed to have a density profile decreasing as 1/r^2 from its centre
    for l=0,2*irr do begin
      for j=0,2*irr do begin
        for i=0,2*irr do begin
          dist=sqrt((i-irr-0.5)^2.+(j-irr-0.5)^2.+(l-irr-0.5)^2.)
          dclump(xc(c)-irr+i,yc(c)-irr+j,zc(c)-irr+l)=ddc(c)/float(1+dist)^2.
        endfor
      endfor
    endfor
  endfor
end





  ;function that computes the profile of data in the file dat and its mean. If r is
  ;supplied, the profile and mean is computed within a radius r_max (given in number of cells)
function prof_radial,dat ;,r_max=r_max
  nn=size(dat)
  xc=[nn(1)*0.5-1,nn(2)*0.5-1,nn(3)*0.5-1]
  ng=nn(1)
  r_max=ng*0.5
  i3=where(dat eq max(dat))

  ic=array_indices(dat,i3)

  rx=indgen(ng)-ic(0)
  ry=indgen(ng)-ic(1)
  if N_elements(ic) EQ 3 then rz=indgen(ng)-ic(2)
  rx=reform(rx,ng,1,1)
  ry=reform(ry,1,ng,1)
  if N_elements(ic) EQ 3 then rz=reform(rz,1,1,ng)
  ra=rebin(rx,ng,ng,ng)^2
  ra=temporary(ra)+rebin(ry,ng,ng,ng)^2                 ; conserve memory
  if N_elements(ic) EQ 3 then ra=temporary(ra)+rebin(rz,ng,ng,ng)^2
  ra=fix(temporary(sqrt(ra)))

  prof0=fltarr(max(ra))
  rms0=fltarr(max(ra))
  if not is_def(r_max) then r_max=max(ra)

  for r=0,r_max-1 do begin

    rk=where(ra ge r and ra lt r+1,nr)

    prof0(r)=mean(dat(rk))     ;...radial profile of field "dat"
    rms0(r)=stddev(dat(rk))
    ;print, 'valori medi e stddev al raggio, ',r,' ',prof0(r),rms0(r)
  endfor

  rm=where(ra le r_max)
  mean_field=mean(dat(rm))

  print, 'mean of field within ', r_max,' ',mean_field

  return,[[prof0],[rms0]]
end

