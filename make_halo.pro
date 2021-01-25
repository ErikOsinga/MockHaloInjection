;modelimage= clean components - just to set header etc...
;pixel=arcsec in one pixel
;scale= arcsec to kpc scale
;I0= central intensity of the halo one whants to inject. The halo is
;modelled as  I(r)=I0exp(-r/re) I0 is given in Jy/arcsec^2, re is given
;in kpc
;xc and yc are the centre of the halo (pixel) default: image centre
;scatter is the fraction of I to add/subtract pixel by pixel to the
;analytical profile default: 0.2

pro go_halo_analytic, modelimage=modelimage,pixel=pixel,I0=I0,re=re,scale=scale,xc=xc,yc=yc,scatter=scatter


a=readfits(modelimage,head_halo)
;pixel=3 ;pixel to arcesc
;scale =1.922 ; arcsec to kpc

ss_a=size(a)
nx=ss_a(1)
ny=ss_a(2)


IF NOT(IS_DEF(xc)) THEN xc=nx/2.
IF NOT(IS_DEF(yc)) THEN yc=ny/2.



random_image=randomn(seed,nx,ny)


random_sign=round(randomu(seed,nx,ny)*10)

I0=I0*pixel^2 ; I0 in Jy/arcsec^2
re=re/scale/pixel ;re in kpc and then / scale/pixel

;re=42/pixel
b=a
b(*,*)=1e-10

;exponential profile 
for yy=0, ny-1 do begin
   for xx=0, nx-1 do begin
      rr=sqrt((xx-xc)^2.+(yy-yc)^2.)
      I=I0*exp(-1.*rr/re)
      
      ;scatter
     IF (IS_DEF(scatter)) THEN begin
        deltaI=random_image(xx,yy)*scatter*I+I
                                ;if deltaI GE I and rr LE 3.*re  then print, deltaI, I
                                ;add or subtract scatter
                                ;if random_sign(xx,yy) mod 2 EQ 0 then b(xx,yy)=I+abs(deltaI)
                                ;if random_sign(xx,yy) mod 2 EQ 1 then b(xx,yy)=I-abs(deltaI)
      
        if deltaI LE 0 and rr LE 5* re then print, 'warning, negative flux within 3 re, change scatter!!'
        if rr le 5*re then b(xx,yy)=deltaI
     endif

     IF NOT(IS_DEF(scatter)) THEN begin
         b(xx,yy)=I
         
        endif

  endfor
endfor



print,'total flux in mJy', total(b)*1e3

outputname=strsplit(modelimage,'.fits',/regex,/extract)+'halo_I0='+string(I0*1e6/pixel^2,format='(f3.1)')+'.fits'
print, 'outputname' ,strsplit(modelimage,'.fits',/regex,/extract)+'halo_I0='+string(I0*1e6/pixel^2,format='(f3.1)')+'.fits'

writefits,outputname,b,head_halo

end








;takes the halo generated above (analytic halo) and a map of kolmogorov flutuation generated with MIRO,
;1)  matches the 2 maps by nesting the miro map into a larger map 
;2) computes the normalization parameters for each radius to have that
;the average profile of the miro map matches the analytic halo map 
;3) multiplies the normalization parameter by the miro map. 
pro fake_halo_ps,analytic=analytic,miro=miro,scale=scale, pixel=pixel, xc=xc,yc=yc,re=re,outfile=outfile


; convert e-folding to radius in pixels
re=re/scale/pixel 

a=readfits(analytic,head_a)
m_r=readfits(miro)

m=m_r(*,*,1) ;slice 1 of the cube made by miro has the magnetic field

s_a=size(a)
s_b=size(m)

nx_a=s_a(1)
ny_a=s_a(2)
nx_b=s_b(1)
ny_b=s_b(2)

;1 nesting of the miro image into a larger one with dimensions equal
;to the analytic image


m_new=fltarr(nx_a,ny_a)
m_new(*,*)=1e-8


;nesting: excluding the 4 pixels on the edges of the miro image
;because they may have numerical artefacts
edge=1e-8
m(0,0)=edge
m(0,ny_b-1)=edge
m(nx_b-1,0)=edge
m(nx_b-1,ny_b-1)=edge

m_new(xc-nx_b/2:xc+nx_b/2-1,yc-ny_b/2:yc+ny_b/2-1)=m(0:nx_b-1,0:ny_b-1)


writefits,'test_nesting.fits',m_new,head_a

tot_m=fltarr(nx_a)
tot_a=fltarr(nx_a)
tot_a(*)=0.
tot_m(*)=0.
count=fltarr(nx_a)
count(*)=1.
m_norm=fltarr(nx_a,ny_a)

;2 Compute normalization constants radius by radius
for yy=0, ny_a-1 do begin
   for xx=0, nx_a-1 do begin
      rr=sqrt((xx-xc)^2.+(yy-yc)^2.)
; ADDED: check whether we don't go outside of boundaries
      IF (rr LE nx_a) THEN BEGIN 
        tot_m(rr)+=m_new(xx,yy)
        tot_a(rr)+=a(xx,yy)
        count(rr)++
      ENDIF
   end
end

cc=where(count GT 0)
tot_a_n=tot_a(cc)
tot_m_n=tot_m(cc)
count_n=count(cc)

avg_a=tot_a_n/count_n
avg_m=tot_m_n/count_n

norm=avg_a/avg_m


;3 multiply by norm constant
for yy=0, ny_a-1 do begin
   for xx=0, nx_a-1 do begin
      rr=sqrt((xx-xc)^2.+(yy-yc)^2.)
      
      if rr LE 5*re then m_norm(xx,yy)=m_new(xx,yy)*norm(rr) else  m_norm(xx,yy)=a(xx,yy)
      ;if rr LE 5*re then print, m_new(xx,yy),m_norm(xx,yy),a(xx,yy)
   end
end

writefits,outfile,m_norm,head_a
print, 'image written in ',outfile
end







pro example


go_halo_analytic,modelimage='~/Documents/PR_relics/UV_data/SIM/ZWCL008/FAKE_HALOS/zwcl0008_uvsub_R1.model_reference.fits',pixel=3.,scale=1.922,I0=0.45e-6,re=094.,xc=1183,yc=1102


;inputs are: analytic= the map with the analytic profile, miro= the
;map with the power spectrum fluctuations, pixel=pixelsize in arcsec,
;scale= kpc/arcsec, re is the r_e parameter of the analytic profile,
;xc and yc are the pixel cooridnates in the map of the position where
;you want the halo to be centered
fake_halo_ps, analytic='~/Documents/PR_relics/UV_data/SIM/ZWCL008/FAKE_HALOS/zwcl0008_uvsub_R1.model_referencehalo_I0=0.4.fits',miro='~/Documents/PR_relics/UV_data/SIM/PROVA_HALO_MIRO/ZWCL0008/maps_256_6_127.fits',pixel=3.,scale=1.922,re=094.,xc=1183,yc=1102,outfile='~/Documents/PR_relics/UV_data/SIM/ZWCL008/FAKE_HALOS/zwcl0008_fake_I00.45_re94_miro2.fits'

end


pro erik

fake_halo_ps,analytic='erik_map.fits',miro='./OUT/map_miro.fits',pixel=1.5,scale=4.726,re=065.,xc=679,yc=502,outfile='./test.fits'
end



pro WHL1

; FOR THE FIRST WHL cluster WHLJ160425.6+543307.fits
; ON CLUSTER COORDS  xc=679,yc=502
; OFF CLUSTER COORDS xc=546,yc=506

anim1='/net/achterrijn/data2/osinga/cluster_pointsource_sub/WHLJ160425.6+543307/WHLJ160425.6+543307_maskROBUST-0.5_TEMP-000'
anim2='-model.fits'
outim1='/net/reusel/data1/osinga/phd/year1/deepfields/power_spectrum_halos/OUT/WHLJ160425.6+543307-000'
outim2='-model.fits'

; ADD PS FLUCT TO MODEL 0 to 6
for modnum=0, 5 do begin
  anim=anim1+STRTRIM(modnum,1)+anim2
  outim=outim1+STRTRIM(modnum,1)+outim2
  print, 'Doing ',outim
  fake_halo_ps, analytic=anim,miro='/net/reusel/data1/osinga/phd/year1/deepfields/power_spectrum_halos/OUT/maps_128_5_63.fits',pixel=1.5,scale=4.086,re=065.,xc=546,yc=506,outfile=outim
end


end


pro WHL2

; FOR THE SECOND WHL cluster WHLJ161135.9+541635.fits
anim1='/net/nederrijn/data1/osinga/cluster_pointsource_sub/WHLJ161135.9+541635/WHLJ161135.9+541635_maskROBUST-0.5_TEMP-000'
anim2='-model.fits'
outim1='/net/reusel/data1/osinga/phd/year1/deepfields/power_spectrum_halos/OUT/WHLJ161135.9+541635-000'
outim2='-model.fits'

; ADD PS FLUCT TO MODEL 0 to 6
for modnum=0, 5 do begin
  anim=anim1+STRTRIM(modnum,1)+anim2
  outim=outim1+STRTRIM(modnum,1)+outim2
  print, 'Doing ',outim
  fake_halo_ps, analytic=anim,miro='/net/reusel/data1/osinga/phd/year1/deepfields/power_spectrum_halos/OUT/maps_128_5_63.fits',pixel=1.5,scale=4.854,re=080.,xc=1194,yc=844,outfile=outim
end


end

pro WHL2off

; FOR THE SECOND WHL cluster WHLJ161135.9+541635.fits OFF CLUSTER
anim1='/net/nederrijn/data1/osinga/cluster_pointsource_sub/WHLJ161135.9+541635/WHLJ161135.9+541635_maskROBUST-0.5_TEMP-000'
anim2='-model.fits'
outim1='/net/reusel/data1/osinga/phd/year1/deepfields/power_spectrum_halos/OUT/WHLJ161135.9+541635-000'
outim2='-model.fits'

; ADD PS FLUCT TO MODEL 0 to 6
for modnum=0, 5 do begin
  anim=anim1+STRTRIM(modnum,1)+anim2
  outim=outim1+STRTRIM(modnum,1)+outim2
  print, 'Doing ',outim
  fake_halo_ps, analytic=anim,miro='/net/reusel/data1/osinga/phd/year1/deepfields/power_spectrum_halos/OUT/maps_128_5_63.fits',pixel=1.5,scale=4.854,re=080.,xc=1049,yc=813,outfile=outim
end


end



pro WHL3

; FOR THE LAST WHL cluster WHLJ161420.1+544254.fits
anim1='/net/nederrijn/data1/osinga/cluster_pointsource_sub/WHLJ161420.1+544254/WHLJ161420.1+544254_TEMP-000'
anim2='-model.fits'
outim1='/net/reusel/data1/osinga/phd/year1/deepfields/power_spectrum_halos/OUT/WHLJ161420.1+544254-000'
outim2='-model.fits'

; ADD PS FLUCT TO MODEL 0 to 6
for modnum=0, 5 do begin
  anim=anim1+STRTRIM(modnum,1)+anim2
  outim=outim1+STRTRIM(modnum,1)+outim2
  print, 'Doing ',outim
  fake_halo_ps, analytic=anim,miro='/net/reusel/data1/osinga/phd/year1/deepfields/power_spectrum_halos/OUT/maps_128_5_63.fits',pixel=1.5,scale=4.726,re=065.,xc=405,yc=1000,outfile=outim
end


end


pro nightinject1

; FOR THE INJECTION GAME. r_e=65kpc
anim1='/net/achterrijn/data2/osinga/inject_halo_recovery/MCXCJ1036.1+5713_1night_maskROBUST-0.5-000'
anim2='-model.fits'
outim1='/net/reusel/data1/osinga/phd/year1/deepfields/power_spectrum_halos/OUT/MCXCJ1036.1+5713_1night-000'
outim2='-model.fits'

; ADD PS FLUCT TO MODEL 0 to 6
for modnum=0, 5 do begin
  anim=anim1+STRTRIM(modnum,1)+anim2
  outim=outim1+STRTRIM(modnum,1)+outim2
  print, 'Doing ',outim
  fake_halo_ps, analytic=anim,miro='/net/reusel/data1/osinga/phd/year1/deepfields/power_spectrum_halos/OUT/maps_128_5_63.fits',pixel=1.5,scale=4.454,re=065.,xc=1194,yc=844,outfile=outim
end


end


pro nightinject2

; FOR THE INJECTION GAME. r_e=65kpc
anim1='/net/achterrijn/data2/osinga/inject_halo_recovery/MCXCJ1036.1+5713_1night_maskROBUST-0.5-000'
anim2='-model.fits'
outim1='/net/reusel/data1/osinga/phd/year1/deepfields/power_spectrum_halos/OUT/MCXCJ1036.1+5713_1night-000'
outim2='-model.fits'

; ADD PS FLUCT TO MODEL 0 to 6
for modnum=0, 5 do begin
  anim=anim1+STRTRIM(modnum,1)+anim2
  outim=outim1+STRTRIM(modnum,1)+outim2
  print, 'Doing ',outim
  fake_halo_ps, analytic=anim,miro='/net/reusel/data1/osinga/phd/year1/deepfields/power_spectrum_halos/OUT/maps_128_5_63.fits',pixel=1.5,scale=4.454,re=065.,xc=1194,yc=844,outfile=outim
end


end

pro nightinject3

; FOR THE INJECTION GAME. r_e=65kpc
anim1='/net/achterrijn/data2/osinga/inject_halo_recovery/MCXCJ1036.1+5713_1night_maskROBUST-0.5-000'
anim2='-model.fits'
outim1='/net/reusel/data1/osinga/phd/year1/deepfields/power_spectrum_halos/OUT/MCXCJ1036.1+5713_1night-000'
outim2='-model.fits'

; ADD PS FLUCT TO MODEL 0 to 6
for modnum=0, 5 do begin
  anim=anim1+STRTRIM(modnum,1)+anim2
  outim=outim1+STRTRIM(modnum,1)+outim2
  print, 'Doing ',outim
  fake_halo_ps, analytic=anim,miro='/net/reusel/data1/osinga/phd/year1/deepfields/power_spectrum_halos/OUT/maps_128_5_63.fits',pixel=1.5,scale=4.454,re=065.,xc=1194,yc=844,outfile=outim
end


end


pro nightinject4

; FOR THE INJECTION GAME. r_e=80kpc
anim1='/net/achterrijn/data2/osinga/inject_halo_recovery/MCXCJ1036.1+5713_1night_maskROBUST-0.5-000'
anim2='-model.fits'
outim1='/net/reusel/data1/osinga/phd/year1/deepfields/power_spectrum_halos/OUT/MCXCJ1036.1+5713_1night-000'
outim2='-model.fits'

; ADD PS FLUCT TO MODEL 0 to 6
for modnum=0, 5 do begin
  anim=anim1+STRTRIM(modnum,1)+anim2
  outim=outim1+STRTRIM(modnum,1)+outim2
  print, 'Doing ',outim
  fake_halo_ps, analytic=anim,miro='/net/reusel/data1/osinga/phd/year1/deepfields/power_spectrum_halos/OUT/maps_128_5_63.fits',pixel=1.5,scale=4.454,re=080.,xc=1194,yc=844,outfile=outim
end


end

pro nightinject5

; FOR THE INJECTION GAME. r_e=100kpc
anim1='/net/achterrijn/data2/osinga/inject_halo_recovery/MCXCJ1036.1+5713_1night_maskROBUST-0.5-000'
anim2='-model.fits'
outim1='/net/reusel/data1/osinga/phd/year1/deepfields/power_spectrum_halos/OUT/MCXCJ1036.1+5713_1night-000'
outim2='-model.fits'

; ADD PS FLUCT TO MODEL 0 to 6
for modnum=0, 5 do begin
  anim=anim1+STRTRIM(modnum,1)+anim2
  outim=outim1+STRTRIM(modnum,1)+outim2
  print, 'Doing ',outim
  fake_halo_ps, analytic=anim,miro='/net/reusel/data1/osinga/phd/year1/deepfields/power_spectrum_halos/OUT/maps_128_5_63.fits',pixel=1.5,scale=4.454,re=100.,xc=1194,yc=844,outfile=outim
end


end

pro nightinject6
; FOR THE INJECTION GAME. r_e=120kpc
anim1='/net/achterrijn/data2/osinga/inject_halo_recovery/MCXCJ1036.1+5713_1night_maskROBUST-0.5-000'
anim2='-model.fits'
outim1='/net/reusel/data1/osinga/phd/year1/deepfields/power_spectrum_halos/OUT/MCXCJ1036.1+5713_1night-000'
outim2='-model.fits'

; ADD PS FLUCT TO MODEL 0 to 6
for modnum=0, 5 do begin
  anim=anim1+STRTRIM(modnum,1)+anim2
  outim=outim1+STRTRIM(modnum,1)+outim2
  print, 'Doing ',outim
  fake_halo_ps, analytic=anim,miro='/net/reusel/data1/osinga/phd/year1/deepfields/power_spectrum_halos/OUT/maps_128_5_63.fits',pixel=1.5,scale=4.454,re=120.,xc=1194,yc=844,outfile=outim
end


end


pro nightinject7
; FOR THE INJECTION GAME. r_e=65kpc
anim1='/net/achterrijn/data2/osinga/inject_halo_recovery/MCXCJ1036.1+5713_1night_maskROBUST-0.5-000'
anim2='-model.fits'
outim1='/net/reusel/data1/osinga/phd/year1/deepfields/power_spectrum_halos/OUT/MCXCJ1036.1+5713_1night-000'
outim2='-model.fits'

; ADD PS FLUCT TO MODEL 0 to 6
for modnum=0, 5 do begin
  anim=anim1+STRTRIM(modnum,1)+anim2
  outim=outim1+STRTRIM(modnum,1)+outim2
  print, 'Doing ',outim
  fake_halo_ps, analytic=anim,miro='/net/reusel/data1/osinga/phd/year1/deepfields/power_spectrum_halos/OUT/maps_128_5_63.fits',pixel=1.5,scale=4.454,re=065.,xc=1194,yc=844,outfile=outim
end


end


pro haloinject0
; anim1 is the model image that was overwritten with the analytical halo
anim1='/net/achterrijn/data2/osinga/forLUCA/PSZ2G059.18+32.91_image_9-000'
anim2='-model.fits'
outim1='/net/reusel/data1/osinga/phd/year1/deepfields/power_spectrum_halos/OUT/forLUCA-000'
outim2='-model.fits'

; scale is calculated with cosmo.kpc_proper_per_arcmin(redshift).value/60 (i.e., proper kpc per arcsec)
; ADD PS FLUCT TO MODEL 0 to 6
for modnum=0, 5 do begin
  anim=anim1+STRTRIM(modnum,1)+anim2
  outim=outim1+STRTRIM(modnum,1)+outim2
  print, 'Doing ',outim
  fake_halo_ps, analytic=anim,miro='/net/reusel/data1/osinga/phd/year1/deepfields/power_spectrum_halos/OUT/maps_256_5_127.fits',pixel=1.5,scale=5.231,re=189.,xc=558,yc=1324,outfile=outim
end


end



pro PSZRX
; FOR THE UPPER LIM PSZRX
anim1='/net/lofar7/data2/osinga/cluster_pointsource_sub/PSZRXG084.01+46.28/PSZRXG084.01+46.28_TEMP-000'
anim2='-model.fits'
outim1='/net/reusel/data1/osinga/phd/year1/deepfields/power_spectrum_halos/OUT/PSZRXG084.01+46.28-000'
outim2='-model.fits'

; ADD PS FLUCT TO MODEL 0 to 6
for modnum=0, 5 do begin
  anim=anim1+STRTRIM(modnum,1)+anim2
  outim=outim1+STRTRIM(modnum,1)+outim2
  print, 'Doing ',outim
  fake_halo_ps, analytic=anim,miro='/net/reusel/data1/osinga/phd/year1/deepfields/power_spectrum_halos/OUT/maps_128_5_63.fits',pixel=1.5,scale=1.292,re=65.,xc=470,yc=730,outfile=outim
end


end