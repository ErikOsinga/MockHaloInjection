pro sim_miro


;Comment for Erik: change input and output folder
  fold='./'  ;input folder
  foldo=fold+'OUT/'  ;output folder, it must already exist

  ;...constants
  mp=1.67e-24
  msol=1.98e33
  kb=1.38e-16 ;..erg/K
  mu=0.59 ;::.
  cmtoMpc=3.08e24


;Comment for Erik: I am 99% confident there's no need to change
;these, only used for gas density profile. 
;...properties of galaxy cluster

  tkT=8;...keV - only used for gas density profile
  tkT*=1.16e7  
  t500=tkT  ;..temperature T500 in kelvin only used for gas density profile
  cc_ncc="ncc"           ;...this cluster is a Non Cool-Core -> different fitting formula from Ghirardini+
   

;Comment for Erik: I think you just have to touch these 5 parameter below 
  z=0.383  ;cluster redshift
  n=256 ; grid size in pixels
  res=10  ;..cell resolution of output model (in kpc)
  kin=5. ;minimum scale for the power-law power-spectrum in Fourier space (so maximum physical scale would be n*res/kmin [kpc])
  kmax=(n*0.5)-1.   ;...maximum k  (n0*0.5-1 is the Nyquist frequency)
  rc=290.  ;cluster core radius in kpc

  print,"Starting Simulation"

;...main function. Check createB_3d.pro for explanations on all possible keywords and default values

    miro,n0=n,res=res,cc_ncc=cc_ncc,t500=t500,temp=temp,kin=kin,kmax=kmax,rc=rc,dir_out=foldo
  
    print,"Simulation Done." 
  
  
  end
