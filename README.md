# MockHaloInjection
Everything related to injection of Mock halos into radio data.


Credits for the IDL scripts (.pro) go to Annalisa Bonafede 
Credits for the imaging script go to Reinout van Weeren



---------------------------------------------------------------------
## Instructions


1. Make a .yaml file to define the parameter values (see e.g., forLUCA_0.yaml)
2. Change directory structures defined in halo_injection_main.py and halo_predict_image.py to your matching singularity image etc. 
3. Edit the file make_halo.pro to add a hardcoded function for the fluctuations (see e.g., haloinject0 in make_halo.pro)
4. Edit run_imaging.py to suit your imaging needs. By default it does normal imaging and compact source subtraction.
5. Run python halo_injection_main.py 
