import sys, os
import numpy as np
import matplotlib.pyplot as plt

# import tqdm

import subprocess
import pyrap.tables as pt

from astropy.io import fits
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

import pidly

sys.path.append('/software/rhel7/lib64/python2.7/site-packages/')
import yaml

"""
Check how well we can recover injected halos where we know their parameters.
Meant to be run on achterrijn in /data2/osinga/inject_halo_recovery/

"""

def power(redshift, total_flux, alpha=-1.2):
    """
    Calculates the power of a source at the reference frequency given the 
    redshift, integrated flux, and assumed average spectral index alpha
    common is -1.2 or -1.3 for halos

    Uses alpha because k correction has to be made for sources at high redshift
    """

    # luminosity distance
    Dl = cosmo.luminosity_distance(redshift)
    Dl = Dl.to(u.kpc)

    # see (https://arxiv.org/pdf/0802.2770.pdf) 9: Appendix (Page 53)
    L = total_flux*4*np.pi*Dl**2 * (1+redshift)**(-1.0 * alpha - 1)
    L = L.to(u.Watt/u.Hz)
    # print (L)

    return L

def flux_given_power(redshift, power, alpha=-1.3):
    """
    Above function rewritten to solve for S_\nu
    Give power with astropy units included (e.g., u.Watt/u.Hz)
    """
    # luminosity distance
    Dl = cosmo.luminosity_distance(redshift)
    Dl = Dl.to(u.kpc)

    Sv = power/(4*np.pi*Dl**2*(1+redshift)**(-1.0*alpha-1))
    Sv = Sv.to(u.mJy)

    return Sv

def predict_expected_flux(mass,redshift, alpha=-1.3):
    """
    Given a Mass and spectral index, we can compute the expected power at 
    1.4 GHz using the correlation from Cassano+13 and then 
    compute the power at 144 MHz and then compute the flux at 144 MHz.
    """

    # BCES Bisector RH only parameters from Cassano+13
    a = 0.125 
    b = 3.77
    linrel = lambda x: 10**(b*np.log10(x/(10**14.9))+a)*10**24.5 # output y in lin space as function of x in linspace

    power14 = linrel(mass)*(u.Watt/u.Hz) # power at 1.4 GHz in W/Hz
    # Power at 144 MHz
    power144 = power14 * (0.144/1.4)**(alpha) 
    # Total flux density at 144 MHz
    flux144 = flux_given_power(redshift, power144, alpha=-1.3)

    return flux144, power14

def power_size_relation(r_e):
    """ 
    Power size relation from Murgia+2009 (https://arxiv.org/abs/0901.1943)

    Given r_e in kpc, return expected power at 1.4GHz
    """

    return 10**(23.52 + 3*np.log10(r_e/100))

def size_power_relation(p14):
    """ 
    Power size relation from Murgia+2009 (https://arxiv.org/abs/0901.1943)

    Given power at 1.4GHz in W/Hz, return expected size r_e in kpc
    """
    return 100*(10**((np.log10(p14)-23.52)/3))


def circle_model(I0, x0, y0, r_e, shape):
    """
    Adapted from Jort Boxelaar. Simple circular exponential profile

    x0, y0, r_e all in pixel coordinates
    """

    x = np.arange(0, shape[1], step=1, dtype='float')
    y = np.arange(0, shape[0], step=1, dtype='float')
    x_pix, y_pix = np.meshgrid(x,y)

    r   = np.sqrt((x_pix-x0)**2+(y_pix-y0)**2)
    Ir  =  I0 * np.exp(-(r/r_e))

    return Ir

def read_fits_freq(fitsfile):
    """Return the central frequency of the image from the header"""
    with fits.open(fitsfile) as hdul:
        head = hdul[0].header
        centralfreq = head['CRVAL3']
        # double check that CUNIT3 is indeed freq
        if not (head['CUNIT3'].strip(' ') == 'Hz' or head['CUNIT3'].strip(' ') == 'HZ'): 
            raise ValueError("Could not check the header of image %s for the freq info"%fitsfile)
    return centralfreq

def model_spectral_index(data, originalmodel, channelsout=6, alpha=-1.5):
    """
    Scale the circle model with an assumed spectral index.
    Reads the frequency information from the modelimage headers.


    INPUT
    data            -- 2D numpy array -- model data at central (MFS) frequency
    originalmodel   -- string -- prefix name to model iamges
    channelsout     -- int -- number of chanout used originally.
    apha            -- float -- assumed spectral index

    RETURNS
    newmodel -- 3D numpy array, shape (channelsout,data.shape[0],data.shape[1])
                model data for every channel image. 
    """

    # First open the MFS model to get the central frequency
    omodel = originalmodel+'-MFS-model.fits'
    centralfreq = read_fits_freq(omodel)

    newmodel = np.zeros((channelsout,data.shape[0],data.shape[1]))

    for i in range(0,channelsout):
        omodel = originalmodel+'-{0:04}-model.fits'.format(i)
        chanfreq = read_fits_freq(omodel) # freq of this channel
        # Model image of this channel correctly scaled with spectral index
        newmodel[i] = data*( (chanfreq/centralfreq)**alpha)
    return newmodel

def make_MFS_image(data, model, channelsout=6):
    """
    Once the channel model images have been made, average them to 
    create the MFS image. (Just simple averaging)
    """
    omodel = model+'-MFS-model.fits'
    with fits.open(omodel, mode='update') as hdul:
        averaged = np.zeros(hdul[0].data.shape)
        # Open all the channel images. Average them
        for i in range(0,channelsout):
            omodel = model+'-{0:04}-model.fits'.format(i)
            with fits.open(omodel) as hdu2:
                averaged += hdu2[0].data
        averaged /= channelsout
        # Write the averaged image to the MFS image
        hdul[0].data = averaged
        hdul.flush()

def kpc_to_pixcoord(r_e, scale, redshift):
    """
    Given r_e in kpc, return r_e in pix coords. 

    r_e      -- in kpc
    scale    -- conversion 1 pixel = ? arcsec
    redshift -- 

    """
    kpcamin = cosmo.kpc_proper_per_arcmin(redshift)
    r_e_in_asec = (r_e*(u.kpc)/kpcamin).to(u.arcsec)
    r_e_pix = r_e_in_asec.value/scale

    return r_e_pix

def pixcoord_to_kpc(r_e_pix, scale, redshift):
    """
    Given r_e in pixcoords, return r_e in kpc

    r_e      -- in kpc
    scale    -- conversion 1 pixel = ? arcsec
    redshift -- 

    """
    kpcamin = cosmo.kpc_proper_per_arcmin(redshift)
    r_e_in_asec = r_e_pix*scale*u.arcsec
    r_e = (kpcamin*r_e_in_asec).to(u.kpc)

    return r_e_pix

def create_model_image(data, originalmodel,channelsout=6):
    """
    Create a number of model images equal to 'channelsout'. Opens the original
    model images and replaces their data by 'data'. 

    """

    backupmodel = originalmodel +'_backup'


    # Original model will end with '-{0000...0005}-model.fits'
    for i in range(0,channelsout):
        bmodel = backupmodel+'-{0:04}-model.fits'.format(i)
        omodel = originalmodel+'-{0:04}-model.fits'.format(i)

        if not os.path.isfile(bmodel):
            print ("Updating original model image, writing backup model %s"%bmodel)
            os.system("cp ./%s ./%s"%(omodel,bmodel))

        else:
            print ("Updating original model image, backup model already exists.")

        with fits.open(omodel, mode='update') as hdul:
            hdul[0].data = data[i]
            hdul.flush() # write changes

def triple_plot(modeldata, model):
    fig, axes = plt.subplots(1,3)
    ax = axes[0]
    im = ax.imshow(modeldata,origin='lower')
    cbar = fig.colorbar(im)
    ax.set_title("Model data")
    ax = axes[1]
    im = ax.imshow(model,origin='lower')
    cbar = fig.colorbar(im)
    ax.set_title("New model")
    ax = axes[2]
    im = ax.imshow(model+modeldata,origin='lower')
    cbar = fig.colorbar(im)
    ax.set_title("Model + new model")

    plt.show()

def unwrap_mslist(mslist):
    """ Turn the list of strings of MSes into one big string"""
    unwrapped = mslist[0]
    for ms in mslist[1:]:
        unwrapped += ' '+ ms
    return unwrapped

def wsclean_predict(modelimage, mslist, channelsout=6):
    """
    image      -- image to predict from
    msfiles    -- ms files to predict into
    chanelsout -- number of channels out

    """
    
    sexec = 'singularity exec -B /tmp,/dev/shm,/disks/paradata,/data1,/net/lofar1,/net/rijn,/net/nederrijn/,/net/bovenrijn,/net/voorrijn/,/net/para10,/net/lofar2,/net/lofar3,/net/lofar4,/net/lofar5,/net/lofar6,/net/lofar7,/net/reusel/data1/osinga/ /net/lofar1/data1/sweijen/software/LOFAR/singularity/lofarddf.simg'
    sexec = '' # make sure python is run inside the singularity image

    print ("Predicting model image into MODEL_DATA...")

    for ms in mslist:
        wsclean = 'wsclean -channels-out %i'%(channelsout)
        wsclean += ' '+ '-predict -name %s'%(modelimage) # the image to predict from
        # wsclean += ' '+ '-weight briggs -0.5'
        wsclean += ' '+ '%s'%(ms) # the MS to predict into
        print (wsclean)
        subprocess.call(sexec+" "+wsclean, shell=True)

def add_predict_to_data(mslist, outname_end):
    """
    Add the predicted MODEL_DATA column to the data column and save the 
    result in a new column called UPPER_LIM

    """
    outcolumn = 'UPPER_LIM'

    # Check if 'outcolumn' exists, if it doesn't: create it.
    for ms in mslist:
        ts  = pt.table(ms, readonly=False)
        colnames = ts.colnames()
        if outcolumn not in colnames:
            print ('Creating column ' + outcolumn + ' in %s..'%ms)
            desc = ts.getcoldesc('DATA')
            desc['name'] = outcolumn
            ts.addcols(desc)
            ts.close() # to write results

        else:
            print (outcolumn + ' already exists in %s. Overwriting..'%ms)
            ts.close()

    # Add CORRECTED_DATA or DATA to MODEL_DATA and save it in 'outcolumn'
    # for ms in tqdm.tqdm(mslist,desc="Doing %s. Adding model data and saving in %s"%(outname_end,outcolumn)):
    print ("Doing %s. Adding model data and saving in %s"%(outname_end,outcolumn))
    for ms in (mslist):
        ts  = pt.table(ms, readonly=False)
        colnames = ts.colnames()
        if 'CORRECTED_DATA' in colnames:
            print ("Opening CORRECTED_DATA")
            data = ts.getcol('CORRECTED_DATA')
        else:
            print ("Opening DATA")
            data = ts.getcol('DATA') 
        print ("OPENING MODEL_DATA")
        model = ts.getcol('MODEL_DATA') 
        print ("PUTTING DATA+MODEL_DATA in %s"%outcolumn)
        ts.putcol(outcolumn,data+model)
        print ("CLOSING TABLE")
        ts.close()

    print ("ADDED SYNTHETIC MODEL TO COLUMN %s"%outcolumn)
    return

def clean(mslist, outname):
    """
    Pretty manually defined WSclean command for now.

    """
    outcolumn = 'UPPER_LIM'

    # Old singularity image
    # sexec = 'singularity exec -B /tmp,/dev/shm,/disks/paradata,/data1,/net/lofar1,/net/rijn,/net/nederrijn/,/net/bovenrijn,/net/voorrijn/,/net/para10,/net/lofar2,/net/lofar3,/net/lofar4,/net/lofar5,/net/lofar6,/net/lofar7,/net/reusel/data1/osinga/ /net/lofar1/data1/sweijen/software/LOFAR/singularity/lofarddf.simg'

    # New achterrijn singularity image.
    # sexec = 'singularity exec -B /data1,/net/lofar1,/net/rijn,/net/nederrijn/,/net/bovenrijn,/net/voorrijn/,/net/para10,/net/lofar2,/net/lofar3,/net/lofar4,/net/lofar5,/net/lofar6,/net/lofar7,/net/reusel/data1/osinga/ /net/lofar1/data1/sweijen/software/LOFAR/singularity/lofar_sksp_ddf.simg'
    sexec = '' # don't use singularity exec, just start the env before script.

    print ("Make sure python or ipython is being ran inside singularity")

    wsclean = 'wsclean -no-update-model-required -minuv-l 80.0 -size 1538 1538'
    wsclean += ' ' + '-reorder -weight briggs -0.5 -weighting-rank-filter 3 -clean-border 1'
    # wsclean += ' ' + '-parallel-reordering 4'
    wsclean += ' ' + '-mgain 0.8 -fit-beam'
    wsclean += ' ' + '-data-column %s'%outcolumn # Mind the UPPER_LIM data col
    wsclean += ' ' + '-join-channels -channels-out 6'
    wsclean += ' ' + '-padding 1.4 -auto-mask 2.5 -auto-threshold 1.0'
    # wsclean += ' ' + '-fits-mask PSZRXG084.01+46.28_image_7-MFS-image.fits.mask.fits'
    wsclean += ' ' + '-fit-spectral-pol 3 -pol i -baseline-averaging 10.6387917669'
    wsclean += ' ' + '-multiscale -multiscale-scales 0,4,8,16,32,64'
    wsclean += ' ' + '-name %s -scale 1.5arcsec'%outname
    wsclean += ' ' + '-niter 20000'
    wsclean += ' ' + unwrap_mslist(mslist)
    
    print (wsclean)
    subprocess.call(sexec+" "+wsclean, shell=True)

def thismslist():
    # Only 1 night of data
    mslist = ['/data2/osinga/inject_halo_recovery/LOCKMAN_MCXCJ1036.1+5713_try2.dysco.sub.shift.avg.weights.ms.archive0.calibrated']
    return mslist

def write_IDL_script(params):
    """
    Write the IDL script from within python.

    Actually only appends the function we need to run to test.pro
    unless it already exists.
    """

    # Name of function in make_halo.pro that is executed.
    # Dots are not allowed and neither are +
    torun = params['name'].replace('.','_')
    torun = torun.replace('+','_')
    torun += '_it%i'%(params['i'])

    with open('./make_halo.pro', 'a+') as f:
        lines = f.readlines()
        # Check if the function already exists
        exists = False
        for line in lines:
            if torun in line:
                exists = True

        if exists:
            print ("WARNING: Function %s already exists in IDL script"%torun)
            print ("!!! Not updating parameters !!!")

        else:
            print ("Writing function %s to IDL script"%torun)

            # template power spectrum fluctuations
            miro = '/net/reusel/data1/osinga/phd/year1/deepfields/power_spectrum_halos/OUT/maps_256_5_127.fits'

            # conversion kpc to arcsec
            scale = cosmo.kpc_proper_per_arcmin(redshift).value/60
            # # Need r_e in kpc for the IDL script
            if params['r_e_pixel']:
                r_e_pix = params['r_e']
                r_e = pixcoord_to_kpc(r_e_pix, params['scale'], params['redshift'])
            else:
                r_e = params['r_e']

            # Use on-cluster coordinates or off-cluster coordinates
            if params['oncluster']:
                outname_end = 'inject_%i'%(params['i'])
                x0, y0 = params['x0_on'], params['y0_on']
            else:
                outname_end = 'inject_%i_off'%(params['i'])
                x0, y0 = params['x0_off'], params['y0_off']


            write  = "\n\npro %s\n"%torun
            write += "; anim1 is the model image that was overwritten with the analytical halo\n"
            write += "; same image as 'modelimage' in the .yaml file\n"
            write += "anim1='%s'\n"%params['modelimage']
            write += "anim2='-model.fits'\n"
            # Final output MODEL image of halo with PS fluctuations
            outimage = params['modelimage'].split('/')
            outimage.insert(-1,'OUT')
            outimage = '/'.join(outimage)
            outimage += '_it%i'%(params['i'])
            write += "outim1='%s'\n"%outimage
            write += "outim2='-model.fits'\n"
            write += "; scale is calculated with cosmo.kpc_proper_per_arcmin(redshift).value/60 (i.e., proper kpc per arcsec)\n"
            write += "; ADD PS FLUCT TO MODEL CHANNEL IMAGES\n"
            write += "for modnum=0, %i do begin\n"%(params['channelsout']+1)
            write += "  anim=anim1+STRTRIM(modnum,1)+anim2\n"
            write += "  outim=outim1+STRTRIM(modnum,1)+outim2\n"
            write += "  print, 'Doing ',outim\n"
            write += "  fake_halo_ps, analytic=anim,miro='%s',pixel=%.1f,scale=%.3f,re=%.1f,xc=%i,yc=%i,outfile=outim\n"%(miro,params['scale'],scale,r_e,x0, y0)
            write += "end\n\n"
            write += "end"

            f.write(write)

            # e.g.,:
            """
            pro haloinject0
            ; anim1 is the model image that was overwritten with the analytical halo
            ; same image as 'modelimage' in the .yaml file
            anim1='/net/bovenrijn/data1/digennaro/HighRedshiftClusters/uGMRT/HaloInjection/PSZ2G160.83+81.66/PSZ2G160.83+81.66_maskROBUST-0.5uvmin80'
            anim2='-model.fits'
            outim1='/net/bovenrijn/data1/digennaro/HighRedshiftClusters/uGMRT/HaloInjection/PSZ2G160.83+81.66/OUT/upperlimit_PSZ2G160.83+81.66_maskROBUST-0.5uvmin80'
            outim2='-model.fits'

            ; scale is calculated with cosmo.kpc_proper_per_arcmin(redshift).value/60 (i.e., proper kpc per arcsec)
            ; ADD PS FLUCT TO MODEL 0 to 6
            for modnum=0, 5 do begin
              anim=anim1+STRTRIM(modnum,1)+anim2
              outim=outim1+STRTRIM(modnum,1)+outim2
              print, 'Doing ',outim
              fake_halo_ps, analytic=anim,miro='/net/bovenrijn/data1/digennaro/HighRedshiftClusters/uGMRT/HaloInjection/PSZ2G160.83+81.66/OUT/maps_256_5_127.fits',pixel=1.,scale=7.865,re=110.,xc=2003,yc=1902,outfile=outim
            end


            end
            """

def run_IDL_script(params):
    """
    Run the IDL script from within python using pIDLy
    """

    # Name of function in make_halo.pro that is executed. Hardcoded function.
    # Need to edit:
    # the pixel scale, the redshift scale, r_e in kpc, xc, yc and the outim1.
    torun = params['name'].replace('.','_')
    torun = torun.replace('+','_')
    torun += '_it%i'%(params['i'])

    print ("Running %s in IDL"%torun)
    # Open idl
    idl = pidly.IDL()
    # compile the script
    idl('.r /net/reusel/data1/osinga/phd/MockHaloInjection/make_halo.pro')
    # run the correct function in the IDL script
    idl(torun)
    # Close idl
    idl.close()

def inject_fake_halo(x0, y0, I0_units, r_e_pixel, r_e, redshift, scale, shape
    ,modelimage, mslist, outname_end,channelsout, alpha, psfluct):
    """
    Creates a theoretical halo to add to the CORRECTED_DATA or DATA column
    A new column is created called 'UPPER_LIM'.

    x0,y0       -- int       -- center of halo in image coordinates
    I0_units    -- float     -- in Jy/arcsec**2; peak surface brightness
    r_e_pixel   -- bool      -- whether r_e is given in pixels or in kpc
    r_e         -- float     -- in kpc or pixels, e-folding radius
    redshift    -- float     -- redshift of cluster
    scale       -- float     -- size of a pixel in arcsec
    shape       -- (int,int) -- shape of the image in pixel coordinates
    modelimage  -- str       -- which model image to use as template model image
    mslist      -- list      -- all measurement sets as a list of strings
    channelsout -- int       -- amount of channels 
    alpha       -- float     -- average spectral index of the halo
    """

    if r_e_pixel:
        r_e_pix = r_e
    else:
        print ("Using manually defined r_e = %.3f kpc"%r_e)
        # # Need r_e in pixel coordinates
        r_e_pix = kpc_to_pixcoord(r_e, scale, redshift)  
    print ("Using r_e = %.3f pixels"%r_e_pix)

    print ("Using manually defined I0=%f muJy arcsec^{-2}"%(I0_units*1e6))
    # arcsec**-2 to pix**-2 is a factor scale**2
    I0 = I0_units*(scale**2) # Translate to Jy/pix**2
    print ("Which is I0=%f muJy/pix"%(I0*1e6))

    # The theoretical halo model    
    modeldata = circle_model(I0, x0, y0, r_e_pix, shape)
    # Scale this model corresponding to a spectral index alpha for all channels
    modeldata = model_spectral_index(modeldata, modelimage, channelsout=channelsout, alpha=alpha)
    # Overwrite original model images. Backups are saved.
    create_model_image(modeldata, modelimage)
    # Overwrite the MFS model image as well
    make_MFS_image(modeldata, modelimage, channelsout=channelsout)
    if psfluct:
        # Run the IDL script to add the powerspectrum fluctuations.
        # Make sure a function is hardcoded into this script
        write_IDL_script(params)
        run_IDL_script(params)

    # Now we need to call another python script since it has to be run inside the
    # singularity image. But cannot run entire script inside sing image because
    # IDL cannot be run inside singularity
    tocall ="%s python /net/reusel/data1/osinga/phd/MockHaloInjection/halo_predict_image.py %s"%(sexec,paramfile)
    print (tocall)
    os.system(tocall)

def initial_clean(mslist, outname):
    """
    Make the initial image with only one night of observations.
    WSclean command with the same settings.

    LOCKMAN_MCXCJ1036.1+5713_try2.dysco.sub.shift.avg.weights.ms.archive0.calibrated

    """
    datacolumn = 'DATA'

    # Old singularity image
    # sexec = 'singularity exec -B /tmp,/dev/shm,/disks/paradata,/data1,/net/lofar1,/net/rijn,/net/nederrijn/,/net/bovenrijn,/net/voorrijn/,/net/para10,/net/lofar2,/net/lofar3,/net/lofar4,/net/lofar5,/net/lofar6,/net/lofar7,/net/reusel/data1/osinga/ /net/lofar1/data1/sweijen/software/LOFAR/singularity/lofarddf.simg'

    # New achterrijn singularity image.
    # sexec = 'singularity exec -B /data1,/net/lofar1,/net/rijn,/net/nederrijn/,/net/bovenrijn,/net/voorrijn/,/net/para10,/net/lofar2,/net/lofar3,/net/lofar4,/net/lofar5,/net/lofar6,/net/lofar7,/net/reusel/data1/osinga/ /net/lofar1/data1/sweijen/software/LOFAR/singularity/lofar_sksp_ddf.simg'
    sexec = '' # don't use singularity exec, just start the env before script.

    print ("Make sure python or ipython is being ran inside singularity")

    wsclean = 'wsclean -no-update-model-required -minuv-l 80.0 -size 1538 1538'
    wsclean += ' ' + '-reorder -weight briggs -0.5 -weighting-rank-filter 3 -clean-border 1'
    # wsclean += ' ' + '-parallel-reordering 4'
    wsclean += ' ' + '-mgain 0.8 -fit-beam'
    wsclean += ' ' + '-data-column %s'%datacolumn 
    wsclean += ' ' + '-join-channels -channels-out 6'
    wsclean += ' ' + '-padding 1.4 -auto-mask 2.5 -auto-threshold 1.0'
    # wsclean += ' ' + '-fits-mask PSZRXG084.01+46.28_image_7-MFS-image.fits.mask.fits'
    wsclean += ' ' + '-fit-spectral-pol 3 -pol i -baseline-averaging 10.0'
    wsclean += ' ' + '-multiscale -multiscale-scales 0,4,8,16,32,64'
    wsclean += ' ' + '-name %s -scale 1.5arcsec'%outname
    wsclean += ' ' + '-niter 20000'
    wsclean += ' ' + unwrap_mslist(mslist)
    
    print (wsclean)
    subprocess.call(sexec+" "+wsclean, shell=True)

def error_fitted_flux(fluxhalo, sigma_fit, sigma_sub, fluxscale=0.0):
    """
    Return the error on the flux (Modified Eq 1 in https://arxiv.org/pdf/1306.4379.pdf)
    Eq 4 in my paper. 
    Given the fitted flux of the halo, the error due to uncertainty in source 
    subtraction and the error in fluxscale.

    Normally, there'd be also an rms*sqrt(Nbeams) term, but this is incorporated
    by the MCMC fitting procedure. 

    Doesnt account for missing flux due to missing short uv spacings.
    """
    sigma_fsquared = (fluxscale*fluxhalo)**2 + sigma_fit**2 + sigma_sub**2
    return np.sqrt(sigma_fsquared)

def read_parameter_file(paramfile):
    """
    Read parameters from a yaml file
    """

    print ("Reading parameters from %s"%paramfile)
    file = open(paramfile)
    params = yaml.load(file)
    file.close()
    return params

if __name__ == '__main__':
    # Singularity command
    sexec = 'singularity exec -B /data1,/net/lofar1,/net/rijn,/net/nederrijn/,/net/bovenrijn,/net/voorrijn/,/net/para10,/net/lofar2,/net/lofar3,/net/lofar4,/net/lofar5,/net/lofar6,/net/lofar7,/net/reusel/data1/osinga/,/net/voorrijn/,/software,/data2/osinga/,/net/achterrijn/ /net/lofar1/data1/sweijen/software/LOFAR/singularity/lofar_sksp_fedora27_ddf.sif'

    # on /net/achterrijn/. Just a check so we don't mess up unwanted directories
    assert os.getcwd() == '/data2/osinga/forLUCA'
    
    # First argument this script is called with. YAML file.
    paramfile = str(sys.argv[1]) 
    params = read_parameter_file(paramfile)
    print ("Doing inject i=%i"%(params['i']))

    # The Measurement Sets where the halo will be predicted into.
    mslist = params['mslist']

    # outname = 'MCXCJ1036.1+5713_1night_maskROBUST-0.5'
    # # Make the initial image (only have to do once)
    # initial_clean(mslist, outname)
    # sys.exit("Done with initial clean")

    scale = params['scale'] # One pixel is this many arcsec
    shape = params['shape'] # Dimensions of the image

    # Use on-cluster coordinates or off-cluster coordinates
    if params['oncluster']:
        outname_end = 'inject_%i'%(params['i'])
        x0, y0 = params['x0_on'], params['y0_on']
    else:
        outname_end = 'inject_%i_off'%(params['i'])
        x0, y0 = params['x0_off'], params['y0_off']

    redshift = params['redshift'] # for conversion r_e (kpc) to r_e (pixels)
    modelimage = params['modelimage'] # which model image to overwrite with halo

    I0_units = params['I0_units'] # in Jy/arcsec2
    r_e_pixel = params['r_e_pixel']
    r_e = params['r_e']      
    channelsout = params['channelsout']
    alpha = params['alpha'] 
    psfluct = params['psfluct'] 

    print ("Doing cluster %s, i=%i"%(params['name'],params['i']))
    if psfluct:
        print ("Adding power spectrum fluctuations")
    else:
        print ("Not adding power spectrum fluctuations")

    inject_fake_halo(x0, y0, I0_units, r_e_pixel, r_e, redshift, scale, shape, modelimage, mslist
        , outname_end, channelsout, alpha, psfluct)







