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

sys.path.append('/software/rhel7/lib64/python2.7/site-packages/')
import yaml

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
    
    # sexec = 'singularity exec -B /tmp,/dev/shm,/disks/paradata,/data1,/net/lofar1,/net/rijn,/net/nederrijn/,/net/bovenrijn,/net/voorrijn/,/net/para10,/net/lofar2,/net/lofar3,/net/lofar4,/net/lofar5,/net/lofar6,/net/lofar7,/net/reusel/data1/osinga/ /net/lofar1/data1/sweijen/software/LOFAR/singularity/lofarddf.simg'
    sexec = '' # make sure python is run inside the singularity image

    print ("Predicting model image into MODEL_DATA...")

    for ms in mslist:
        wsclean = 'wsclean -channels-out %i'%(channelsout)
        wsclean += ' '+ '-predict -name %s'%(modelimage) # the image to predict from
        # wsclean += ' '+ '-weight briggs -0.5'
        wsclean += ' '+ '%s'%(ms) # the MS to predict into
        print (wsclean)
        subprocess.call(sexec+" "+wsclean, shell=True)

def add_predict_to_data(mslist):
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

    # Add DIFFUSE_SUB to MODEL_DATA and save it in 'outcolumn'
    for ms in mslist:#,desc="Adding model data and saving in %s"%outcolumn):
        ts  = pt.table(ms, readonly=False)
        colnames = ts.colnames()
        # """
        if 'CORRECTED_DATA' in colnames:
            print ("Opening CORRECTED_DATA")
            data = ts.getcol('CORRECTED_DATA')
        else:
            print ("Opening DATA")
            data = ts.getcol('DATA') 
        # """
        print ("OPENING MODEL_DATA")
        model = ts.getcol('MODEL_DATA') 
        print ("PUTTING DIFFUSE_SUB+MODEL_DATA in %s"%outcolumn)
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

    # New singularity image.
    # sexec = 'singularity exec -B /data1,/net/lofar1,/net/rijn,/net/nederrijn/,/net/bovenrijn,/net/voorrijn/,/net/para10,/net/lofar2,/net/lofar3,/net/lofar4,/net/lofar5,/net/lofar6,/net/lofar7,/net/reusel/data1/osinga/,/net/voorrijn/ /net/lofar1/data1/sweijen/software/LOFAR/singularity/lofar_sksp_fedora27_ddf.sif'
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
    wsclean += ' ' + '-fit-spectral-pol 3 -pol i -baseline-averaging 10.0'
    wsclean += ' ' + '-multiscale -multiscale-scales 0,4,8,16,32,64'
    wsclean += ' ' + '-name %s -scale 1.5arcsec'%outname
    wsclean += ' ' + '-niter 20000'
    wsclean += ' ' + unwrap_mslist(mslist)
    
    print (wsclean)
    subprocess.call(sexec+" "+wsclean, shell=True)

def save_params_header(outname):
    """
    After imaging, save the params used to the header of the MFS-image.
    """
    imname = outname + '-MFS-image.fits'
    with fits.open(imname, mode='update') as hdul:
        header = hdul[0].header

        for key in params:
            comment = '%s: %s'%(key,str(params[key]))
            header.add_comment(comment) 

        hdul.flush()

def image_and_pssub(params):
    imsize = params['shape'][0]
    niter = 25000 # default
    channelsout = params['channelsout']
    minuv = 80 # default
    pixelscale = params['scale']
    redshift = params['redshift']
    sourceLLS = 0.5 # Mpc
    outname = params['name']
    mses = ' '.join(params['mslist'])


    cmd = 'python'
    cmd += ' ./run_imaging.py'
    cmd += ' --imsize %i -n %i --channelsout %i --minuv %.5f'%(imsize,niter,channelsout,minuv)
    cmd += ' --pixelscale %.1f --sourceLLS %.2f --z %.3f'%(pixelscale,sourceLLS,redshift)
    cmd += ' --weighting ROBUST-0.5'
    cmd += ' --array uGMRT'
    cmd += ' -i %s_%i'%(outname,params['i'])
    cmd += ' --maskthreshold 3 %s'%mses

    print (cmd)
    os.system(cmd)

def copy_injected_mses(params, outdir):
    """
    After we're done, save the ms with the injected halo column
    to directory 'outdir'
    
    Make sure outdir is a function of iteration number, otherwise they will be
    overwitten again anyways
    """

    if (not os.path.exists(outdir)):
        print ("Creating directories where MSes are saved")
        print (outdir)
        os.system('mkdir -p %s'%outdir)

    print ("Copying injected MSes to output directory %s"%outdir)

    for ms in params['mslist']:
        cmd = 'cp -r'
        cmd += ' %s'%ms
        cmd += ' %s'%(outdir)
        print(cmd)
        os.system(cmd)

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

    
    paramfile = str(sys.argv[1])
    params = read_parameter_file(paramfile)

    psfluct = params['psfluct'] 

    if psfluct:
        PFmodel = '/net/reusel/data1/osinga/phd/MockHaloInjection/OUT/'
        # PFmodel = '/net/reusel/data1/osinga/phd/year1/deepfields/power_spectrum_halos/OUT/'
        PFmodel += params['name'] 

    else:
        print ("Assuming no power spectrum fluctuations")
        modelimage = params['modelimage'] # which model image to overwrite with halo
        # No power spectrum fluctuations, just use regular model.
        PFmodel = modelimage

    i = params['i']

    if params['oncluster']:
        outname_end = 'inject_%i'%(i)
        x0, y0 = params['x0_on'], params['y0_on']
    else:
        outname_end = 'inject_%i_off'%(i)
        x0, y0 = params['x0_off'], params['y0_off']

    print ("Doing cluster %s, i=%i"%(params['name'],params['i']))
    print ("Doing wsclean predict, add and imaging step. Should be called from halo_injection_main.py")
    
    # Predict the new model image WITH or WITHOUT Power spectrum Fluctuations 
    # to the MODEL_DATA column of the mslist mses
    wsclean_predict(PFmodel, params['mslist'])
    # Add the model data column to the CORRECTED_DATA or DATA column 
    add_predict_to_data(params['mslist'])
    print ("Done adding predict to CORRECTED DATA. DOING IMAGING AND PS SUB NOW.")
    image_and_pssub(params)
    outdir = './%s/'%params['name'] +  'inject_%i/'%params['i']
    copy_injected_mses(params, outdir)    
    print ("Done with i=%i. Can run Halo-FDCA now."%params['i'])