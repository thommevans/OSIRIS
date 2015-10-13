import fitsio
import pdb, sys, os
import numpy as np
import matplotlib.pyplot as plt
from spectroscopy_dev.spectroscopy import stellar as stellar_module

# Attributes:
#    stellar.ddir
#    stellar.science_images_list = []
#    stellar.badpix_map = None

# Stellar attributes that could be lists, with each entry corresponding
# to a different image extension:
#    stellar.star_names = [ [], [], ... ]
#    stellar.disp_axis = [ ]
#    stellar.crossdisp_axis = [ ]



def median_combine_frames( list_filepath, ddir_output, median_filename, frame_type=None ):
    """
    Median combines the bias frames into a master bias.
    """

    frames = np.loadtxt( list_filepath, dtype=str )
    nframes = len( frames )
    nexts = 2
    exptimes = np.zeros( nframes )
    print '\nReading in individual frames...'
    cube1 = []
    cube2 = []
    for i in range( nframes ):
        filepath = os.path.join( ddir_output, frames[i] )
        print '{0} of {1}. {2}'.format( i+1, nframes, filepath )
        hdu = fitsio.FITS( filepath, 'r' )
        h0 = hdu[0].read_header()
        exptimes[i] = h0['EXPTIME']
        cube1 += [ hdu[1].read() ] # 1st chip
        cube2 += [ hdu[2].read() ] # 2nd chip
        hdu.close()
    if len( np.unique( exptimes ) )!=1:
        pdb.set_trace() # all exptimes should be the same
    else:
        exptime = exptimes[0]
    cube1 = np.dstack( cube1 )
    cube2 = np.dstack( cube2 )
    median1 = np.median( cube1, axis=2 )
    median2 = np.median( cube2, axis=2 )
    # TRIM THE EDGES LIKE FOR THE NOT??
    header = { 'OBJECT':'master_{0}'.format( frame_type ), 'EXPTIME':exptime }
    save_filepath = os.path.join( ddir_output, median_filename )
    if os.path.isfile( save_filepath ):
        os.remove( save_filepath )
    hdu = fitsio.FITS( save_filepath, 'rw' )
    hdu.write( None, header=None )
    hdu.write( median1, header=header )
    hdu.write( median2, header=header )
    hdu.close()
    print '\nSaved master {0}:\n{1}'.format( frame_type, save_filepath )
    
    return None



def prep_stellar_obj():

    arc_ext = 'arc'
    bias_ext = 'bias'
    flat_ext = 'flat'
    science_ext = 'object'

    n_first = 332896
    n_last = 333614
    framen = range( n_first, n_last+1 )
    nframes = len( framen )
    ddir_science = os.path.join( ddir, 'object' )
    science_frames = []
    for i in range( nframes ):
        science_frames += [ '0000{0:.0f}-20130520-OSIRIS-OsirisLongSlitSpectroscopy.fits'.format( framen[i] ) ]
    science_frames = np.array( science_frames, dtype=str )


    stellar = stellar_module.stellar()
    stellar.ddir = ddir_science
    stellar.adir = '/home/tevans/analysis/gtc'
    stellar.star_names = [ [ 'hatp18' ], [ 'reference' ] ]
    stellar.science_images_list = os.path.join( stellar.adir, 'temp.lst' )
    np.savetxt( stellar.science_images_list, science_frames, fmt='%s' )
    stellar.badpix_map = None
    stellar.disp_axis = 0
    stellar.crossdisp_axis = 1
    stellar.next = 2

    # SHOULD THESE IMAGES BE CALIBRATED (DARK, FLAT ETC) BEFORE THIS STEP?
    image_filenames = np.loadtxt( stellar.science_images_list, dtype=str )
    nimages = len( image_filenames )
    stellar.jds = np.ones( nimages )
    stellar.exptime_secs = np.ones( nimages )
    for i in range( nimages ):
        print i+1, nimages
        image = os.path.join( stellar.ddir, image_filenames[i] )
        hdu = fitsio.FITS( image )
        h0 = hdu[0].read_header()
        stellar.jds[i] = h0['MJD-OBS'] + 2400000.5
        stellar.exptime_secs[i] = h0['EXPTIME']
    stellar.header_kws['gain'] = 'GAIN'    

    stellar.wcal_kws = {}
    marc_filename = '' # todo = make a master arc
    stellar.wcal_kws['marc_file'] = os.path.join( adir_full, marc_filename )
    stellar.wcal_kws['fiducial_lines'] = [] # todo = fiducial wavelengths in units of nm

    return stellar
