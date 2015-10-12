import fitsio
import pdb, sys, os
import numpy as np
import matplotlib.pyplot as plt
from spectroscopy_dev.spectroscopy import stellar as stellar_module

# HAT-P-18
ddir = '/media/hdisk1/data1/gtc/hatp18/GTC2018-09ESO/OB0051'

# Attributes:
#    stellar.ddir
#    stellar.science_images_list = []
#    stellar.badpix_map = None

# Stellar attributes that could be lists, with each entry corresponding
# to a different image extension:
#    stellar.star_names = [ [], [], ... ]
#    stellar.disp_axis = [ ]
#    stellar.crossdisp_axis = [ ]


def make_marc():

    n_first = 333646
    n_last = 333648
    framen = range( n_first, n_last+1 )
    nframes = len( framen )
    ddir_bias = os.path.join( ddir, 'arc' )
    science_frames = []
    for i in range( nframes ):
        science_frames += [ '0000{0:.0f}-20130520-OSIRIS-OsirisLongSlitSpectroscopy.fits'.format( framen[i] ) ]
    science_frames = np.array( science_frames, dtype=str )
    
    return None


def make_mbias():

    n_first = 332777
    n_last = 332877
    framen = range( n_first, n_last+1 )
    nframes = len( framen )
    ddir_bias = os.path.join( ddir, 'bias' )
    bias_frames = []
    for i in range( nframes ):
        bias_frames += [ '0000{0:.0f}-20130520-OSIRIS-OsirisBias.fits'.format( framen[i] ) ]
    nframes = len( bias_frames )
    nexts = 2
    mbias_filename = 'mbias.fits'
    mbias_path = os.path.join( ddir_bias, mbias_filename )
    exptimes = np.zeros( nframes )
    print '\nReading in bias arrays...'
    bias_cube1 = []
    bias_cube2 = []
    for i in range( nframes ):
        bias_filepath = os.path.join( ddir_bias, bias_frames[i] )
        print '{0} of {1}. {2}'.format( i+1, nframes, bias_filepath )
        hdu = fitsio.FITS( bias_filepath, 'r' )
        h0 = hdu[0].read_header()
        exptimes[i] = h0['EXPTIME']
        bias_cube1 += [ hdu[1].read() ] # 1st chip
        bias_cube2 += [ hdu[2].read() ] # 2nd chip
        hdu.close()
    if len( np.unique( exptimes ) )!=1:
        pdb.set_trace() # all exptimes should be the same
    else:
        exptime = exptimes[0]
    bias_cube1 = np.dstack( bias_cube1 )
    bias_cube2 = np.dstack( bias_cube2 )
    mbias1 = np.median( bias_cube1, axis=2 )
    mbias2 = np.median( bias_cube2, axis=2 )
    # TRIM THE EDGES LIKE FOR THE NOT??
    header = { 'OBJECT':'master_bias', 'EXPTIME':exptime }
    if os.path.isfile( mbias_path ):
        os.remove( mbias_path )
    hdu = fitsio.FITS( mbias_path, 'rw' )
    hdu.write( None, header=None )
    hdu.write( mbias1, header=header )
    hdu.write( mbias2, header=header )
    hdu.close()
    print '\nSaved master bias:\n{0}'.format( mbias_path )
    pdb.set_trace()
    
    return None


def make_mflat():

    n_first = 333651
    n_last = 333751
    framen = range( n_first, n_last+1 )
    nframes = len( framen )
    ddir_bias = os.path.join( ddir, 'bias' )
    science_frames = []
    for i in range( nframes ):
        science_frames += [ '0000{0:.0f}-20130520-OSIRIS-OsirisLongSlitSpectroscopy.fits'.format( framen[i] ) ]
    science_frames = np.array( science_frames, dtype=str )
    
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
