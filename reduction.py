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



def median_combine_frames( list_filepath, ddir_output, median_filename, frame_type=None, n_exts=2 ):
    """
    Median combines frames into a master frame.
    """

    frames = np.loadtxt( list_filepath, dtype=str )
    nframes = len( frames )
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
        print 'Different exposure times used:'
        print np.unique( exptimes )
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


def calibrate_raw_science( raw_ddir='', cal_ddir='', mbias_filepath='', mflat_filepath='', raw_science_list_filepath='', n_exts=2 ):
    """
    For each raw science image, this routine trims the array edges and
    performs dark subtraction, before saving the processed image to file.
    Flatfielding is optional. The mflat_filename argument controls whether
    or not flatfielding is done, and if so, which file contains the master
    flat. Only the windows defined by the cross-dispersion and dispersion
    bounds will be flatfielded (after removal of the lamp continuum emission).
    """

    if ( mflat_filepath=='' )+( mflat_filepath==None ):
        cal_science_list_filename='science_noflatfield.lst'
        mflat_filepath = None
    else:
        cal_science_list_filename='science_flatfielded.lst'
        
    mbias_hdu = fitsio.FITS( mbias_filepath, 'r' )
    mbias1 = mbias_hdu[1].read()
    mbias2 = mbias_hdu[2].read()
    mbias_arr = [ mbias1, mbias2 ]
    ## READ OUT BOTH EXTENSIONS!!

    if mflat_filepath!=None:
        # Read in the master flatfield:
        mflat_path = os.path.join( night_dir, mflat_filename )
        mflat_hdu = fitsio.FITS( mflat_path, 'r' )
        mflat1_raw = mflat_hdu[1].read()
        mflat2_raw = mflat_hdu[2].read()
        #mflat_corr = remove_lamp_spectrum( mflat_raw, night )
        #mflat_arr = [ mflat1_corr, mflat2_corr ]
        pdb.set_trace() # todo = remove_lamp_spectrum()
    else:
        mflat_arr = None
    raw_science_filenames = np.loadtxt( raw_science_list_filepath, dtype='str' )
    nframes = len( raw_science_filenames )
    list_dir = os.path.dirname( raw_science_list_filepath )
    cal_science_list_path = os.path.join( list_dir, cal_science_list_filename )
    cal_science_list_ofile = open( cal_science_list_path, 'w' )
    for i in range( nframes ):

        # Read the raw image:
        ifilename = os.path.join( raw_ddir, raw_science_filenames[i] )
        raw_hdu = fitsio.FITS( ifilename, 'r' )

        # Prepare the HDU for writing:
        calext = 'c{0}'.format( raw_science_filenames[i] )
        if mflat_filepath==None:
            calext = calext.replace( '.', '_noflatfielded.' )
        else:
            calext = calext.replace( '.', '_flatfielded.' )
        ofilepath = os.path.join( cal_ddir, calext )
        if os.path.isfile( ofilepath ):
            os.remove( ofilepath )
        cal_hdu = fitsio.FITS( ofilepath, 'rw' )
        # Copy the header for ext 0:
        cal_hdu.write( None, header=raw_hdu[0].read_header() )
        for j in range( n_exts ):

            raw_arr_j = raw_hdu[j+1].read()

            # Trim edges: # DONT DO THIS FOR GTC?
            #raw_arr = raw_arr[ VTRIM[0]:-VTRIM[1], HTRIM[0]:-HTRIM[1] ]

            # Bias subtraction:
            cal_arr_j = raw_arr_j - mbias_arr[j]

            # Flatfield:
            if mflat_arr!=None:
                cal_arr_j = cal_arr_j/mflat_arr[j]

            # Write calibrated data for current extension:
            cal_hdu.write( cal_arr_j, header=raw_hdu[j].read_header() )
            #cal_hdu[0].write_key( 'XOVERSC', 0, comment='Overscan has been trimmed' )
            #cal_hdu[0].write_key( 'YOVERSC', 0, comment='Overscan has been trimmed' )

            #if i>20:
            #    plt.figure()
            #    plt.imshow( cal_arr_j, interpolation='nearest', aspect='auto' )
            #    pdb.set_trace()
        raw_hdu.close()
        cal_hdu.close()
        print 'Image {0} of {1}. Saved {2}'.format( i+1, nframes, ofilepath )
        cal_science_list_ofile.write( '{0}\n'.format( calext ) )
    cal_science_list_ofile.close()

    mbias_hdu.close()
    if mflat_filepath!=None:
        mflat_hdu.close()
    cal_hdu.close()

    return None

def generate_spectra( nstars=2, crossdisp_bounds=[], disp_bounds=[] ):
    
    stellar = prep_stellar_obj()

    # Bounding boxes for the stars in pixel coordinates:
    stellar.nstars = nstars
    stellar.crossdisp_bounds = crossdisp_bounds
    stellar.disp_bounds = disp_bounds

    # Generate the master science subarrays: (???)
    stellar.identify_badpixels()

    # Spectral trace fitting:
    
    

def prep_stellar_obj():
    """
    THIS SHOULD BE GENERALISED TO APPLY TO ARBITRARY GTC DATASETS
    """

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
    stellar.n_exts = 2

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
