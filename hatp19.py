import numpy as np
import pdb, sys, os
import reduction
import fitsio
import matplotlib.pyplot as plt


# Number of FITS extensions for the images in this dataste:
N_EXTS = 1

ddir = '/media/hdisk1/data1/gtc/hatp19/GTC2018-09ESO/OB0058'
adir = '/home/tevans/analysis/gtc/hatp19'


def make_master_cals():

    # Create the list of biases:
    n_first = 407479
    n_last = 407528
    framen = range( n_first, n_last+1 )
    nframes = len( framen )
    ddir_bias = os.path.join( ddir, 'bias' )
    bias_frames = []
    for i in range( nframes ):
        bias_frames += [ '0000{0:.0f}-20130728-OSIRIS-OsirisBias.fits'.format( framen[i] ) ]
    bias_list_filename = 'bias.lst'
    bias_list_filepath = os.path.join( adir, bias_list_filename )
    np.savetxt( bias_list_filepath, np.array( bias_frames, dtype=str ), fmt='%s' )
    # Generate the master bias:
    ddir_bias = os.path.join( ddir, 'bias' )
    reduction.median_combine_frames( bias_list_filepath, ddir_bias, 'mbias.fits', frame_type='bias', n_exts=N_EXTS )
    
    # Create the list of flats:
    n_first = 407531
    n_last = 407630
    framen = range( n_first, n_last+1 )
    nframes = len( framen )
    ddir_flat = os.path.join( ddir, 'flat' )
    flat_frames = []
    for i in range( nframes ):
        flat_frames += [ '0000{0:.0f}-20130728-OSIRIS-OsirisSpectralFlat.fits'.format( framen[i] ) ]
    flat_list_filename = 'flat.lst'
    flat_list_filepath = os.path.join( adir, flat_list_filename )
    np.savetxt( flat_list_filepath, np.array( flat_frames, dtype=str ), fmt='%s' )
    # Generate the master flat:
    ddir_flat = os.path.join( ddir, 'flat' )
    reduction.median_combine_frames( flat_list_filepath, ddir_flat, 'mflat.fits', frame_type='flat', n_exts=N_EXTS )

    # Create the list of arcs:
    framens = [ 407634, 407636 ]
    ddir_arc = os.path.join( ddir, 'arc' )
    arc_frames = []
    for i in framens:
        arc_frames += [ '0000{0:.0f}-20130728-OSIRIS-OsirisCalibrationLamp.fits'.format( i ) ]
    arc_list_filename = 'arc.lst'
    arc_list_filepath = os.path.join( adir, arc_list_filename )
    np.savetxt( arc_list_filepath, np.array( arc_frames, dtype=str ), fmt='%s' )
    # Generate the master arc:
    ddir_arc = os.path.join( ddir, 'arc' )
    reduction.median_combine_frames( arc_list_filepath, ddir_arc, 'marc.fits', frame_type='arc', n_exts=N_EXTS )

    return None


def calibrate_raw_science():

    # Create the list of raw science frames:
    n_first = 407201
    n_last = 407475
    framen = range( n_first, n_last+1 )
    nframes = len( framen )
    ddir_science = os.path.join( ddir, 'object' )
    raw_science_frames = []
    for i in range( nframes ):
        raw_science_frames += [ '0000{0:.0f}-20130728-OSIRIS-OsirisLongSlitSpectroscopy.fits'.format( framen[i] ) ]
    raw_science_list_filename = 'raw_science.lst'
    raw_science_list_filepath = os.path.join( adir, raw_science_list_filename )
    np.savetxt( raw_science_list_filepath, np.array( raw_science_frames, dtype=str ), fmt='%s' )

    # Call the calibration routine:
    bias_dir = os.path.join( ddir, 'bias' )
    mbias_filepath = os.path.join( bias_dir, 'mbias.fits' )
    mflat_filepath = None
    reduction.calibrate_raw_science( n_exts=N_EXTS, raw_ddir=ddir_science, cal_ddir=ddir_science, \
                                     mbias_filepath=mbias_filepath, mflat_filepath=None, \
                                     raw_science_list_filepath=raw_science_list_filepath )

    return None


def prep_stellar_obj():
    #todo + include call to combine_spectra() in reduction.py
    adir = '/home/tevans/analysis/gtc/hatp19'
    science_images_full_list_filename = 'cal_science_noflatfield.lst'
    badpix_maps_full_list_filename = 'badpix_full.lst'
    science_images_list_filename = 'science_images.lst' # should apply to all stars...
    badpix_maps_list_filename = 'badpix_maps.lst' # should apply to all stars...
    science_traces_list_filename = [ 'science_traces_reference.lst', 'science_traces_hatp19.lst' ]
    science_spectra_list_filename = [ 'science_spectra_reference.lst', 'science_spectra_hatp19.lst' ]
    ddir_science = os.path.join( ddir, 'object' )
    ddir_arc = os.path.join( ddir, 'arc' )
    star_names = [ 'reference',  'hatp19' ]
    disp_axis = 0
    crossdisp_axis = 0
    crossdisp_bounds = [ [ 150, 350 ], [ 510, 710 ] ]
    disp_bounds = [ [ 400, 2050 ], [ 400, 2050 ] ]
    stellar = reduction.prep_stellar_obj( adir=adir, \
                                          science_images_full_list_filename=science_images_full_list_filename, \
                                          badpix_maps_full_list_filename=badpix_maps_full_list_filename, \
                                          science_images_list_filename=science_images_list_filename, \
                                          badpix_maps_list_filename=badpix_maps_list_filename, \
                                          science_traces_list_filename=science_traces_list_filename, \
                                          science_spectra_list_filename=science_spectra_list_filename, \
                                          ddir_science=ddir_science, ddir_arc=ddir_arc, n_exts=N_EXTS, \
                                          star_names=star_names, \
                                          disp_axis=disp_axis, crossdisp_axis=crossdisp_axis, \
                                          crossdisp_bounds=crossdisp_bounds, disp_bounds=disp_bounds )
    return stellar


def run_all():
    #identify_bad_pixels()
    fit_traces()
    extract_spectra( spectral_ap_radius=15, sky_inner_radius=25, sky_band_width=5 )
    return None

def identify_bad_pixels():
    stellar = prep_stellar_obj()
    stellar.identify_bad_pixels()
    return None

def fit_traces():
    stellar = prep_stellar_obj()
    stellar.fit_traces( make_plots=True )
    return None

def extract_spectra( spectral_ap_radius=30, sky_inner_radius=45, sky_band_width=5 ):
    stellar = prep_stellar_obj()
    stellar.spectral_ap_radius = spectral_ap_radius
    stellar.sky_inner_radius = sky_inner_radius
    stellar.sky_band_width = sky_band_width
    stellar.extract_spectra()
    return None

def combine_spectra():
    # todo = add something to this routines that extracts
    # auxiliary variables, including opening each image and
    # reading stuff out of the headers like time stamp.
    stellar = prep_stellar_obj()
    image_list_filepath = os.path.join( stellar.adir, stellar.science_images_list )
    image_list = np.loadtxt( image_list_filepath, dtype=str )
    nimages = len( image_list )
    hatp19_spectra_list_filepath = os.path.join( stellar.adir, stellar.science_spectra_list[0] )
    hatp19_spectra_list = np.loadtxt( hatp19_spectra_list_filepath, dtype=str )
    if len( hatp19_spectra_list )!=nimages:
        pdb.set_trace()
    ref_spectra_list_filepath = os.path.join( stellar.adir, stellar.science_spectra_list[1] )
    ref_spectra_list = np.loadtxt( ref_spectra_list_filepath, dtype=str )    
    if len( ref_spectra_list )!=nimages:
        pdb.set_trace()
    mjds = []
    hatp19_spectra = []
    hatp19_disp_pixs = []
    hatp19_skyppix = []
    hatp19_nappixs = []
    hatp19_fwhm = []
    ref_spectra = []
    ref_disp_pixs = []
    ref_skyppix = []
    ref_nappixs = []
    ref_fwhm = []    
    for i in range( nimages ):
        print '... reading spectrum {0} of {1}'.format( i+1, nimages )
        image_filepath = os.path.join( stellar.ddir, image_list[i] )
        hdu = fitsio.FITS( image_filepath )
        h0 = hdu[0].read_header()
        mjds += [ h0['MJD-OBS'] ]
        hatp19_spectrum_filepath = os.path.join( stellar.adir, hatp19_spectra_list[i] )
        hatp19_spectrum = fitsio.FITS( hatp19_spectrum_filepath )
        hatp19_spectra += [ hatp19_spectrum[1].read_column( 'apflux' ) ]
        hatp19_disp_pixs += [ hatp19_spectrum[1].read_column( 'disp_pixs' ) ]
        hatp19_skyppix += [ hatp19_spectrum[1].read_column( 'skyppix' ) ]
        hatp19_nappixs += [ hatp19_spectrum[1].read_column( 'nappixs' ) ]
        hatp19_fwhm += [ hatp19_spectrum[1].read_header()['FWHM'] ]
        hatp19_spectrum.close()
        ref_spectrum_filepath = os.path.join( stellar.adir, ref_spectra_list[i] )
        ref_spectrum = fitsio.FITS( ref_spectrum_filepath )
        ref_spectra += [ ref_spectrum[1].read_column( 'apflux' ) ]
        ref_disp_pixs += [ ref_spectrum[1].read_column( 'disp_pixs' ) ]
        ref_skyppix += [ ref_spectrum[1].read_column( 'skyppix' ) ]
        ref_nappixs += [ ref_spectrum[1].read_column( 'nappixs' ) ]
        ref_fwhm += [ ref_spectrum[1].read_header()['FWHM'] ]
        ref_spectrum.close()
        #print 'aaa'
        #plt.figure()
        #plt.plot(hatp19_spectra[-1],'-b')
        #plt.plot(ref_spectra[-1],'-r')
        #pdb.set_trace()
        #print 'bbb'
    mjds = np.array( mjds )
    jds = mjds + 2400000.5
    tmins = ( jds-jds.min() )*24.*60.
    hatp19_spectra = np.row_stack( hatp19_spectra )
    hatp19_disp_pixs = np.row_stack( hatp19_disp_pixs )
    hatp19_skyppix = np.row_stack( hatp19_skyppix )
    hatp19_nappixs = np.row_stack( hatp19_nappixs )
    hatp19_fwhm = np.array( hatp19_fwhm )
    
    ref_spectra = np.row_stack( ref_spectra )
    ref_disp_pixs = np.row_stack( ref_disp_pixs )
    ref_skyppix = np.row_stack( ref_skyppix )
    ref_nappixs = np.row_stack( ref_nappixs )
    ref_fwhm = np.array( ref_fwhm )
    
    y1 = np.sum( hatp19_spectra, axis=1 )
    y2 = np.sum( ref_spectra, axis=1 )
    plt.figure()
    plt.subplot( 311 )
    ax1 = plt.gca()
    plt.plot( tmins, y1, '.r' )
    plt.plot( tmins, y2, '.b' )
    plt.subplot( 312, sharex=ax1 )
    plt.plot( tmins, y1/y2, '.k' )
    plt.subplot( 313, sharex=ax1 )
    plt.plot( tmins, hatp19_fwhm, '-r' )
    plt.plot( ref_fwhm, '-b' )
    plt.figure()
    plt.plot( tmins, y2/y1, '.k' )
    plt.title('HAT-P-19')
    plt.ylabel('Relative Flux')
    plt.xlabel('Time (minutes)')
    pdb.set_trace()
    return None
    
    
    
