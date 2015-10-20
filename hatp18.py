import numpy as np
import pdb, sys, os
import reduction
import fitsio
import matplotlib.pyplot as plt

# Number of FITS extensions for the images in this dataste:
N_EXTS = 1

ddir = '/media/hdisk1/data1/gtc/hatp18/GTC2018-09ESO/OB0051'
adir = '/home/tevans/analysis/gtc/hatp18'


def make_master_cals():

    # Create the list of biases:
    n_first = 332777
    n_last = 332877
    framen = range( n_first, n_last+1 )
    nframes = len( framen )
    ddir_bias = os.path.join( ddir, 'bias' )
    bias_frames = []
    for i in range( nframes ):
        bias_frames += [ '0000{0:.0f}-20130520-OSIRIS-OsirisBias.fits'.format( framen[i] ) ]
    bias_list_filename = 'bias.lst'
    bias_list_filepath = os.path.join( adir, bias_list_filename )
    np.savetxt( bias_list_filepath, np.array( bias_frames, dtype=str ), fmt='%s' )
    # Generate the master bias:
    ddir_bias = os.path.join( ddir, 'bias' )
    reduction.median_combine_frames( bias_list_filepath, ddir_bias, 'mbias.fits', frame_type='bias', n_exts=N_EXTS )
    
    # Create the list of flats:
    n_first = 333651
    n_last = 333751
    framen = range( n_first, n_last+1 )
    nframes = len( framen )
    ddir_flat = os.path.join( ddir, 'flat' )
    flat_frames = []
    for i in range( nframes ):
        flat_frames += [ '0000{0:.0f}-20130520-OSIRIS-OsirisDomeFlat.fits'.format( framen[i] ) ]
    flat_list_filename = 'flat.lst'
    flat_list_filepath = os.path.join( adir, flat_list_filename )
    np.savetxt( flat_list_filepath, np.array( flat_frames, dtype=str ), fmt='%s' )
    # Generate the master flat:
    ddir_flat = os.path.join( ddir, 'flat' )
    reduction.median_combine_frames( flat_list_filepath, ddir_flat, 'mflat.fits', frame_type='flat', n_exts=N_EXTS )

    # Create the list of arcs:
    n_first = 333646
    n_last = 333647
    framen = range( n_first, n_last+1 )
    nframes = len( framen )
    ddir_arc = os.path.join( ddir, 'arc' )
    arc_frames = []
    for i in range( nframes ):
        arc_frames += [ '0000{0:.0f}-20130520-OSIRIS-OsirisCalibrationLamp.fits'.format( framen[i] ) ]
    arc_list_filename = 'arc.lst'
    arc_list_filepath = os.path.join( adir, arc_list_filename )
    np.savetxt( arc_list_filepath, np.array( arc_frames, dtype=str ), fmt='%s' )
    # Generate the master arc:
    ddir_arc = os.path.join( ddir, 'arc' )
    reduction.median_combine_frames( arc_list_filepath, ddir_arc, 'marc.fits', frame_type='arc', n_exts=N_EXTS )

    return None


def calibrate_raw_science():

    # Create the list of raw science frames:
    n_first = 332896
    n_last = 333614
    framen = range( n_first, n_last+1 )
    nframes = len( framen )
    ddir_science = os.path.join( ddir, 'object' )
    raw_science_frames = []
    for i in range( nframes ):
        raw_science_frames += [ '0000{0:.0f}-20130520-OSIRIS-OsirisLongSlitSpectroscopy.fits'.format( framen[i] ) ]
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
    adir = '/home/tevans/analysis/gtc/hatp18'
    science_images_full_list_filename = 'cal_science_noflatfield.lst'
    badpix_maps_full_list_filename = 'badpix_full.lst'
    science_images_list_filename = 'science_images.lst' # should apply to all stars...
    badpix_maps_list_filename = 'badpix_maps.lst' # should apply to all stars...
    science_traces_list_filename = [ [ 'science_traces_hatp18.lst' ], [ 'science_traces_reference.lst' ] ]
    science_spectra_list_filename = [ [ 'science_spectra_hatp18.lst' ], [ 'science_spectra_reference.lst' ] ]
    ddir_science = os.path.join( ddir, 'object' )
    ddir_arc = os.path.join( ddir, 'arc' )
    n_exts = 2
    star_names_chip1 = [ 'hatp18' ]
    star_names_chip2 = [ 'reference' ]
    star_names = [ star_names_chip1, star_names_chip2 ]
    disp_axis = 0
    crossdisp_axis = 1
    crossdisp_bounds_chip1 = [ [ 300, 900 ] ]
    disp_bounds_chip1 = [ [ 600, 2050 ] ]
    crossdisp_bounds_chip2 = [ [ 200, 400 ] ]
    disp_bounds_chip2 = [ [ 600, 2050 ] ]
    crossdisp_bounds = [ crossdisp_bounds_chip1, crossdisp_bounds_chip2 ]
    disp_bounds = [ disp_bounds_chip1, disp_bounds_chip2 ]
    stellar = reduction.prep_stellar_obj( adir=adir, \
                                          science_images_full_list_filename=science_images_full_list_filename, \
                                          badpix_maps_full_list_filename=badpix_maps_full_list_filename, \
                                          science_images_list_filename=science_images_list_filename, \
                                          badpix_maps_list_filename=badpix_maps_list_filename, \
                                          science_traces_list_filename=science_traces_list_filename, \
                                          science_spectra_list_filename=science_spectra_list_filename, \
                                          ddir_science=ddir_science, ddir_arc=ddir_arc, n_exts=n_exts, \
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

def extract_spectra( spectral_ap_radius=15, sky_inner_radius=25, sky_band_width=5 ):
    stellar = prep_stellar_obj()
    #stellar.science_traces_lists = [ [ 'science_traces_hatp18.lst' ], [ 'science_traces_reference.lst' ] ]
    #stellar.science_images_lists = [ [ 'science_images_hatp18.lst' ], [ 'science_images_reference.lst' ] ]
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
    hatp18_spectra_list_filepath = os.path.join( stellar.adir, stellar.science_spectra_list[0][0] )
    hatp18_spectra_list = np.loadtxt( hatp18_spectra_list_filepath, dtype=str )
    if len( hatp18_spectra_list )!=nimages:
        pdb.set_trace()
    ref_spectra_list_filepath = os.path.join( stellar.adir, stellar.science_spectra_list[1][0] )
    ref_spectra_list = np.loadtxt( ref_spectra_list_filepath, dtype=str )    
    if len( ref_spectra_list )!=nimages:
        pdb.set_trace()
    mjds = []
    hatp18_spectra = []
    hatp18_disp_pixs = []
    hatp18_skyppix = []
    hatp18_nappixs = []
    hatp18_fwhm = []
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
        hatp18_spectrum_filepath = os.path.join( stellar.adir, hatp18_spectra_list[i] )
        hatp18_spectrum = fitsio.FITS( hatp18_spectrum_filepath )
        hatp18_spectra += [ hatp18_spectrum[1].read_column( 'apflux' ) ]
        hatp18_disp_pixs += [ hatp18_spectrum[1].read_column( 'disp_pixs' ) ]
        hatp18_skyppix += [ hatp18_spectrum[1].read_column( 'skyppix' ) ]
        hatp18_nappixs += [ hatp18_spectrum[1].read_column( 'nappixs' ) ]
        hatp18_fwhm += [ hatp18_spectrum[1].read_header()['FWHM'] ]
        hatp18_spectrum.close()
        ref_spectrum_filepath = os.path.join( stellar.adir, ref_spectra_list[i] )
        ref_spectrum = fitsio.FITS( ref_spectrum_filepath )
        ref_spectra += [ ref_spectrum[1].read_column( 'apflux' ) ]
        ref_disp_pixs += [ ref_spectrum[1].read_column( 'disp_pixs' ) ]
        ref_skyppix += [ ref_spectrum[1].read_column( 'skyppix' ) ]
        ref_nappixs += [ ref_spectrum[1].read_column( 'nappixs' ) ]
        ref_fwhm += [ ref_spectrum[1].read_header()['FWHM'] ]
        ref_spectrum.close()
    mjds = np.array( mjds )
    jds = mjds + 2400000.5
    tmins = ( jds-jds.min() )*24.*60.
    hatp18_spectra = np.row_stack( hatp18_spectra )
    hatp18_disp_pixs = np.row_stack( hatp18_disp_pixs )
    hatp18_skyppix = np.row_stack( hatp18_skyppix )
    hatp18_nappixs = np.row_stack( hatp18_nappixs )
    hatp18_fwhm = np.array( hatp18_fwhm )
    
    ref_spectra = np.row_stack( ref_spectra )
    ref_disp_pixs = np.row_stack( ref_disp_pixs )
    ref_skyppix = np.row_stack( ref_skyppix )
    ref_nappixs = np.row_stack( ref_nappixs )
    ref_fwhm = np.array( ref_fwhm )
    
    y1 = np.sum( hatp18_spectra, axis=1 )
    y2 = np.sum( ref_spectra, axis=1 )
    plt.figure()
    plt.subplot( 311 )
    ax1 = plt.gca()
    plt.plot( y1, '.r' )
    plt.plot( y2, '.b' )
    plt.subplot( 312, sharex=ax1 )
    plt.plot( y1/y2, '.k' )
    plt.subplot( 313, sharex=ax1 )
    plt.plot( hatp18_fwhm, '-r' )
    plt.plot( ref_fwhm, '-b' )
    plt.figure()
    plt.plot( tmins, y1/y2, '.k' )
    plt.title('HAT-P-18')
    plt.ylabel('Relative Flux')
    plt.xlabel('Time (minutes)')
    pdb.set_trace()
    return None
    
    
    
