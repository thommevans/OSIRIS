import numpy as np
import pdb, sys, os
import reduction
import fitsio
import matplotlib.pyplot as plt
import cPickle
import utils


# Number of FITS extensions for the images in this dataste:
N_EXTS = 1

ddir = '/media/hdisk1/data1/gtc/hatp19/GTC2018-09ESO/OB0058'
adir = '/home/tevans/analysis/gtc/hatp19'


def make_master_cals( bias=True, flat=True, arc=True ):

    # Create the list of biases:
    if bias==True:
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
    if flat==True:
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

    # Create the list of HgAr arcs:
    if arc==True:
        framens = [ 407634 ]
        ddir_arc = os.path.join( ddir, 'arc' )
        arc_frames = []
        for i in framens:
            arc_frames += [ '0000{0:.0f}-20130728-OSIRIS-OsirisCalibrationLamp.fits'.format( i ) ]
        arc_list_filename = 'arc_HgAr.lst'
        arc_list_filepath = os.path.join( adir, arc_list_filename )
        np.savetxt( arc_list_filepath, np.array( arc_frames, dtype=str ), fmt='%s' )
        # Generate the master HgAr arc:
        ddir_arc = os.path.join( ddir, 'arc' )
        reduction.median_combine_frames( arc_list_filepath, ddir_arc, 'marc_HgAr.fits', frame_type='arc', n_exts=N_EXTS )

    # Create the list of Ne arcs:
    if arc==True:
        framens = [ 407636 ]
        ddir_arc = os.path.join( ddir, 'arc' )
        arc_frames = []
        for i in framens:
            arc_frames += [ '0000{0:.0f}-20130728-OSIRIS-OsirisCalibrationLamp.fits'.format( i ) ]
        arc_list_filename = 'arc_Ne.lst'
        arc_list_filepath = os.path.join( adir, arc_list_filename )
        np.savetxt( arc_list_filepath, np.array( arc_frames, dtype=str ), fmt='%s' )
        # Generate the master Ne arc:
        ddir_arc = os.path.join( ddir, 'arc' )
        reduction.median_combine_frames( arc_list_filepath, ddir_arc, 'marc_Ne.fits', frame_type='arc', n_exts=N_EXTS )

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


def prep_stellar_obj( spectral_ap_radius=None, sky_inner_radius=None, flatfield=False ):
    #todo + include call to combine_spectra() in reduction.py
    adir = '/home/tevans/analysis/gtc/hatp19'
    science_images_full_list_filename = 'cal_science_noflatfield.lst'
    badpix_maps_full_list_filename = 'badpix_full.lst'
    science_images_list_filename = 'science_images.lst' # should apply to all stars...
    badpix_maps_list_filename = 'badpix_maps.lst' # should apply to all stars...
    science_traces_list_filename = [ 'science_traces_reference.lst', 'science_traces_hatp19.lst' ]
    ext = '_ap{0:.2f}bg{1:.2f}'.format( spectral_ap_radius, sky_inner_radius )
    science_spectra_list_filename = [ 'science_spectra_reference{0}.lst'.format( ext ), \
                                      'science_spectra_hatp19{0}.lst'.format( ext ) ]
    ddir_science = os.path.join( ddir, 'object' )
    ddir_arc = os.path.join( ddir, 'arc' )
    star_names = [ 'reference',  'hatp19' ]
    disp_axis = 0
    crossdisp_axis = 1
    crossdisp_bounds = [ [ 150, 350 ], [ 510, 710 ] ]
    disp_bounds = [ [ 600, 2050 ], [ 600, 2050 ] ]
    marc_filename = 'marc_HgAr.fits'
    star0_wavcal_crossdisp_bounds = crossdisp_bounds[0]
    star1_wavcal_crossdisp_bounds = crossdisp_bounds[1]
    wavcal_fiducial_lines = [ 4046.563, 5460.735, 7067.21, 7272.936, 7383.981, 7635.106, \
                              7723.98, 7948.176, 8264.522, 8521.442, 8667.944 ]
    star0_wavcal_disp_bounds = [ [ 785, 797 ], \
                                 [ 1186, 1204 ], \
                                 [ 1608, 1620 ], \
                                 [ 1659, 1675 ], \
                                 [ 1687, 1705 ], \
                                 [ 1752, 1767 ], \
                                 [ 1775, 1795 ], \
                                 [ 1830, 1847 ], \
                                 [ 1910, 1930, ], \
                                 [ 1967, 1979], \
                                 [ 1979, 1990] ]
    star1_wavcal_disp_bounds = star0_wavcal_disp_bounds
    wavcal_crossdisp_bounds = [ star0_wavcal_crossdisp_bounds, star1_wavcal_crossdisp_bounds ]
    wavcal_disp_bounds = [ star0_wavcal_disp_bounds, star1_wavcal_disp_bounds ]

    stellar = reduction.prep_stellar_obj( adir=adir, \
                                          science_images_full_list_filename=science_images_full_list_filename, \
                                          badpix_maps_full_list_filename=badpix_maps_full_list_filename, \
                                          science_images_list_filename=science_images_list_filename, \
                                          badpix_maps_list_filename=badpix_maps_list_filename, \
                                          science_traces_list_filename=science_traces_list_filename, \
                                          science_spectra_list_filename=science_spectra_list_filename, \
                                          ddir_science=ddir_science, ddir_arc=ddir_arc, n_exts=N_EXTS, \
                                          star_names=star_names, marc_filename=marc_filename, \
                                          disp_axis=disp_axis, crossdisp_axis=crossdisp_axis, \
                                          wavcal_fiducial_lines=wavcal_fiducial_lines, \
                                          wavcal_crossdisp_bounds=wavcal_crossdisp_bounds, \
                                          wavcal_disp_bounds=wavcal_disp_bounds, \
                                          crossdisp_bounds=crossdisp_bounds, disp_bounds=disp_bounds  )

    return stellar


def calibrate_wavelength_scale( spectral_ap_radius=None, sky_inner_radius=None, \
                                flatfield=False, make_plots=True, poly_order=3 ):
    stellar = prep_stellar_obj( spectral_ap_radius=spectral_ap_radius, sky_inner_radius=sky_inner_radius, \
                                flatfield=flatfield )
    stellar.calibrate_wavelength_scale( make_plots=make_plots, poly_order=poly_order )

def run_all():
    #identify_bad_pixels()
    fit_traces()
    extract_spectra( spectral_ap_radius=15, sky_inner_radius=25, sky_band_width=5 )
    return None


def run_extraction_grid():
    flatfield = False
    wavcal_poly_order = 3
    #spectral_ap_radii = [ 20, 30, 40, 50 ]
    #sky_inner_radii_ext = [ 0, 15, 30 ]
    spectral_ap_radii = [ 30 ]
    sky_inner_radii_ext = [ 15 ]
    bw = 30
    n = len( spectral_ap_radii )
    m = len( sky_inner_radii_ext )
    for i in range( n ):
        print '\n\n{0}\nRunning apradius {1} of {2}:\n'.format( 50*'#', i+1, n )
        r = spectral_ap_radii[i]
        for j in range( m ):
            s = r + sky_inner_radii_ext[j]
            extract_spectra( spectral_ap_radius=r, sky_inner_radius=s, sky_band_width=bw, flatfield=flatfield )
            calibrate_wavelength_scale( spectral_ap_radius=r, sky_inner_radius=s, sky_band_width=bw, \
                                        make_plots=True, flatfield=flatfield, poly_order=wavcal_poly_order )
            pdb.set_trace()
    return None


def identify_bad_pixels():
    stellar = prep_stellar_obj()
    stellar.identify_bad_pixels()
    return None

def fit_traces():
    stellar = prep_stellar_obj()
    stellar.fit_traces( make_plots=True )
    return None

def extract_spectra( spectral_ap_radius=30, sky_inner_radius=45, sky_band_width=5, flatfield=False ):
    stellar = prep_stellar_obj( spectral_ap_radius=spectral_ap_radius, sky_inner_radius=sky_inner_radius, \
                                flatfield=flatfield )
    image_filenames = np.loadtxt( os.path.join( stellar.adir, stellar.science_images_list ), dtype=str )
    stellar.science_spectra_filenames = make_spectra_filenames( image_filenames, \
                                                                spectral_ap_radius=spectral_ap_radius, \
                                                                sky_inner_radius=sky_inner_radius )
    stellar.spectral_ap_radius = spectral_ap_radius
    stellar.sky_inner_radius = sky_inner_radius
    stellar.sky_band_width = sky_band_width
    stellar.extract_spectra()
    return None

def make_spectra_filenames( image_filenames, spectral_ap_radius=30, sky_inner_radius=45, shifted=False ):
    spectra_filenames = []
    nimages = len( image_filenames )
    for i in range( nimages ):
        f = image_filenames[i].replace( '-OSIRIS-OsirisLongSlitSpectroscopy', '' )
        ext = '_ap{0:.2f}bg{1:.2f}'.format( spectral_ap_radius, sky_inner_radius )
        spectra_filenames += [ f.replace( '.fits', '{0}.fits'.format( ext ) ) ]
        print spectra_filenames[-1]
    if shifted==True:
        spectra_filenames = create_shifted_spectra_filenames( spectra_filenames )
    return spectra_filenames

def make_shifted_spectra_filenames( spectra_filenames ):
    n = len( spectra_filenames )
    new_filenames = []
    for i in range( n ):
        new_filenames += [ spectra_filenames[i].replace( '.fits', '_wavshifted.fits' ) ]
    return new_filenames
    

def crosscor_spectra( spectral_ap_radius=30, sky_inner_radius=45, remove_continuum=False, spectra_smoothing_sig=5, continuum_smoothing_sig=20 ):

    stellar = prep_stellar_obj()
    image_filenames = np.loadtxt( os.path.join( stellar.adir, stellar.science_images_list ), dtype=str )
    stellar.spectra_filenames = make_spectra_filenames( image_filenames, spectral_ap_radius=spectral_ap_radius, \
                                                        sky_inner_radius=sky_inner_radius, shifted=False )
    nframes = len( stellar.spectra_filenames )
    max_shift_pix = 2
    resample_factor = 1000
    
    star_names = [ 'hatp19', 'reference' ]
    disp_bound_ixs = [ 400, 1405 ]
    m = len( star_names )
    dwavs = []
    ref_spectrum = []
    spectra = []
    spectra_shifted = []
    dspec = []
    #plt.figure()
    for k in range( m ):
        print 'Star {0} of {1}'.format( k+1, m )
        spec_dir = os.path.join( stellar.adir, 'spectra/{0}'.format( star_names[k] ) )
        spectra_k = []
        for i in range( nframes ):
            print i+1, nframes
            ipath = os.path.join( spec_dir, stellar.spectra_filenames[i] )
            hdu = fitsio.FITS( ipath )
            spectra_k += [ hdu[1].read_column( 'apflux' ) ]
            hdu.close()
        spectra_k = np.row_stack( spectra_k )
        ref_spectrum_k = np.median( spectra_k[:12,:], axis=0 )
        dwavs_k, spectra_shifted_k, dspec_k = reduction.crosscor( spectra_k, ref_spectrum_k, \
                                                                  max_shift=max_shift_pix, \
                                                                  spectra_smoothing_sig=spectra_smoothing_sig, \
                                                                  continuum_smoothing_sig=continuum_smoothing_sig, \
                                                                  remove_continuum=remove_continuum, \
                                                                  disp_bound_ixs=disp_bound_ixs, \
                                                                  resample_factor=resample_factor )
        dwavs += [ dwavs_k ]
        ref_spectrum +=  [ ref_spectrum_k ]
        spectra += [ spectra_k ]        
        spectra_shifted += [ spectra_shifted_k ]
        dspec += [ dspec_k ]
        
    shifted_spectra_filenames = make_shifted_spectra_filenames( stellar.spectra_filenames )

    plt.figure()
    plt.subplot( 211 )
    plt.plot( spectra[0][50,:], '-k' )
    plt.plot( spectra[0][120,:], '-c' )
    ax = plt.gca()
    plt.subplot( 212, sharex=ax, sharey=ax )
    plt.plot( spectra_shifted[0][50,:], '-g' )
    plt.plot( spectra_shifted[0][120,:], '-r' )

    plt.figure()
    plt.subplot( 211 )
    plt.plot( spectra[1][50,:], '-m' )
    plt.plot( spectra[1][120,:], '-k' )
    ax = plt.gca()
    plt.subplot( 212, sharex=ax, sharey=ax )
    plt.plot( spectra_shifted[1][50,:], '-c' )
    plt.plot( spectra_shifted[1][120,:], '-b' )

    
    ix0 = disp_bound_ixs[0]
    ix1 = disp_bound_ixs[1]
    n1 = 100
    n2 = 250
    for g in range(2):
        r0 = spectra[g][n1,:]-spectra[g][n2,:]
        rms0 = np.sqrt( np.mean( r0[ix0:ix1+1]**2. ) )
        r1 = spectra_shifted[g][n1,:]-spectra_shifted[g][n2,:]
        rms1 = np.sqrt( np.mean( r1[ix0:ix1+1]**2. ) )
        print '\nthis ', rms1
        print 'should be lower than this ', rms0
    plt.figure()
    plt.subplot(211)
    w0 = np.sum( spectra_shifted[0], axis=1 )
    w1 = np.sum( spectra_shifted[1], axis=1 )
    w0 /= w0[0]
    w1 /= w1[0]
    plt.plot( w0/w1, '-k' )
    #plt.plot( w1, '-b' )
    ax = plt.gca()
    plt.subplot(212,sharex=ax)
    plt.plot( dwavs[0] )
    plt.plot( dwavs[1] )
    
    pdb.set_trace()
        
    return None

def combine_spectra( spectral_ap_radius=None, sky_inner_radius=None, flatfield=False ):
    # todo = add something to this routines that extracts
    # auxiliary variables, including opening each image and
    # reading stuff out of the headers like time stamp.
    stellar = prep_stellar_obj( spectral_ap_radius=spectral_ap_radius, sky_inner_radius=sky_inner_radius, \
                                flatfield=flatfield )
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
    mjd = []
    airmass = []
    hatp19_spectra = []
    hatp19_disp_pixs = []
    hatp19_wav = []
    hatp19_skyppix = []
    hatp19_nappixs = []
    hatp19_fwhm = []
    ref_spectra = []
    ref_disp_pixs = []
    ref_wav = []
    ref_skyppix = []
    ref_nappixs = []
    ref_fwhm = []    
    for i in range( nimages ):
        print '... reading spectrum {0} of {1}'.format( i+1, nimages )
        image_filepath = os.path.join( stellar.ddir, image_list[i] )
        hdu = fitsio.FITS( image_filepath )
        h0 = hdu[0].read_header()
        mjd += [ h0['MJD-OBS'] ]
        airmass += [ h0['AIRMASS'] ]
        hatp19_spectrum_filepath = os.path.join( stellar.adir, hatp19_spectra_list[i] )
        hatp19_spectrum = fitsio.FITS( hatp19_spectrum_filepath )
        hatp19_spectra += [ hatp19_spectrum[1].read_column( 'apflux' ) ]
        hatp19_disp_pixs += [ hatp19_spectrum[1].read_column( 'disp_pixs' ) ]
        hatp19_wav += [ hatp19_spectrum[1].read_column( 'wav' ) ]
        hatp19_skyppix += [ hatp19_spectrum[1].read_column( 'skyppix' ) ]
        hatp19_nappixs += [ hatp19_spectrum[1].read_column( 'nappixs' ) ]
        hatp19_fwhm += [ hatp19_spectrum[1].read_header()['FWHM'] ]
        hatp19_spectrum.close()
        ref_spectrum_filepath = os.path.join( stellar.adir, ref_spectra_list[i] )
        ref_spectrum = fitsio.FITS( ref_spectrum_filepath )
        ref_spectra += [ ref_spectrum[1].read_column( 'apflux' ) ]
        ref_disp_pixs += [ ref_spectrum[1].read_column( 'disp_pixs' ) ]
        ref_wav += [ ref_spectrum[1].read_column( 'wav' ) ]
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
    mjd = np.array( mjd )
    airmass = np.array( airmass )
    jd = mjd + 2400000.5
    tmins = ( jd-jd.min() )*24.*60.
    hatp19_spectra = np.row_stack( hatp19_spectra )
    hatp19_disp_pixs = np.row_stack( hatp19_disp_pixs )[0,:]
    hatp19_wav = np.row_stack( hatp19_wav )[0,:]
    hatp19_skyppix = np.row_stack( hatp19_skyppix )
    hatp19_nappixs = np.row_stack( hatp19_nappixs )
    hatp19_fwhm = np.array( hatp19_fwhm )
    
    ref_spectra = np.row_stack( ref_spectra )
    ref_disp_pixs = np.row_stack( ref_disp_pixs )[0,:]
    ref_wav = np.row_stack( ref_wav )[0,:]
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

    hatp19 = { 'spectra':hatp19_spectra, 'disp_pixs':hatp19_disp_pixs, 'wav':hatp19_wav, 'skyppix':hatp19_skyppix, 'nappixs':hatp19_nappixs, 'fwhm':hatp19_fwhm  }
    ref = { 'spectra':ref_spectra, 'disp_pixs':ref_disp_pixs, 'wav':ref_wav, 'skyppix':ref_skyppix, 'nappixs':ref_nappixs, 'fwhm':ref_fwhm  }
    output = { 'jd':jd, 'airmass':airmass, 'hatp19':hatp19, 'ref':ref }

    opath = '/home/tevans/analysis/gtc/hatp19/temp_spectra.pkl'
    ofile = open( opath, 'w' )
    cPickle.dump( output, ofile )
    ofile.close()
    print '\nSaved: {0}'.format( opath )

    return None

def Na_lightcurve():
    ifile = open( '/home/tevans/analysis/gtc/hatp19/temp_spectra.pkl' )
    z = cPickle.load( ifile )
    ifile.close()
    jd = z['jd']
    spectra1 = z['hatp19']['spectra']
    spectra2 = z['ref']['spectra']
    wav1 = z['hatp19']['wav']
    wav2 = z['ref']['wav']
    Na_center = 5892
    Na_bw = 30
    wide_bw = 300
    Na_cuton = Na_center-Na_bw
    Na_cutoff = Na_center+Na_bw
    wide_cuton = Na_center-wide_bw
    wide_cutoff = Na_center+wide_bw

    nav = 30 
    ixsNa1 = (wav1>=Na_cuton)*(wav1<=Na_cutoff)
    ixsNa2 = (wav2>=Na_cuton)*(wav2<=Na_cutoff)
    fNa1 = np.sum( spectra1[:,ixsNa1], axis=1 )
    fNa2 = np.sum( spectra2[:,ixsNa2], axis=1 )
    fNa = fNa2/fNa1
    fNa /= np.mean( fNa[:nav] )
    ixsw1 = (wav1>=wide_cuton)*(wav1<=wide_cutoff)
    ixsw2 = (wav2>=wide_cuton)*(wav2<=wide_cutoff)
    fw1 = np.sum( spectra1[:,ixsw1], axis=1 )
    fw2 = np.sum( spectra2[:,ixsw2], axis=1 )
    fw = fw2/fw1
    fw /= np.mean( fw[:nav] )
    plt.figure()
    plt.plot( jd, fNa, '-b' )
    plt.plot( jd, fw, '-r' )
    plt.figure()
    plt.plot(wav1,spectra1[0,:],'-r')
    plt.plot(wav1,spectra1[0,:],'.r')
    plt.plot(wav1[ixsNa1],spectra1[0,:][ixsNa1],'.k')
    plt.axvline(Na_cuton)
    plt.axvline(Na_cutoff)
    pdb.set_trace()

def quick_test():

    ifile = open( '/home/tevans/analysis/gtc/hatp19/temp_spectra.pkl' )
    z = cPickle.load( ifile )
    ifile.close()
    jd = z['jd']
    airmass = z['airmass']
    x2 = z['hatp19']['disp_pixs']
    fwhm2 = z['hatp19']['fwhm']
    sky2 = np.sum( z['hatp19']['skyppix'], axis=1 )
    x1 = z['ref']['disp_pixs']
    fwhm1 = z['ref']['fwhm']
    sky1 = np.sum( z['ref']['skyppix'], axis=1 )
    s2 = z['hatp19']['spectra']
    s1 = z['ref']['spectra']

    ix = 400

    w1 = np.sum( s1[:,ix:], axis=1 )
    w2 = np.sum( s2[:,ix:], axis=1 )
    wa = w2/w1

    wb = np.sum( s2[:,ix:]/s1[:,ix:], axis=1 )
    
    wa /= wa[0]
    wb /= wb[0]

    plt.figure()
    plt.plot( jd, wa, '.r' )
    plt.plot( jd, wb, '.b' )


    n = len( jd )
    offset = np.ones( n )

    dxs = []

    dpix = 15
    ix = 600
    chB1 = np.sum( s1[:,ix:ix+dpix+1], axis=1 )
    chB2 = np.sum( s2[:,ix:ix+dpix+1], axis=1 )
    dxs += [ [ ix, ix+dpix ] ]
    chB = chB2/chB1
    chB /= chB[0]
    chBdiff = chB/wa
    phi = np.column_stack( [ offset, jd ] )
    c = np.linalg.lstsq( phi, chBdiff )[0]
    ttrend = np.dot( phi, c )
    chBdiff /= ttrend
    med = np.median( chBdiff )
    chBdiff -= med
    chBdiff *= 1e6
    sig = np.std( chBdiff )
    nsig = np.abs( chBdiff/sig )
    ixsB = nsig<4

    ix = 1115
    chV1 = np.sum( s1[:,ix:ix+dpix+1], axis=1 )
    chV2 = np.sum( s2[:,ix:ix+dpix+1], axis=1 )
    dxs += [ [ ix, ix+dpix ] ]
    chV = chV2/chV1
    chV /= chV[0]
    chVdiff = chV/wa
    phi = np.column_stack( [ offset, jd ] )
    c = np.linalg.lstsq( phi, chVdiff )[0]
    ttrend = np.dot( phi, c )
    chVdiff /= ttrend
    med = np.median( chVdiff )
    chVdiff -= med
    chVdiff *= 1e6
    sig = np.std( chVdiff )
    nsig = np.abs( chVdiff/sig )
    ixsV = nsig<4

    ix = 1445
    chR1 = np.sum( s1[:,ix:ix+dpix+1], axis=1 )
    chR2 = np.sum( s2[:,ix:ix+dpix+1], axis=1 )
    dxs += [ [ ix, ix+dpix ] ]
    chR = chR2/chR1
    chR /= chR[0]
    chRdiff = chR/wa
    chRdiff /= chRdiff[0]
    phi = np.column_stack( [ offset, jd ] )
    c = np.linalg.lstsq( phi, chRdiff )[0]
    ttrend = np.dot( phi, c )
    chRdiff /= ttrend
    med = np.median( chRdiff )
    chRdiff -= med
    chRdiff *= 1e6
    sig = np.std( chRdiff )
    nsig = np.abs( chRdiff/sig )
    ixsR = nsig<4
    
    plt.figure()
    plt.plot(s2[0,:],'-b')
    plt.plot(s1[0,:],'-r')
    for i in range( 3 ):
        print dxs[i]
        plt.axvspan( dxs[i][0], dxs[i][1], color=0.6*np.ones(3), alpha=0.5 )
    plt.xlabel( 'Dispersion axis (pixel)' )
    plt.ylabel( 'Flux' )
    plt.title( 'HAT-P-19' )

    plt.figure()
    plt.subplot( 411 )
    ax = plt.gca()
    plt.plot( jd, wa, '-b' )
    plt.plot( jd, wb, '-r' )
    plt.ylabel( 'white' )
    plt.subplot( 412, sharex=ax )
    plt.plot( jd, chB, '-g' )
    plt.ylabel( 'B' )
    plt.subplot( 413, sharex=ax )
    plt.plot( jd, chV, '-g' )
    plt.ylabel( 'V' )
    plt.subplot( 414, sharex=ax )
    plt.plot( jd, chR, '-g' )
    plt.ylabel( 'R' )
    plt.xlabel( 'JD' )
    plt.suptitle( 'HAT-P-19' )

    nbins = int( np.round( n/25. ) )
    t = ( jd-jd.min() )*24*60

    plt.figure()

    ylow = -10000
    yupp = 10000

    plt.subplot( 611 )
    ax = plt.gca()
    stdB_1 = np.std( chBdiff[ixsB] )
    plt.plot( t[ixsB], chBdiff[ixsB], '-g', label='std_1={0:.0f}ppm'.format( stdB_1 ) )
    ax0 = plt.gca()
    tb, yb, stdvs, npb = utils.bin_1d( t[ixsB], chBdiff[ixsB], nbins=nbins )
    stdB_25 = np.std( yb[npb>20] )
    plt.plot( tb[npb>20], yb[npb>20], '-or', label='std_25={0:.0f}ppm'.format( stdB_25 ) )
    plt.legend( numpoints=1 )
    plt.ylabel( 'B resids' )
    plt.ylim( [ ylow, yupp ] )

    plt.subplot( 612, sharex=ax )
    stdV_1 = np.std( chVdiff[ixsV] )
    plt.plot( t[ixsV], chVdiff[ixsV], '-g', label='std_1={0:.0f}ppm'.format( stdV_1 ) )
    ax0 = plt.gca()
    tb, yb, stdvs, npb = utils.bin_1d( t[ixsV], chVdiff[ixsV], nbins=nbins )
    stdV_25 = np.std( yb[npb>20] )
    plt.plot( tb[npb>20], yb[npb>20], '-or', label='std_25={0:.0f}ppm'.format( stdV_25 ) )
    plt.legend( numpoints=1 )
    plt.ylabel( 'V resids' )
    plt.ylim( [ ylow, yupp ] )

    plt.subplot( 613, sharex=ax )
    stdR_1 = np.std( chRdiff[ixsR] )
    plt.plot( t[ixsR], chRdiff[ixsR], '-g', label='std_1={0:.0f}ppm'.format( stdR_1 ) )
    ax0 = plt.gca()
    tb, yb, stdvs, npb = utils.bin_1d( t[ixsR], chRdiff[ixsR], nbins=nbins )
    stdR_25 = np.std( yb[npb>20] )
    plt.plot( tb[npb>20], yb[npb>20], '-or', label='std_25={0:.0f}ppm'.format( stdR_25 ) )
    plt.legend( numpoints=1 )
    plt.ylabel( 'R resids' )
    plt.ylim( [ ylow, yupp ] )

    plt.subplot( 614, sharex=ax )
    plt.plot( t, fwhm1, '-r' )
    plt.plot( t, fwhm2, '-b' )
    plt.ylabel( 'FWHM' )
    plt.subplot( 615, sharex=ax )
    plt.plot( t, sky1/np.median( sky1 ), '-r' )
    plt.plot( t, sky2/np.median( sky1 ), '-b' )
    plt.ylabel( 'Sky' )
    plt.subplot( 616, sharex=ax )
    plt.plot( t, airmass, '-g' )
    plt.ylabel( 'Airm' )
    plt.xlabel( 'time (minutes)' )
    
    plt.suptitle( 'HAT-P-19' )
    ax.set_xlim( [ 0, 370 ] )

    pdb.set_trace()
    return None



    
    
