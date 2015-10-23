import fitsio
import pdb, sys, os, copy
import numpy as np
import matplotlib.pyplot as plt
from spectroscopy_dev.spectroscopy import stellar as stellar_module
import scipy.interpolate
import scipy.ndimage

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
    if np.ndim( frames )==0:
        nframes = 1
        frames = [ str( frames ) ]
    else:
        nframes = len( frames )
    exptimes = np.zeros( nframes )
    print '\nReading in individual frames...'
    cube1 = []
    if n_exts==2:
        cube2 = []
    for i in range( nframes ):
        filepath = os.path.join( ddir_output, frames[i] )
        print '{0} of {1}. {2}'.format( i+1, nframes, filepath )
        hdu = fitsio.FITS( filepath, 'r' )
        h0 = hdu[0].read_header()
        exptimes[i] = h0['EXPTIME']
        cube1 += [ hdu[1].read() ] # 1st chip
        if n_exts==2:
            cube2 += [ hdu[2].read() ] # 2nd chip
        hdu.close()
    if len( np.unique( exptimes ) )!=1:
        print 'Different exposure times used:'
        print np.unique( exptimes )
        pdb.set_trace() # all exptimes should be the same
    else:
        exptime = exptimes[0]
    cube1 = np.dstack( cube1 )
    median1 = np.median( cube1, axis=2 )    
    if n_exts==2:
        cube2 = np.dstack( cube2 )
        median2 = np.median( cube2, axis=2 )
    # TRIM THE EDGES LIKE FOR THE NOT??
    header = { 'OBJECT':'master_{0}'.format( frame_type ), 'EXPTIME':exptime }
    save_filepath = os.path.join( ddir_output, median_filename )
    if os.path.isfile( save_filepath ):
        os.remove( save_filepath )
    hdu = fitsio.FITS( save_filepath, 'rw' )
    hdu.write( None, header=None )
    hdu.write( median1, header=header )
    if n_exts==2:
        hdu.write( median2, header=header )
    hdu.close()
    print '\nSaved master {0}:\n{1}'.format( frame_type, save_filepath )
    
    return None


def calibrate_raw_science( n_exts=1, raw_ddir='', cal_ddir='', mbias_filepath='', mflat_filepath='', raw_science_list_filepath='' ):
    """
    For each raw science image, this routine trims the array edges and
    performs dark subtraction, before saving the processed image to file.
    Flatfielding is optional. The mflat_filename argument controls whether
    or not flatfielding is done, and if so, which file contains the master
    flat. Only the windows defined by the cross-dispersion and dispersion
    bounds will be flatfielded (after removal of the lamp continuum emission).
    """

    if ( mflat_filepath=='' )+( mflat_filepath==None ):
        cal_science_list_filename = 'cal_science_noflatfield.lst'
        mflat_filepath = None
    else:
        cal_science_list_filename = 'cal_science_flatfielded.lst'
        
    mbias_hdu = fitsio.FITS( mbias_filepath, 'r' )
    mbias1 = mbias_hdu[1].read()
    if n_exts==1:
        mbias_arr = [ mbias1 ]
    elif n_exts==2:
        mbias2 = mbias_hdu[2].read()
        mbias_arr = [ mbias1, mbias2 ]
    else:
        pdb.set_trace()

    if mflat_filepath!=None:
        # Read in the master flatfield:
        mflat_path = os.path.join( night_dir, mflat_filename )
        mflat_hdu = fitsio.FITS( mflat_path, 'r' )
        mflat1_raw = mflat_hdu[1].read()
        if n_exts==1:
            mflat_arr = [ mflat1_raw ]
        elif n_exts==2:
            mflat2_raw = mflat_hdu[2].read()
            mflat_arr = [ mflat1_raw, mflat2_raw ]
        #mflat_corr = remove_lamp_spectrum( mflat_raw, night )
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
            calext = calext.replace( '.', '_noflatfield.' )
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

#def generate_spectra( adir='', science_images_list_filename='', badpix_maps_list_filename='', ddir_science='', ddir_arc='', n_exts=2, star_names=[], disp_axis=0, crossdisp_axis=1, crossdisp_bounds=[], disp_bounds=[] ):
def identify_bad_pixels( stellar ): # THIS IS A REDUNDANT ROUTINE PROBABLY

    #stellar = prep_stellar_obj( adir=adir, science_images_list_filename=science_images_list_filename, badpix_maps_list_filename=badpix_maps_list_filename, ddir_science=ddir_science, n_exts=n_exts, star_names=star_names, disp_axis=disp_axis, crossdisp_axis=crossdisp_axis )

    #stellar.wcal_kws = {}
    #marc_filename = '' # todo = make a master arc
    #stellar.wcal_kws['marc_file'] = os.path.join( ddir_arc, marc_filename )
    #stellar.wcal_kws['fiducial_lines'] = [] # todo = fiducial wavelengths in units of nm

    # Bounding boxes for the stars in pixel coordinates:
    #if n_exts>1:
    #    stellar.nstars = []
    #    for i in range( n_exts ):
    #        stellar.nstars += [ len( star_names[i] ) ]
    #else:
    #    stellar.nstars = len( star_names )
                                
    #stellar.crossdisp_bounds = crossdisp_bounds
    #stellar.disp_bounds = disp_bounds

    # Generate the list of filenames that the bad pixel maps
    # will be saved to:
    #science_image_filenames = np.loadtxt( stellar.science_images_list, dtype=str )
    #nimages = len( science_image_filenames )
    #badpix_map_filenames = []
    #for i in range( nimages ):
    #    badpix_map_filenames += [ 'bp{0}'.format( science_image_filenames[i][1:] ) ]
    #stellar.badpix_maps_list = os.path.join( os.path.dirname( stellar.science_images_list ), 'badpix_maps.lst' )
    #np.savetxt( stellar.badpix_maps_list, badpix_map_filenames, fmt='%s' )
    stellar.identify_badpixels()

    return stellar

def fit_traces( stellar ): # THIS IS A REDUNDANT ROUTINE PROBABLY
    # Spectral trace fitting:
    #stellar.fit_traces( make_plots=True )
    stellar.fit_traces( 3333 )
    pdb.set_trace()
    
    return None
    

def prep_stellar_obj( adir='', ddir_science='', ddir_arc='', n_exts=2, star_names=[], science_images_full_list_filename='', badpix_maps_full_list_filename='', science_images_list_filename='', badpix_maps_list_filename='', science_traces_list_filename='', marc_filename='', science_spectra_list_filename='', disp_axis=0, crossdisp_axis=1, crossdisp_bounds=[], disp_bounds=[], wavcal_fiducial_lines=[], wavcal_crossdisp_bounds=[], wavcal_disp_bounds=[] ):
    """
    THIS SHOULD BE GENERALISED TO APPLY TO ARBITRARY GTC DATASETS
    """

    #arc_ext = 'arc'
    #bias_ext = 'bias'
    #flat_ext = 'flat'
    #science_ext = 'object'

    #n_first = 332896
    #n_last = 333614
    #framen = range( n_first, n_last+1 )
    #nframes = len( framen )
    #ddir_science = os.path.join( ddir, 'object' )
    #science_frames = []
    #for i in range( nframes ):
    #    science_frames += [ '0000{0:.0f}-20130520-OSIRIS-OsirisLongSlitSpectroscopy.fits'.format( framen[i] ) ]
    #science_frames = np.array( science_frames, dtype=str )


    stellar = stellar_module.stellar()
    stellar.ddir = ddir_science
    stellar.adir = adir
    stellar.star_names = star_names
    stellar.science_images_full_list = science_images_full_list_filename
    stellar.science_images_list = science_images_list_filename
    stellar.badpix_maps_full_list = badpix_maps_full_list_filename
    stellar.badpix_maps_list = badpix_maps_list_filename
    stellar.science_traces_list = science_traces_list_filename
    stellar.science_spectra_list = science_spectra_list_filename

    # Generate the list of filenames that the bad pixel maps
    # will be saved to:
    science_images_list_filepath = os.path.join( stellar.adir, stellar.science_images_full_list )
    science_image_filenames = np.loadtxt( science_images_list_filepath, dtype=str )
    nimages = len( science_image_filenames )
    badpix_map_filenames = []
    for i in range( nimages ):
        badpix_map_filenames += [ 'bp{0}'.format( science_image_filenames[i][1:] ) ]
    stellar.badpix_maps_full_list = badpix_maps_full_list_filename
    badpix_opath = os.path.join( stellar.adir, stellar.badpix_maps_full_list )
    np.savetxt( badpix_opath, badpix_map_filenames, fmt='%s' )
    #np.savetxt( stellar.science_images_list, science_frames, fmt='%s' ) # this should be deleted?
    stellar.badpix_map = None
    stellar.disp_axis = disp_axis
    stellar.crossdisp_axis = crossdisp_axis
    stellar.n_exts = n_exts

    stellar.tracefit_kwargs = { 'method':'linear_interpolation', 'binw_disp':30 }

    stellar.wcal_kws = {}
    stellar.wcal_kws['marc_file'] = os.path.join( ddir_arc, marc_filename )
    stellar.wcal_kws['fiducial_lines'] = wavcal_fiducial_lines
    stellar.wcal_kws['crossdisp_pixbounds'] = wavcal_crossdisp_bounds
    stellar.wcal_kws['disp_pixbounds'] = wavcal_disp_bounds

    # Bounding boxes for the stars in pixel coordinates:
    if n_exts>1:
        stellar.nstars = []
        for i in range( n_exts ):
            stellar.nstars += [ len( star_names[i] ) ]
    else:
        stellar.nstars = len( star_names )
                                
    stellar.crossdisp_bounds = crossdisp_bounds
    stellar.disp_bounds = disp_bounds

    # SHOULD THESE IMAGES BE CALIBRATED (DARK, FLAT ETC) BEFORE THIS STEP?
    #image_filenames = np.loadtxt( science_images_list, dtype=str )
    #nimages = len( image_filenames )
    stellar.jds = np.ones( nimages )
    stellar.exptime_secs = np.ones( nimages )
    for i in range( nimages ):
        print i+1, nimages
        image = os.path.join( stellar.ddir, science_image_filenames[i] )
        hdu = fitsio.FITS( image )
        h0 = hdu[0].read_header()
        stellar.jds[i] = h0['MJD-OBS'] + 2400000.5
        stellar.exptime_secs[i] = h0['EXPTIME']
    stellar.header_kws['gain'] = 'GAIN'    

    return stellar


def crosscor( spectra_in, ref_spectrum, max_shift=1, resample_factor=10, disp_bound_ixs=[], spectra_smoothing_sig=None, continuum_smoothing_sig=None, remove_continuum=False ):
    spectra = spectra_in.copy()
    # Max shift has to be an integer number of pixels:
    max_shift = int( np.round( max_shift ) )
    if spectra_smoothing_sig!=None:
        ref_spectrum_smooth = scipy.ndimage.filters.gaussian_filter1d( ref_spectrum, spectra_smoothing_sig )
    else:
        ref_spectrum_smooth = ref_spectrum
    if remove_continuum==True:
        ref_spectrum_smooth -= scipy.ndimage.filters.gaussian_filter1d( ref_spectrum_smooth, continuum_smoothing_sig )
    ix0 = disp_bound_ixs[0]
    ix1 = disp_bound_ixs[1]
    nframes, ndisp = np.shape( spectra )
    x = np.arange( ndisp )
    dshift = 1./resample_factor
    nshifts = 0
    while nshifts*dshift<max_shift:
        nshifts += 1
    max_shift = int( nshifts*dshift )
    dshifts = np.r_[ -max_shift : max_shift+dshift : dshift ]
    nshifts = len( dshifts )
    npad = max_shift# + 1
    xi = np.arange( -npad, ndisp+npad )
    spectra_padded = np.zeros( [ nframes, npad+ndisp+npad ] )
    if spectra_smoothing_sig!=None:
        for j in range( nframes ):
            spectra[j,:] = scipy.ndimage.filters.gaussian_filter1d( spectra[j,:], spectra_smoothing_sig )
    spectra_padded[:,npad:npad+ndisp] = spectra
    pad = np.zeros( npad )
    ref_spectrumi = np.concatenate( [ pad, ref_spectrum_smooth, pad ] )
    print '\nInterpolating reference spectrum onto super-sampled grid...'
    interpolation = 'linear'
    interpf = scipy.interpolate.interp1d( xi, ref_spectrumi, kind=interpolation )
    print 'Done.'
    shifted = np.zeros( [ nshifts, ndisp ] )
    for j in range( nshifts ):
        shifted[j,:] = interpf( x+dshifts[j] )
        if spectra_smoothing_sig!=None:
            shifted[j,:] = scipy.ndimage.filters.gaussian_filter1d( shifted[j,:], spectra_smoothing_sig )
        if remove_continuum==True:
            shifted[j,:] -= scipy.ndimage.filters.gaussian_filter1d( shifted[j,:], continuum_smoothing_sig )

    # Now loop over the individual spectra and determine which of the 
    # shifted reference spectra gives the best match:
    wavshifts = np.zeros( nframes )
    vstretches = np.zeros( nframes )
    dspec = np.zeros( [ nframes, ndisp ] )
    spectra_shifted = np.zeros( [ nframes, ndisp ] )
    A = np.ones( [ ndisp, 2 ] )
    print '\nDetermining wavelength shifts:'
    for i in range( nframes ):
        rms_i = np.zeros( nshifts )
        diffs = np.zeros( [ nshifts, ndisp ] )
        vstretches_i = np.zeros( nshifts )
        spectrum_smooth_i = spectra[i,:]
        target = spectrum_smooth_i
        if remove_continuum==True:
            target_continuum = scipy.ndimage.filters.gaussian_filter1d( target, continuum_smoothing_sig )
            target -= target_continuum
        fits = []
        for j in range( nshifts ):
            basis = shifted[j,:]
            if dshifts[j]==0:
                plt.figure()
                plt.plot(target,'-r')
                plt.plot(basis,'-k')
                pdb.set_trace()
                
            A[:,1] = basis
            b = np.reshape( target, [ ndisp, 1 ] )
            res = np.linalg.lstsq( A, b )
            c = res[0].flatten()
            fit = np.dot( A, c )
            vstretches_i[j] = c[1]
            diffs[j,:] = target - fit
            rms_i[j] = np.sqrt( np.mean( diffs[j,:][ix0:ix1+1]**2. ) )
            fits += [ fit ]
        ix = np.argmin( rms_i )
        dspec[i,:] = diffs[ix,:]/ref_spectrum_smooth # this needs updating for remove_continuum==True
        wavshifts[i] = dshifts[ix]
        vstretches[i] = vstretches_i[ix]
        interpf_i = scipy.interpolate.interp1d( xi, spectra_padded[i,:], kind=interpolation )
        spectra_shifted[i,:] = interpf_i( x-dshifts[ix] )
        print '... frame {0} of {1} --> {2:.3f} pix'.format( i+1, nframes, dshifts[ix] )
        if 0*(np.abs(dshifts[ix])>0.2):
            plt.figure()
            plt.plot(target,'-k',label='target science spectrum')
            plt.plot(fits[ix],'-r',label='shifted')
            plt.plot(ref_spectrum_smooth,'-c',label='unshifted')
            plt.legend()
            pdb.set_trace()
        if 0*np.abs(dshifts[ix])>0.2:
            plt.figure()
            A[:,1] = shifted[ix,:]
            c = np.linalg.lstsq( A, b )[0].flatten()
            fitbest = np.dot( A, c )
            plt.plot( spectrum_smooth_i, '-k', label='target to fit' )
            plt.plot( fitbest, '-r', label='should be good' )
            plt.plot( fits[ix], '--c' )
            A[:,1] = ref_spectrum_smooth 
            c = np.linalg.lstsq( A[ix0:ix1+1,:], b[ix0:ix1+1,:] )[0].flatten()
            fitnothing = np.dot( A, c )
            plt.plot( fitnothing, '-m', label='should be worse' )
            #plt.plot( ref_spectrum_smooth, '-g' )
            plt.legend()
            plt.axvline(ix0)
            plt.axvline(ix1)
            pdb.set_trace()
                  
        #z_spectrum_smooth = spectrum_smooth_i # this one, should align with this one:
        #z_ref_spectrum_shifted_smooth = scipy.ndimage.filters.gaussian_filter1d( spectra[i,:], smoothing_sig )
        # but this one should be offset:
        #z_ref_spectrum_smooth = ref_spectrum_smooth
        #plt.plot( z_spectrum_smooth/vstretches[i], '-k' )
        #plt.plot( z_ref_spectrum_shifted_smooth, '--r' )
        ##plt.plot( spectra_shifted[i,:], '-g' )
        #plt.plot( z_ref_spectrum_smooth, '-g' )
        #pdb.set_trace()
    return wavshifts, spectra_shifted, dspec
    
    
    
