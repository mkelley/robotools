######################################################################
def copy_files(source_dir, target_dir, config, logger):
    import os
    import re
    import shutil
    from glob import glob
    import astropy.units as u
    from astropy.io import fits
    from ccdproc import CCDData, gain_correct

    pat = os.sep.join((source_dir.strip('/'), '*'))
    in_list = [f for f in glob(pat)
               if len(re.findall(config.file_template, f)) == 1]

    logger.info('Found {} raw files in {}.'.format(len(in_list), source_dir))
    assert len(in_list) > 0, 'No files matched: {}'.format(
        config.file_template)

    files = [os.path.basename(fn) + '.fits' for fn in in_list]
    copied = 0
    for ifn, ofn in zip(in_list, files):
        ofn = os.sep.join((target_dir, ofn))
        if os.path.exists(ofn) and not config.reprocess_data:
            continue

        shutil.copy(ifn, ofn)
        os.chmod(ofn, 0o644)

        logger.debug('{} -> {}'.format(ifn, ofn))
        copied += 1

        # header fix
        ccd = CCDData.read(ofn, unit=u.adu)
        if ccd.meta['OBJECT'] == '2016R2 Pan-STARRS':
            ccd.meta['OBJECT'] = 'C/2016 R2 (PanSTARRS)'
        elif ccd.meta['OBJECT'] == '2016M1 Pan-STARRS':
            ccd.meta['OBJECT'] = 'C/2016 M1 (PanSTARRS)'
        elif ccd.meta['OBJECT'] == '2017T1 Heinze':
            ccd.meta['OBJECT'] = 'C/2017 T1 (Heinze)'
        elif ccd.meta['OBJECT'] == '2019Y1 ATLAS':
            ccd.meta['OBJECT'] = 'C/2019 Y1 (ATLAS)'
        elif ccd.meta['OBJECT'] == '21P Gia-Zin':
            ccd.meta['OBJECT'] = '21P/Giacobini-Zinner'
        elif ccd.meta['OBJECT'] == '29P Schwas-Wach':
            ccd.meta['OBJECT'] = '29P/Schwassmann-Wachmann 1'
        elif ccd.meta['OBJECT'] == '29P Schwas-Wach (Low':
            ccd.meta['OBJECT'] = '29P/Schwassmann-Wachmann 1'
        elif ccd.meta['OBJECT'] == '38P Ste-Ote':
            ccd.meta['OBJECT'] = '38P/Stephan-Oterma'
        elif ccd.meta['OBJECT'] == '46P Wirtanen':
            ccd.meta['OBJECT'] = '46P/Wirtanen'
        elif ccd.meta['OBJECT'] == '64P Swift-Gehrels':
            ccd.meta['OBJECT'] = '64P/Swift-Gehrels'
        elif ccd.meta['OBJECT'] == '123P Wes-Har':
            ccd.meta['OBJECT'] = '123P/West-Hartley'
        elif ccd.meta['OBJECT'] == '123P West-Hartley':
            ccd.meta['OBJECT'] = '123P/West-Hartley'
        elif ccd.meta['OBJECT'] == '2018Y1 Iwamoto':
            ccd.meta['OBJECT'] = 'C/2018 Y1 (Iwamoto)'
        elif ccd.meta['OBJECT'] == '2019Q4 Borisov (Lowe':
            ccd.meta['OBJECT'] = 'C/2019 Q4 (Borisov)'
        elif ccd.meta['OBJECT'] == '2019Q4 Borisov':
            ccd.meta['OBJECT'] = 'C/2019 Q4 (Borisov)'
        elif ccd.meta['OBJECT'] == '2I Borisov':
            ccd.meta['OBJECT'] = 'C/2019 Q4 (Borisov)'
        elif ccd.meta['OBJECT'] == '155P Shoemaker 3':
            ccd.meta['OBJECT'] = '155P/Shoemaker 3'

        ccd.meta['FILTER'] = (ccd.meta['FILTER1'] +
                              ccd.meta['FILTER2']).replace('Open', '').strip()
        ccd = gain_correct(ccd, ccd.meta['gain'], gain_unit=u.electron / u.adu)
        ccd.write(ofn, overwrite=True)

    logger.info('{} files copied and gain corrected.'.format(copied))
    return files

######################################################################


def load_files(target_dir, config, files, logger):
    from ccdproc import ImageFileCollection

    keywords = ['object', 'imagetyp', 'date-obs', 'telra', 'teldec',
                'airmass', 'filter', 'oairtemp', 'relhum',
                'subbias', 'flatcor']
    ic = ImageFileCollection(location=target_dir, filenames=files,
                             keywords=keywords)
    nbias = len(ic.files_filtered(imagetyp='bias'))
    nflat = len(ic.files_filtered(imagetyp='flat'))
    ndata = len(ic.files_filtered(imagetyp='object'))

    filters = []
    for h in ic.headers():
        if h['IMAGETYP'] in ['object', 'flat']:
            if h['FILTER'] not in filters:
                filters.append(h['FILTER'])

    filters.sort()
#    flat_breakdown = []
#    data_breakdown = []
#    for filt in filters:
#        for k in config.flats.values():
#            flat_breakdown.append('{} {} {}'.format(
#                len(ic.files_filtered(
#                    imagetyp='flat', object=k, filter=filt)),
#                filt, k))
#        data_breakdown.append('{} {}'.format(
#            len(ic.files_filtered(imagetyp='OBJECT', filter=filt)), filt))

    logger.info('{} files: {} bias, {} flats, {} object ({} filters)'.format(
        len(ic.files), nbias, nflat, ndata, len(filters)))

    return ic

######################################################################


def find_last(target_dir, config, fn):
    '''fn may be a tuple of file names in order of precedence.'''

    import os
    from glob import glob

    if isinstance(fn, str):
        fn = (fn, )
    else:
        assert isinstance(fn, (tuple, list))

    pat = os.sep.join((target_dir.strip('/'), '..', '..',
                       '20[12][0-9][01][0-9][0123][0-9]',
                       'ppp'))
    src_directories = sorted([os.path.abspath(p) for p in glob(pat)])
    i = src_directories.index(os.path.abspath(target_dir))
    assert i > 0, 'No older directories found.'

    # step through in reverse order to find last file
    for src in src_directories[:i][::-1]:
        for f in (os.sep.join((src, f)) for f in fn):
            if os.path.exists(f):
                break
            else:
                f = None

        if f is not None:
            break

    assert f is not None, '{} not found in older directories.'.format(fn)

    return f

######################################################################


def bias_correction(target_dir, ic, config, logger):
    import os
    import numpy as np
    import astropy.units as u
    from ccdproc import CCDData, combine
    from ccdproc import subtract_bias, subtract_overscan, trim_image

    logger.debug('Bias correction.')
    fn = os.sep.join((target_dir, 'bias.fits'))
    if os.path.exists(fn) and not config.reprocess_all:
        logger.info('Read bias = {}.'.format(fn))
        bias = CCDData.read(fn)
    elif len(ic.files_filtered(imagetyp='bias')) == 0:
        logger.warning('No bias files provided.')
        logger.warning('{} not found'.format(fn))
        try:
            last_bias = find_last(target_dir, config, 'bias.fits')
            bias = CCDData.read(last_bias)
            logger.info('Using {}.'.format(last_bias))
        except AssertionError:
            logger.warning('No previous bias found.  Not subtracting bias.')
            bias = 0 * u.electron
    else:
        logger.info('Create bias = {}.'.format(fn))
        files = ic.files_filtered(include_path=True, imagetyp='bias')
        bias = combine(files, method='average', clip_extrema=True, nlow=1,
                       nhigh=1)
        bias.meta['FILENAME'] = os.path.basename(fn)

        n = str([int(f.split('.')[-2]) for f in files])
        bias.meta.add_history(
            'Created from file numbers: {}'.format(n))

        bias.write(fn, overwrite=True)

        logger.info('Bias subtract and trim data.')

    i = ic.summary['subbias'].mask
    logger.debug('Bias subtract {} files.'.format(sum(i)))
    meanbias = np.mean(bias)
    for fn in ic.summary['file'][i]:
        ccd = CCDData.read(os.sep.join([ic.location, fn]), unit='electron')
        logger.debug(fn)

        ccd = subtract_bias(ccd, bias)
        ccd.meta['BIASFILE'] = (bias.meta['FILENAME'],
                                'Name of the bias frame used.')
        ccd.meta['MEANBIAS'] = (meanbias, 'Mean bias level, electrons')

        if config.overscan_correct:
            overscan = CCDData(np.hstack([ccd[:, :50], ccd[:, 2115:]]),
                               unit=ccd.unit)
            ccd = subtract_overscan(ccd, overscan, median=True)

        for s in config.trim_sections:
            ccd = trim_image(ccd[s])

        ccd.write(os.sep.join([ic.location, fn]), overwrite=True)

    ic.refresh()
    return ic

######################################################################


def mode_scaler(im):
    import numpy as np
    return 1 / (3 * np.ma.median(im) - 2 * np.ma.mean(im))

######################################################################


def flatfield_correction(target_dir, ic, config, logger):
    import os
    import numpy as np
    import scipy.ndimage as nd
    from ccdproc import CCDData, combine, cosmicray_lacosmic, flat_correct

    logger.debug('Flat fields.')

    i = ((ic.summary['imagetyp'] == 'object')
         + (ic.summary['imagetyp'] == 'flat'))
    filters = np.unique(ic.summary['filter'][i].data.data)

    flats = dict()
    for filt in filters:
        flats[filt] = 1
        for flat_key, flat_name in config.flats.items():
            fn = '{}-{}.fits'.format(flat_key, filt)
            fn = os.sep.join((target_dir, fn))
            if os.path.exists(fn) and not config.reprocess_all:
                logger.info('Reading {}.'.format(fn))
                flats[filt] = CCDData.read(fn)
            elif len(ic.files_filtered(object=flat_name, filter=filt)) == 0:
                logger.warning(
                    'No {} files provided for {} and {} not found.'.format(flat_name, filt, fn))
            else:
                logger.info('Generating {}.'.format(fn))
                files = ic.files_filtered(include_path=True,
                                          object=flat_name,
                                          filter=filt)
                flat = combine(files, method='median', scale=mode_scaler)
                flat.mask = (flat.data > 1.2) + (flat.data < 0.8)

                n = str([int(f.split('.')[-2]) for f in files])
                flat.meta.add_history(
                    'Created from file numbers: {}'.format(n))
                flat.meta['FILENAME'] = os.path.basename(fn)[1]

                flat.write(fn, overwrite=True)

                if flats[filt] == 1:
                    logger.info('Using {} for {}'.format(
                        flat_name.lower(), filt))
                    flats[filt] = flat

        if flats[filt] == 1:
            # find last flat
            fn = tuple(('{}-{}.fits'.format(k, filt) for k in config.flats))
            try:
                last_flat = find_last(target_dir, config, fn)
                flats[filt] = CCDData.read(last_flat)
                logger.info('Using {}.'.format(last_flat))
            except AssertionError:
                logger.warning(
                    'No previous flat found.  Not flat correcting {} data.'.format(filt))

    i = ((ic.summary['imagetyp'] != 'bias')
         & ic.summary['flatcor'].mask
         & ~ic.summary['subbias'].mask)
    logger.info('{} files to flat correct.'.format(sum(i)))
    for fn in ic.summary['file'][i]:
        ccd = CCDData.read(os.sep.join([ic.location, fn]))

        filt = ccd.meta['FILTER']
        if flats[filt] == 1:
            logger.debug(
                '{} skipped (no {} flat field provided).'.format(fn, filt))
            continue

        logger.debug(fn)
        ccd = flat_correct(ccd, flats[filt])

        if config.lacosmic:
            cleaned = cosmicray_lacosmic(
                ccd, pssl=ccd.meta['meanbias'],
                sigclip=3.5,
                satlevel=config.saturation,
                readnoise=ccd.meta['rdnoise'])
            ccd.mask += nd.binary_dilation(cleaned.mask)
            ccd.header['LACOSMIC'] = 1, 'L.A.Cosmic processing flag.'
        else:
            ccd.header['LACOSMIC'] = 0, 'L.A.Cosmic processing flag.'

        ccd.meta['FLATFILE'] = (flats[filt].meta['FILENAME'],
                                'Name of the flat field correction used.')
        ccd.write(os.sep.join([ic.location, fn]), overwrite=True)

    ic.refresh()
    return ic

######################################################################


def array_corrections(ic, config, logger):
    import numpy as np
    import astropy.units as u
    from ccdproc import CCDData
    from mskpy import rebin

    logger.debug('Array corrections.')

    nmodified = 0
    for f in ic.files_filtered(include_path=True, imagetyp='object'):
        ccd = CCDData.read(f)

        modified = False

        if ccd.meta.get('SATLEVEL') is None:
            modified = True
            ccd.meta['SATLEVEL'] = config.saturation, 'Saturation level, unbinned pixels, e-'
            binning = ccd.meta.get('REBIN', 1)
            if binning > 1:
                ccd.meta.add_comment(
                    'Warning: Saturation masking binned pixels.')
            i = ccd > config.saturation * binning**2 * u.electron
            if np.any(i):
                ccd.mask[i] = True

        if ccd.meta.get('XFLIP') is None:
            modified = True
            ccd.meta['XFLIP'] = 0, 'Data was not x-flipped'
            if config.xflip:
                ccd.data = ccd.data[:, ::-1]
                if ccd.mask is not None:
                    ccd.mask = ccd.mask[:, ::-1]
                ccd.meta['XFLIP'] = 1, 'Data was x-flipped'

        if ccd.meta.get('ROTATE') is None:
            modified = True
            ccd.meta['ROTATE'] = config.rotate, 'Data was not rotated CCW'
            if config.rotate:
                ccd.data = np.rot90(ccd.data, config.rotate // 90, axes=(1, 0))
                if ccd.mask is not None:
                    ccd.mask = np.rot90(ccd.mask, config.rotate // 90,
                                        axes=(1, 0))
                ccd.meta['ROTATE'] = config.rotate, 'Data was rotated CCW'

        if ccd.meta.get('REBIN') is None:
            modified = True
            ccd.meta['REBIN'] = config.rebin, 'Data was not rebinned.'
            if config.rebin != 1:
                ccd.data = rebin(ccd.data, -config.rebin, flux=True, trim=True)
                if ccd.mask is not None:
                    m = rebin(ccd.mask, -config.rebin,
                              flux=True, trim=True) > 0
                    ccd.mask = m

                ccd.meta['REBIN'] = config.rebin, 'Data was rebinned.'

        if modified:
            ccd.write(f, overwrite=True)
            logger.debug('  {} updated.'.format(f))
            nmodified += 1

    logger.debug(
        '{} files updated (includes header additions).'.format(nmodified))

######################################################################


def solve_wcs(ic, config, logger):
    import os
    import subprocess
    from astropy.io import fits
    from astropy.wcs import WCS
    from astropy.coordinates import SkyCoord

    logger.debug('Solve WCS.')

    os.system('mkdir -p /tmp/robo-rx')
    cmd = 'solve-field --dir /tmp/robo-rx --overwrite -L 0.8 -H 1.0 -u arcsecperpix --ra {ra} --dec {dec} --radius 10 --no-verify --no-tweak {filename}'

    files = ic.files_filtered(include_path=True, imagetyp='object')
    nsolved = 0
    for f in files:
        with fits.open(f, 'update') as hdu:
            logger.debug(f)
            if hdu[0].header.get('WCSSOLVE', False) and not config.reprocess_data:
                continue

            # use masked data
            masked = '/tmp/robo-rx/{}'.format(os.path.basename(f))
            try:
                mask = hdu[1].data.astype(bool)
            except IndexError:
                # missing mask, probably not flat-fielded
                mask = False

            fits.writeto(masked, hdu[0].data * (~mask), hdu[0].header,
                         overwrite=True)

            try:
                subprocess.check_call(
                    cmd.format(ra=hdu[0].header['telra'],
                               dec=hdu[0].header['teldec'],
                               filename=masked).split())
            except subprocess.CalledProcessError as e:
                logger.error(
                    'astrometry.net solve-field error:'.format(str(e)))
                continue

            wcs = WCS(masked.replace('.fits', '.wcs'))
            hdu[0].header.update(wcs.to_header())
            hdu[0].header['WCSSOLVE'] = 1, 'Updated WCS solution with astrometry.net'
            nsolved += 1

            c = SkyCoord(*wcs.wcs.crval, unit='deg')
            radec = c.to_string('hmsdms', sep=':')
            logger.info('{} new WCS solution: {}'.format(f, radec))

    logger.debug('{} files updated'.format(nsolved))
