__all__ = [
    'RoboPhotException',
    'ObjectFilesError',
    'BackgroundError',
    'CalibrationError',
    'CatalogError',
    'SourceListError',
    'ImageError',
    'photometry',
    'find_files',
    'plot'
]

import os
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import pixel_to_skycoord
from astropy.modeling import models, fitting
from astropy.stats import sigma_clip, sigma_clipped_stats
from astropy.table import Table
from ccdproc import CCDData, ImageFileCollection
import sep
from calviacat import PanSTARRS1, SkyMapper
from mskpy import bgphot
from . import background


class RoboPhotException(Exception):
    pass


class ObjectFilesError(RoboPhotException):
    pass


class BackgroundError(RoboPhotException):
    pass


class CalibrationError(RoboPhotException):
    pass


class CatalogError(RoboPhotException):
    pass


class SourceListError(RoboPhotException):
    pass


class ImageError(RoboPhotException):
    pass

######################################################################


def find_files(path, db, config, logger):
    files = glob(os.sep.join((path, config.file_template)))
    files = sorted([os.path.basename(f) for f in files])

    keywords = ['object', 'imagetyp', 'date-obs', 'ut', 'telra', 'teldec',
                'airmass', 'filter', 'exptime', 'oairtemp', 'relhum',
                'seeing', 'subbias', 'flatcor']
    ic = ImageFileCollection(location=path, filenames=files,
                             keywords=keywords)

    files = sorted(list(ic.files_filtered(imagetyp='object')))
    logger.debug('Found {} object files in {}.'.format(len(files), path))
    if len(files) == 0:
        raise ObjectFilesError('No object files found in {}.'.format(path))

    ic = ImageFileCollection(location=path,
                             filenames=files,
                             keywords=keywords)

    rows = []
    for row in ic.summary:
        fn = row['file'].replace('.fits', '')
        c = db.execute('SELECT * FROM obs WHERE filename=?', [fn])
        if c.fetchone() is None or config.update_obs:
            with fits.open(os.sep.join((path, row['file']))) as hdu:
                wcs = WCS(hdu[0])
                shape = hdu[0].data.shape

            if ((wcs.wcs.crpix[0] == 0)
                or (wcs.wcs.crpix[0] == 1)
                    or ('PIXEL' in wcs.wcs.ctype)):
                # missing WCS solution
                ra, dec = '0', '0'
            else:
                ra, dec = wcs.all_pix2world(shape[1] / 2, shape[0] / 2, 0)
                c = SkyCoord(ra, dec, unit='deg')
                ra, dec = c.to_string('hmsdms', sep=':').split()

            rows.append((fn, row['object'],
                         row['date-obs'] + ' ' + row['ut'],
                         row['airmass'], row['filter'], row['exptime'],
                         ra, dec, row['oairtemp'], row['relhum'],
                         row['seeing'] / 2))

    if len(rows) > 0:
        if config.update_obs:
            cols = list(zip(*rows))
            update = zip(*(cols + [cols[0]]))
            logger.info('Updating with FITS header meta data.')
            db.executemany('''
            UPDATE OR IGNORE obs SET
              filename=?,object=?,obsdate=?,airmass=?,filter=?,exptime=?,
              ra=?,dec=?,oairtemp=?,relhum=?,seeing=?
            WHERE filename=?
            ''', update)

        # insert whatever wasn't updated
        db.executemany('''
        INSERT OR IGNORE INTO obs
          (filename,object,obsdate,airmass,filter,exptime,ra,dec,
          oairtemp,relhum,seeing)
        VALUES
          (?,?,?,?,?,?,?,?,?,?,?)''', rows)

    return ic

########################################################################


def border_mask(ccd, config):
    import numpy as np
    if config.border_mask > 0:
        if ccd.mask is None:
            ccd.mask = np.zeros_like(ccd.data, bool)

        ccd.mask[:config.border_mask] = True
        ccd.mask[-config.border_mask:] = True
        ccd.mask[:, :config.border_mask] = True
        ccd.mask[:, -config.border_mask:] = True

    return ccd

########################################################################


def photometry(path, basename, db, config, logger):
    """Master Photometry Program"""
    imf = os.sep.join((path, basename)) + '.fits'
    bgf = imf.replace('.fits', '.bg.fits')
    catf = imf.replace('.fits', '.cat.fits')

    # Background check
    row = db.execute('''
      SELECT count() FROM obs
      WHERE bg IS NOT null AND filename=?
    ''', [basename]).fetchone()
    count = row[0]
    run_bg = any((count == 0, config.reprocess_bg, not os.path.exists(bgf)))

    # Photometry check
    run_phot = any((run_bg, not os.path.exists(catf), config.reprocess_phot))

    run_cal = any((run_phot, config.reprocess_cal))

    if not any((run_bg, run_phot, run_cal)):
        return

    logger.info(basename)
    ccd = border_mask(CCDData.read(imf), config)
    ccd.mask[np.isnan(ccd.data)] = True

    if run_bg:
        background.background(ccd, bgf, db, config, logger)

    if run_phot:
        catalog(ccd, bgf, catf, db, config, logger)

    if run_cal:
        calibrate(ccd, basename, imf, catf, db, config, logger)

########################################################################


def catalog(ccd, bgf, catf, db, config, logger):
    bg = fits.open(bgf)
    im = ccd.data - bg['background'].data
    ps = ccd.meta['SCALE'] * ccd.meta.get('REBIN', 1)

    bgrms = bg['background'].header['bgrms']
    objects = sep.extract(im, 2, err=bgrms, mask=ccd.mask)
    logger.info('Found {} sources.'.format(len(objects)))

    rap = max(ccd.meta['SEEING'] * 2, 5 / ps)
    flux, fluxerr, flag = sep.sum_circle(
        im, objects['x'], objects['y'], rap, err=bgrms)

    # avoid theta rounding error
    theta = np.maximum(np.minimum(objects['theta'], np.pi / 2.00001),
                       -np.pi / 2.00001)
    kronrad, krflag = sep.kron_radius(
        im, objects['x'], objects['y'], objects['a'], objects['b'],
        theta, 6.0)
    krflux, krfluxerr, _flag = sep.sum_ellipse(
        im, objects['x'], objects['y'], objects['a'], objects['b'],
        theta, 2.5 * kronrad, subpix=1, err=bgrms)
    krflag |= _flag

    # an additional background estimate, which should help when there are
    # large extended sources in scene: IN TESTS, THIS DID NOT AFFECT RESULTS
    # for i in range(len(objects)):
    #    krflux[i], krfluxerr[i] = bg_subtract2(im, objects[i], krflux[i],
    #                                           krfluxerr[i])
    #    flux[i], fluxerr[i] = bg_subtract2(im, objects[i], flux[i],
    #                                       fluxerr[i], r=rap)

    if ccd.wcs.wcs.crval[0] == ccd.wcs.wcs.crval[1]:
        ra, dec = np.zeros((2, len(objects)))
    else:
        ra, dec = ccd.wcs.all_pix2world(objects['x'], objects['y'], 0)

    tab = Table((objects['x'], objects['y'], ra, dec, flux, fluxerr,
                 flag, objects['a'], objects['b'], theta,
                 kronrad, krflux, krfluxerr, krflag),
                names=('x', 'y', 'ra', 'dec', 'flux', 'fluxerr', 'flag',
                       'a', 'b', 'theta', 'kronrad', 'krflux',
                       'krfluxerr', 'krflag'))

    hdu = fits.HDUList()
    hdu.append(fits.BinTableHDU(tab, name='cat'))
    hdu['cat'].header['RADIUS'] = (
        rap * ps, 'aperture photometry radius, arcsec')
    hdu.writeto(catf, overwrite=True)

########################################################################


def bg_subtract2(im, obj, flux, fluxerr, r=None):
    if r is None:
        a = obj['a']
        b = obj['b']
    else:
        a = r
        b = r

    theta = max(min(obj['theta'], np.pi / 2.00001), -np.pi / 2.00001)

    outer = np.zeros(im.shape, bool)
    sep.mask_ellipse(outer, obj['x'], obj['y'], a, b, theta, r=5)

    inner = np.zeros(im.shape, bool)
    sep.mask_ellipse(inner, obj['x'], obj['y'], a, b, theta, r=3)

    area = np.zeros(im.shape, bool)
    sep.mask_ellipse(area, obj['x'], obj['y'], a, b, theta)
    area = area.sum()

    mask = outer * ~inner
    bgap = im * mask
    bgap = bgap[bgap != 0]
    mms = sigma_clipped_stats(bgap)

    flux -= mms[1] * area
    fluxerr = np.sqrt(fluxerr**2 + area * mms[2] * (1 + area / len(bgap)))
    return flux, fluxerr

########################################################################


def calibrate(ccd, basename, imf, catf, db, config, logger):
    # only calibrate R frames
    filt = fits.getval(imf, 'FILTER')
    if filt != 'R':
        return

    if ccd.wcs.wcs.crval[0] == ccd.wcs.wcs.crval[1]:
        logger.info('Calibration: None (no WCS)')
        return

    phot = Table(fits.getdata(catf))
    phot = phot[phot['krflag'] == 0]
    sources = SkyCoord(phot['ra'], phot['dec'], unit='deg')

    # default: calibrate with PS1, otherwise SkyMapper
    for Catalog in (PanSTARRS1, SkyMapper):
        try:
            cat = Catalog('cat.db', match_limit=config.match_limit,
                          logger=logger)
            check_coverage(cat, sources, ccd.shape, ccd.wcs)
            break
        except CatalogError:
            cov = False

    try:
        objids, distances = cat.xmatch(sources)
    except TypeError:
        raise CalibrationError('Catalog cross-matching failed.')

    m_inst = -2.5 * np.log10(phot['krflux'])
    m_err = phot['krfluxerr'] / phot['krflux'] * 1.0857

    i = m_err < 0.2
    if sum(i) == 0:
        raise CalibrationError('No sources passed SNR cut.')

    color = 'g-r'
    try:
        zp, C, unc, m, mcolor, gmi = cat.cal_color(
            objids[i], m_inst[i], filt.lower(), color, C=config.color_cor_r,
            mlim=config.mlim, gmi_lim=config.gmi_lim)
    except Exception as e:
        raise CalibrationError(str(e))

    j = ~m.mask
    ncal = len(m[j])
    if ncal == 0:
        raise CalibrationError('All catalog sources masked.')

    with fits.open(catf, mode='update') as hdu:
        tab = Table((objids[i], m_inst[i], m, mcolor, gmi),
                    names=('objid', 'm inst', 'm', color, 'g-i'))
        tab = tab[j]
        h = fits.Header({
            'catalog': str(cat).split()[0].split('.')[-1],
            'magzp': zp,
            'colorcor': C,
            'zpunc': unc,
            'color': color,
            'ncal': ncal
        })
        hdu.append(fits.BinTableHDU(tab, h, name='cal'))

    logger.info(
        "  Calibration: r = m_inst + {:.3f} + (g-r) * {:.3f}  ({:.3f})"
        .format(zp, C, unc))

    db.execute('''
      UPDATE obs SET zp=?,colorcor=?,ncal=?,calunc=?
      WHERE filename=?
    ''', (zp, C, ncal, unc, basename))


######################################################################


def check_coverage(cat, sources, shape, wcs):
    # check for objects in all corners
    test_cat = cat.search(sources)[1]
    if len(test_cat) == 0:
        n = np.r_[0]
    else:
        x, y = wcs.all_world2pix(test_cat.ra, test_cat.dec, 0)
        n = np.r_[
            sum((x < shape[1] / 3) * (y < shape[0] / 3)),
            sum((x > shape[1] * 2 / 3) * (y < shape[0] / 3)),
            sum((x > shape[1] * 2 / 3) * (y > shape[0] * 2 / 3)),
            sum((x < shape[1] / 3) * (y > shape[0] * 2 / 3))
        ]

    if any(n < cat.max_records / 50):
        cat.fetch_field(sources)
        if len(cat.search(sources)[0]) == 0:
            raise CatalogError('Not enough sources.')

######################################################################


def plot(path, basename, logger):
    hdu = fits.open(os.path.join(path, basename) + '.cat.fits')
    try:
        tab = Table(hdu['cal'].data)
    except KeyError:
        logger.error('Not plotting {} (file not calibrated)'
                     .format(basename))
        return

    zp = hdu['cal'].header['magzp']
    C = hdu['cal'].header['colorcor']
    zpunc = hdu['cal'].header['zpunc']
    color = hdu['cal'].header['color']
    ncal = hdu['cal'].header['ncal']

    plt.clf()
    ax = plt.gca()
    ax.scatter(tab[color].data, (tab['m'] - tab['m inst']).data,
               marker='.', color='k')

    x = np.linspace(0, 1.5)
    y = C * x + zp
    label = '{:.3f} + {:.3f}({}), (σ={:.3f}, n={})'.format(
        zp, C, color, zpunc, ncal)
    ax.plot(x, y, 'r-', label=label)
    ax.fill_between(x, y - zpunc, y + zpunc, color='r', alpha=0.33)

    plt.setp(ax, xlabel='${}$ (mag)'.format(color),
             ylabel=r'$m-m_{\rm inst}$ (mag)')
    plt.legend()
    plt.tight_layout()
    plt.show()
