import os

import numpy as np
from sbpy.data import Ephem
from astropy.coordinates import SkyCoord, Angle
from astropy.modeling import models, fitting
from astropy.stats import sigma_clip, sigma_clipped_stats
from astropy.time import Time
from astropy.io import fits
from astropy.nddata import StdDevUncertainty, NoOverlapError
import astropy.units as u
from astropy.wcs.utils import proj_plane_pixel_scales, pixel_to_skycoord
from ccdproc import CCDData
from photutils import CircularAperture, CircularAnnulus, aperture_photometry, centroid_2dg
from photutils.utils import calc_total_error, cutout_footprint
from sbpy.activity import Afrho

from ..database import fetchall


class RoboCometException(Exception):
    pass


class DataDiscoveryError(RoboCometException):
    pass


class EphemerisError(RoboCometException):
    pass


class PhotometryError(RoboCometException):
    pass


class CentroidingError(RoboCometException):
    pass

######################################################################


def find_data(db, config, logger):
    if config.reprocess_all:
        c = db.execute('''
        SELECT rowid,object,filename FROM obs WHERE object LIKE ? || "%"
        ORDER BY obsdate
        ''', [config.comet])
    elif config.reprocess is not None:
        c = db.execute('''
        SELECT rowid,object,filename FROM obs
        WHERE object LIKE ? || "%"
          AND obsdate LIKE ? || "%"
        ORDER BY obsdate
        ''', [config.comet, config.reprocess])
    else:
        c = db.execute('''
        SELECT obs.rowid,obs.object,filename FROM obs
        LEFT JOIN comet ON obs.rowid=comet.obsid
        WHERE obs.object LIKE ? || "%"
          AND comet.obsid IS null
        ORDER BY obsdate
        ''', [config.comet])

    rows = c.fetchall()
    if len(rows) == 0:
        raise DataDiscoveryError(
            'No data found for object {}'.format(config.comet))

    logger.info('Processing {} files for {}.'.format(len(rows), config.comet))

    obsid, obj, filename = [list(a) for a in list(zip(*rows))]
    data_dirs = set(['20' + f[:6] for f in filename])
    for data_dir in data_dirs:
        if os.path.exists(os.sep.join((data_dir, 'ppp', 'processing'))):
            for i in range(len(obsid) - 1, -1, -1):
                if data_dir in filename[i]:
                    del obsid[i]
                    del obj[i]
                    del filename[i]

    if len(obsid) == 0:
        raise DataDiscoveryError(
            'All data in directories currently being processed')

    logger.debug(
        '{} files after pruing directories currently being processed.'.format(len(obsid)))

    data = dict(zip(obsid, filename))
    db.executemany('''
      INSERT OR IGNORE INTO comet (object,obsid) VALUES (?,?)
    ''', zip(obj, obsid))

    return data

########################################################################


def ephemeris(db, config, logger):
    c = db.execute('''
      SELECT comet.rowid,obs.obsdate
      FROM comet
      INNER JOIN obs ON comet.obsid=obs.rowid
      WHERE comet.ra_jpl IS null
    ''')
    rows = c.fetchall()
    if len(rows) == 0:
        return

    logger.info('Updating ephemeris for {} observations.'.format(len(rows)))
    rowid, obsdate = list(zip(*rows))
    #q = callhorizons.query(config.horizons, closest_apparition=True)
    jd = Time(obsdate).jd
    # q.set_discreteepochs(jd)
    #n = q.get_ephemerides(688)
    eph = Ephem.from_horizons(config.horizons, id_type='designation',
                              location='688', epochs=jd,
                              closest_apparition=True)
    if len(eph) != len(jd):
        raise EphemerisError(
            'Expected {} rows from JPL/HORIZONS, but {} were returned.'.format(len(jd), n))

    sangle = Angle(eph['sunTargetPA'] - 180 * u.deg).wrap_at(360 * u.deg)
    vangle = Angle(eph['velocityPA'] - 180 * u.deg).wrap_at(360 * u.deg)

    c = db.executemany('''
      UPDATE comet SET
        ra_jpl=?,dec_jpl=?,ra_rate=?,dec_rate=?,
        rh=?,delta=?,phase=?,sangle=?,vangle=?
      WHERE rowid=?
    ''', zip(eph['RA'].to('deg').value, eph['DEC'].to('deg').value,
             eph['RA_rate'].to('arcsec/hr').value,
             eph['DEC_rate'].to('arcsec/hr').value,
             eph['r'].value, eph['delta'].value,
             eph['alpha'].to('deg').value,
             sangle.deg, vangle.deg,
             rowid))

########################################################################


def photometry(obsid, basename, db, config, logger):
    obs = db.execute('SELECT * FROM obs WHERE rowid=?', [obsid]).fetchone()

    if obs['zp'] is None:
        raise PhotometryError('Frame not calibrated.')

    obj = db.execute('SELECT rowid,* FROM comet WHERE obsid=?',
                     [obsid]).fetchone()

    sources = fetchall(db, ['rowid', 'xcentroid', 'ycentroid', 'ra', 'dec',
                            'nbg', 'bg', 'flux', 'func'], 'phot',
                       'filename=? AND flux IS NOT null',
                       [basename])

    f = '20{}/ppp/{}.fits'.format(basename[:6], basename)
    ccd = CCDData.read(f)
    bg = fits.open(f.replace('.fits', '.bg.fits'))
    ccd.data = (ccd.data - bg['background'].data) / ccd.meta['exptime']
    ccd.uncertainty = StdDevUncertainty(
        np.ones_like(ccd) * bg[0].header['STDEV']
        / ccd.meta['exptime'] * u.electron)

    # find comet
    xy = locate_comet(ccd, obj, obs, logger)
    logger.info('  Comet at {:.1f}, {:.1f}.'.format(*xy))
    c = pixel_to_skycoord(xy[0], xy[1], ccd.wcs)
    db.execute('''
    UPDATE comet SET x=?,y=?,ra=?,dec=?
      WHERE rowid=?
    ''', (xy[0], xy[1], c.ra.deg, c.dec.deg, obj['rowid']))

    obj = db.execute('SELECT rowid,* FROM comet WHERE obsid=?',
                     [obsid]).fetchone()

    date = obs['obsdate'][:10]
    cal = db.execute('SELECT rowid,* FROM cal WHERE obsdate=?',
                     [date]).fetchone()
    if cal['zp0'] is None:
        raise PhotometryError('Night not calibrated.')

    # measure comet
    comet_phot(ccd, sources, bg, obj, obs, cal, db, config, logger)

########################################################################


def locate_comet(ccd, obj, obs, logger):
    xyi = ccd.wcs.all_world2pix(obj['ra_jpl'], obj['dec_jpl'], 0)
    xy = xyi
    for scale in [1, 0.5, 0]:
        box = int(obs['seeing'] * 2**scale)
        if box < 3:
            break
        try:
            subim, submask, suberr, s = cutout_footprint(
                ccd.data, xy, box_size=box, mask=ccd.mask,
                error=ccd.uncertainty)
        except NoOverlapError:
            raise CentroidingError(
                'Cannot centroid comet.  Expected at (x, y) {}, last iteration {}.'.format(xyi, xy))

        nxy = centroid_2dg(subim, submask) + np.r_[s[1].start, s[0].start]
        d = xy - nxy
        d = np.sqrt(np.sum(d**2))

        xy = nxy
        if d < 0.25:
            break

    return [float(z) for z in xy]

########################################################################


def comet_phot(ccd, sources, bg, obj, obs, cal, db, config, logger):
    zp = cal['zp0'] - cal['extcor'] * obs['airmass'] - \
        cal['colorcor'] * config.gmr_sun
    zp_unc = max(cal['zp0_std'], cal['zp0_unc'])

    # Revise background estimation with nearby stars.  Which stars?
    # Select based on rough photometry of the comet: we want stars
    # outside the radius at which the surface brightness is 0.5 sigma
    # per pixel.
    ap = CircularAperture((obj['x'], obj['y']), r=obs['seeing'])
    rough_flux = aperture_photometry(ccd.data, ap)
    C = rough_flux['aperture_sum'].data[0] / 2 / np.pi / ap.r
    d = C / 0.5 / obs['bgstdev'] * obs['exptime']
    d = min(max((100, C / 0.5 / obs['bgstdev'] * obs['exptime'])), 500)

    delta_bg = 0
    if len(sources) > 0:
        r = np.sqrt((sources['xcentroid'] - obj['x'])**2
                    + (sources['ycentroid'] - obj['y'])**2)
        i = r > d
    else:
        i = []

    if i.sum() < 10:
        logger.info('  Using large aperture background.')
        apmask = CircularAnnulus((x, y), d, d + 100).to_mask()
        bgap = ((apmask.to_image(ccd.shape) == 1) * ~ccd.mask
                * ~bg['source mask'])

        mms = sigma_clipped_stats(im[bgap])

        delta_bg = mms[1]
        bgstd = mms[2]

        db.execute('''
        UPDATE comet SET nbg=?,bg=?,bgsig=? WHERE rowid=?
        ''', (int(np.sum(bgap)), obj['bg'] + delta_bg, bgstd, obj['rowid']))
    else:
        logger.info('  Using source photometry table background.')
        model = models.Polynomial2D(1, c0_0=1, c1_0=1e-4, c0_1=1e-4)
        fitter = fitting.FittingWithOutlierRemoval(
            fitting.LinearLSQFitter(), sigma_clip, niter=3, sigma=3.)
        x = sources['xcentroid'][i]
        y = sources['ycentroid'][i]
        filtered, fit = fitter(model, x, y, sources['bg'][i])
        y, x = np.indices(ccd.shape)

        nbg = np.sum(sources['nbg'][i])
        delta_bg = fit(x, y)
        bgstd = filtered.std()

        f = '20{}/ppp/{}.dbg.fits'.format(
            obs['filename'][:6], obs['filename'])
        hdu = fits.HDUList()
        hdu.append(fits.PrimaryHDU(delta_bg))
        hdu[0].data += delta_bg
        hdu.writeto(f, overwrite=True)

        db.execute(
            'UPDATE comet SET nbg=?,bg=0,bgsig=? WHERE rowid=?', (nbg, bgstd, obj['rowid']))

    im = ccd.data - delta_bg

    ps = np.mean(proj_plane_pixel_scales(ccd.wcs)) * 3600
    r10k = 1e4 / 725 / obj['delta'] / 0.92
    ap10k = CircularAperture((obj['x'], obj['y']), r=r10k)
    area = ap10k.area()
    flux10k = aperture_photometry(im, ap10k)['aperture_sum'].data[0]
    calunc = zp_unc * flux10k / 1.0857
    err10k = np.sqrt(flux10k / obs['exptime']
                     + calunc**2
                     + area * (obs['bgstdev'] / obs['exptime'])**2)

    m10k = -2.5 * np.log10(flux10k) + zp
    merr10k = 1.0857 * err10k / flux10k

    fnu = 3631 * 10**(-0.4 * m10k) * u.Jy
    eph = {'rh': obj['rh'] * u.au, 'delta': obj['delta'] * u.au}
    afrho10k = Afrho.from_fluxd(0.617 * u.um, fnu, 1e4 * u.km, eph,
                                S=216555717567966.12 * u.Jy)
    afrho10k_unc = afrho10k * err10k / flux10k
    db.execute('''
    UPDATE comet SET calid=?,m=?,munc=?,afrho10k=?,afrho10k_unc=?
    WHERE rowid=?
    ''', (cal['rowid'], m10k, merr10k, afrho10k.cm, afrho10k_unc.cm,
          obj['rowid']))
