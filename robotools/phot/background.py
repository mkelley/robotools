import numpy as np
import scipy.ndimage as nd
from astropy.io import fits
from astropy.modeling import models, fitting
from astropy.stats import sigma_clip, sigma_clipped_stats
import sep

########################################################################


def background(ccd, bgf, db, config, logger):
    # first pass
    bkg, source_map = estimate_bg(ccd, db, box=256, thresh=3)

    # second pass
    mask = source_map > 0
    bkg, source_map = estimate_bg(ccd, db, box=64, thresh=1.5, mask=mask)

    db.execute('UPDATE obs SET bg=?,bgrms=? WHERE filename=?',
               (bkg.globalback, bkg.globalrms, ccd.meta['ccdfname']))

    hdu = fits.HDUList()
    hdu.append(fits.PrimaryHDU(bkg.back()))
    hdu[0].name = 'background'
    hdu[0].header['BG'] = bkg.globalback, 'Global background'
    hdu[0].header['BGRMS'] = bkg.globalrms, 'Global background RMS'
    hdu.append(fits.ImageHDU(source_map, name='source map'))

    hdu.writeto(bgf, overwrite=True)
    logger.debug('  2D background and source map saved to {}.'.format(bgf))

########################################################################


def estimate_bg(ccd, db, box=128, thresh=3, mask=None):
    im = ccd.data.astype('<f4')
    if mask is None:
        mask = ccd.mask
    else:
        mask = ccd.mask + mask

    bkg = sep.Background(im, mask=mask, bw=box, bh=box, fw=3, fh=3)
    objects, source_map = sep.extract(im - bkg, thresh, err=bkg.globalrms,
                                      segmentation_map=True)
    return bkg, source_map
