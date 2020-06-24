from .core import *

def configure():
    import os
    import argparse
    from collections import OrderedDict
    import numpy as np
    import astropy.units as u

    def Array(a):
        return np.fromstring(a, sep=',')

    parser = argparse.ArgumentParser(
        description='Measure photometry in Lowell 31-in Robo telescope data.')
    parser.add_argument(
        'path',
        help='directory to process, data are in a subdirectory labeled ppp')
    parser.add_argument(
        '-f', action='store_true',
        help='force processing of directory')
    parser.add_argument(
        '--update-obs', action='store_true',
        help=('update observations database with meta data '
              'from FITS headers'))
    parser.add_argument(
        '--plot', action='store_true',
        help='plot catalog calibration')
    parser.add_argument(
        '--reprocess-all', action='store_true',
        help='rerun all steps')
    parser.add_argument(
        '--reprocess-bg', action='store_true',
        help=('rerun background, source detection, photometry, and '
              'calibration'))
    parser.add_argument(
        '--reprocess-phot', action='store_true',
        help='rerun photometry, and frame calibration')
    parser.add_argument(
        '--reprocess-cal', action='store_true',
        help='rerun frame calibration')
    parser.add_argument(
        '--match-limit', type=float, default=1.4,
        help='maximum distance to PS1 catalog for a match. [arcsec]')
    parser.add_argument(
        '--color-cor-r', type=float,
        help='g-r color correction')
    parser.add_argument(
        '--mlim', type=Array, default=np.r_[12, 17],
        help='magnitude limits for catalog calibration')
    parser.add_argument(
        '--gmi-lim', type=Array, default=np.r_[0.2, 3.0],
        help='g-i color limits for catalog calibration')
    args = parser.parse_args()

    assert os.path.isdir(
        args.path), '"{}" must be an existing directory.'.format(args.path)

    class config:
        file_template = '[12][0-9][01][0-9][0-3][0-9]*.[0-9][0-9][0-9].fits'
        border_mask = 5  # number of pixels to mask along border
        database = 'robo.db'

    config.force_processing = args.f
    config.update_obs = args.update_obs
    config.plot = args.plot

    config.reprocess_bg = any((args.reprocess_all, args.reprocess_bg))
    config.reprocess_phot = any((args.reprocess_all,
                                 args.reprocess_bg,
                                 args.reprocess_phot))
    config.reprocess_cal = any((args.reprocess_all,
                                args.reprocess_bg,
                                args.reprocess_phot,
                                args.reprocess_cal))

    config.source_dir = args.path.strip('/')
    config.match_limit = args.match_limit * u.arcsec
    config.color_cor_r = args.color_cor_r
    config.mlim = args.mlim
    config.gmi_lim = args.gmi_lim

    return config


def setup_tables(db):
    db.execute('''
CREATE TABLE IF NOT EXISTS
  obs(
    obsid INTEGER PRIMARY KEY,
    filename TEXT(32) UNIQUE,

    object TEXT(128),
    obsdate TEXT(32),
    airmass FLOAT,
    filter TEXT(8),
    exptime FLOAT,
    ra TEXT(16),
    dec TEXT(16),
    oairtemp FLOAT,
    relhum FLOAT,
    seeing FLOAT,

    bg FLOAT,
    bgrms FLOAT,

    zp FLOAT,
    colorcor FLOAT,
    ncal INTEGER,
    calunc FLOAT
)''')

#    db.execute('''
#    CREATE TABLE IF NOT EXISTS
#    cal(
#      obsid INTGER PRIMARY KEY,
#      obsdate TEXT(8) UNIQUE,
#
#      nframes INTEGER,
#      ngood INTEGER,
#
#      zp0 FLOAT,
#      zp0_std FLOAT,
#      zp0_unc FLOAT,
#      extcor FLOAT,
#      colorcor FLOAT,
#
#      FOREIGN KEY (obsid) REFERENCES obs(obsid)
#    )''')

########################################################################


def main():
    import re
    import os
    import sys
    from glob import glob
    import warnings
    from ..logging import setup_logging
    from ..utils import file_lock
    from ..database import open_database

    config = configure()
    logger = setup_logging('Robo Photometry', 'robo-phot.log')

    warnings.simplefilter("ignore")

    data_dir = os.sep.join((config.source_dir, 'ppp'))
    db_dir = './'

    with open_database(db_dir, config, logger) as db:
        setup_tables(db)

        if not os.path.isdir(data_dir):
            logger.error('{} is not a directory'.format(data_dir))
            return

        lock = os.sep.join((data_dir, 'processing'))
        if os.path.exists(lock) and not config.force_processing:
            logger.debug('{} has a status of "processing."  Skip.'.format(
                data_dir))
            return

        with file_lock(lock):
            try:
                ic = find_files(data_dir, db, config, logger)
            except RoboPhotException as e:
                logger.error(str(e))
                return

            db.commit()

            for f in ic.files:
                try:
                    basename = f.replace('.fits', '')
                    photometry(data_dir, basename, db, config, logger)

                    if config.plot:
                        plot(data_dir, basename, logger)

                except RoboPhotException as e:
                    logger.error(str(e))
                finally:
                    db.commit()
