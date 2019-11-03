from .core import *

######################################################################
def configure():
    import os
    import re
    import argparse
    from collections import OrderedDict
    import numpy as np
    import astropy.units as u
    
    parser = argparse.ArgumentParser(description='Measure photometry for comets in Lowell 31-in Robo telescope data.')
    parser.add_argument('desg', nargs='+', help='Designation of the comet to process, and the basis for file outputs.  Match the beginning of the string in the FITS header object fields.')
    parser.add_argument('--reprocess', help='Reprocess specific date.')
    parser.add_argument('--reprocess-all', action='store_true', help='Reprocess all data.')
    parser.add_argument('--horizons', help='Optional name to use in HORIZONS ephemerides queries.  Otherwise <desg> will be used.')
    args = parser.parse_args()

    class config:
        database = 'robo.db'
        gmr_sun = 0.405  # Colina et al. 1996 + Tonry et al. 2012

    config.comet = ' '.join(args.desg).strip()

    date = args.reprocess
    if date is not None:
        assert len(date) in [6, 8, 12], "Invalid date {}.  Use: YYMMDD, YYYYMMDD, or YYYY-MM-DD".format(date)    
        if len(date) == 6:
            date = '20{}-{}-{}'.format(date[:2], date[2:4], date[4:])
        elif len(date) == 8:
            date = '{}-{}-{}'.format(date[:4], date[4:6], date[6:])
        elif len(date) == 10:
            assert date.count('-') == 2, "Invalid date {}.  Use: YYMMDD, YYYYMMDD, or YYYY-MM-DD".format(date)

    config.reprocess = date
    
    config.reprocess_all = args.reprocess_all
    if args.horizons is None:
        config.horizons = config.comet
    else:
        config.horizons = args.horizons

    config.abbrv = config.comet.replace(' ', '').replace('/', '').lower()
    try:
        config.abbrv = config.abbrv[:config.abbrv.index('(')]
    except ValueError:
        pass

    return config

######################################################################
def create_tables(db):
    db.execute('''
      CREATE TABLE IF NOT EXISTS
        comet(
          object TEXT(128),
          obsid INTEGER,
          calid INTEGER,

          x FLOAT,
          y FLOAT,
          ra FLOAT,
          dec FLOAT,

          ra_jpl FLOAT,
          dec_jpl FLOAT,
          ra_rate FLOAT,
          dec_rate FLOAT,
          rh FLOAT,
          delta FLOAT,
          phase FLOAT,
          sangle FLOAT,
          vangle FLOAT,

          nbg FLOAT,
          bg FLOAT,
          bgsig FLOAT,

          rap  BLOB,
          area BLOB,
          flux BLOB,
          func BLOB,

          m FLOAT,
          munc FLOAT,

          afrho10k FLOAT,
          afrho10k_unc FLOAT,

          nearest_photid INTEGER,
          nearest_d FLOAT,
          nearest_m FLOAT,
          ap_sources BLOB
        )''')

    db.execute('''
    CREATE UNIQUE INDEX IF NOT EXISTS objobs
    ON comet(object,obsid)
    ''')

    db.execute('''
    CREATE VIEW IF NOT EXISTS cometobs AS
    SELECT * FROM comet
    INNER JOIN obs ON comet.obsid=obs.rowid
    ''')
    
########################################################################
def main():
    import sys
    import warnings
    from ..logging import setup_logging
    from ..database import open_database

    config = configure()
    logger = setup_logging('Robo Comet', 'robo-comet.log')

    with open_database('.', config, logger) as db:
        create_tables(db)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            try:
                data = find_data(db, config, logger)
                db.commit()

                ephemeris(db, config, logger)
                db.commit()
            except RoboCometException as e:
                logger.error(str(e))
                sys.exit(0)

            for rowid, basename in data.items():
                logger.info(basename)

                try:
                    photometry(rowid, basename, db, config, logger)
                except RoboCometException as e:
                    logger.error(str(e))
                finally:
                    db.commit()
