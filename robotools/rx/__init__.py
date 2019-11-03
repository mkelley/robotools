from .core import *

def configure():
    import os
    import argparse
    from collections import OrderedDict
    import numpy as np

    parser = argparse.ArgumentParser(
        description='Reduce Lowell 31-in Robo telescope data.',
        epilog='''
robo-rx expects the following directory tree for raw data:

.
├── raw
│   ├── 170827
│   │   ├── 170827.123
│   │   ├── 170827.124
│   │   ├── 170827.125
│   │   └── ...
│   └── 170830
│       ├── 170830.310
│       ├── 170830.311
│       ├── 170830.312
│       └── ...
└── robo-rx.log

where `raw` is in the current directory.
''', formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('date', default=None, nargs='?', help='Date to process: YYYYMMDD, YYYY-MM-DD, or YYMMDD.  The default is to process all.')
    parser.add_argument('--reprocess-data', action='store_true', help='Reprocess data files.')
    parser.add_argument('--reprocess-all', action='store_true', help='Recreate bias, flat, and all data files.')
    args = parser.parse_args()
    
    class config:
        file_template = '.*/[12][0-9][01][0-9][0-3][0-9][a-z]?.[0-9][0-9][0-9]$'

        # in order of precedence
        flats = OrderedDict()
        flats['sky'] = 'Twilight flat'

        overscan_correct = False
        
        # in order of operations
        trim_sections = (np.s_[:, :2105], np.s_[:, 53:])

        saturation = 50000 * 2.1  # electrons
        
        lacosmic = True

        rebin = 2
        xflip = False
        rotate = 90  # degrees CCW, must be multiple of 90

    config.reprocess_all = args.reprocess_all
    config.reprocess_data = True if args.reprocess_all else args.reprocess_data

    if args.date is not None:
        date = args.date.replace('-', '')
        assert (len(date) in [6, 8]) and date.isnumeric(), 'Date format: YYYYMMDD, YYYY-MM-DD, or YYMMDD'
        if len(date) == 6:
            date = '20' + date
    
        config.source_dir = os.sep.join(('raw', date[2:]))
    else:
        config.source_dir = 'raw'

    assert os.path.isdir(config.source_dir), '{} does not exist'.format(config.source_dir)
    
    return config

def main():
    import os
    from glob import glob
    from ..logging import setup_logging
    from ..utils import file_lock

    config = configure()
    logger = setup_logging('Robo Reduction', 'robo-rx.log')
    
    if config.source_dir == 'raw':
        source_dirs = sorted(glob('raw/[12][0-9][01][0-9][0-3][0-9]'))
    else:
        source_dirs = [config.source_dir]

    for source_dir in source_dirs:
        if os.path.exists(os.sep.join((source_dir, 'downloading'))):
            logger.debug('{} has a status of "downloading."  Skip.'.format(
                source_dir))
            continue

        target_dir = os.sep.join(('20' + source_dir[-6:], 'ppp'))
        os.system('mkdir -p {}'.format(target_dir))
        assert os.path.isdir(target_dir), 'Error creating target directory {}'.format(target_dir)

        lock = os.sep.join((target_dir, 'processing'))
        if os.path.exists(lock):
            logger.debug('{} has a status of "processing."  Skip.'.format(
                target_dir))
            continue

        with file_lock(lock):
            files = copy_files(source_dir, target_dir, config, logger)
            ic = load_files(target_dir, config, files, logger)
            ic = bias_correction(target_dir, ic, config, logger)
            ic = flatfield_correction(target_dir, ic, config, logger)
            array_corrections(ic, config, logger)
            solve_wcs(ic, config, logger)


