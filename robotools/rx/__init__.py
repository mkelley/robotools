from .core import *

def configure():
    import os
    import argparse
    from collections import OrderedDict
    import numpy as np

    parser = argparse.ArgumentParser(
        description='Reduce Lowell 31-in Robo telescope data.',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('source', type=str,
                        help='directory to process')
    parser.add_argument('target', type=str,
                        help=('directory for output (date is automatically'
                              ' appended)'))
    parser.add_argument('--reprocess-data', action='store_true',
                        help='reprocess data files')
    parser.add_argument('--reprocess-all', action='store_true',
                        help='recreate bias, flat, and all data files')
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

    date = args.source.split(os.sep)[-1]
    assert len(date) == 6 and date.isnumeric(), 'Directory format: path/to/YYMMDD'
    date = '20' + date
    
    config.source_dir = args.source
    config.target_dir = os.path.join(args.target, date, 'ppp')

    return config

def main():
    import os
    from ..logging import setup_logging
    from ..utils import file_lock

    config = configure()
    logger = setup_logging('Robo Reduction', 'robo-rx.log')

    os.system('mkdir -p {}'.format(config.target_dir))
    assert os.path.isdir(config.target_dir), 'Error creating target directory {}'.format(target_dir)

    lock = os.sep.join((config.target_dir, 'processing'))
    if os.path.exists(lock):
        logger.debug('{} has a status of "processing"'.format(
            config.target_dir))
        return

    with file_lock(lock):
        files = copy_files(config.source_dir, config.target_dir, config, logger)
        ic = load_files(config.target_dir, config, files, logger)
        ic = bias_correction(config.target_dir, ic, config, logger)
        ic = flatfield_correction(config.target_dir, ic, config, logger)
        array_corrections(ic, config, logger)
        solve_wcs(ic, config, logger)


if __name__ == '__main__':
    main()
