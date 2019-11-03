def open_database(path, config, logger):
    import os
    import sqlite3
    import numpy as np

    sqlite3.register_adapter(np.int64, int)
    sqlite3.register_adapter(np.int32, int)
    sqlite3.register_adapter(np.float64, float)
    sqlite3.register_adapter(np.float32, float)
    
    fn = os.sep.join((path, config.database))
    db = sqlite3.connect(fn)
    db.row_factory = sqlite3.Row

    logger.info('Opened database: {}'.format(fn))

    return db

########################################################################
def fetchall(db, columns, table, where, parameters):
    import numpy as np
    
    c = db.execute('SELECT {} FROM {} WHERE {}'.format(
        ','.join(columns), table, where), parameters)
    rows = c.fetchall()

    if len(rows) > 0:
        d = dict(((k, np.array(v)) for k, v in zip(columns, zip(*rows))))
    else:
        d = {}
        
    return d

