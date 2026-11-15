# coding: utf-8

"""
.. module::posada.Common.logger

POSADA logging

:copyright: Copyright 2011-2016 SOLab.

.. sectionauthor:: Denis Spiridonov <sdi@rshu.ru>, Alexander Myasoedov <mag@rshu.ru>
.. codeauthor:: Denis Spiridonov <sdi@rshu.ru>, Alexander Myasoedov <mag@rshu.ru>
.. author:: SOLab Dev Team

"""

import logging
import logging.handlers
import os

__author__ = 'Alexander Myasoedov, Denis Spiridonov'

loggers = {}


def get_logger(log_filepath=None):

    if not log_filepath:
        log_dirname = '/home3/homedir7/perso/atarasen/codes/python/posada/posadalog/'
        if not os.path.isdir(log_dirname):
	    print(log_dirname)
            os.makedirs(log_dirname)

        log_filename = 'posada.out'

        log_filepath = os.path.join(log_dirname, log_filename)

    if loggers.get('posada'):
        return loggers.get('posada')
    else:
        logger = logging.getLogger('posada')
        logger.setLevel(logging.DEBUG)

        # Add the log message handler to the logger
        handler = logging.handlers.RotatingFileHandler(log_filepath,
                                                       maxBytes=10000000,
                                                       backupCount=5)
        formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s")

        handler.setFormatter(formatter)
        logger.addHandler(handler)
        loggers.update(dict(posada=logger))
        return logger

# log_dirname = '/var/log/posada'
# if not os.path.isdir(log_dirname):
#     os.makedirs(log_dirname)

# log_filename = 'posada.out'

# log_filepath = os.path.join(log_dirname, log_filename)

# logger = logging.getLogger('posada')
# logger.setLevel(logging.DEBUG)

# # Add the log message handler to the logger
# handler = logging.handlers.RotatingFileHandler(log_filepath,
#                                                maxBytes=10000000,
#                                                backupCount=5)
# formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s")

# handler.setFormatter(formatter)
# logger.addHandler(handler)
