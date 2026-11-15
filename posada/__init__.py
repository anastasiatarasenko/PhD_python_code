__import__('pkg_resources').declare_namespace(__name__)

# encoding: utf-8
__all__ = ['Common', 'Extractors', 'Handlers', 'Models', 'Processors', 'Readers', 'Tests']

# Constants, taken from BaseReader to resolve cirular imports.
READ_ONLY = 'r'
WRITE_NEW = 'w'
READ_WRITE = 'r+'

# Status
SAVED = 0
NOTSAVED = 1

# Module level default values
DEFAULT_TIME_UNITS = 'seconds since 1979-01-01 00:00:00'
CF_AUTHORITY = 'CF-1.7'
