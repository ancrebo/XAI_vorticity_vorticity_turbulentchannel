from __future__ import absolute_import


# ----------------------------------------------------------------------------------------------------------------------
# Modification for reading locally slicer
# from cloudpickle.cloudpickle import *  # noqa
# from cloudpickle.cloudpickle_fast import CloudPickler, dumps, dump  # noqa
# ----------------------------------------------------------------------------------------------------------------------
from py_bin.py_packages.cloudpickle.cloudpickle import *  # noqa
from py_bin.py_packages.cloudpickle.cloudpickle_fast import CloudPickler, dumps, dump  # noqa


# Conform to the convention used by python serialization libraries, which
# expose their Pickler subclass at top-level under the  "Pickler" name.
Pickler = CloudPickler

__version__ = '2.0.0'
