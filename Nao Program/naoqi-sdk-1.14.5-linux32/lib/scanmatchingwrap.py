# This file was automatically generated by SWIG (http://www.swig.org).
# Version 1.3.40
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.
# This file is compatible with both classic and new-style classes.

"""
scanmatching docstring.
"""

from sys import version_info
if version_info >= (2,6,0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_scanmatchingwrap', [dirname(__file__)])
        except ImportError:
            import _scanmatchingwrap
            return _scanmatchingwrap
        if fp is not None:
            try:
                _mod = imp.load_module('_scanmatchingwrap', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _scanmatchingwrap = swig_import_helper()
    del swig_import_helper
else:
    import _scanmatchingwrap
del version_info
try:
    _swig_property = property
except NameError:
    pass # Python < 2.2 doesn't have 'property'.
def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "thisown"): return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static) or hasattr(self,name):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    if (name == "thisown"): return self.this.own()
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError(name)

def _swig_repr(self):
    try: strthis = "proxy of " + self.this.__repr__()
    except: strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0


class ScanMatchResultWrapper(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, ScanMatchResultWrapper, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, ScanMatchResultWrapper, name)
    __repr__ = _swig_repr
    def __init__(self): 
        this = _scanmatchingwrap.new_ScanMatchResultWrapper()
        try: self.this.append(this)
        except: self.this = this
    __swig_setmethods__["transform2D"] = _scanmatchingwrap.ScanMatchResultWrapper_transform2D_set
    __swig_getmethods__["transform2D"] = _scanmatchingwrap.ScanMatchResultWrapper_transform2D_get
    if _newclass:transform2D = _swig_property(_scanmatchingwrap.ScanMatchResultWrapper_transform2D_get, _scanmatchingwrap.ScanMatchResultWrapper_transform2D_set)
    __swig_setmethods__["maxNumInliers"] = _scanmatchingwrap.ScanMatchResultWrapper_maxNumInliers_set
    __swig_getmethods__["maxNumInliers"] = _scanmatchingwrap.ScanMatchResultWrapper_maxNumInliers_get
    if _newclass:maxNumInliers = _swig_property(_scanmatchingwrap.ScanMatchResultWrapper_maxNumInliers_get, _scanmatchingwrap.ScanMatchResultWrapper_maxNumInliers_set)
    __swig_setmethods__["numInliers"] = _scanmatchingwrap.ScanMatchResultWrapper_numInliers_set
    __swig_getmethods__["numInliers"] = _scanmatchingwrap.ScanMatchResultWrapper_numInliers_get
    if _newclass:numInliers = _swig_property(_scanmatchingwrap.ScanMatchResultWrapper_numInliers_get, _scanmatchingwrap.ScanMatchResultWrapper_numInliers_set)
    __swig_setmethods__["inlierError"] = _scanmatchingwrap.ScanMatchResultWrapper_inlierError_set
    __swig_getmethods__["inlierError"] = _scanmatchingwrap.ScanMatchResultWrapper_inlierError_get
    if _newclass:inlierError = _swig_property(_scanmatchingwrap.ScanMatchResultWrapper_inlierError_get, _scanmatchingwrap.ScanMatchResultWrapper_inlierError_set)
    __swig_setmethods__["valid"] = _scanmatchingwrap.ScanMatchResultWrapper_valid_set
    __swig_getmethods__["valid"] = _scanmatchingwrap.ScanMatchResultWrapper_valid_get
    if _newclass:valid = _swig_property(_scanmatchingwrap.ScanMatchResultWrapper_valid_get, _scanmatchingwrap.ScanMatchResultWrapper_valid_set)
    __swig_destroy__ = _scanmatchingwrap.delete_ScanMatchResultWrapper
    __del__ = lambda self : None;
ScanMatchResultWrapper_swigregister = _scanmatchingwrap.ScanMatchResultWrapper_swigregister
ScanMatchResultWrapper_swigregister(ScanMatchResultWrapper)

class RansacIcpMatcher2D(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, RansacIcpMatcher2D, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, RansacIcpMatcher2D, name)
    __repr__ = _swig_repr
    def __init__(self): 
        this = _scanmatchingwrap.new_RansacIcpMatcher2D()
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _scanmatchingwrap.delete_RansacIcpMatcher2D
    __del__ = lambda self : None;
    def getMaxNumSampledPoints(self): return _scanmatchingwrap.RansacIcpMatcher2D_getMaxNumSampledPoints(self)
    def setMaxNumSampledPoints(self, *args): return _scanmatchingwrap.RansacIcpMatcher2D_setMaxNumSampledPoints(self, *args)
    def getRequiredConfidence(self): return _scanmatchingwrap.RansacIcpMatcher2D_getRequiredConfidence(self)
    def setRequiredConfidence(self, *args): return _scanmatchingwrap.RansacIcpMatcher2D_setRequiredConfidence(self, *args)
    def getOutlierRatio(self): return _scanmatchingwrap.RansacIcpMatcher2D_getOutlierRatio(self)
    def setOutlierRatio(self, *args): return _scanmatchingwrap.RansacIcpMatcher2D_setOutlierRatio(self, *args)
    def getInlierThreshold(self): return _scanmatchingwrap.RansacIcpMatcher2D_getInlierThreshold(self)
    def setInlierThreshold(self, *args): return _scanmatchingwrap.RansacIcpMatcher2D_setInlierThreshold(self, *args)
    def getFinalIcpNumIterations(self): return _scanmatchingwrap.RansacIcpMatcher2D_getFinalIcpNumIterations(self)
    def setFinalIcpNumIterations(self, *args): return _scanmatchingwrap.RansacIcpMatcher2D_setFinalIcpNumIterations(self, *args)
    def getNumInliers(self): return _scanmatchingwrap.RansacIcpMatcher2D_getNumInliers(self)
    def getInlierIndices(self): return _scanmatchingwrap.RansacIcpMatcher2D_getInlierIndices(self)
    def getErrors(self): return _scanmatchingwrap.RansacIcpMatcher2D_getErrors(self)
    def getNumTries(self): return _scanmatchingwrap.RansacIcpMatcher2D_getNumTries(self)
    def setNumTries(self, *args): return _scanmatchingwrap.RansacIcpMatcher2D_setNumTries(self, *args)
    def resetNumTries(self): return _scanmatchingwrap.RansacIcpMatcher2D_resetNumTries(self)
    def isFixedNumTries(self): return _scanmatchingwrap.RansacIcpMatcher2D_isFixedNumTries(self)
    def scanMatch(self, *args): return _scanmatchingwrap.RansacIcpMatcher2D_scanMatch(self, *args)
    def particlesScanMatch(self, *args): return _scanmatchingwrap.RansacIcpMatcher2D_particlesScanMatch(self, *args)
RansacIcpMatcher2D_swigregister = _scanmatchingwrap.RansacIcpMatcher2D_swigregister
RansacIcpMatcher2D_swigregister(RansacIcpMatcher2D)

class SwigPyIterator(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, SwigPyIterator, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, SwigPyIterator, name)
    def __init__(self, *args, **kwargs): raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    __swig_destroy__ = _scanmatchingwrap.delete_SwigPyIterator
    __del__ = lambda self : None;
    def value(self): return _scanmatchingwrap.SwigPyIterator_value(self)
    def incr(self, n = 1): return _scanmatchingwrap.SwigPyIterator_incr(self, n)
    def decr(self, n = 1): return _scanmatchingwrap.SwigPyIterator_decr(self, n)
    def distance(self, *args): return _scanmatchingwrap.SwigPyIterator_distance(self, *args)
    def equal(self, *args): return _scanmatchingwrap.SwigPyIterator_equal(self, *args)
    def copy(self): return _scanmatchingwrap.SwigPyIterator_copy(self)
    def next(self): return _scanmatchingwrap.SwigPyIterator_next(self)
    def __next__(self): return _scanmatchingwrap.SwigPyIterator___next__(self)
    def previous(self): return _scanmatchingwrap.SwigPyIterator_previous(self)
    def advance(self, *args): return _scanmatchingwrap.SwigPyIterator_advance(self, *args)
    def __eq__(self, *args): return _scanmatchingwrap.SwigPyIterator___eq__(self, *args)
    def __ne__(self, *args): return _scanmatchingwrap.SwigPyIterator___ne__(self, *args)
    def __iadd__(self, *args): return _scanmatchingwrap.SwigPyIterator___iadd__(self, *args)
    def __isub__(self, *args): return _scanmatchingwrap.SwigPyIterator___isub__(self, *args)
    def __add__(self, *args): return _scanmatchingwrap.SwigPyIterator___add__(self, *args)
    def __sub__(self, *args): return _scanmatchingwrap.SwigPyIterator___sub__(self, *args)
    def __iter__(self): return self
SwigPyIterator_swigregister = _scanmatchingwrap.SwigPyIterator_swigregister
SwigPyIterator_swigregister(SwigPyIterator)

class vectorFloat(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, vectorFloat, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, vectorFloat, name)
    __repr__ = _swig_repr
    def iterator(self): return _scanmatchingwrap.vectorFloat_iterator(self)
    def __iter__(self): return self.iterator()
    def __nonzero__(self): return _scanmatchingwrap.vectorFloat___nonzero__(self)
    def __bool__(self): return _scanmatchingwrap.vectorFloat___bool__(self)
    def __len__(self): return _scanmatchingwrap.vectorFloat___len__(self)
    def pop(self): return _scanmatchingwrap.vectorFloat_pop(self)
    def __getslice__(self, *args): return _scanmatchingwrap.vectorFloat___getslice__(self, *args)
    def __setslice__(self, *args): return _scanmatchingwrap.vectorFloat___setslice__(self, *args)
    def __delslice__(self, *args): return _scanmatchingwrap.vectorFloat___delslice__(self, *args)
    def __delitem__(self, *args): return _scanmatchingwrap.vectorFloat___delitem__(self, *args)
    def __getitem__(self, *args): return _scanmatchingwrap.vectorFloat___getitem__(self, *args)
    def __setitem__(self, *args): return _scanmatchingwrap.vectorFloat___setitem__(self, *args)
    def append(self, *args): return _scanmatchingwrap.vectorFloat_append(self, *args)
    def empty(self): return _scanmatchingwrap.vectorFloat_empty(self)
    def size(self): return _scanmatchingwrap.vectorFloat_size(self)
    def clear(self): return _scanmatchingwrap.vectorFloat_clear(self)
    def swap(self, *args): return _scanmatchingwrap.vectorFloat_swap(self, *args)
    def get_allocator(self): return _scanmatchingwrap.vectorFloat_get_allocator(self)
    def begin(self): return _scanmatchingwrap.vectorFloat_begin(self)
    def end(self): return _scanmatchingwrap.vectorFloat_end(self)
    def rbegin(self): return _scanmatchingwrap.vectorFloat_rbegin(self)
    def rend(self): return _scanmatchingwrap.vectorFloat_rend(self)
    def pop_back(self): return _scanmatchingwrap.vectorFloat_pop_back(self)
    def erase(self, *args): return _scanmatchingwrap.vectorFloat_erase(self, *args)
    def __init__(self, *args): 
        this = _scanmatchingwrap.new_vectorFloat(*args)
        try: self.this.append(this)
        except: self.this = this
    def push_back(self, *args): return _scanmatchingwrap.vectorFloat_push_back(self, *args)
    def front(self): return _scanmatchingwrap.vectorFloat_front(self)
    def back(self): return _scanmatchingwrap.vectorFloat_back(self)
    def assign(self, *args): return _scanmatchingwrap.vectorFloat_assign(self, *args)
    def resize(self, *args): return _scanmatchingwrap.vectorFloat_resize(self, *args)
    def insert(self, *args): return _scanmatchingwrap.vectorFloat_insert(self, *args)
    def reserve(self, *args): return _scanmatchingwrap.vectorFloat_reserve(self, *args)
    def capacity(self): return _scanmatchingwrap.vectorFloat_capacity(self)
    __swig_destroy__ = _scanmatchingwrap.delete_vectorFloat
    __del__ = lambda self : None;
vectorFloat_swigregister = _scanmatchingwrap.vectorFloat_swigregister
vectorFloat_swigregister(vectorFloat)

class pairInt(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, pairInt, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, pairInt, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _scanmatchingwrap.new_pairInt(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_setmethods__["first"] = _scanmatchingwrap.pairInt_first_set
    __swig_getmethods__["first"] = _scanmatchingwrap.pairInt_first_get
    if _newclass:first = _swig_property(_scanmatchingwrap.pairInt_first_get, _scanmatchingwrap.pairInt_first_set)
    __swig_setmethods__["second"] = _scanmatchingwrap.pairInt_second_set
    __swig_getmethods__["second"] = _scanmatchingwrap.pairInt_second_get
    if _newclass:second = _swig_property(_scanmatchingwrap.pairInt_second_get, _scanmatchingwrap.pairInt_second_set)
    def __len__(self): return 2
    def __repr__(self): return str((self.first, self.second))
    def __getitem__(self, index): 
      if not (index % 2): 
        return self.first
      else:
        return self.second
    def __setitem__(self, index, val):
      if not (index % 2): 
        self.first = val
      else:
        self.second = val
    __swig_destroy__ = _scanmatchingwrap.delete_pairInt
    __del__ = lambda self : None;
pairInt_swigregister = _scanmatchingwrap.pairInt_swigregister
pairInt_swigregister(pairInt)

class pairFloat(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, pairFloat, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, pairFloat, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _scanmatchingwrap.new_pairFloat(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_setmethods__["first"] = _scanmatchingwrap.pairFloat_first_set
    __swig_getmethods__["first"] = _scanmatchingwrap.pairFloat_first_get
    if _newclass:first = _swig_property(_scanmatchingwrap.pairFloat_first_get, _scanmatchingwrap.pairFloat_first_set)
    __swig_setmethods__["second"] = _scanmatchingwrap.pairFloat_second_set
    __swig_getmethods__["second"] = _scanmatchingwrap.pairFloat_second_get
    if _newclass:second = _swig_property(_scanmatchingwrap.pairFloat_second_get, _scanmatchingwrap.pairFloat_second_set)
    def __len__(self): return 2
    def __repr__(self): return str((self.first, self.second))
    def __getitem__(self, index): 
      if not (index % 2): 
        return self.first
      else:
        return self.second
    def __setitem__(self, index, val):
      if not (index % 2): 
        self.first = val
      else:
        self.second = val
    __swig_destroy__ = _scanmatchingwrap.delete_pairFloat
    __del__ = lambda self : None;
pairFloat_swigregister = _scanmatchingwrap.pairFloat_swigregister
pairFloat_swigregister(pairFloat)

class vectorPairInt(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, vectorPairInt, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, vectorPairInt, name)
    __repr__ = _swig_repr
    def iterator(self): return _scanmatchingwrap.vectorPairInt_iterator(self)
    def __iter__(self): return self.iterator()
    def __nonzero__(self): return _scanmatchingwrap.vectorPairInt___nonzero__(self)
    def __bool__(self): return _scanmatchingwrap.vectorPairInt___bool__(self)
    def __len__(self): return _scanmatchingwrap.vectorPairInt___len__(self)
    def pop(self): return _scanmatchingwrap.vectorPairInt_pop(self)
    def __getslice__(self, *args): return _scanmatchingwrap.vectorPairInt___getslice__(self, *args)
    def __setslice__(self, *args): return _scanmatchingwrap.vectorPairInt___setslice__(self, *args)
    def __delslice__(self, *args): return _scanmatchingwrap.vectorPairInt___delslice__(self, *args)
    def __delitem__(self, *args): return _scanmatchingwrap.vectorPairInt___delitem__(self, *args)
    def __getitem__(self, *args): return _scanmatchingwrap.vectorPairInt___getitem__(self, *args)
    def __setitem__(self, *args): return _scanmatchingwrap.vectorPairInt___setitem__(self, *args)
    def append(self, *args): return _scanmatchingwrap.vectorPairInt_append(self, *args)
    def empty(self): return _scanmatchingwrap.vectorPairInt_empty(self)
    def size(self): return _scanmatchingwrap.vectorPairInt_size(self)
    def clear(self): return _scanmatchingwrap.vectorPairInt_clear(self)
    def swap(self, *args): return _scanmatchingwrap.vectorPairInt_swap(self, *args)
    def get_allocator(self): return _scanmatchingwrap.vectorPairInt_get_allocator(self)
    def begin(self): return _scanmatchingwrap.vectorPairInt_begin(self)
    def end(self): return _scanmatchingwrap.vectorPairInt_end(self)
    def rbegin(self): return _scanmatchingwrap.vectorPairInt_rbegin(self)
    def rend(self): return _scanmatchingwrap.vectorPairInt_rend(self)
    def pop_back(self): return _scanmatchingwrap.vectorPairInt_pop_back(self)
    def erase(self, *args): return _scanmatchingwrap.vectorPairInt_erase(self, *args)
    def __init__(self, *args): 
        this = _scanmatchingwrap.new_vectorPairInt(*args)
        try: self.this.append(this)
        except: self.this = this
    def push_back(self, *args): return _scanmatchingwrap.vectorPairInt_push_back(self, *args)
    def front(self): return _scanmatchingwrap.vectorPairInt_front(self)
    def back(self): return _scanmatchingwrap.vectorPairInt_back(self)
    def assign(self, *args): return _scanmatchingwrap.vectorPairInt_assign(self, *args)
    def resize(self, *args): return _scanmatchingwrap.vectorPairInt_resize(self, *args)
    def insert(self, *args): return _scanmatchingwrap.vectorPairInt_insert(self, *args)
    def reserve(self, *args): return _scanmatchingwrap.vectorPairInt_reserve(self, *args)
    def capacity(self): return _scanmatchingwrap.vectorPairInt_capacity(self)
    __swig_destroy__ = _scanmatchingwrap.delete_vectorPairInt
    __del__ = lambda self : None;
vectorPairInt_swigregister = _scanmatchingwrap.vectorPairInt_swigregister
vectorPairInt_swigregister(vectorPairInt)

class vectorPairFloat(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, vectorPairFloat, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, vectorPairFloat, name)
    __repr__ = _swig_repr
    def iterator(self): return _scanmatchingwrap.vectorPairFloat_iterator(self)
    def __iter__(self): return self.iterator()
    def __nonzero__(self): return _scanmatchingwrap.vectorPairFloat___nonzero__(self)
    def __bool__(self): return _scanmatchingwrap.vectorPairFloat___bool__(self)
    def __len__(self): return _scanmatchingwrap.vectorPairFloat___len__(self)
    def pop(self): return _scanmatchingwrap.vectorPairFloat_pop(self)
    def __getslice__(self, *args): return _scanmatchingwrap.vectorPairFloat___getslice__(self, *args)
    def __setslice__(self, *args): return _scanmatchingwrap.vectorPairFloat___setslice__(self, *args)
    def __delslice__(self, *args): return _scanmatchingwrap.vectorPairFloat___delslice__(self, *args)
    def __delitem__(self, *args): return _scanmatchingwrap.vectorPairFloat___delitem__(self, *args)
    def __getitem__(self, *args): return _scanmatchingwrap.vectorPairFloat___getitem__(self, *args)
    def __setitem__(self, *args): return _scanmatchingwrap.vectorPairFloat___setitem__(self, *args)
    def append(self, *args): return _scanmatchingwrap.vectorPairFloat_append(self, *args)
    def empty(self): return _scanmatchingwrap.vectorPairFloat_empty(self)
    def size(self): return _scanmatchingwrap.vectorPairFloat_size(self)
    def clear(self): return _scanmatchingwrap.vectorPairFloat_clear(self)
    def swap(self, *args): return _scanmatchingwrap.vectorPairFloat_swap(self, *args)
    def get_allocator(self): return _scanmatchingwrap.vectorPairFloat_get_allocator(self)
    def begin(self): return _scanmatchingwrap.vectorPairFloat_begin(self)
    def end(self): return _scanmatchingwrap.vectorPairFloat_end(self)
    def rbegin(self): return _scanmatchingwrap.vectorPairFloat_rbegin(self)
    def rend(self): return _scanmatchingwrap.vectorPairFloat_rend(self)
    def pop_back(self): return _scanmatchingwrap.vectorPairFloat_pop_back(self)
    def erase(self, *args): return _scanmatchingwrap.vectorPairFloat_erase(self, *args)
    def __init__(self, *args): 
        this = _scanmatchingwrap.new_vectorPairFloat(*args)
        try: self.this.append(this)
        except: self.this = this
    def push_back(self, *args): return _scanmatchingwrap.vectorPairFloat_push_back(self, *args)
    def front(self): return _scanmatchingwrap.vectorPairFloat_front(self)
    def back(self): return _scanmatchingwrap.vectorPairFloat_back(self)
    def assign(self, *args): return _scanmatchingwrap.vectorPairFloat_assign(self, *args)
    def resize(self, *args): return _scanmatchingwrap.vectorPairFloat_resize(self, *args)
    def insert(self, *args): return _scanmatchingwrap.vectorPairFloat_insert(self, *args)
    def reserve(self, *args): return _scanmatchingwrap.vectorPairFloat_reserve(self, *args)
    def capacity(self): return _scanmatchingwrap.vectorPairFloat_capacity(self)
    __swig_destroy__ = _scanmatchingwrap.delete_vectorPairFloat
    __del__ = lambda self : None;
vectorPairFloat_swigregister = _scanmatchingwrap.vectorPairFloat_swigregister
vectorPairFloat_swigregister(vectorPairFloat)



