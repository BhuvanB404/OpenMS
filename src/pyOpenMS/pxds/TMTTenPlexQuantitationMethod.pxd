from Types cimport *
from IsobaricQuantitationMethod cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/TMTTenPlexQuantitationMethod.h>" namespace "OpenMS":

    cdef cppclass TMTTenPlexQuantitationMethod(IsobaricQuantitationMethod) :
        # wrap-inherits:
        #  IsobaricQuantitationMethod
        TMTTenPlexQuantitationMethod() except + nogil 
        TMTTenPlexQuantitationMethod(TMTTenPlexQuantitationMethod &) except + nogil 
