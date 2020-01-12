import numpy as np
cimport numpy as np


DTYPE = np.float
ctypedef np.float_t DTYPE_t

cpdef void add_fields(np.ndarray age_vals,
                     np.ndarray yields,
                     np.ndarray mass_p):


    mass_p = np.sum(np.transpose(age_vals) * yields, axis=1)

    return
