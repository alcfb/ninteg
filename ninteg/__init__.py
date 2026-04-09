#!/usr/bin/env python3
import copy
import numpy
import ctypes
import pathlib

c_char_p   = ctypes.c_char_p
c_double   = ctypes.c_double
c_int      = ctypes.c_int
c_short    = ctypes.c_short
c_long     = ctypes.c_long
c_void_p   = ctypes.c_void_p
p_c_int    = ctypes.POINTER(ctypes.c_int)
p_c_double = ctypes.POINTER(ctypes.c_double)
p_c_char   = ctypes.POINTER(ctypes.c_char)

def ptr2value (ptr):
    return numpy.ctypeslib.as_array(ptr, shape=(1,))[0]

def ptr2array (n, pointer_index, dtype=c_double):
    pointer_to_array = ctypes.cast (pointer_index, ctypes.POINTER(dtype))
    values = numpy.ctypeslib.as_array (pointer_to_array, shape=(n,))
    return values

lib_path = next(pathlib.Path(__file__).parent.glob("libninteg.*"))
lib = ctypes.CDLL(str(lib_path))

class Info:
    def __init__(self):
        self.output = numpy.zeros (15, dtype='float64')
        self.values = []
        self.labels = ['total_iterations', 'successful_steps', 'rejected_iterations', 'rejected_steps', 'step_attempts', 'function_calls', 'previous_order','suggested_order', 'time', 'last_step_size', 'suggested_step_size', 'relative_error']
    def update (self):
        self.values.append (copy.copy(self.output))
    def __call__ (self):
        for label, value in zip(self.labels, self.values[-1]):
            print (f"{label:20s}{value:12.6f}")

class IVP:
    def __init__(self):

        self.info = Info()
        self.lib = lib

    def setup (self, t0, x0, dynamics):

        self.t0 = t0
        self.x = numpy.array (x0, dtype='float64')
        self.dynamics = dynamics

        self.lib.make.argtypes = [p_c_int]
        self.lib.make (ctypes.byref (c_int (self.x.size)))

        self.p_solv = ctypes.CFUNCTYPE (None, p_c_int, p_c_double, p_c_double, c_void_p, c_void_p, p_c_double)(self._solv)

        self.lib.init_solv (ctypes.byref (self.p_solv))

        self.lib.init_t0.argtypes = [p_c_double]
        self.lib.init_t1.argtypes = [p_c_double]
        self.lib.init_h0.argtypes = [p_c_double]
        self.lib.init_rtol.argtypes = [p_c_double]
        self.lib.init_atol.argtypes = [p_c_double]
        self.lib.init_im.argtypes = [p_c_int]
        self.lib.init_q_max.argtypes = [p_c_int]
        self.lib.init_h_min.argtypes = [p_c_double]
        self.lib.init_h_max.argtypes = [p_c_double]

        self.lib.ivp_init.argtypes = []
        self.lib.ivp_next.argtypes = [p_c_int]

        self.lib.get_t.argtypes = [p_c_double]
        self.lib.get_x.argtypes = [p_c_int, c_void_p]
        self.lib.get_info.argtypes = [c_void_p]


    def solve (self, t=1, tol=1.E-6, h=1.E-6, method='anderson', hmin=1.E-10, hmax=1.E+10, qmax=5):

        imethod = {'picard':12, 'anderson':13}[method]

        self.lib.ivp_init ()

        self.lib.init_rtol (ctypes.byref (c_double (tol)))
        self.lib.init_h0 (ctypes.byref (c_double (h)))
        self.lib.init_t0 (ctypes.byref (c_double (self.t0)))
        self.lib.init_t1 (ctypes.byref (c_double (t)))
        self.lib.init_im (ctypes.byref (c_int (imethod)))
        self.lib.init_q_max (ctypes.byref (c_int (qmax)))
        self.lib.init_h_min (ctypes.byref (c_double (hmin)))
        self.lib.init_h_max (ctypes.byref (c_double (hmax)))

        self.lib.init_x0.argtypes = [p_c_int, c_void_p]
        self.lib.init_x0 (ctypes.byref (c_int (self.x.size)), self.x.ctypes)

        time = c_double()

        istatus = c_int ()

        while True:

            self.lib.ivp_next (ctypes.byref (istatus))
            self.lib.get_info (self.info.output.ctypes)
            self.lib.get_x (c_int (self.x.size), self.x.ctypes)
            self.lib.get_t (ctypes.byref (time))

            self.info.update ()

            yield time.value, self.x

            if istatus.value == 3: break

    def __del__(self):
        if hasattr(self, 'lib'): self.lib.clean()

    def close(self):
        self.lib.clean()

    def _solv (self, pn, ph, pt, pb, px, pe):
        n = ptr2value (pn)
        t = ptr2value (pt)
        h = ptr2value (ph)
        e = ptr2value (pe)
        b = ptr2array (n, pb)
        x = ptr2array (n, px)
        self.dynamics (h, t, b, x, e)
