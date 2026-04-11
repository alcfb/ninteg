#!/usr/bin/env python3
import copy
import numpy
import ctypes
import pathlib

c_double   = ctypes.c_double
c_int      = ctypes.c_int
c_void_p   = ctypes.c_void_p
c_int_p    = ctypes.POINTER(ctypes.c_int)
c_double_p = ctypes.POINTER(ctypes.c_double)
c_function = ctypes.CFUNCTYPE(None, c_int_p, c_double_p, c_double_p, c_void_p, c_void_p, c_double_p)

def as_array(ptr, n):
    return numpy.ctypeslib.as_array(ptr, shape=(n,))

lib_path = next(pathlib.Path(__file__).parent.glob("libninteg.*"))
lib = ctypes.CDLL(str(lib_path))

class Info:

    labels = [
        'total_iterations', 'successful_steps', 'rejected_iterations',
        'rejected_steps', 'step_attempts', 'function_calls',
        'previous_order', 'suggested_order', 'time',
        'last_step_size', 'suggested_step_size', 'relative_error',
        'dummy13', 'dummy14', 'dummy15'
    ]

    def __init__(self):
        self.data = numpy.zeros(len(self.labels), dtype=numpy.float64)
        self._index = {k: i for i, k in enumerate(self.labels)}

    def __getattr__(self, name):
        try:
            return self.data[self._index[name]]
        except KeyError:
            raise AttributeError(name)

class IVP:

    def __init__(self, x):

        self.x = numpy.ascontiguousarray(x, dtype=numpy.float64)

        self.lib = lib
        self.info = Info()

        self.lib.make.argtypes = [c_int_p]
        self.lib.make (ctypes.byref (c_int (self.x.size)))

        self._callback = c_function (self._solv)

        self.lib.init_solv (ctypes.byref (self._callback))

        self.lib.init_t0.argtypes = [c_double_p]
        self.lib.init_t1.argtypes = [c_double_p]
        self.lib.init_h0.argtypes = [c_double_p]
        self.lib.init_rtol.argtypes = [c_double_p]
        self.lib.init_atol.argtypes = [c_double_p]
        self.lib.init_q_max.argtypes = [c_int_p]
        self.lib.init_h_min.argtypes = [c_double_p]
        self.lib.init_h_max.argtypes = [c_double_p]

        self.lib.ivp_init.argtypes = []
        self.lib.ivp_next.argtypes = [c_int_p]

        self.lib.get_t.argtypes = [c_double_p]
        self.lib.get_x.argtypes = [c_int_p, c_void_p]
        self.lib.get_info.argtypes = [c_void_p]

        self.lib.init_x0.argtypes = [c_int_p, c_void_p]

        self.t = c_double()
        self.t_ref = ctypes.byref (self.t)

        self.istatus = c_int()
        self.istatus_ref = ctypes.byref (self.istatus)

    def _solv (self, pn, ph, pt, pb, px, pe):

        n = pn[0]
        t = pt[0]
        h = ph[0]
        e = pe[0]

        b = as_array (ctypes.cast(pb, c_double_p), n)
        x = as_array (ctypes.cast(px, c_double_p), n)

        self.solv (h, t, b, x, e)


    def solve (self, time, solv, rtol, atol, h0, hmin, hmax, qmax):

        self.solv = solv

        self.lib.ivp_init ()

        self.lib.init_x0 (ctypes.byref (c_int (self.x.size)), self.x.ctypes)
        self.lib.init_rtol (ctypes.byref (c_double (rtol)))
        self.lib.init_atol (ctypes.byref (c_double (atol)))
        self.lib.init_h0 (ctypes.byref (c_double (h0)))
        self.lib.init_t0 (ctypes.byref (c_double (time[0])))
        self.lib.init_t1 (ctypes.byref (c_double (time[1])))
        self.lib.init_q_max (ctypes.byref (c_int (qmax)))
        self.lib.init_h_min (ctypes.byref (c_double (hmin)))
        self.lib.init_h_max (ctypes.byref (c_double (hmax)))

        while True:

            self.lib.ivp_next (self.istatus_ref)
            self.lib.get_info (self.info.data.ctypes)
            self.lib.get_x (c_int (self.x.size), self.x.ctypes)
            self.lib.get_t (self.t_ref)

            yield self.t.value, self.x, self.info

            if self.istatus.value == 3: break

    def close(self):
        self.lib.clean()




def integrate (time, x0, solv, rtol=1.E-6, atol=1.E-10, h0=1.E-6, hmin=1.E-10, hmax=1.E+10, qmax=5):

    ivp = IVP (x0)

    solution = ivp.solve (time, solv, rtol, atol, h0, hmin, hmax, qmax)

    yield from solution

    ivp.close()
