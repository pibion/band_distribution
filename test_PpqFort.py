import numpy as np
import ctypes
import sys
import os

# Path to the directory containing PpqDist_v6.py
module_dir = './python'
sys.path.append(module_dir)

# File PpqDist_v6 contains our existing (tested!) python functions
# For P(Ep, Eq)
import pq_dist_v6 as ppq

folderpath = '' 
if os.name == 'posix': #Linux/Mac
    DLLname = 'lib/libband_distribution.so'
else:
    print ('OS other than Linux/Mac are not supported')

api = np.ctypeslib.load_library(DLLname,folderpath)

# Define all our parameters and variables
a = 0.16
b = 0.18
F0 = 0.122
s = 0.0
eps = 3E-3
V=3.0
p0=0.06421907
p10=0.48998486
q0=0.06421907
q10=0.48998486

Er = 234.0
Ep = 337.5
Eq = 100

# Define argument and return types for the function Y
api.Y.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_double]
api.Y.restype = ctypes.c_double

result_fort = api.Y(ctypes.c_double(Er), ctypes.c_double(a), ctypes.c_double(b))  
result_py = a * Er ** b
print('back in python after running function Y')
print('The fortran result is ', result_fort)
print('The python result is ', result_py)

# Define argument and return types for the function F
api.F.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_double]
api.F.restype = ctypes.c_double

result_fort = api.F(ctypes.c_double(Er), ctypes.c_double(F0), ctypes.c_double(s))
result_py = F0 + s * Er
print('back in python after running function F')
print('The fortran result is ', result_fort)
print('The python result is ', result_py)

# Define argument and return types for the function Nbar
api.Nbar.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double]
api.Nbar.restype = ctypes.c_double

result_fort = api.Nbar(ctypes.c_double(Er), ctypes.c_double(a), ctypes.c_double(b), ctypes.c_double(eps))
result_py = ppq.Nbar(Er, a=a, b=b, eps=eps)
print('back in python after running function Nbar')
print('The fortran result is ', result_fort)
print('The python result is ', result_py)

# Define argument and return types for the function sigp
api.sigp.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double]
api.sigp.restype = ctypes.c_double

result_fort = api.sigp(ctypes.c_double(Ep), ctypes.c_double(eps), ctypes.c_double(V), 
                       ctypes.c_double(p0), ctypes.c_double(p10))
result_py = ppq.sigp(Ep, eps=eps, V=V, p0=p0, p10=p10)
print('back in python after running function sigp')
print('The fortran result is ', result_fort)
print('The python result is ', result_py)

# Define argument and return types for the function sigq
api.sigq.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_double]
api.sigq.restype = ctypes.c_double

result_fort = api.sigq(ctypes.c_double(Eq), ctypes.c_double(q0), ctypes.c_double(q10))
result_py = ppq.sigq(Eq, q0=q0, q10=q10)
print('back in python after running function sigq')
print('The fortran result is ', result_fort)
print('The python result is ', result_py)

# Define argument and return types for the function PErN
api.PErN.argtypes = [ctypes.c_double]
api.PErN.restype = ctypes.c_double

result_fort = api.PErN(ctypes.c_double(Er))
result_py = ppq.PErN(Er)
print('back in python after running function PErN')
print('The fortran result is ', result_fort)
print('The python result is ', result_py)

# Define argument and return types for the function PErG
api.PErG.argtypes = [ctypes.c_double]
api.PErG.restype = ctypes.c_double

result_fort = api.PErG(ctypes.c_double(Er))
result_py = ppq.PErG(Er)
print('back in python after running function PErG')
print('The fortran result is ', result_fort)
print('The python result is ', result_py)

# Define argument and return types for the function aN
api.aN.argtypes = [ctypes.c_double] * 11
api.aN.restype = ctypes.c_double

result_fort = api.aN(ctypes.c_double(Er), ctypes.c_double(Ep), ctypes.c_double(Eq), 
                     ctypes.c_double(F0), ctypes.c_double(s), ctypes.c_double(eps),
                     ctypes.c_double(V), ctypes.c_double(p0), ctypes.c_double(p10),
                     ctypes.c_double(q0), ctypes.c_double(q10))
result_py = ppq.aN(Er, Ep, Eq, F0=F0, s=s, eps=eps, V=V, p0=p0, p10=p10, q0=q0, q10=q10)
print('back in python after running function aN')
print('The fortran result is ', result_fort)
print('The python result is ', result_py)

# Define argument and return types for the function bN
api.bN.argtypes = [ctypes.c_double] * 11
api.bN.restype = ctypes.c_double

result_fort = api.bN(ctypes.c_double(Er), ctypes.c_double(Ep), ctypes.c_double(Eq), 
                     ctypes.c_double(a), ctypes.c_double(b),
                     ctypes.c_double(F0), ctypes.c_double(s), ctypes.c_double(eps),
                     ctypes.c_double(V), ctypes.c_double(p0), ctypes.c_double(p10),
                     ctypes.c_double(q0), ctypes.c_double(q10))
result_py = ppq.bN(Er, Ep, Eq, a=a, b=b, F0=F0, s=s, eps=eps, V=V, p0=p0, p10=p10, q0=q0, q10=q10)
print('back in python after running function bN')
print('The fortran result is ', result_fort)
print('The python result is ', result_py)

# Define argument and return types for the function cN1
api.cN1.argtypes = [ctypes.c_double] * 11
api.cN1.restype = ctypes.c_double

result_fort = api.cN1(ctypes.c_double(Er), ctypes.c_double(Ep), ctypes.c_double(Eq), 
                      ctypes.c_double(a), ctypes.c_double(b),
                      ctypes.c_double(F0), ctypes.c_double(s), ctypes.c_double(eps),
                      ctypes.c_double(V), ctypes.c_double(p0), ctypes.c_double(p10),
                      ctypes.c_double(q0), ctypes.c_double(q10))
result_py = ppq.cN(Er, Ep, Eq, a=a, b=b, F0=F0, s=s, eps=eps, V=V, p0=p0, p10=p10, q0=q0, q10=q10)
print('back in python after running function cN1')
print('The fortran result is ', result_fort)
print('The python result is ', result_py)

# Define argument and return types for the function cN2
api.cN2.argtypes = [ctypes.c_double] * 13
api.cN2.restype = ctypes.c_double

result_fort = api.cN2(ctypes.c_double(Er), ctypes.c_double(Ep), ctypes.c_double(Eq), 
                      ctypes.c_double(a), ctypes.c_double(b),
                      ctypes.c_double(F0), ctypes.c_double(s), ctypes.c_double(eps),
                      ctypes.c_double(V), ctypes.c_double(p0), ctypes.c_double(p10),
                      ctypes.c_double(q0), ctypes.c_double(q10))
result_py = ppq.CN(Er, Ep, Eq, a=a, b=b, F0=F0, s=s, eps=eps, V=V, p0=p0, p10=p10, q0=q0, q10=q10)
print('back in python after running function cN2')
print('The fortran result is ', result_fort)
print('The python result is ', result_py)

# Define argument and return types for the function PpqExp
api.PpqExp.argtypes = [ctypes.c_double] * 13
api.PpqExp.restype = ctypes.c_double

result_fort = api.PpqExp(ctypes.c_double(Er), ctypes.c_double(Ep), ctypes.c_double(Eq), 
                         ctypes.c_double(a), ctypes.c_double(b),
                         ctypes.c_double(F0), ctypes.c_double(s), ctypes.c_double(eps),
                         ctypes.c_double(V), ctypes.c_double(p0), ctypes.c_double(p10),
                         ctypes.c_double(q0), ctypes.c_double(q10))
result_py = ppq.PpqExp(Er, Ep, Eq, a=a, b=b, F0=F0, s=s, eps=eps, V=V, p0=p0, p10=p10, q0=q0, q10=q10)
print('back in python after running function PpqExp')
print('The fortran result is ', result_fort)
print('The python result is ', result_py)

# Define argument and return types for the function PpqPreExp
api.PpqPreExp.argtypes = [ctypes.c_double] * 13
api.PpqPreExp.restype = ctypes.c_double

result_fort = api.PpqPreExp(ctypes.c_double(Er), ctypes.c_double(Ep), ctypes.c_double(Eq), 
                            ctypes.c_double(a), ctypes.c_double(b),
                            ctypes.c_double(F0), ctypes.c_double(s), ctypes.c_double(eps),
                            ctypes.c_double(V), ctypes.c_double(p0), ctypes.c_double(p10),
                            ctypes.c_double(q0), ctypes.c_double(q10))
result_py = ppq.PpqPreExp(Er, Ep, Eq, a=a, b=b, F0=F0, s=s, eps=eps, V=V, p0=p0, p10=p10, q0=q0, q10=q10)
print('back in python after running function PpqPreExp')
print('The fortran result is ', result_fort)
print('The python result is ', result_py)

# Define argument and return types for the function PpqFullN
api.PpqFullN.argtypes = [ctypes.c_double] * 13
api.PpqFullN.restype = ctypes.c_double

result_fort = api.PpqFullN(ctypes.c_double(Er), ctypes.c_double(Ep), ctypes.c_double(Eq), 
                           ctypes.c_double(a), ctypes.c_double(b),
                           ctypes.c_double(F0), ctypes.c_double(s), ctypes.c_double(eps),
                           ctypes.c_double(V), ctypes.c_double(p0), ctypes.c_double(p10),
                           ctypes.c_double(q0), ctypes.c_double(q10))
result_py = ppq.PpqFullN(Er, Ep, Eq, a=a, b=b, F0=F0, s=s, eps=eps, V=V, p0=p0, p10=p10, q0=q0, q10=q10)
print('back in python after running function PpqFullN')
print('The fortran result is ', result_fort)
print('The python result is ', result_py)

# Define argument and return types for the function PpqFullG
api.PpqFullG.argtypes = [ctypes.c_double] * 11
api.PpqFullG.restype = ctypes.c_double

result_fort = api.PpqFullG(ctypes.c_double(Er), ctypes.c_double(Ep), ctypes.c_double(Eq), 
                           ctypes.c_double(F0), ctypes.c_double(s), ctypes.c_double(eps),
                           ctypes.c_double(V), ctypes.c_double(p0), ctypes.c_double(p10),
                           ctypes.c_double(q0), ctypes.c_double(q10))
result_py = ppq.PpqFullG(Er, Ep, Eq, F0=F0, s=s, eps=eps, V=V, p0=p0, p10=p10, q0=q0, q10=q10)
print('back in python after running function PpqFullG')
print('The fortran result is ', result_fort)
print('The python result is ', result_py)

# Define argument and return types for the function PpqN
api.PpqN.argtypes = [ctypes.c_double] * 12
api.PpqN.restype = ctypes.c_double

result_fort = api.PpqN(ctypes.c_double(Ep), ctypes.c_double(Eq), 
                       ctypes.c_double(a), ctypes.c_double(b),
                       ctypes.c_double(F0), ctypes.c_double(s), ctypes.c_double(eps),
                       ctypes.c_double(V), ctypes.c_double(p0), ctypes.c_double(p10),
                       ctypes.c_double(q0), ctypes.c_double(q10))
(result_py, _), _, _ = ppq.PpqN_safe_inspect(Ep, Eq, a=a, b=b, F0=F0, s=s, eps=eps, V=V, p0=p0, p10=p10, q0=q0, q10=q10)
print('back in python after running function PpqN')
print('The fortran result is ', result_fort)
print('The python result is ', result_py)

# Define argument and return types for the function PpqG
api.PpqG.argtypes = [ctypes.c_double] * 10
api.PpqG.restype = ctypes.c_double

result_fort = api.PpqG(ctypes.c_double(Ep), ctypes.c_double(Eq), 
                       ctypes.c_double(F0), ctypes.c_double(s), ctypes.c_double(eps),
                       ctypes.c_double(V), ctypes.c_double(p0), ctypes.c_double(p10),
                       ctypes.c_double(q0), ctypes.c_double(q10))
(result_py, _), _, _ = ppq.PpqG_safe_inspect(Ep, Eq, F0=F0, s=s, eps=eps, V=V, p0=p0, p10=p10, q0=q0, q10=q10)
print('back in python after running function PpqG')
print('The fortran result is ', result_fort)
print('The python result is ', result_py)

# Set argtypes and restype for the function PpqN_vector
api.PpqN_vector.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # Ep_arr
    np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # Eq_arr
    ctypes.c_int,              # n
    ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
    ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
    ctypes.c_double, ctypes.c_double,  # scalar params
    np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")   # res_arr (output)
]
api.PpqN_vector.restype = None

# Inputs
Eq_arr = np.array([100.0, 100.0, 100.0], dtype=np.float64)
Ep_arr = np.array([347.0, 346.0, 348.0], dtype=np.float64)
n = Ep_arr.size

# Output array
res_arr = np.empty(n, dtype=np.float64)

# Call the Fortran vectorized function
api.PpqN_vector(Ep_arr, Eq_arr, n, a, b, F0, s, eps, V, p0, p10, q0, q10, res_arr)

# Set argtypes and restype for the function PpqG_vector
api.PpqG_vector.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),      # Ep_arr
    np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),      # Eq_arr
    ctypes.c_int,                                                        # n
    ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
    ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,  # scalar params
    np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")       # res_arr (output)
]
api.PpqG_vector.restype = None

# Inputs
Eq_arr = np.array([100.0, 100.0, 100.0], dtype=np.float64)
Ep_arr = np.array([347.0, 346.0, 348.0], dtype=np.float64)
n = Ep_arr.size

# Output array
res_arr = np.empty(n, dtype=np.float64)

# Call the Fortran vectorized function
api.PpqG_vector(Ep_arr, Eq_arr, n, F0, s, eps, V, p0, p10, q0, q10, res_arr)

print("Vectorized result:", res_arr)

print("Parameters Are ####################")
print("Er: ", Er)
print("Ep: ", Ep)
print("Eq: ", Eq)
print("a: ", a)
print("b: ", b)
print("F0: ", F0)
print("s: ", s)
print("eps: ", eps)
print("V: ", V)
print("p0: ", p0)
print("p10: ", p10)
print("q0: ", q0)
print("q10: ", q10)
