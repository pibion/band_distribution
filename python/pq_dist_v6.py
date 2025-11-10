import numpy as np
from scipy import special
import scipy.integrate as integrate
import scipy.optimize as so
from scipy.optimize import curve_fit
import os
try:
    import ML_estimation_functions as ml
except:
    pass

# order of parameters: specified by this dictionary 
# When building functions: leave out un-needed parameters, always give default values
def get_par_idx():
    return {'a':0, 'b':1, 'F':2, 'F0':3, 's':4, 'eps':5, 'V':6, 'p0':7, 'p10':8,'q0':9,'q10':10,\
           'PNa':11, 'PNb':12, 'PNd':13, 'PGa':14 , 'PGb':15 , 'PGd':16}

#get default
def get_par_default():
    return {'a':0.16, 'b':0.18, 'F':0.122, 'F0':0.0, 's':0.0, 'eps':3.0e-3, 'V':3.0, 'p0':0.06421907, 'p10':0.48998486,'q0':0.23718488,'q10':0.27093151,\
           'PNa':0.53693208, 'PNb':6.41515782, 'PNd':23.71789286, 'PGa':5.73211975e-01 , 'PGb':1.69520023e-01 , 'PGd':2.79552394e+02}

def get_par_old_default():
    return {'a':0.16, 'b':0.18, 'F':0.0, 'F0':0.12, 's':0.0, 'eps':3.0e-3, 'V':3.0, 'p0':0.5, 'p10':0.956,'q0':0.1,'q10':0.306,\
           'PNa':0.9629, 'PNb':5.468, 'PNd':29.4900, 'PGa':1.0 , 'PGb':100.0 , 'PGd':10000000000000000000000.0}

#get default ERs
def get_par_defaultER():
    #only difference is yield parameters must be a=1, b=0; also F0=s=0
    pars = get_par_default()
    pars['a'] = 1.0
    pars['b'] = 0.0
    return pars

#get old default ERs
def get_par_old_defaultER():
    #only difference is yield parameters must be a=1, b=0; also F0=s=0
    pars = get_par_old_default()
    pars['a'] = 1.0
    pars['b'] = 0.0
    return pars

#the ionization yield. Could use fancier models later.
def Y(Er,
	*,
	a=0.16,b=0.18):
    #ER should be in keV
    #a=0.16; b=0.18 from previous NRFano paper
    #return a*np.absolute(Er)**0.6
    if isinstance(Er, (list, tuple, np.ndarray)):
        print ("Er in Y is an array")

    if isinstance(a*np.absolute(Er)**b, (list, tuple, np.ndarray)):
        print ("Somehow the yield is evaluating to an array")
        print ("a, Er, and b are ", a, Er, b)

    return a*np.absolute(Er)**b

#define the Fano Factors for ERs and NRs
#default for ER is average of fano range from: https://arxiv.org/pdf/2206.13639.pdf
def F(Er,
	*,
	F0=0.122,s=0.0):
    #s is slope in fano per keV
    #F is Fano at zero energy
    return F0 + s*Er

#average numbers of e/h pairs
def Nbar(Er,
	  *,
	  a=0.16,b=0.18,
	  eps=3e-3):
    
    if isinstance(Y(Er,a=a,b=b), (list, tuple, np.ndarray)):
        print ("Y in Nbar is an array")
    return Y(Er,a=a,b=b)*Er/eps

#phonon and ionization resolutions
def sigp(Ep,
	  *,
	  eps=3.0e-3,
	  V=3.0,
	  p0=0.06421907, p10=0.48998486):
    e2pre = (p10**2 - p0**2)
    e2e = (Ep/(10*(1+(V/eps/1000))))**2
    #e2e = (Ep/(10*((V/eps/1000))))**2 #THIS ERROR WAS IN pq_dist v1!! it makes a difference!
    return np.sqrt(p0**2+e2pre*e2e)

def sigq(Eq,
	  *,
	  q0=0.23718488, q10=0.27093151):
    e2pre = (q10**2 - q0**2)
    e2e = (Eq/10)**2
    return np.sqrt(q0**2+e2pre*e2e)

#probability distribution of _true_ recoil energy, separate for gammas, Cfgammas or neutrons
#fit parameters A=0.9629
def PErN(Er,
              *,
              PNa=0.53693208,
              PNb=6.41515782, PNd=23.71789286):
    Er = np.asarray(Er)  # Convert input to a NumPy array
    result = np.where(Er < 0, 0, PNa * (1 / PNb) * np.exp(-Er / PNb) + (1 - PNa) * (1 / PNd) * np.exp(-Er / PNd))
    return result 
#...TBD gammas

# Cfgammas

def PErG(Er, *, PGa=5.73211975e-01, PGb=1.69520023e-01, PGd=2.79552394e+02):
    Er = np.asarray(Er)  # Convert Er to a NumPy array
    result = np.where(Er < 0, 0, PGa * (1 / PGb) * np.exp(-Er / PGb) + (1 - PGa) * (1 / PGd) * np.exp(-Er / PGd))
    return result   

def aN(Er,Ep,Eq,
	      *,
	      F0=0.122,s=0.0,
	      eps=3.0e-3,
	      V=3.0,
	      p0=0.06421907, p10=0.48998486,
	      q0=0.23718488, q10=0.27093151):
    t1 = (V/1e3)*(Ep-Er)/sigp(Ep,eps=eps,V=V,p0=p0,p10=p10)**2
    t2 = eps*Eq/sigq(Eq,q0=q0,q10=q10)**2
    t3 = 1/F(Er,F0=F0,s=s)
    return t1+t2+t3

def bN(Er,Ep,Eq,
	      *,
	      a=0.16,b=0.18,
	      F0=0.122,s=0.0,
	      eps=3.0e-3,
	      V=3.0,
	      p0=0.06421907, p10=0.48998486,
	      q0=0.23718488, q10=0.27093151):
    t1 = 1/(2*Nbar(Er,a=a,b=b,eps=eps)*F(Er,F0=F0,s=s))
    t2 = eps**2/(2*sigq(Eq,q0=q0,q10=q10)**2)
    t3 = V**2/(2*(sigp(Ep,eps=eps,V=V,p0=p0,p10=p10)*1e3)**2)

    if isinstance(t1, (list, tuple, np.ndarray)):
        print ("t1 in bN is an array")

    if isinstance(t2, (list, tuple, np.ndarray)):
        print ("t2 in bN is an array")

    if isinstance(t3, (list, tuple, np.ndarray)):
        print ("t3 in bN is an array")

    return t1+t2+t3

def cN(Er,Ep,Eq,
	      *,
	      a=0.16,b=0.18,
	      F0=0.122,s=0.0,
	      eps=3.0e-3,
	      V=3.0,
	      p0=0.06421907, p10=0.48998486,
	      q0=0.23718488, q10=0.27093151):
    t1 = -(Ep-Er)**2/(2*sigp(Ep,eps=eps,V=V,p0=p0,p10=p10)**2)
    t2 = -Eq**2/(2*sigq(Eq,q0=q0,q10=q10)**2)
    t3 = -Nbar(Er,a=a,b=b,eps=eps)/(2*F(Er,F0=F0,s=s))
    return t1+t2+t3

def CN(Er,Ep,Eq,
	      *,
	      a=0.16,b=0.18,
	      F0=0.122,s=0.0,
	      eps=3.0e-3,
	      V=3.0,
	      p0=0.06421907, p10=0.48998486,
	      q0=0.23718488, q10=0.27093151):
    return 1/(2*np.pi*np.sqrt(2*np.pi)*sigp(Ep,eps=eps,V=V,p0=p0,p10=p10)*sigq(Eq,q0=q0,q10=q10)*np.sqrt(Nbar(Er,a=a,b=b,eps=eps)*F(Er,F0=F0,s=s)))

def PpqExp(Er,Ep,Eq,
	          *,
	          a=0.16,b=0.18,
	          F0=0.122,s=0.0,
	          eps=3.0e-3,
	          V=3.0,
	          p0=0.06421907, p10=0.48998486,
	          q0=0.23718488, q10=0.27093151):
    cN_val = cN(Er,Ep,Eq,a=a,b=b,F0=F0,s=s,eps=eps,V=V,p0=p0,p10=p10,q0=q0,q10=q10)
    aN_val = aN(Er,Ep,Eq,F0=F0,s=s,eps=eps,V=V,p0=p0,p10=p10,q0=q0,q10=q10)
    bN_val = bN(Er,Ep,Eq,a=a,b=b,F0=F0,s=s,eps=eps,V=V,p0=p0,p10=p10,q0=q0,q10=q10)

    if isinstance(cN_val, (list, tuple, np.ndarray)):
        print ("cN is returning an array")

    if isinstance(aN_val, (list, tuple, np.ndarray)):
        print ("aN is returning an array")

    if isinstance(bN_val, (list, tuple, np.ndarray)):
        print ("bN is returning an array")

    exponent = cN_val + (aN_val**2 / (4*bN_val))
    return exponent

def PpqPreExp(Er,Ep,Eq,
		     *,
		     a=0.16,b=0.18,
		     F0=0.122,s=0.0,
		     eps=3.0e-3,
		     V=3.0,
		     p0=0.06421907, p10=0.48998486,
		     q0=0.23718488, q10=0.27093151):
    exponent = PpqExp(Er,Ep,Eq,a=a,b=b,F0=F0,s=s,eps=eps,V=V,p0=p0,p10=p10,q0=q0,q10=q10)
    return (1/2)*np.sqrt(np.pi)*CN(Er,Ep,Eq,a=a,b=b,F0=F0,s=s,eps=eps,V=V,p0=p0,p10=p10,q0=q0,q10=q10)*np.exp(exponent)*(1/np.sqrt(bN(Er,Ep,Eq,a=a,b=b,F0=F0,s=s,eps=eps,V=V,p0=p0,p10=p10,q0=q0,q10=q10)))

def PpqFullN(Er,Ep,Eq,
		    *,
		    a=0.16,b=0.18,
		    F0=0.122,s=0.0,
		    eps=3.0e-3,
		    V=3.0,
		    p0=0.06421907, p10=0.48998486,
		    q0=0.23718488, q10=0.27093151):
    prefac=PpqPreExp(Er,Ep,Eq,a=a,b=b,F0=F0,s=s,eps=eps,V=V,p0=p0,p10=p10,q0=q0,q10=q10)
    erffac=special.erf(aN(Er,Ep,Eq,F0=F0,s=s,eps=eps,V=V,p0=p0,p10=p10,q0=q0,q10=q10)/(2*np.sqrt(bN(Er,Ep,Eq,a=a,b=b,F0=F0,s=s,eps=eps,V=V,p0=p0,p10=p10,q0=q0,q10=q10))))+1
    return prefac*erffac*PErN(Er)

def PpqN(Ep,Eq,
	     *,
	     a=0.16,b=0.18,
	     F0=0.122,s=0.0,
	     eps=3.0e-3,
	     V=3.0,
	     p0=0.06421907, p10=0.48998486,
	     q0=0.23718488, q10=0.27093151):
    Er0 = Ep-(V*Eq/eps)*1e-3
    Er0Y = a*(np.absolute(Er0)**(b))*Er0
    f = lambda er,ep,eq: -1*PpqExp(er,ep,eq,a=a,b=b,F0=F0,s=s,eps=eps,V=V,p0=p0,p10=p10,q0=q0,q10=q10)
    loc = so.fmin(f,x0 = np.absolute(Er0Y), args=(Ep,Eq), disp = False)
    g = lambda er,ep,eq: PpqFullN(er,ep,eq,a=a,b=b,F0=F0,s=s,eps=eps,V=V,p0=p0,p10=p10,q0=q0,q10=q10)
    sca = 10 #np.sqrt(1/(2*W0))
    NoDevs = 5
    minx = np.max([(loc[0]-NoDevs*sca),0])
    maxx = np.max([(loc[0]+NoDevs*sca),0])
    ans = integrate.quad(g, minx, maxx, args = (Ep,Eq,))
    # (integral, abs_error), (er_arr, er_integrand_arr), (Er_low, Er_high, Er_max)
    return ans, (None, None), (minx, maxx, loc)

def PpqNfast(Ep,Eq,
		 *,
		 a=0.16,b=0.18,
		 F0=0.122,s=0.0,
		 eps=3.0e-3,
		 V=3.0,
		 p0=0.06421907, p10=0.48998486,
		 q0=0.23718488, q10=0.27093151,
		 imn=5,imx=50):
    g = lambda er,ep,eq: PpqFullN(er,ep,eq,a=a,b=b,F0=F0,s=s,eps=eps,V=V,p0=p0,p10=p10,q0=q0,q10=q10)
    ans = integrate.quad(g, imn, imx, args = (Ep,Eq,))[0]
    return ans

def PpqNfast2_inspect(Ep,Eq,
		  *,
		  a=0.16,b=0.18,
		  F0=0.122,s=0.0,
		  eps=3.0e-3,
		  V=3.0,
		  p0=0.06421907, p10=0.48998486,
		  q0=0.23718488, q10=0.27093151,
		  imn=0.01,imx=400,N=1000):
    g = lambda er,ep,eq: PpqFullN(er,ep,eq,a=a,b=b,F0=F0,s=s,eps=eps,V=V,p0=p0,p10=p10,q0=q0,q10=q10)

    bins=np.linspace(imn,imx,N)
    de=np.mean(np.diff(bins))
    er_integrand_arr = np.zeros(np.size(bins))
    er_arr = np.zeros(np.size(bins))
    ans=0
    for idx, e in enumerate(bins):
      ans+=g(e,Ep,Eq)*de
      er_arr[idx] = e
      #er_integrand_arr[idx] = g(e,Ep,Eq)

    return ans, (er_arr, er_integrand_arr)

def PpqNfast2a_inspect(Ep, Eq,
              *,
              a=0.16, b=0.18,
              F0=0.122, s=0.0,
              eps=3.0e-3,
              V=3.0,
              p0=0.06421907, p10=0.48998486,
              q0=0.23718488, q10=0.27093151,
              imn=5,imx=50,N=1000, tol=1e-4,
              log_dir="logs"):
    # Ensure the log directory exists
    os.makedirs(log_dir, exist_ok=True)
    do_logging = False

    # Define the log file path
    log_file_path = os.path.join(log_dir, "PpqNfast2a_log.txt")
    with open(log_file_path, "w") as log_file:
        
        # Define the function (replace PpqFullN with your actual implementation)
        g = lambda er, ep, eq: PpqFullN(er, ep, eq, a=a, b=b, F0=F0, s=s, eps=eps, V=V, p0=p0, p10=p10, q0=q0, q10=q10)
        de = (imx-imn)/N
        # Get the most likely Er
        Ermx = Ep - (V / 1e3) * (Eq / eps) * 0.3
        if Ermx < 0:
            log_file.write(f"Error encountered in PpqNfast2a: {Ermx} is less than zero")
            Ermx =0

        # Initialize variables
        ans = 0
        ansprev = 1.0
        highidx = 1
        lowidx = 0
        er_arr = []
        er_integrand_arr = []
        log_file.write("Starting PpqNfast2a calculations\n")
        log_file.write(f"Initial Ermx: {Ermx}\n")
        log_file.write(f"Initial highidx: {highidx}, lowidx: {lowidx}\n\n")
        # Loop from the highest contribution
        while (ansprev / (ansprev + ans)) > tol:
            if do_logging:
                log_file.write(f"Iteration start: ans={ans}, ansprev={ansprev},Ep and Eq={Ep,Eq},comparison with tol={ansprev / (ansprev + ans)}\n" )

            g_val = g(Ermx + de*highidx, Ep, Eq)
            high_contribution = g_val * de
            er_integrand_arr.append(g_val)
            er_arr.append(Ermx + de*highidx)
            ansprev = high_contribution
            highidx += 1
            if do_logging:
                log_file.write(f"High index {highidx - 1}: contribution={high_contribution}\n")

            if Ermx - de*lowidx > 0:
                g_val = g(Ermx - de*lowidx, Ep, Eq)
                low_contribution = g_val * de
                er_integrand_arr.insert(0, g_val)
                er_arr.insert(0, Ermx - de*lowidx)
                ansprev += low_contribution
                lowidx += 1
                if do_logging:
                    log_file.write(f"Low index {lowidx + 1}: contribution={low_contribution}\n")

            ans += ansprev
            if do_logging:
                log_file.write(f"Updated ans={ans}, ansprev={ansprev}\n\n")

    return ans, (er_arr, er_integrand_arr)

def PpqN_fast3_inspect(Ep, Eq,
              *,
              a=0.16, b=0.18,
              F0=0.122, s=0.0,
              eps=3.0e-3,
              V=3.0,
              p0=0.06421907, p10=0.48998486,
              q0=0.23718488, q10=0.27093151,
              tol=1e-6):
    Ermx, width = ml.estimateN_ML(Ep, Eq, a, b, F0, s, eps, p0, p10, q0, q10, model_suffix="600params")
    g = lambda er, ep, eq: PpqFullN(er, ep, eq, a=a, b=b, F0=F0, s=s, eps=eps, V=V, p0=p0, p10=p10, q0=q0, q10=q10)
    return integrate_g_fast_inspect(g, Ep, Eq, Ermx, width, tol)

def PpqN_fast4_inspect(Ep, Eq, interpolator,
              *,
              a=0.16, b=0.18,
              F0=0.122, s=0.0,
              eps=3.0e-3,
              V=3.0,
              p0=0.06421907, p10=0.48998486,
              q0=0.23718488, q10=0.27093151,
              tol=1e-6, threshold=1e-20):
    # function call is interpolate(self, ep, eq, q10, F0, a, b, p0, p10, method="nearest"):
    Ermx, width, max_integrand_val = interpolator.interpolate(Ep, Eq, q10, F0, a, b, p0, p10, method="nearest")
    print ("max integrand value is ", max_integrand_val)
    if max_integrand_val * width * 12 > threshold:
        # only integrate if there's a chance the value
        # is over the user-specificied threshold
        g = lambda er, ep, eq: PpqFullN(er, ep, eq, a=a, b=b, F0=F0, s=s, eps=eps, V=V, p0=p0, p10=p10, q0=q0, q10=q10)
        return integrate_g_fast_inspect(g, Ep, Eq, Ermx, width, tol)
    else:
        # in most cases the function evaluates to zero 
        return (0, threshold), (None, None), None

def integrate_g_fast_inspect(g, Ep, Eq, Ermx, width, tol):
    # how many sigmas to include in the interval
    sigma = 5
    fraction_reduced = 1e-5
    peak_info_array = [{"Ermx": Ermx, "width": width}]

    n_start = 1000
    minx = max(0, Ermx - sigma*width)
    maxx = Ermx + sigma*width
    er_arr = np.linspace(minx, maxx, n_start)
    er_integrand_arr = np.zeros(np.size(er_arr))    
    for idx, er in enumerate(er_arr):
        val = g(er, Ep, Eq)
        er_integrand_arr[idx] = val

    # sometimes the value at 0 keV is NaN
    # don't want to include that when getting the max value
    # max(any array that contains a NaN) = NaN
    max_val = max(er_integrand_arr[~np.isnan(er_integrand_arr)])
    print("max val is ", max_val)

    # check that the first value is small enough
    # the outcome of this loop is that minx is set correctly
    # er_integrand_arr[0] is sometimes NaN
    first_val = er_integrand_arr[1]
    print("first val is ", first_val)
    i = 0
    adjust_bounds_bool = False
    while True:
        if not first_val == 0:
            if minx > 0 and first_val/max_val > fraction_reduced:
                i += 1
                minx = max(0, Ermx - (sigma + 2**i)*width)
                first_val = g(minx, Ep, Eq)
                if first_val > max_val:
                    max_val = first_val
                adjust_bounds_bool = True
                print ("Need to adjust lower bound")
            else:
                break
        else:
            break

    # check that the last value is small enough
    # the outcome of this loop is that maxx is set correctly
    last_val = er_integrand_arr[-1]
    i = 0
    while True:
        if not last_val == 0:
            if last_val/max_val > fraction_reduced:
                i += 1
                maxx = Ermx + (sigma + 2**i)*width
                last_val = g(maxx, Ep, Eq)
                if last_val > max_val:
                    max_val = last_val
                adjust_bounds_bool = True
                print ("Need to adjust upper bound")
            else:
                break
        else:
            break

    peak_info_array[0].update({"minx": minx, "maxx": maxx})

    # see formula at https://en.wikipedia.org/wiki/Riemann_sum
    # for error of integration estimate due to midpoint rule
    # also: for a guassian,
    # the absolute value of the maximum value of the second derivative
    # is twice the max amplitude
    n = np.sqrt(2*max_val*np.power(maxx - minx,3) / (24 * tol))
    # make n an integer so we can use it in the linspace call
    n = int(np.ceil(n))

    if n > n_start or adjust_bounds_bool:
        if n > n_start:
            print ("have to resample to achieve requested tolerance, required n is ", n)
        else:
            # don't want to decrease the sampling 
            n = n_start

        if adjust_bounds_bool:
            print ("Bounds were adjusted")
        er_arr = np.linspace(minx, maxx, n)
        er_integrand_arr = np.zeros(np.size(er_arr))
        for idx, er in enumerate(er_arr):
            val = g(er, Ep, Eq)
            er_integrand_arr[idx] = val
        err = 2*max_val*np.power(maxx - minx,3) / (24 * n**2)
        peak_info_array[0].update({"n": n})
    else:
        err = 2*max_val*np.power(maxx - minx,3) / (24 * n_start**2)
        peak_info_array[0].update({"n": n_start})

    #print("er_arr is ", er_arr)
    #print("minx: ", minx)
    #print("maxx: ", maxx)
    #print("n: ", n)

    # now do integration using the midpoint rule
    # sometimes the value at 0 keV is NaN
    # don't want to include that in the calculation
    delta = np.diff(er_arr)[0]
    integral = np.sum(er_integrand_arr[~np.isnan(er_integrand_arr)]) * delta

    return (integral, err), (er_arr, er_integrand_arr), peak_info_array

def PpqN_safe_inspect_vec(Ep, Eq,
              *,
              a=0.16, b=0.18,
              F0=0.122, s=0.0,
              eps=3.0e-3,
              V=3.0,
              p0=0.06421907, p10=0.48998486,
              q0=0.23718488, q10=0.27093151,
              res=0.1,
              log_dir="logs"):
    ppqNArr = []
    for this_Ep, this_Eq in zip(Ep, Eq):
        PpqNval, _, _ = PpqN_safe_inspect(this_Ep, this_Eq, a=a, b=b, F0=F0, s=s, eps=eps, V=V, p0=p0, p10=p10, q0=q0, q10=q10, res=res, log_dir=log_dir)
        ppqNArr.append(PpqNval[0])
    return ppqNArr

"""
Usage to just get the value of PpqN:
(PpqN, _), _, _ = PpqN_safe_inspect(args)
"""
def PpqN_safe_inspect(Ep, Eq,
              *,
              a=0.16, b=0.18,
              F0=0.122, s=0.0,
              eps=3.0e-3,
              V=3.0,
              p0=0.06421907, p10=0.48998486,
              q0=0.23718488, q10=0.27093151,
              res=0.1,
              log_dir="logs"):
    # find Er_max
    f = lambda er,ep,eq: -1*PpqExp(er,ep,eq,a=a,b=b,F0=F0,s=s,eps=eps,V=V,p0=p0,p10=p10,q0=q0,q10=q10)

    # WARNING
    # this will only work if your peaks have a width
    # greater than 0.1 keV
    # which seems to be true for all reasonable CDMS parameters
   # the advantage is that the minimizer sometimes fails and this never does
    maxx = max(Ep, Eq) + 5
    er_arr = np.arange(0, maxx, res)
    f_arr = np.zeros(len(er_arr), dtype=np.float64)
    for idx, er in enumerate(er_arr):
        #print (f(er, Ep, Eq), type(f(er, Ep, Eq)))
        try:
            f_arr[idx] = -1*PpqExp(er,Ep,Eq,a=a,b=b,F0=F0,s=s,eps=eps,V=V,p0=p0,p10=p10,q0=q0,q10=q10)
        except Exception as e:
            pass
            # print (er, Ep, Eq)
            # print (-1*PpqExp(er,Ep,Eq,a=a,b=b,F0=F0,s=s,eps=eps,V=V,p0=p0,p10=p10,q0=q0,q10=q10))
    real_val_idx = np.where(f_arr > 0)[0]
    if len(real_val_idx) == 0:
        if np.all(f_arr == 0):
            # print ("we really couldn't find any real values, there is likely no peak")
            return (0, 0), (er_arr, f_arr), []
        elif np.all(np.isnan(f_arr)):
            # print ("all the evaluated values are NaN!")
            return (np.nan, np.nan), (er_arr, f_arr), []
    else:
        min_idx = np.argmin(f_arr[~np.isnan(f_arr)])
        Ermx = er_arr[min_idx]
        # print ("We found a peak at ", Ermx, " in the exponent term")

    # Define the function (replace PpqFullG with your actual implementation)
    g = lambda er, ep, eq: PpqFullN(er, ep, eq, a=a, b=b, F0=F0, s=s, eps=eps, V=V, p0=p0, p10=p10, q0=q0, q10=q10)

    return integrate_g_safe_inspect(g, Ep, Eq, Ermx, res)

"""
Return just the answer of PpqNfast2
"""
def PpqNfast2(Ep,Eq,
		  *,
		  a=0.16,b=0.18,
		  F0=0.122,s=0.0,
		  eps=3.0e-3,
		  V=3.0,
		  p0=0.06421907, p10=0.48998486,
		  q0=0.23718488, q10=0.27093151,
		  imn=5,imx=50,N=1000):
    
    return PpqNfast2_inspect(Ep, Eq, a=a, b=b, F0=F0, s=s, eps=eps, V=V, p0=p0, p10=p10, q0=q0, q10=q10, imn=imn, imx=imx, N=N)[0]

"""
Return just the answer of PpqNfast2a
"""
def PpqNfast2a(Ep,Eq,
		  *,
		  a=0.16,b=0.18,
		  F0=0.122,s=0.0,
		  eps=3.0e-3,
		  V=3.0,
		  p0=0.06421907, p10=0.48998486,
		  q0=0.23718488, q10=0.27093151,
		  imn=5,imx=50,N=1000):
    
    return PpqNfast2a_inspect(Ep, Eq, a=a, b=b, F0=F0, s=s, eps=eps, V=V, p0=p0, p10=p10, q0=q0, q10=q10, imn=imn, imx=imx, N=N)[0]

def PpqNfast3(Ep,Eq,
		  *,
		  a=0.16,b=0.18,
		  F0=0.122,s=0.0,
		  eps=3.0e-3,
		  V=3.0,
		  p0=0.06421907, p10=0.48998486,
		  q0=0.23718488, q10=0.27093151,
		  imn=5,imx=50,algo='gk21'):
    g = lambda er,ep,eq: PpqFullN(er,ep,eq,a=a,b=b,F0=F0,s=s,eps=eps,V=V,p0=p0,p10=p10,q0=q0,q10=q10)
    ans = integrate.quad_vec(g, imn, imx, args = (Ep,Eq,),quadrature=algo)[0]
    return ans

def PpqFullG(Er,Ep,Eq,
		    *,
		    F0=0.122,s=0.0,
		    eps=3.0e-3,
		    V=3.0,
		    p0=0.06421907, p10=0.48998486,
		    q0=0.23718488, q10=0.27093151):
    prefac=PpqPreExp(Er,Ep,Eq,a=1.0,b=0.0,F0=F0,s=s,eps=eps,V=V,p0=p0,p10=p10,q0=q0,q10=q10)
    erffac=special.erf(aN(Er,Ep,Eq,F0=F0,s=s,eps=eps,V=V,p0=p0,p10=p10,q0=q0,q10=q10)/(2*np.sqrt(bN(Er,Ep,Eq,a=1.0,b=0.0,F0=F0,s=s,eps=eps,V=V,p0=p0,p10=p10,q0=q0,q10=q10))))+1
    return prefac*erffac*PErG(Er)

def PpqG(Ep,Eq,
	     *,
	     F0=0.122,s=0.0,
	     eps=3.0e-3,
	     V=3.0,
	     p0=0.06421907, p10=0.48998486,
	     q0=0.23718488, q10=0.27093151):
    Er0 = Ep-(V*Eq/eps)*1e-3
    f = lambda er,ep,eq: -1*PpqExp(er,ep,eq,a=1.0,b=0.0,F0=F0,s=s,eps=eps,V=V,p0=p0,p10=p10,q0=q0,q10=q10)
    loc = so.fmin(f,x0 = np.absolute(Er0), args=(Ep,Eq), disp = False)
    g = lambda er,ep,eq: PpqFullG(er,ep,eq,F0=F0,s=s,eps=eps,V=V,p0=p0,p10=p10,q0=q0,q10=q10)
    sca = 10 #np.sqrt(1/(2*W0))
    NoDevs = 5
    minx = np.max([(loc[0]-NoDevs*sca),0])
    maxx = np.max([(loc[0]+NoDevs*sca),0])
    ans = integrate.quad(g, minx, maxx, args = (Ep,Eq,))[0]
    return ans

def PpqGfast(Ep,Eq,
		 *,
		 F0=0.122,s=0.0,
		 eps=3.0e-3,
		 V=3.0,
		 p0=0.06421907, p10=0.48998486,
		 q0=0.23718488, q10=0.27093151,
		 imn=5,imx=50):
    g = lambda er,ep,eq: PpqFullG(er,ep,eq,F0=F0,s=s,eps=eps,V=V,p0=p0,p10=p10,q0=q0,q10=q10)
    ans = integrate.quad(g, imn, imx, args = (Ep,Eq,))[0]
    return ans

def PpqGfast2_inspect(Ep,Eq,
		  *,
		  F0=0.122,s=0.0,
		  eps=3.0e-3,
		  V=3.0,
		  p0=0.06421907, p10=0.48998486,
		  q0=0.23718488, q10=0.27093151,
		  imn=0.01,imx=400,N=1000):
    g = lambda er,ep,eq: PpqFullG(er,ep,eq,F0=F0,s=s,eps=eps,V=V,p0=p0,p10=p10,q0=q0,q10=q10)

    bins=np.linspace(imn,imx,N)
    de=np.mean(np.diff(bins))
    er_integrand_arr = np.zeros(np.size(bins))
    er_arr = np.zeros(np.size(bins))
    ans=0
    for idx, e in enumerate(bins):
      g_val = g(e,Ep,Eq)
      ans+=g_val*de
      er_arr[idx] = e
      er_integrand_arr[idx] =g_val

    return ans, (er_arr, er_integrand_arr)

#***********************
def PpqGfast2a_inspect(Ep, Eq,
              *,
              F0=0.122, s=0.0,
              eps=3.0e-3,
              V=3.0,
              p0=0.06421907, p10=0.48998486,
              q0=0.23718488, q10=0.27093151,
              imn=5, imx=50, N=1000, tol=1e-4,
              log_dir="logs"):
    # Ensure the log directory exists
    os.makedirs(log_dir, exist_ok=True)

    # Define the log file path
    log_file_path = os.path.join(log_dir, "PpqGfast2a_log.txt")
    with open(log_file_path, "w") as log_file:
        log_file.write("Starting PpqGfast2a calculations\n")
        # Define the function (replace PpqFullG with your actual implementation)
        g = lambda er, ep, eq: PpqFullG(er, ep, eq, F0=F0, s=s, eps=eps, V=V, p0=p0, p10=p10, q0=q0, q10=q10)
        de = (imx-imn)/N
        # Get the most likely Er
        Ermx = Ep - (V / 1e3) * (Eq / eps) 
        if Ermx < 0:
            log_file.write("Error encountered in PpqNfast2a: {Ermx} is less than zero")
            Ermx=0
            #print ("Error encountered in PpqNfast2a: Ermx is set to zero")
        
        # Initialize variables
        ans = 0
        ansprev = 1.0
        highidx = 1
        lowidx =0
        er_arr = []
        er_integrand_arr = []

        # Open the log file in write mode
    
        log_file.write(f"Initial Ermx: {Ermx}\n")
        log_file.write(f"Initial highidx: {highidx}, lowidx: {lowidx}\n\n")

        # Loop from the highest contribution
        while (ansprev / (ansprev + ans)) > tol:
            log_file.write(f"Iteration start: ans={ans}, ansprev={ansprev}\n")
            g_val = g(Ermx + de*highidx, Ep, Eq)
            high_contribution = g_val * de
            er_integrand_arr.append(g_val)
            er_arr.append(Ermx + de*highidx)
            ansprev = high_contribution
            highidx += 1
            log_file.write(f"High index {highidx - 1}: contribution={high_contribution}\n")

            if (Ermx - de*lowidx) > 0:
                g_val = g(Ermx - de*lowidx, Ep, Eq)
                low_contribution = g_val * de
                er_integrand_arr.insert(0, g_val)
                er_arr.insert(0, Ermx - de*lowidx)
                ansprev += low_contribution
                lowidx += 1
                log_file.write(f"Low index {lowidx + 1}: contribution={low_contribution}\n")

            ans += ansprev
            log_file.write(f"Updated ans={ans}, ansprev={ansprev}\n\n")
            
    return ans, (er_arr, er_integrand_arr)

def consecutive(data, stepsize=1):
    return np.split(data, np.where(np.diff(data) != stepsize)[0]+1)

# Define Gaussian function
def gaussian(x, A, mu, sigma):
    return A * np.exp(-((x - mu) ** 2) / (2 * sigma ** 2))

def integrate_g_safe_inspect(g, Ep, Eq, Ermx, resolution): 
    minx = 0
    maxx = max(Ep, Eq) + 5
    # print("minx and maxx are ", minx, maxx)

    # Amy thinks the following block might be dead code that could be deleted
    # it came from when we were trying to estimate Ermx,
    # but were doing so poorly
    # so we were seeing all zeros because we weren't near the maximum
    # Test for all zeros
    if g(Ermx, Ep, Eq) == 0 or np.isnan(g(Ermx, Ep, Eq)):
        # print ("The maximum value is smaller than machine precision")
        # print ("The value of the integrand is effectively zero")
        return (0, 0), (None, None), []
    else:
        pass
        # print ("The maximum value is ", g(Ermx, Ep, Eq))

    er_arr = np.arange(minx, maxx, resolution / 2)
    er_integrand_arr = np.zeros(np.size(er_arr))
    for idx, er in enumerate(er_arr):
        val = g(er, Ep, Eq)
        er_integrand_arr[idx] = val
    #print ("are any of the values non-zero? ", np.where(er_integrand_arr > 0))

    # we expect a structure that is peaks, surrounded by zeros
    nonZero_idx = np.where(er_integrand_arr > 0)[0]
    nonZero_idx_groups = consecutive(nonZero_idx)

    peak_info_array = []

    for peak_idx in nonZero_idx_groups:
        while True:
            min_idx = max(0, min(peak_idx) - 5)
            max_idx = min(np.size(er_integrand_arr) - 1, max(peak_idx) + 5)

            minx = er_arr[min_idx]
            maxx = er_arr[max_idx]

            # if the value array looks like 0 0 0 4.5 6.7 8 9 10.5 5 1 0 0 0 0 
            # then the peak_idx array will look like 3 4 5 6 7 8 9
            # and the max_idx will be less than the size of er_integrand_arr
            # if this is the case then we don't to expand the limits
            # because the upper limit already goes past the place where the integrand
            # drops to zero
            if max_idx < np.size(er_integrand_arr) - 1:
                break
            else:
                #print ("Expanding the upper limit")
                maxx = 1.5*maxx

                er_arr = np.arange(minx, maxx, resolution / 2)
                er_integrand_arr = np.zeros(np.size(er_arr))
                for idx, er in enumerate(er_arr):
                    val = g(er, Ep, Eq)
                    er_integrand_arr[idx] = val
            
                peak_idx = np.where(er_integrand_arr > 0)[0]

        # if we've broken out of the loop, we have our minx and maxx
        # for this peak
        peak_info_array.append({"minx": minx, "maxx": maxx})

    # print(peak_info_array)
    # if there are two peak ranges, check that they are distinct
    if len(peak_info_array) > 1:
        if peak_info_array[0]['maxx'] > peak_info_array[1]['minx']:
            # the peak ranges are overlapping, you should remove the first one
            peak_info_array.pop(0)

    # Integrate using quad now that you have accurate limits
    for jdx, peak_info in enumerate(peak_info_array):
        minx = peak_info["minx"]
        maxx = peak_info["maxx"]
        # print ("final minx and maxx for peak ", jdx, " are ", minx, maxx)
        # integrate and record the value
        ans = integrate.quad(g, minx, maxx, args = (Ep,Eq,))
        peak_info_array[jdx]["integral"] = ans[0]
        peak_info_array[jdx]["error"] = ans[1]
        
        # this array is for the fitting
        er_arr = np.linspace(minx, maxx, 2000)
        er_integrand_arr = np.zeros(np.size(er_arr))
        for idx, er in enumerate(er_arr):
            er_integrand_arr[idx] = g(er, Ep, Eq)

        # find the width
        # Find maximum function value and corresponding energy 
        er_arr = er_arr.flatten()

        # Filter data to only include points where the value is not NaN
        mask_nan = np.isnan(er_integrand_arr)
        mask_zeros = er_integrand_arr == 0
        filtered_energies = er_arr[~mask_nan & ~mask_zeros]
        filtered_values = er_integrand_arr[~mask_nan & ~mask_zeros]

        max_index = np.argmax(filtered_values)
        max_value = filtered_values[max_index]
        max_energy = filtered_energies[max_index]
        width_estimate = np.nan
        #print (filtered_energies, filtered_values)

        try:
            # estimate the width of the function
            curvature = np.diff(np.diff(filtered_values))

            if np.all(curvature == 0):
                print ("this is a flat line!!")
                peak_info_array[jdx]["width"] = width_estimate
                peak_info_array[jdx]["Ermx"] = max_energy
                peak_info_array[jdx]['integrand_max'] = max_value

            else:
                # Walk to the right of max_index
                for i in range(max_index + 1, len(curvature) - 1):
                    if curvature[i] > 0 and curvature[i-1] <= 0:  
                        inflection_index = i + 2 # Return the index where it switches from negative to positive
                        break
                    
                width_estimate = filtered_energies[inflection_index] - filtered_energies[max_index]

                # Fit Gaussian to filtered data
                # note that we fit the values filtered_values * 1/max_value
                # this is so that when the values are small the fit still works!
                popt, _ = curve_fit(gaussian, filtered_energies, filtered_values * 1/max_value, p0=[1, max_energy, width_estimate])
                # print(f"Fitted Gaussian width (sigma): {width}")

                # Er_max is the value of Er for which g is maximum
                peak_info_array[jdx]['integrand_max'] = popt[0] * max_value
                peak_info_array[jdx]["width"] = popt[2]
                peak_info_array[jdx]["Ermx"] = popt[1]
        except Exception as e:
            print ("Could not fit width")
            peak_info_array[jdx]["width"] = width_estimate
            # not setting an Ermx will mean this call will fail to fill any values
            # I kind of want this to fail in a way that's noticeable
            #peak_info_array[jdx]["Ermx"] = max_energy
            print ("Starting parameters are, ", max_value, max_energy, width_estimate)
            print ("Ep and Eq are ", Ep, Eq)
            print (e)

    # sum the integrals from all the peaks to get the final total integral and error
    integral = 0
    error = 0
    for peak_info in peak_info_array:
        integral += peak_info["integral"]
        error += peak_info["error"]

    # find the overall minimum and maximum and evaluate er_integrand_arr over that range
    min_minx = min(d["minx"] for d in peak_info_array)
    max_maxx = max(d["maxx"] for d in peak_info_array)
    er_arr = np.arange(min_minx, max_maxx, resolution / 2)
    er_integrand_arr = np.zeros(np.size(er_arr))
    for idx, er in enumerate(er_arr):
        val = g(er, Ep, Eq)
        er_integrand_arr[idx] = val

    # now we're ready to return everything to the user
    return (integral, error), (er_arr, er_integrand_arr), peak_info_array

def PpqG_safe_inspect_vec(Ep, Eq,
              *,
              F0=0.122, s=0.0,
              eps=3.0e-3,
              V=3.0,
              p0=0.06421907, p10=0.48998486,
              q0=0.23718488, q10=0.27093151,
              res=0.1,
              log_dir="logs"):
    ppqGArr = []
    for this_Ep, this_Eq in zip(Ep, Eq):
        PpqGval, _, _ = PpqG_safe_inspect(this_Ep, this_Eq, F0, s, eps, V, p0, p10, q0, q10, res, log_dir)
        ppqGArr.append(PpqGval[0])
    return ppqGArr

"""
Usage to just get the value of PpqN:
(PpqG, _), _, _ = PpqG_safe_inspect(args)
"""
def PpqG_safe_inspect(Ep, Eq,
              *,
              F0=0.122, s=0.0,
              eps=3.0e-3,
              V=3.0,
              p0=0.06421907, p10=0.48998486,
              q0=0.23718488, q10=0.27093151,
              res=0.1,
              log_dir="logs"):
    # Set parameters for integration
    # with CDMS parameters, peaks are
    # never less than 1 keV wide
    # so a resolution of 0.1 will "catch" them reliably

    # Define the function (replace PpqFullG with your actual implementation)
    g = lambda er, ep, eq: PpqFullG(er, ep, eq, F0=F0, s=s, eps=eps, V=V, p0=p0, p10=p10, q0=q0, q10=q10)

    f = lambda er,ep,eq: -1*PpqExp(er,ep,eq,a=1.0,b=0.0,F0=F0,s=s,eps=eps,V=V,p0=p0,p10=p10,q0=q0,q10=q10)

    # WARNING
    # this will only work if your peaks have a width
    # greater than 0.1 keV
    # which seems to be true for all reasonable CDMS parameters
    maxx = max(Ep, Eq) + 5
    er_arr = np.arange(0, maxx, res)
    f_arr = np.zeros(len(er_arr))
    for idx, er in enumerate(er_arr):
        f_arr[idx] = f(er, Ep, Eq)

    if not np.any(f_arr > 0):
        if np.all(f_arr == 0):
            # print ("we really couldn't find any real values, there is likely no peak")
            return (0, 0), (er_arr, f_arr), []
        elif np.all(np.isnan(f_arr)):
           # print ("all the evaluated values are NaN!")
            return (np.nan, np.nan), (er_arr, f_arr), []
    else:
        min_idx = np.argmin(f_arr[~np.isnan(f_arr)])
        Ermx = er_arr[min_idx]
        # print ("We found a peak at ", Ermx, " in the exponent term")

    return integrate_g_safe_inspect(g, Ep, Eq, Ermx, res)

def PpqG_fast2b(Ep, Eq,
              *,
              F0=0.122, s=0.0,
              eps=3.0e-3,
              V=3.0,
              p0=0.06421907, p10=0.48998486,
              q0=0.23718488, q10=0.27093151,
              imn=5, imx=50, N=1000, tol=1e-6,
              log_dir="logs"):
    # Set parameters for integration
    change_tol = 0.1
    de_limit = 0.001

    # Define the function (replace PpqFullG with your actual implementation)
    g = lambda er, ep, eq: PpqFullG(er, ep, eq, F0=F0, s=s, eps=eps, V=V, p0=p0, p10=p10, q0=q0, q10=q10)
    de_high = (imx-imn)/N
    de_low = (imx-imn)/N
    # Get the most likely Er
    Ermx = Ep - (V / 1e3) * (Eq / eps) 
    if Ermx < 0:
        Ermx=0
        #print ("Error encountered in PpqNfast2a: Ermx is set to zero")
    print("Max Er for gamma is calculated as ", Ermx)

    # Initialize variables
    ans_high = 0
    ansprev_high = 1.0
    ans_low = 0
    ansprev_low = 1.0
    Er_val = Ermx
    er_arr = []
    er_integrand_arr = []
    first_pass = True
    de_high_first = 0

    # Do an initial check of the de 
    # Initially de_high and de_low are the same value
    de = de_high
    g_val = g(Er_val, Ep, Eq)
    g_val_advance =  g(Er_val + 0.5*de, Ep, Eq)
    if Er_val - de > 0:
        g_val_behind = g(Er_val - de, Ep, Eq)
    else:
        g_val_behind = g_val

    de_high_is_large = False
    de_low_is_large = False
    if np.abs((g_val_advance - g_val) / g_val) > change_tol:
        de_high_is_large = True
        print ("have to change de right away for upwards integral")
    
    if np.abs((g_val_behind - g_val) / g_val) > change_tol:
        de_low_is_large = True
        print ("have to change de right away for downwards integral")

    # Loop from the highest contribution
    while (ansprev_high / (ansprev_high + ans_high)) > tol/2:
        g_val = g(Er_val + 0.5*de_high, Ep, Eq)
        high_contribution = g_val * de_high

        if((np.abs((ansprev_high-high_contribution) / ansprev_high) < change_tol or de_high < de_limit) and not de_high_is_large):
            er_integrand_arr.append(g_val)
            er_arr.append(Er_val)
            ansprev_high = high_contribution
            Er_val += de_high
            ans_high += ansprev_high
            if first_pass:
                de_high_first = de_high
                first_pass = False
        elif (de_high > de_limit):
            de_high = de_high/5
            print ("changing de_high to ", de_high)
            de_high_is_large = False
        else:
            print ("WARNING: reached de_limit without reaching function change tolerance")

    # Reset Er_val so we can integrate downwards
    Er_val = Ermx - 0.5*de_low - 9.5*de_high_first

    while ((True if ans_low <= 0 else ansprev_low / ans_low > tol/2)) and (Er_val > 0):
        
        g_val = g(Er_val, Ep, Eq)
        low_contribution = g_val * de_low

        if((np.abs((ansprev_low-low_contribution) / ansprev_low) < change_tol or de_low < de_limit) and not de_low_is_large):
            er_integrand_arr.insert(0, g_val)
            er_arr.insert(0, Er_val)
            ansprev_low = low_contribution
            Er_val -= de_low
            ans_low += ansprev_low
        elif(de_low > de_limit):
            de_low = de_low/5
            de_low_is_large = False
            print ("changing de_low to ", de_low)
        else:
            print ("WARNING: reached de_limit without reaching function change tolerance")

    return (ans_low + ans_high, error), (er_arr, er_integrand_arr)

"""
Return just the answer of PpqGfast2a
"""
def PpqGfast2a(Ep,Eq,
		  *,
		  F0=0.122,s=0.0,
		  eps=3.0e-3,
		  V=3.0,
		  p0=0.06421907, p10=0.48998486,
		  q0=0.23718488, q10=0.27093151,
		  imn=5,imx=50,N=1000):
    
    return PpqGfast2a_inspect(Ep, Eq, F0=F0, s=s, eps=eps, V=V, p0=p0, p10=p10, q0=q0, q10=q10, imn=imn, imx=imx, N=N)[0]
    
"""
Return just the answer of PpqGfast2
"""
def PpqGfast2(Ep,Eq,
		  *,
		  F0=0.122,s=0.0,
		  eps=3.0e-3,
		  V=3.0,
		  p0=0.06421907, p10=0.48998486,
		  q0=0.23718488, q10=0.27093151,
		  imn=5,imx=50,N=1000):
    
    return PpqGfast2_inspect(Ep, Eq, F0=F0, s=s, eps=eps, V=V, p0=p0, p10=p10, q0=q0, q10=q10, imn=imn, imx=imx, N=N)[0]
    
