import numpy as np
import scipy.integrate as integrate
from scipy.optimize import curve_fit

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
    exponent = cN_val + (aN_val**2 / (4*bN_val))
    return exponent

def PpqFullN(Er, Ep, Eq,
             *,
             a=0.16, b=0.18,
             F0=0.122, s=0.0,
             eps=3.0e-3,
             V=3.0,
             p0=0.06421907, p10=0.48998486,
             q0=0.23718488, q10=0.27093151):
    F_val    = F(Er, F0=F0, s=s)
    Nbar_val = Nbar(Er, a=a, b=b, eps=eps)
    if Nbar_val <= 0.0 or F_val <= 0.0:
        return 0.0
    sigma_N    = np.sqrt(Nbar_val * F_val)
    norm_const = 1.0 / (2.0 * np.pi) ** 1.5
    sp0 = p0 ** 2
    sp1 = (p10 ** 2 - p0 ** 2) / (10.0 * (1.0 + V / (eps * 1.0e3))) ** 2
    sq0 = q0 ** 2
    sq1 = (q10 ** 2 - q0 ** 2) / 100.0
    EP_mean      = Er + (V * 1.0e-3) * Nbar_val
    EQ_mean      = eps * Nbar_val
    sigp_mean_sq = sp0 + sp1 * EP_mean ** 2
    sigq_mean_sq = sq0 + sq1 * EQ_mean ** 2
    a_coeff = ((V * 1.0e-3) * (Ep - Er) / sigp_mean_sq
               + eps * Eq / sigq_mean_sq
               + 1.0 / F_val)
    b_coeff = (1.0 / (2.0 * Nbar_val * F_val)
               + eps ** 2 / (2.0 * sigq_mean_sq)
               + V ** 2 / (2.0 * 1.0e6 * sigp_mean_sq))
    N_star     = a_coeff / (2.0 * b_coeff)
    sigma_Neff = 1.0 / np.sqrt(2.0 * b_coeff)
    if N_star - 5.0 * sigma_Neff < 0.0:
        N_lo, N_hi = 0.0, 10.0 * sigma_Neff
    else:
        N_lo, N_hi = N_star - 5.0 * sigma_Neff, N_star + 5.0 * sigma_Neff
    N_pts = np.linspace(N_lo, N_hi, 21)
    h_N   = (N_hi - N_lo) / 20.0
    EP_nl = Er + (V * 1.0e-3) * N_pts
    EQ_nl = eps * N_pts
    sig_p = np.sqrt(sp0 + sp1 * EP_nl ** 2)
    sig_q = np.sqrt(sq0 + sq1 * EQ_nl ** 2)
    exp_arg = (-0.5 * ((Ep - EP_nl) / sig_p) ** 2
               - 0.5 * ((Eq - EQ_nl) / sig_q) ** 2
               - 0.5 * ((N_pts - Nbar_val) / sigma_N) ** 2)
    integrand = np.zeros(21)
    mask = exp_arg >= -700.0
    integrand[mask] = (norm_const / (sigma_N * sig_p[mask] * sig_q[mask])
                       * np.exp(exp_arg[mask]))
    simps_w = np.array([1,4,2,4,2,4,2,4,2,4,2,4,2,4,2,4,2,4,2,4,1], dtype=float)
    return float((h_N / 3.0) * np.dot(simps_w, integrand) * PErN(Er))

def PpqN_safe_inspect_vec(Ep, Eq,
              *,
              a=0.16, b=0.18,
              F0=0.122, s=0.0,
              eps=3.0e-3,
              V=3.0,
              p0=0.06421907, p10=0.48998486,
              q0=0.23718488, q10=0.27093151,
              res=0.1):
    ppqNArr = []
    for this_Ep, this_Eq in zip(Ep, Eq):
        PpqNval, _, _ = PpqN_safe_inspect(this_Ep, this_Eq, a=a, b=b, F0=F0, s=s, eps=eps, V=V, p0=p0, p10=p10, q0=q0, q10=q10, res=res)
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
              res=0.1):
    if F0 == 0 and s == 0:
        raise ValueError("Fano factor F(Er) = F0 + s*Er is zero for all Er; F0 and s cannot both be zero")

    # find Er_max
    f = lambda er,ep,eq: -1*PpqExp(er,ep,eq,a=a,b=b,F0=F0,s=s,eps=eps,V=V,p0=p0,p10=p10,q0=q0,q10=q10)

    # WARNING
    # this will only work if your peaks have a width
    # greater than 0.1 keV
    # which seems to be true for all reasonable CDMS parameters
    # the advantage is that the minimizer sometimes fails and this never does
    maxx = max(Ep, Eq) + 5
    er_arr = np.arange(5e-21, maxx, res)
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

def PpqFullG(Er, Ep, Eq,
             *,
             F0=0.122, s=0.0,
             eps=3.0e-3,
             V=3.0,
             p0=0.06421907, p10=0.48998486,
             q0=0.23718488, q10=0.27093151):
    F_val    = F(Er, F0=F0, s=s)
    Nbar_val = Er / eps          # Y=1 for electron recoils
    if Nbar_val <= 0.0 or F_val <= 0.0:
        return 0.0
    sigma_N    = np.sqrt(Nbar_val * F_val)
    norm_const = 1.0 / (2.0 * np.pi) ** 1.5
    sp0 = p0 ** 2
    sp1 = (p10 ** 2 - p0 ** 2) / (10.0 * (1.0 + V / (eps * 1.0e3))) ** 2
    sq0 = q0 ** 2
    sq1 = (q10 ** 2 - q0 ** 2) / 100.0
    EP_mean      = Er + (V * 1.0e-3) * Nbar_val
    EQ_mean      = eps * Nbar_val
    sigp_mean_sq = sp0 + sp1 * EP_mean ** 2
    sigq_mean_sq = sq0 + sq1 * EQ_mean ** 2
    a_coeff = ((V * 1.0e-3) * (Ep - Er) / sigp_mean_sq
               + eps * Eq / sigq_mean_sq
               + 1.0 / F_val)
    b_coeff = (1.0 / (2.0 * Nbar_val * F_val)
               + eps ** 2 / (2.0 * sigq_mean_sq)
               + V ** 2 / (2.0 * 1.0e6 * sigp_mean_sq))
    N_star     = a_coeff / (2.0 * b_coeff)
    sigma_Neff = 1.0 / np.sqrt(2.0 * b_coeff)
    if N_star - 5.0 * sigma_Neff < 0.0:
        N_lo, N_hi = 0.0, 10.0 * sigma_Neff
    else:
        N_lo, N_hi = N_star - 5.0 * sigma_Neff, N_star + 5.0 * sigma_Neff
    N_pts = np.linspace(N_lo, N_hi, 21)
    h_N   = (N_hi - N_lo) / 20.0
    EP_nl = Er + (V * 1.0e-3) * N_pts
    EQ_nl = eps * N_pts
    sig_p = np.sqrt(sp0 + sp1 * EP_nl ** 2)
    sig_q = np.sqrt(sq0 + sq1 * EQ_nl ** 2)
    exp_arg = (-0.5 * ((Ep - EP_nl) / sig_p) ** 2
               - 0.5 * ((Eq - EQ_nl) / sig_q) ** 2
               - 0.5 * ((N_pts - Nbar_val) / sigma_N) ** 2)
    integrand = np.zeros(21)
    mask = exp_arg >= -700.0
    integrand[mask] = (norm_const / (sigma_N * sig_p[mask] * sig_q[mask])
                       * np.exp(exp_arg[mask]))
    simps_w = np.array([1,4,2,4,2,4,2,4,2,4,2,4,2,4,2,4,2,4,2,4,1], dtype=float)
    return float((h_N / 3.0) * np.dot(simps_w, integrand) * PErG(Er))

def consecutive(data, stepsize=1):
    return np.split(data, np.where(np.diff(data) != stepsize)[0]+1)

# Define Gaussian function
def gaussian(x, A, mu, sigma):
    return A * np.exp(-((x - mu) ** 2) / (2 * sigma ** 2))

def integrate_g_safe_inspect(g, Ep, Eq, Ermx, resolution):
    minx = 5e-21
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
        # epsabs=0 forces quad to meet the *relative* tolerance: these
        # integrands can be ~1e-20, far below the default epsabs=1.49e-8,
        # so with the default quad accepts its first coarse estimate over
        # the wide window without ever subdividing around the narrow peak
        ans = integrate.quad(g, minx, maxx, args = (Ep,Eq,), epsabs=0, limit=200)
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
              res=0.01):
    ppqGArr = []
    for this_Ep, this_Eq in zip(Ep, Eq):
        PpqGval, _, _ = PpqG_safe_inspect(this_Ep, this_Eq, F0=F0, s=s, eps=eps, V=V, p0=p0, p10=p10, q0=q0, q10=q10, res=res)
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
              res=0.01):
    if F0 == 0 and s == 0:
        raise ValueError("Fano factor F(Er) = F0 + s*Er is zero for all Er; F0 and s cannot both be zero")

    # Set parameters for integration
    # ER peaks can be as narrow as the zero-energy ionization
    # resolution q0 (~0.06 keV for test parameters), so the scan
    # grid must be finer than that to set accurate integration
    # limits; res=0.01 scans at 0.005 keV, finer than the 0.01 keV
    # grid the Fortran PpqG integrates on

    # Define the function (replace PpqFullG with your actual implementation)
    g = lambda er, ep, eq: PpqFullG(er, ep, eq, F0=F0, s=s, eps=eps, V=V, p0=p0, p10=p10, q0=q0, q10=q10)

    f = lambda er,ep,eq: -1*PpqExp(er,ep,eq,a=1.0,b=0.0,F0=F0,s=s,eps=eps,V=V,p0=p0,p10=p10,q0=q0,q10=q10)

    # WARNING
    # this will only work if your peaks have a width
    # greater than 0.1 keV
    # which seems to be true for all reasonable CDMS parameters
    maxx = max(Ep, Eq) + 5
    er_arr = np.arange(5e-21, maxx, res)
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

