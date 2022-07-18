from constants import *
import exoplanet as xo


def calc_aor(P, rho):
    '''
    Compute semimajor axis in units of R_star
    
    Args:
        P: period in units of days
        rho: stellar density in units of g/cc
    Out:
        aor: semimajor axis in units of R_star
    '''
    per = P * 86400.
    aor = ( (per**2 * G * rho) / (3*np.pi) )**(1/3)
    return aor


def calc_T14(P, rho, b, ror, e, w):
    '''
    T14 equation from Winn 2010 
    
    Args:
        P: period in units of days
        rho: stellar density in units of g/cc
        b: impact parameter
        ror: radius ratio
        ecc: eccentricity
        omega: argument of periastron in radians
    Out:
        T14: duration in units of days
    '''
    aor = calc_aor(P, rho)
    f_ew = (1 - e**2)**(1/2) / (1 + e*np.sin(w))

    per = P * 86400.
    T14 = per/np.pi * np.arcsin( ( (aor**2 - b**2) / ((1+ror)**2 - b**2) )**(-1/2) ) * f_ew 
    return T14 / 86400.


def calc_rho_star(P, T14, b, ror, ecc, omega):
    '''
    Inverting T14 equation from Winn 2010 
    (rather than using Kipping 2010a which introduced issues at high ecc)
    
    Args:
        P: period in units of days
        T14: duration in units of days
        b: impact parameter
        ror: radius ratio
        ecc: eccentricity
        omega: argument of periastron in radians
    Out:
        rho_star: stellar density in units of g/cc
    '''
    per = P * 86400.
    dur = T14 * 86400.
    rho_star = 3*np.pi / (per**2 * G) * (  (((1+ror)**2-b**2) / np.sin((np.pi*dur/per) * (1+ecc*np.sin(omega))/(1-ecc**2))**2) + b**2)**(3/2)
    return rho_star



def calc_noise(T14, N, ror, SNR, cadence=2):
    '''
    Find sigma_noise given an SNR, duration, radius, and N_transits
    Args:
        T14: duration in units of days
        N: # of transits
        ror: radius ratio
        SNR: desired signal-to-noise
        cadence: desired cadence in minutes (optional; default 2-minute)
    Out:
        noise: sigma_noise in normalized flux units
    '''
    
    texp = (60.*cadence)
    
    noise = ( (T14*86400.)/texp * N )**(1/2) * ror**2 / SNR
    return noise


def calc_snr(T14, N, ror, noise, cadence=2):
    '''
    Find SNR given a sigma_noise, duration, radius, and N_transits
    Args:
        T14: duration in units of days
        N: # of transits
        ror: radius ratio
        noise: desired sigma_noise in normalized flux units
        cadence: desired cadence in minutes (optional; default 2-minute)
    Out:
        SNR: signal-to-noise
    '''    
    texp = (60.*cadence)
    
    SNR = ( (T14*86400.)/texp * N )**(1/2) * ror**2 / noise
    return SNR


def simple_model(x, params, synth=False, oversample=7):
    '''
    Fit a light curve model across x values using exoplanet code by DFM
    Args:
        x: x values to evaluate light curve across
        params: dictionary containing info for constructing the light curve,
                i.e. RHOSTAR, PERIOD, T0, IMPACT, ROR, ECC, OMEGA, LD_U1, LD_U2, ROR
        synth: if you just want a clean synthetic curve to plot, regardless of exposure time;
               sets exposure time to a very small value (default = False)
        oversample: determines the assumed oversampling rate of the data (default = 7)

    Out:
        light_curve: flux values of the light curve model
    '''   
    texp = np.min(np.diff(x))
    if synth == True:
        texp /= 100

    _keys = list(params.keys())
    if 'ECC' in _keys and 'OMEGA' in _keys:
        orbit = xo.orbits.KeplerianOrbit(rho_star=params['RHOSTAR'], period=params['PERIOD'], 
            t0=params['T0'], b=params['IMPACT'], ror=params['ROR'],
            ecc=params['ECC'], omega=params['OMEGA'])
    elif 'DUR14' in _keys:
        orbit = xo.orbits.KeplerianOrbit(rho_star=params['RHOSTAR'], period=params['PERIOD'], 
            t0=params['T0'], b=params['IMPACT'], ror=params['ROR'],
            duration=params['DUR14'])
    else:
        print('WARNING: dictionary of inputs requires either ECC and OMEGA or DUR14!')
        return np.zeros(len(x))

    # Compute a limb-darkened light curve using starry
    light_curve = (
        xo.LimbDarkLightCurve(np.array([params['LD_U1'],params['LD_U2']]))
        .get_light_curve(orbit=orbit, r=params['ROR'], t=x, texp=texp, oversample=oversample)
        .eval()
    )

    return light_curve