import rasterio, math, numpy as np, os, glob
from rasterio import windows
from pyproj import Transformer
from scipy.optimize import brentq
from scipy.interpolate import interp1d
import xarray as xr
from ice_refractive_interpol import get_nk, get_alpha, alpha, MERSI2_BANDS


# Path and station

# geo_folder = r"D:\KAN_U\geo"
# lon,lat=-47.03,67 #KAN_U h=1840m

# geo_folder=r"D:\EGP\geo"
# lon,lat=-35.97,75.63 #EGP h=2660m

geo_folder = r"D:\QAS_L\geo"
lon,lat=-46.85,61.03 #QAS_L h=280m

# geo_folder=r"D:\SCO_U\MERSI2\geo"
# lon, lat = -27.23, 72.39 #SCO_U h=970m

# geo_folder = r"D:\DomeC\geo"
# lon,lat=123.3,-75.1 #Dome C h=3233m

# geo_folder = r"D:\DomeA\geo"
# lon,lat=77.124,-80.428 #Dome A h=4093m

window_size =1 #Window size is modifiable

transformer = Transformer.from_crs("EPSG:4326", "EPSG:3995", always_xy=True)
x, y = transformer.transform(lon, lat)


def read_window(tif_path, win):
    with rasterio.open(tif_path) as ds:
        return ds.read(1, window=win)

def safe_mean(arr):
    arr = np.where((arr == -999) | np.isnan(arr), np.nan, arr)
    return np.nanmean(arr)

def u(mu):
    return (3.0/5.0)*mu + (1.0 + np.sqrt(mu))/3.0

def rinff(cos_sza, cos_vza, theta):
    a = 1.247
    b = 1.186
    c = 5.157
    p = 11.1 * np.exp(-0.087 * theta) + 1.1 * np.exp(-0.014 * theta)
    return ((a + b * (cos_sza + cos_vza) + c * cos_sza * cos_vza + p) 
            / 4.0 / (cos_sza + cos_vza))


# Solar irradiance flux calculation

f0 = 32.38
f1 = -160140.33
f2 = 7959.53
bet = 11.72
gam = 2.49

def sol(x):
    sol1a = f0 * x
    sol1b = -f1 * np.exp(-bet * x) / bet
    sol1c = -f2 * np.exp(-gam * x) / gam
    return sol1a + sol1b + sol1c

sol_vis = sol(0.7) - sol(0.33)
sol_nir = sol(2.4) - sol(0.7)
sol_sw = sol_vis + sol_nir
asol = sol(0.865) - sol(0.7)

def analyt_func(z1, z2):
    ak1 = (z2**2 - z1**2)/2
    ak2 = (z2/bet + 1/bet**2)*np.exp(-bet*z2) - (z1/bet + 1/bet**2)*np.exp(-bet*z1)
    ak3 = (z2/gam + 1/gam**2)*np.exp(-gam*z2) - (z1/gam + 1/gam**2)*np.exp(-gam*z1)
    am1 = (z2**3 - z1**3)/3
    am2 = (z2**2/bet + 2*z2/bet**2 + 2/bet**3)*np.exp(-bet*z2) - (z1**2/bet + 2*z1/bet**2 + 2/bet**3)*np.exp(-bet*z1)
    am3 = (z2**2/gam + 2*z2/gam**2 + 2/gam**3)*np.exp(-gam*z2) - (z1**2/gam + 2*z1/gam**2 + 2/gam**3)*np.exp(-gam*z1)
    return (f0*ak1 - f1*ak2 - f2*ak3), (f0*am1 - f1*am2 - f2*am3)

def quad_func(x0, x1, x2, y0, y1, y2):
    d1 = (x0-x1)*(x0-x2)
    d2 = (x1-x0)*(x1-x2)
    d3 = (x2-x0)*(x2-x1)
    a1 = x1*x2*y0/d1 + x0*x2*y1/d2 + x0*x1*y2/d3
    b1 = -(x1+x2)*y0/d1 - (x0+x2)*y1/d2 - (x0+x1)*y2/d3
    c1 = y0/d1 + y1/d2 + y2/d3
    sa = a1 + b1*x1 + c1*x1*x1
    return sa, a1, b1, c1

def build_alb_vector(rs_412, rs_555, rs_709, rs_746, rs_865, rs_1030):
    alb = np.zeros(19)
    alb[7]  = rs_412   # 412nm
    alb[10] = rs_555   # 555nm
    alb[12] = rs_709   # 709nm
    alb[13] = rs_746   # 746nm
    alb[14] = rs_865   # 865nm
    alb[18] = rs_1030  # 1030nm
    return alb

def BBA_calc_pol(alb):
    alam2, alam3, alam5 = 0.412, 0.555, 0.709
    alam6, alam7, alam8 = 0.746, 0.865, 1.03
    r2, r3, r5 = alb[7], alb[10], alb[12]
    r6, r7, r8 = alb[13], alb[14], alb[18]
    
    # 0.33-0.7μm
    _, a1, b1, c1 = quad_func(alam2, alam3, alam5, r2, r3, r5)
    coef1, coef2 = analyt_func(0.33, 0.70)
    aj1 = a1*sol_vis + b1*coef1 + c1*coef2
    
    # 0.7-0.865μm
    _, a2, b2, c2 = quad_func(alam5, alam6, alam7, r5, r6, r7)
    coef3, coef4 = analyt_func(0.70, 0.865)
    aj2 = a2*asol + b2*coef3 + c2*coef4
    
    # 0.865-2.4μm
    rati = r7 / r8
    alasta = (alam8 - alam7) / np.log(rati)
    an = 1.0 / alasta
    p = r7 * np.exp(alam7 / alasta)
    z1, z2 = 0.865, 2.4
    aj31 = (1.0 / an) * (np.exp(-an * z2) - np.exp(-an * z1))
    aj32 = (1.0 / (bet + an)) * (np.exp(-(bet + an) * z2) - np.exp(-(an + bet) * z1))
    aj33 = (1.0 / (gam + an)) * (np.exp(-(gam + an) * z2) - np.exp(-(an + gam) * z1))
    aj3 = (-f0 * aj31 - f1 * aj32 - f2 * aj33) * p
    
    BBA_sw = (aj1 + aj2 + aj3) / sol_sw
    return BBA_sw

# Snow albedo retrieval 
def snow_albedo_iterative(TOA, angles, snow, atmosphere):
    rs_idx = [0, 1, 2, 3, 4, 5]
    toa_all = TOA['toa']
    Ta = atmosphere['Ta']
    ra = atmosphere['ra']
    Ra = atmosphere['Ra']
    R0 = snow['R0']
    u1 = angles['u1']
    u2 = angles['u2']
    nb, ny, nx = toa_all.shape
    rs = np.full_like(toa_all, np.nan)
    
    for ib in rs_idx:
        for j in range(ny):
            for i in range(nx):
                def root(rs_):
                    return (Ta[ib, j, i] * R0 * rs_**(u1[j, i]*u2[j, i]/R0[j, i])
                            / (1.0 - rs_ * ra[ib, j, i]) + Ra[ib, j, i] - toa_all[ib, j, i])
                
                try:
                    rs[ib, j, i] = brentq(root, 0.1, 1.5, xtol=1e-10, maxiter=100)
                except ValueError:
                    pass 
    
    # 
    rs = np.clip(rs, 0.0, 1.5)
    return TOA, snow, rs


# Atmospheric correction

def atmos_corr_mersi2(wls, sza, vza, vaa, saa, height=280, aot500=0.07, ang=1.3):  #Different sites have different elevations
    cos_sza = np.cos(np.deg2rad(sza))
    cos_vza = np.cos(np.deg2rad(vza))
    cos_raa = np.cos(np.deg2rad(vaa - saa))
    inv_cos_za = 1/cos_sza + 1/cos_vza
    cos_sa = -cos_sza*cos_vza + np.sin(np.deg2rad(sza))*np.sin(np.deg2rad(vza))*cos_raa
    
    tau_aer = aot500 * (wls/0.5)**(-ang)
    g_aer = 0.5263 + 0.4627*np.exp(-wls/0.4685)
    tau_mol = 0.008735 * (wls/1)**(-4.08) * np.exp(-height/6000)
    tau = tau_aer + tau_mol
    
    g11, g22 = 0.80, -0.45
    pa1 = (1-g11**2)/(1-2*g11*cos_sa+g11**2)**1.5
    pa2 = (1-g22**2)/(1-2*g22*cos_sa+g22**2)**1.5
    cp = (g_aer - g11)/(g11 - g22)
    pa = cp*pa1 + (1-cp)*pa2
    pr = 0.75*(1 + cos_sa**2)
    p = (tau_mol*pr + tau_aer*pa)/tau
    
    b1 = 1 + 1.5*cos_sza + (1-1.5*cos_sza)*np.exp(-tau/cos_sza)
    b2 = 1 + 1.5*cos_vza + (1-1.5*cos_vza)*np.exp(-tau/cos_vza)
    astra = (1 - np.exp(-tau*inv_cos_za))/4/(cos_sza+cos_vza)
    oskar = 4 + 3*(1-g_aer)*tau
    rms = 1 - b1*b2/oskar + (3*(1+g_aer)*cos_sza*cos_vza - 2*(cos_sza+cos_vza))*astra
    Ra = p*astra + rms
    
    tz = 1 + g_aer**2 + (1-g_aer**2)*np.sqrt(1+g_aer**2)
    Baer = 0.5 + g_aer*(g_aer**2-3)/tz/2
    Baer = np.where(g_aer>=0.001, (1-g_aer)*((1+g_aer)/np.sqrt(1+g_aer**2)-1)/2/g_aer, Baer)
    B = (0.5*tau_mol + Baer*tau_aer)/tau
    Ta = np.exp(-B*tau/cos_sza)*np.exp(-B*tau/cos_vza)
    
    gasa = 0.5772157
    y = (1+tau)*tau*np.exp(-tau)/4
    z3 = tau**2*(-np.log(tau)-gasa)
    z4 = tau**2*(tau - tau**2/4 + tau**3/18)
    Z = z3 + z4
    W1 = 1 + (1+0.5*tau)*Z/2 - y
    W2 = 1 + 0.75*tau*(1-g_aer)
    ra = 1 - W1/W2
    
    return {'Ta':Ta, 'ra':ra, 'Ra':Ra}


# Main loop

B15_files = sorted(glob.glob(os.path.join(geo_folder, "*B15_TOA_1km.tif")))

# 
R08_thv0 = 0.75  # 412nm TOA threshold

for B15_path in B15_files:
    base = os.path.basename(B15_path)
    parts = base.split("_")
    timestamp = f"{parts[4]}_{parts[5]}"
    geo_tag = "_".join(parts[0:6])
    
    # 
    B08_path = B15_path.replace("B15", "B08")
    B11_path = B15_path.replace("B15", "B11")
    B13_path = B15_path.replace("B15", "B13")
    B14_path = B15_path.replace("B15", "B14")
    B19_path = B15_path.replace("B15", "B19")
    SZA_path = f"{geo_folder}\\{geo_tag}_GEO1K_MS_SZA_1km.tif"
    VZA_path = f"{geo_folder}\\{geo_tag}_GEO1K_MS_VZA_1km.tif"
    SAA_path = f"{geo_folder}\\{geo_tag}_GEO1K_MS_SAA_1km.tif"
    VAA_path = f"{geo_folder}\\{geo_tag}_GEO1K_MS_VAA_1km.tif"
    RAA_path = f"{geo_folder}\\{geo_tag}_GEO1K_MS_RAA_1km.tif"
    
    with rasterio.open(B15_path) as src:
        row, col = src.index(x, y)
        win = windows.Window(col, row, window_size, window_size)
    
    # 
    R08 = safe_mean(read_window(B08_path, win))   # 412nm
    R11 = safe_mean(read_window(B11_path, win))   # 555nm
    R13 = safe_mean(read_window(B13_path, win))   # 709nm
    R14 = safe_mean(read_window(B14_path, win))   # 746nm
    R15 = safe_mean(read_window(B15_path, win))   # 865nm
    R19 = safe_mean(read_window(B19_path, win))   # 1030nm
    SZA = safe_mean(read_window(SZA_path, win))
    VZA = safe_mean(read_window(VZA_path, win))
    SAA = safe_mean(read_window(SAA_path, win))
    VAA = safe_mean(read_window(VAA_path, win))
    RAA = safe_mean(read_window(RAA_path, win))
    
    
    if np.isnan(R08) or np.isnan(R19) or R08 < 0.2 or R19 < 0.1 or SZA > 75:
        continue
    
    mu0 = math.cos(math.radians(SZA))
    mu = math.cos(math.radians(VZA))
    
    # Calculate angular parameters
    u1 = u(mu0)
    u2 = u(mu)
    theta = np.degrees(np.arccos(-mu0 * mu + math.sin(math.radians(SZA))*math.sin(math.radians(VZA))*math.cos(math.radians(RAA))))
    
    
    # rinff
    rinff_i = rinff(mu0, mu, theta)
    
    # R08 < threshold
    msk = (R08 < R08_thv0) 
    
    #  calculate factor
    if msk:
        factor = R08 / rinff_i 
    else:
        factor =1.0 
    
    # 
    TOA = np.array([R08, R11, R13, R14, R15, R19])
    if msk:
        TOA_corr = TOA / factor 
    else:
        TOA_corr = TOA
    
    #R0
    alpha1 = alpha(865)
    alpha2 = alpha(1030)
    b = math.sqrt(alpha1/alpha2)
    beta = 1/(1-b)
    
    R15_proc = TOA_corr[4]  # 865nm
    R19_proc = TOA_corr[5]  # 1030nm
    R0 = (R15_proc**beta) * (R19_proc**(1-beta))
    # print(R15_proc,R19_proc,R0,L)
    
    # L
    W = 1/alpha2
    kesei = u1 * u2 / R0
    L = W * (math.log(R19_proc/R0)**2)* (kesei**(-2)) * 10

    
    # 
    wls = np.array([0.412, 0.555, 0.709, 0.746, 0.865, 1.030])
    atm = atmos_corr_mersi2(wls, SZA, VZA, VAA, SAA)
    
    # 
    angles_dict = {'u1': np.array([[u1]]), 'u2': np.array([[u2]])}
    snow_dict = {'R0': np.array([[R0]])}
    atmos_dict = {
        'Ta': atm['Ta'].reshape(6, 1, 1),
        'ra': atm['ra'].reshape(6, 1, 1),
        'Ra': atm['Ra'].reshape(6, 1, 1)
    }
    TOA_dict = {'toa': TOA_corr.reshape(6, 1, 1)}
    
    _, _, rs_6band = snow_albedo_iterative(TOA_dict, angles_dict, snow_dict, atmos_dict)
    rs_6band = rs_6band[:, 0, 0]
    
    # 
    if np.any(np.isnan(rs_6band)):
        print(f"{timestamp} 反演失败(NaN)，跳过")
        continue    

    rs_6band = np.clip(rs_6band, 0.0, 1.4)    
    alb_sph_412 = factor * rs_6band[0]
    
    # classify    
    if alb_sph_412 >0.98:
        isnow = 1  # clean snow
    elif factor > 0.99:
        isnow = 2  # snow with impurities
    else:
        isnow = 3  # patchy snow
    
    #  BBA
    rp_pure = rs_6band ** u1

    if isnow == 1:
        BBA =0.5271 + 0.3612 * math.exp(-u1 * math.sqrt(0.0235 * L))

    else:
            alb_rp = build_alb_vector(
        rp_pure[0],
        rp_pure[1],
        rp_pure[2],
        rp_pure[3],
        rp_pure[4],
        rp_pure[5]
    )

            BBA = factor*BBA_calc_pol(alb_rp)


    print(f"==> {timestamp}  isnow={isnow}  BBA={BBA:.4f}   "
          f"factor={factor:.4f}  alb_sph_412={alb_sph_412:.4f}  "
          f"msk={msk}  R08={R08:.3f}")
    
    # save             
    import pandas as pd
    import os
    excel_path = os.path.join(os.path.dirname(geo_folder), 'BBA.xlsx')
    new_row = pd.DataFrame({'timestamp': [timestamp], 'BBA': [BBA]})

    if os.path.isfile(excel_path):
        old = pd.read_excel(excel_path)
        updated = pd.concat([old, new_row], ignore_index=True)
    else:
        updated = new_row

    updated.to_excel(excel_path, index=False)


