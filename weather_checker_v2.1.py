“””
Weather Risk Checker V2.1
Severe Weather Sounding Analyzer — WFO Louisville Tip Sheet Integration
Uses IEM model soundings (GFS, NAM, HRRR, RAP)

Parameters:
MLCAPE, MLCIN, MLLCL, LFC, EL
SRH 0-1 / 0-3 km (Bunkers RM/LM)
Bulk shear 0-1 / 0-3 / 0-6 / 0-8 km
Effective inflow layer, Eff SRH, Eff Bulk Shear
BRN, BRN Shear
Lapse rates 0-3 km / 700-500 mb
Boundary layer RH, Precipitable Water
Hodograph shape (straight/curved + tornadic spike signature)
DCAPE (400-500 hPa)
SCP, STP (effective layer)
Storm mode classifier (tornadic supercell / supercell / bow / MCS / multicell / non-meso)
“””

import tkinter as tk
from tkinter import ttk, scrolledtext
import requests
import math
import threading
from datetime import datetime, timezone

# ── Physical constants ────────────────────────────────────────────────────────

Rd = 287.05
Cp = 1005.7
g  = 9.81
Lv = 2.501e6

# ── Unit helpers ──────────────────────────────────────────────────────────────

def c_to_k(c):      return c + 273.15
def k_to_c(k):      return k - 273.15

# ══════════════════════════════════════════════════════════════════════════════

# THERMODYNAMICS

# ══════════════════════════════════════════════════════════════════════════════

def sat_vapor_pressure(T_c):
“”“Saturation vapor pressure (hPa), Bolton 1980.”””
return 6.112 * math.exp(17.67 * T_c / (T_c + 243.5))

def mixing_ratio_from_dewpoint(Td_c, pres_hpa):
“”“Saturation mixing ratio (kg/kg) from dewpoint and pressure.”””
e = sat_vapor_pressure(Td_c)
return 0.622 * e / max(pres_hpa - e, 0.001)

def dewpoint_from_mixing_ratio(w_kgkg, pres_hpa):
“”“Dewpoint (°C) from mixing ratio (kg/kg) and pressure (hPa).”””
e = w_kgkg * pres_hpa / (0.622 + w_kgkg)
e = max(e, 1e-6)
return (243.5 * math.log(e / 6.112)) / (17.67 - math.log(e / 6.112))

def theta_e(T_c, Td_c, pres_hpa):
“”“Equivalent potential temperature (K), Bolton 1980.”””
T_k  = c_to_k(T_c)
Td_k = c_to_k(Td_c)
e    = sat_vapor_pressure(Td_c)
w    = 0.622 * e / max(pres_hpa - e, 0.001)
TL   = 56 + 1.0 / (1.0 / (Td_k - 56) + math.log(T_k / Td_k) / 800.0)
return (T_k * (1000.0 / pres_hpa) ** (0.2854 * (1 - 0.28 * w))
* math.exp((3376.0 / TL - 2.54) * w * (1 + 0.81 * w)))

def lcl_temperature(T_c, Td_c):
“”“LCL temperature (K), Bolton 1980.”””
T_k  = c_to_k(T_c)
Td_k = c_to_k(Td_c)
return 1.0 / (1.0 / (Td_k - 56) + math.log(T_k / Td_k) / 800.0) + 56

def lcl_pressure(T_c, Td_c, pres_hpa):
“”“LCL pressure (hPa).”””
T_lcl = lcl_temperature(T_c, Td_c)
return pres_hpa * (T_lcl / c_to_k(T_c)) ** (Cp / Rd)

def moist_adiabat_temp(te_K, pres_hpa, T_first_guess=None):
“”“Temperature on moist adiabat at pres_hpa conserving theta-e (Newton iteration).”””
if T_first_guess is None:
T_first_guess = k_to_c(te_K * (pres_hpa / 1000.0) ** (Rd / Cp))
T = T_first_guess
for _ in range(100):
ws  = mixing_ratio_from_dewpoint(T, pres_hpa)
Td  = dewpoint_from_mixing_ratio(ws, pres_hpa)
te  = theta_e(T, Td, pres_hpa)
dT  = (te_K - te) * 0.25
T  += dT
if abs(dT) < 0.001:
break
return T

def pressure_to_height_interp(levels, target_pres):
“”“Interpolate height AGL (m) for a target pressure (hPa).”””
for i in range(1, len(levels)):
p0 = levels[i-1][“pres”]
p1 = levels[i][“pres”]
if p1 <= target_pres <= p0:
frac = (p0 - target_pres) / max(p0 - p1, 0.001)
return levels[i-1][“hght”] + frac * (levels[i][“hght”] - levels[i-1][“hght”])
return 0.0

def height_to_temp_interp(levels, target_hght):
“”“Interpolate temperature (°C) at a target height AGL (m).”””
for i in range(1, len(levels)):
h0 = levels[i-1][“hght”]
h1 = levels[i][“hght”]
if h0 <= target_hght <= h1:
frac = (target_hght - h0) / max(h1 - h0, 0.001)
return levels[i-1][“temp”] + frac * (levels[i][“temp”] - levels[i-1][“temp”])
return None

# ══════════════════════════════════════════════════════════════════════════════

# MIXED-LAYER PARCEL

# ══════════════════════════════════════════════════════════════════════════════

def mixed_layer_parcel(levels):
“””
Average theta-e and mixing ratio over lowest 100 hPa.
Returns (ml_temp_c, ml_dwpt_c, ml_pres_hpa).
“””
sfc_pres = levels[0][“pres”]
top_pres = sfc_pres - 100.0

```
te_sum = w_sum = 0.0
count  = 0

for lv in levels:
    if lv["pres"] < top_pres:
        break
    te_sum += theta_e(lv["temp"], lv["dwpt"], lv["pres"])
    w_sum  += mixing_ratio_from_dewpoint(lv["dwpt"], lv["pres"])
    count  += 1

if count == 0:
    return levels[0]["temp"], levels[0]["dwpt"], levels[0]["pres"]

avg_te = te_sum / count
avg_w  = w_sum  / count
ml_p   = sfc_pres - 50.0

T = levels[0]["temp"]
for _ in range(60):
    Td  = dewpoint_from_mixing_ratio(avg_w, ml_p)
    te  = theta_e(T, Td, ml_p)
    dT  = (avg_te - te) * 0.1
    T  += dT
    if abs(dT) < 0.001:
        break

Td_ml = dewpoint_from_mixing_ratio(avg_w, ml_p)
return T, Td_ml, ml_p
```

# ══════════════════════════════════════════════════════════════════════════════

# CAPE / CIN

# ══════════════════════════════════════════════════════════════════════════════

def cape_cin(levels, parcel_temp, parcel_dwpt, parcel_pres):
“””
Integrate CAPE and CIN from parcel trace.
Returns (cape, cin, lcl_hght_m, lfc_pres, el_pres).
“””
p_lcl = lcl_pressure(parcel_temp, parcel_dwpt, parcel_pres)
te_p  = theta_e(parcel_temp, parcel_dwpt, parcel_pres)

```
trace      = []
lfc_pres   = None
el_pres    = None
prev_buoy  = None

for lv in levels:
    p = lv["pres"]
    if p > parcel_pres:
        continue

    if p >= p_lcl:
        T_parcel = k_to_c(c_to_k(parcel_temp) * (p / parcel_pres) ** (Rd / Cp))
    else:
        T_parcel = moist_adiabat_temp(te_p, p)

    buoy = (c_to_k(T_parcel) - c_to_k(lv["temp"])) / c_to_k(lv["temp"])

    if p <= p_lcl and prev_buoy is not None:
        if prev_buoy <= 0 < buoy and lfc_pres is None:
            lfc_pres = p
        if prev_buoy > 0 >= buoy and lfc_pres is not None:
            el_pres = p

    trace.append((p, buoy, lv["hght"]))
    prev_buoy = buoy

cape = cin = 0.0
for i in range(1, len(trace)):
    p1, b1, h1 = trace[i-1]
    p2, b2, h2 = trace[i]
    dz    = abs(h2 - h1)
    avg_b = (b1 + b2) / 2.0
    energy = g * avg_b * dz

    above_lfc = lfc_pres and p1 <= lfc_pres
    below_el  = el_pres  and p1 >= el_pres

    if above_lfc and below_el and energy > 0:
        cape += energy
    elif lfc_pres and p1 > lfc_pres and energy < 0:
        cin += energy

lcl_hght = pressure_to_height_interp(levels, p_lcl)
return max(cape, 0), min(cin, 0), lcl_hght, lfc_pres, el_pres
```

# ══════════════════════════════════════════════════════════════════════════════

# WIND / SHEAR HELPERS

# ══════════════════════════════════════════════════════════════════════════════

def wind_to_uv(spd_kt, dir_deg):
r = math.radians(dir_deg)
return -spd_kt * math.sin(r), -spd_kt * math.cos(r)

def uv_to_dir(u, v):
return (270 - math.degrees(math.atan2(v, u))) % 360

def uv_to_spd(u, v):
return math.sqrt(u**2 + v**2)

def layer_wind_interp(levels, target_m):
“”“Interpolate (u, v) in knots at target height AGL.”””
for i in range(1, len(levels)):
h0, h1 = levels[i-1][“hght”], levels[i][“hght”]
if h0 <= target_m <= h1:
frac = (target_m - h0) / max(h1 - h0, 0.001)
u0, v0 = wind_to_uv(levels[i-1][“wspd”], levels[i-1][“wdir”])
u1, v1 = wind_to_uv(levels[i][“wspd”],   levels[i][“wdir”])
return u0 + frac*(u1-u0), v0 + frac*(v1-v0)
# Extrapolate with top level if above sounding
u, v = wind_to_uv(levels[-1][“wspd”], levels[-1][“wdir”])
return u, v

def bulk_shear(levels, top_m):
“”“Bulk shear magnitude (kt) from surface to top_m AGL.”””
u_sfc, v_sfc = wind_to_uv(levels[0][“wspd”], levels[0][“wdir”])
u_top, v_top = layer_wind_interp(levels, top_m)
return uv_to_spd(u_top - u_sfc, v_top - v_sfc)

# ══════════════════════════════════════════════════════════════════════════════

# BUNKERS STORM MOTION

# ══════════════════════════════════════════════════════════════════════════════

def bunkers_storm_motion(levels):
“””
Bunkers (2000) internal dynamics method.
Returns (rm_u, rm_v, lm_u, lm_v, mw_u, mw_v) in knots.
“””
u_sum = v_sum = 0.0
count = 0
for lv in levels:
if lv[“hght”] > 6000:
break
u, v = wind_to_uv(lv[“wspd”], lv[“wdir”])
u_sum += u; v_sum += v; count += 1

```
if count == 0:
    return (0,) * 6

mw_u = u_sum / count
mw_v = v_sum / count

u_lo, v_lo = layer_wind_interp(levels, 500)
u_hi, v_hi = layer_wind_interp(levels, 5500)

su = u_hi - u_lo
sv = v_hi - v_lo
smag = max(uv_to_spd(su, sv), 1e-6)
D = 7.5  # kt

perp_u =  sv / smag * D
perp_v = -su / smag * D

return (mw_u + perp_u, mw_v + perp_v,
        mw_u - perp_u, mw_v - perp_v,
        mw_u, mw_v)
```

# ══════════════════════════════════════════════════════════════════════════════

# STORM-RELATIVE HELICITY

# ══════════════════════════════════════════════════════════════════════════════

def compute_srh(levels, storm_u, storm_v, top_m):
“”“SRH (m²/s²) from surface to top_m AGL, Davies-Jones 1990.”””
srh   = 0.0
layer = [lv for lv in levels if lv[“hght”] <= top_m]
for i in range(1, len(layer)):
u0, v0 = wind_to_uv(layer[i-1][“wspd”], layer[i-1][“wdir”])
u1, v1 = wind_to_uv(layer[i][“wspd”],   layer[i][“wdir”])
sr_u0, sr_v0 = u0 - storm_u, v0 - storm_v
sr_u1, sr_v1 = u1 - storm_u, v1 - storm_v
srh += sr_u0 * sr_v1 - sr_v0 * sr_u1
return srh

# ══════════════════════════════════════════════════════════════════════════════

# EFFECTIVE INFLOW LAYER

# ══════════════════════════════════════════════════════════════════════════════

def effective_inflow_layer(levels):
“””
Effective inflow layer: CAPE >= 100 J/kg and CIN >= -250 J/kg.
Returns (base_idx, top_idx) of first contiguous layer, or (None, None).
“””
base_idx = top_idx = None
for i, lv in enumerate(levels):
c, ci, _, _, _ = cape_cin(levels, lv[“temp”], lv[“dwpt”], lv[“pres”])
if c >= 100 and ci >= -250:
if base_idx is None:
base_idx = i
top_idx = i
else:
if base_idx is not None:
break
return base_idx, top_idx

def effective_bulk_shear(levels, base_idx, top_idx):
“”“Effective bulk shear (kt) through the effective inflow layer.”””
if base_idx is None or top_idx is None:
return 0.0
base_lv = levels[base_idx]
top_lv  = levels[top_idx]
cap_h   = min(top_lv[“hght”] + (top_lv[“hght”] - base_lv[“hght”]) / 2.0, 6000)
u_base, v_base = wind_to_uv(base_lv[“wspd”], base_lv[“wdir”])
u_top,  v_top  = layer_wind_interp(levels, cap_h)
return uv_to_spd(u_top - u_base, v_top - v_base)

# ══════════════════════════════════════════════════════════════════════════════

# BULK RICHARDSON NUMBER

# ══════════════════════════════════════════════════════════════════════════════

def compute_brn(mlcape, levels):
“””
Bulk Richardson Number = MLCAPE / (0.5 * U_mean²)
U_mean = magnitude of mean wind vector over 0-6 km minus 0-500 m mean.
BRN Shear = 0.5 * U_mean² (m²/s²), converted to kt² for display.
“””
if mlcape <= 0:
return None, 0.0

```
# Mean wind 0-6 km
u6 = v6 = 0.0; n6 = 0
for lv in levels:
    if lv["hght"] > 6000:
        break
    u, v = wind_to_uv(lv["wspd"], lv["wdir"])
    u6 += u; v6 += v; n6 += 1

# Mean wind 0-500 m
u5 = v5 = 0.0; n5 = 0
for lv in levels:
    if lv["hght"] > 500:
        break
    u, v = wind_to_uv(lv["wspd"], lv["wdir"])
    u5 += u; v5 += v; n5 += 1

if n6 == 0 or n5 == 0:
    return None, 0.0

du = u6/n6 - u5/n5
dv = v6/n6 - v5/n5
# Convert kt → m/s for energy calculation
U_ms = uv_to_spd(du, dv) * 0.514444
brn_shear = 0.5 * U_ms ** 2   # J/kg
brn = mlcape / brn_shear if brn_shear > 0 else None
return brn, brn_shear
```

# ══════════════════════════════════════════════════════════════════════════════

# LAPSE RATES

# ══════════════════════════════════════════════════════════════════════════════

def lapse_rate_layer(levels, bot_hght_m, top_hght_m):
“””
Environmental lapse rate (°C/km) between two heights AGL.
Positive = unstable.
“””
T_bot = height_to_temp_interp(levels, bot_hght_m)
T_top = height_to_temp_interp(levels, top_hght_m)
if T_bot is None or T_top is None:
return None
dz_km = (top_hght_m - bot_hght_m) / 1000.0
if dz_km <= 0:
return None
return (T_bot - T_top) / dz_km   # positive = temp decreasing with height

def lapse_rate_pressure(levels, bot_hpa, top_hpa):
“”“Lapse rate (°C/km) between two pressure levels.”””
h_bot = pressure_to_height_interp(levels, bot_hpa)
h_top = pressure_to_height_interp(levels, top_hpa)
T_bot = height_to_temp_interp(levels, h_bot)
T_top = height_to_temp_interp(levels, h_top)
if T_bot is None or T_top is None:
return None
dz_km = (h_top - h_bot) / 1000.0
if dz_km <= 0:
return None
return (T_bot - T_top) / dz_km

# ══════════════════════════════════════════════════════════════════════════════

# BOUNDARY LAYER RH AND PRECIPITABLE WATER

# ══════════════════════════════════════════════════════════════════════════════

def boundary_layer_rh(levels):
“”“Mean relative humidity (%) in the lowest 100 hPa (boundary layer).”””
sfc_p    = levels[0][“pres”]
top_p    = sfc_p - 100.0
rh_sum   = 0.0
count    = 0
for lv in levels:
if lv[“pres”] < top_p:
break
e_s = sat_vapor_pressure(lv[“temp”])
e   = sat_vapor_pressure(lv[“dwpt”])
rh  = 100.0 * e / max(e_s, 0.001)
rh_sum += rh
count  += 1
return rh_sum / count if count > 0 else 0.0

def precipitable_water(levels):
“””
Precipitable water (inches) integrated through the sounding.
Uses the virtual temperature method.
“””
pw = 0.0
for i in range(1, len(levels)):
p0 = levels[i-1][“pres”]
p1 = levels[i][“pres”]
w0 = mixing_ratio_from_dewpoint(levels[i-1][“dwpt”], p0)
w1 = mixing_ratio_from_dewpoint(levels[i][“dwpt”],   p1)
dp = (p0 - p1) * 100.0   # Pa
pw += 0.5 * (w0 + w1) * dp / g
return pw / 25.4   # kg/m² ≈ mm → inches

# ══════════════════════════════════════════════════════════════════════════════

# HODOGRAPH SHAPE

# ══════════════════════════════════════════════════════════════════════════════

def hodograph_shape(levels):
“””
Classify hodograph shape.

```
Returns dict with:
  shape        : "Straight" | "Curved-CW" | "Curved-CCW" | "Backed"
  tornadic_sig : True/False — 0-1 km spike + curvature above
  notes        : list of descriptive strings
"""
# Sample hodograph at key levels
pts = {}
for h in [0, 500, 1000, 2000, 3000, 6000]:
    u, v = layer_wind_interp(levels, h)
    pts[h] = (u, v)

def cross(a, b, c):
    """Cross product of vectors AB and AC. Positive = CCW turn."""
    return (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])

# ── Curvature 0-3 km (three points: sfc, 1.5 km, 3 km) ──────────────────
p0 = pts[0]
p1 = layer_wind_interp(levels, 1500)
p2 = pts[3000]
turn_03 = cross(p0, p1, p2)   # positive=CCW(LM favored), negative=CW(RM favored)

# ── Curvature 1-3 km (above the spike) ───────────────────────────────────
p1km = pts[1000]
p2km = pts[2000]
p3km = pts[3000]
turn_13 = cross(p1km, p2km, p3km)

# ── 0-1 km straightness (spike) ──────────────────────────────────────────
u0, v0 = pts[0]
u1, v1 = pts[1000]
# Deviation of 500 m point from straight line 0→1 km
um, vm = layer_wind_interp(levels, 500)
# Vector projection of midpoint onto 0-1 km line
seg_u = u1 - u0
seg_v = v1 - v0
seg_len = max(uv_to_spd(seg_u, seg_v), 0.001)
dev_u = um - u0 - seg_u * 0.5
dev_v = vm - v0 - seg_v * 0.5
spike_dev = uv_to_spd(dev_u, dev_v) / seg_len   # normalized deviation
spike_straight = spike_dev < 0.15   # < 15% deviation = spike is straight

# ── Shear vector direction (backed / veered) ──────────────────────────────
u6, v6 = pts[6000]
shear_dir = uv_to_dir(u6 - u0, v6 - v0)
# Backed low-level wind: surface wind from SE/S/SW backing relative to upper flow
sfc_dir = levels[0]["wdir"]
backed = (sfc_dir < 180) and (shear_dir > 180)

# ── Classify overall shape ────────────────────────────────────────────────
abs_turn = abs(turn_03)
if abs_turn < 30:
    shape = "Straight"
elif turn_03 < 0:
    shape = "Curved-CW"    # clockwise → RM supercell favored
else:
    shape = "Curved-CCW"   # counter-clockwise → LM supercell favored

if backed and shape == "Straight":
    shape = "Backed"

# ── Tornadic signature: straight 0-1 km spike + CW curvature 1-3 km ──────
tornadic_sig = spike_straight and (turn_13 < -20)

notes = []
if shape == "Curved-CW":
    notes.append("Clockwise hodograph curvature → right-mover supercell favored")
elif shape == "Curved-CCW":
    notes.append("Counter-clockwise curvature → left-mover or splitting storm possible")
elif shape == "Straight":
    notes.append("Straight hodograph → storm splitting likely; both movers possible")
elif shape == "Backed":
    notes.append("Backed surface winds with straight shear vector")

if tornadic_sig:
    notes.append("Tornadic hodograph signature: straight 0-1 km spike with curvature above")
if spike_straight and not tornadic_sig:
    notes.append("Straight 0-1 km hodograph spike present")

return {"shape": shape, "tornadic_sig": tornadic_sig, "notes": notes}
```

# ══════════════════════════════════════════════════════════════════════════════

# DCAPE

# ══════════════════════════════════════════════════════════════════════════════

def compute_dcape(levels):
“”“DCAPE (J/kg) from 400-500 hPa layer, most moist parcel.”””
layer = [lv for lv in levels if 400 <= lv[“pres”] <= 500]
if not layer:
return 0.0
parcel = min(layer, key=lambda lv: lv[“temp”] - lv[“dwpt”])

```
T_start  = parcel["temp"]
Td_start = parcel["dwpt"]
p_start  = parcel["pres"]
p_lcl    = lcl_pressure(T_start, Td_start, p_start)
te       = theta_e(T_start, Td_start, p_start)
sfc_p    = levels[0]["pres"]

dcape    = 0.0
prev_lv  = None

for lv in reversed(levels):
    if lv["pres"] < p_start or lv["pres"] > sfc_p:
        continue
    p = lv["pres"]
    if p >= p_lcl:
        T_parcel = moist_adiabat_temp(te, p, T_start)
    else:
        T_parcel = k_to_c(c_to_k(T_start) * (p / p_start) ** (Rd / Cp))

    buoy = g * (c_to_k(T_parcel) - c_to_k(lv["temp"])) / c_to_k(lv["temp"])
    if prev_lv is not None and buoy < 0:
        dz = abs(lv["hght"] - prev_lv["hght"])
        dcape += abs(buoy) * dz
    prev_lv = lv

return dcape
```

# ══════════════════════════════════════════════════════════════════════════════

# COMPOSITE INDICES

# ══════════════════════════════════════════════════════════════════════════════

def supercell_composite(mlcape, eff_srh, eff_shear_kt):
“”“SCP — Thompson 2004 effective layer.”””
if eff_shear_kt < 10:
return 0.0
return (mlcape / 1000.0) * (eff_srh / 50.0) * (min(eff_shear_kt, 60) / 40.0)

def significant_tornado(mlcape, mllcl_m, mlcin, eff_srh, eff_shear_kt):
“”“STP — effective layer, Thompson 2004/2012.”””
if eff_shear_kt < 10 or mlcape < 100:
return 0.0
lcl_t = 1.0 if mllcl_m <= 1000 else max((2000 - mllcl_m) / 1000.0, 0)
cin_t = 1.0 if mlcin >= -50 else max((200 + mlcin) / 150.0, 0)
return ((mlcape / 1500.0) * lcl_t * (eff_srh / 150.0)
* (min(eff_shear_kt, 60) / 45.0) * cin_t)

# ══════════════════════════════════════════════════════════════════════════════

# STORM MODE CLASSIFIER  (WFO Louisville Tip Sheet)

# ══════════════════════════════════════════════════════════════════════════════

def classify_storm_mode(
mlcape, mlcin, mllcl_m,
shear_01, shear_03, shear_06, shear_08,
srh_01, srh_03,
scp, stp,
brn, brn_shear,
lr_03, lr_700_500,
bl_rh, pw_in,
hodo,
dcape
):
“””
Returns (primary_mode, modes_list, flags, issues).
All thresholds from WFO Louisville Tip Sheet.
“””
modes  = []
flags  = []
issues = []

```
# ── Tornadic Supercell ────────────────────────────────────────────────────
tornadic_score = 0
if shear_06 >= 40:   tornadic_score += 1
if srh_01   >= 150:  tornadic_score += 2
if srh_03   >= 300:  tornadic_score += 1
if mllcl_m  <= 1000: tornadic_score += 2
if bl_rh    >= 65:   tornadic_score += 1
if shear_01 >= 20:   tornadic_score += 1
if hodo and hodo["tornadic_sig"]: tornadic_score += 2
if stp >= 1:          tornadic_score += 2
if stp >= 4:          tornadic_score += 2

if tornadic_score >= 7 and scp >= 2:
    modes.append("TORNADIC SUPERCELL")
    if stp >= 4:
        flags.append("⚠ SIGNIFICANT TORNADO THREAT — STP {:.2f}".format(stp))
    elif stp >= 2:
        flags.append("⚠ Elevated Tornado Threat — STP {:.2f}".format(stp))
    else:
        flags.append("Tornado Possible — STP {:.2f}".format(stp))
    if hodo and hodo["tornadic_sig"]:
        flags.append("⚠ Tornadic hodograph signature detected")
    if srh_01 >= 300:
        flags.append("⚠ Extreme 0-1 km SRH ({:.0f} m²/s²)".format(srh_01))
    elif srh_01 >= 150:
        flags.append("Significant 0-1 km SRH ({:.0f} m²/s²)".format(srh_01))

# ── Long-track tornado indicators ────────────────────────────────────────
if shear_01 >= 30 and srh_01 >= 200 and shear_08 >= 60:
    flags.append("⚠ Long-track tornado environment (0-1 shear ≥30 kt, 0-1 SRH ≥200, 0-8 shear ≥60 kt)")

# ── Supercell (non-tornadic or marginal) ─────────────────────────────────
if shear_06 >= 40 and scp >= 2:
    if "TORNADIC SUPERCELL" not in modes:
        modes.append("SUPERCELL")
    flags.append("Supercell Composite ≥ 2 (SCP {:.1f})".format(scp))
elif shear_06 >= 40 and mlcape >= 2500:
    # Marginal SRH but very high CAPE — still possible per tip sheet
    if "SUPERCELL" not in modes and "TORNADIC SUPERCELL" not in modes:
        modes.append("SUPERCELL (marginal shear, high CAPE)")
    flags.append("High CAPE compensating marginal shear — supercell still possible")

# ── BRN / supercell mode discrimination ──────────────────────────────────
if brn is not None:
    if 10 <= brn <= 50:
        flags.append("BRN {:.0f} — supercell range (10-50); higher end if boundary present".format(brn))
    elif brn > 50:
        flags.append("BRN {:.0f} — multicell storms more likely".format(brn))
    if brn_shear >= 80:
        flags.append("BRN Shear {:.0f} m²/s² — long-lived supercell potential".format(brn_shear))
    elif brn_shear >= 40:
        flags.append("BRN Shear {:.0f} m²/s² — supercell supportive".format(brn_shear))

# ── Bow Echo / Derecho ────────────────────────────────────────────────────
bow_score = 0
if shear_03 >= 30: bow_score += 2
if shear_03 >= 40: bow_score += 1
if mlcape   >= 2000: bow_score += 1
if lr_03 and lr_03 >= 7.0: bow_score += 1   # steep low-level lapse rates

if bow_score >= 3 and shear_06 < 50:   # bow echoes favor lower deep-layer shear
    modes.append("BOW ECHO / QLCS")
    flags.append("0-3 km shear {:.0f} kt — bow echo / bowing segment possible".format(shear_03))

# Derecho — high instability + fast unidirectional flow
if mlcape >= 2500 and shear_03 >= 30 and shear_06 >= 40:
    if "BOW ECHO / QLCS" in modes:
        modes[modes.index("BOW ECHO / QLCS")] = "BOW ECHO / DERECHO"
    flags.append("⚠ Derecho-favorable: high CAPE + 0-3 km shear ≥ 30 kt + 0-6 km shear ≥ 40 kt")

# ── MCS / Multicell ───────────────────────────────────────────────────────
if 20 <= shear_06 <= 35:
    modes.append("ORGANIZED MULTICELL / MCS")
    flags.append("0-6 km shear {:.0f} kt — organized multicell / MCS range".format(shear_06))
elif shear_06 < 20:
    modes.append("DISORGANIZED MULTICELL / PULSE")
    issues.append("0-6 km shear {:.0f} kt — disorganized convection likely".format(shear_06))

# MCS maintenance criteria
if mlcape >= 1700 and lr_700_500 and lr_700_500 >= 6.5 and shear_06 >= 40:
    flags.append("MCS maintenance criteria met (MUCAPE, lapse rates, shear)")

# ── Large Hail ────────────────────────────────────────────────────────────
if mlcape >= 2500 and lr_700_500 and lr_700_500 >= 7.0 and shear_06 >= 40:
    flags.append("⚠ Large hail potential: high CAPE + steep 700-500 mb lapse rates + shear")
elif mlcape >= 1500 and lr_700_500 and lr_700_500 >= 6.5:
    flags.append("Hail growth zone favorable (CAPE + lapse rates)")

# ── Mini-supercell ────────────────────────────────────────────────────────
if mlcape <= 1000 and srh_01 >= 100 and shear_06 >= 25:
    modes.append("MINI-SUPERCELL POSSIBLE")
    flags.append("Low CAPE / strong low-level shear — mini-supercell environment")

# ── Non-mesocyclone tornado ────────────────────────────────────────────────
if lr_03 and lr_03 >= 8.0 and mlcin >= -10 and mllcl_m >= 1000:
    modes.append("NON-MESO TORNADO POSSIBLE")
    flags.append("Steep 0-3 km lapse rates + weak CIN — non-mesocyclone tornado environment")

# ── Wet microburst ────────────────────────────────────────────────────────
if dcape >= 800 and bl_rh >= 60 and lr_03 and lr_03 >= 7.0:
    flags.append("⚠ Wet microburst potential: DCAPE {:.0f} J/kg + steep lapse rates".format(dcape))

# ── Limiting factors / fail modes ────────────────────────────────────────
if mlcin < -100:
    issues.append("Strong CIN ({:.0f} J/kg) — cap may prevent convective initiation".format(mlcin))
elif mlcin < -50:
    issues.append("Moderate CIN ({:.0f} J/kg) — initiation may be difficult without strong forcing".format(mlcin))

if mllcl_m > 1500 and "TORNADIC SUPERCELL" in modes:
    issues.append("High LCL ({:.0f} m) — cold RFD possible, reduces tornadogenesis potential".format(mllcl_m))

if shear_06 < 35 and srh_03 >= 150:
    issues.append("0-6 km shear < 35 kt — any rotation acquired may not persist (tip sheet)")

if shear_06 >= 40 and srh_03 < 150:
    issues.append("SRH < 150 m²/s² despite strong shear — supercell possible but rotation may be limited")

if mlcape < 500:
    issues.append("MLCAPE very low ({:.0f} J/kg) — convective initiation and intensity uncertain".format(mlcape))

if bl_rh < 50 and "TORNADIC SUPERCELL" in modes:
    issues.append("Low boundary layer RH ({:.0f}%) — high LCL / cold RFD risk".format(bl_rh))

# ── Overall risk ──────────────────────────────────────────────────────────
if "TORNADIC SUPERCELL" in modes and stp >= 2:
    risk = "HIGH"
elif "TORNADIC SUPERCELL" in modes or (scp >= 4 and srh_01 >= 120):
    risk = "MODERATE"
elif scp >= 2 or srh_01 >= 120 or "BOW ECHO / DERECHO" in modes:
    risk = "MARGINAL"
else:
    risk = "LOW"

return risk, modes, flags, issues
```

# ══════════════════════════════════════════════════════════════════════════════

# GEOCODING / SOUNDING FETCH

# ══════════════════════════════════════════════════════════════════════════════

def geocode_city(city_name):
url = “https://nominatim.openstreetmap.org/search”
r   = requests.get(url,
params={“q”: city_name, “format”: “json”, “limit”: 1},
headers={“User-Agent”: “WeatherRiskCheckerV2/1.0”},
timeout=10)
r.raise_for_status()
res = r.json()
if not res:
raise ValueError(f”City not found: {city_name}”)
return float(res[0][“lat”]), float(res[0][“lon”])

MODEL_MAP = {“GFS”: “gfs”, “NAM”: “nam”, “HRRR”: “hrrr”, “RAP”: “rap”}

def fetch_sounding(lat, lon, model=“gfs”):
now_utc  = datetime.now(timezone.utc)
run_hour = now_utc.hour if model in (“hrrr”, “rap”) else (now_utc.hour // 6) * 6
time_str = now_utc.strftime(f”%Y-%m-%dT{run_hour:02d}:00:00Z”)

```
url = (f"https://mesonet.agron.iastate.edu/json/sounding.py"
       f"?lon={lon}&lat={lat}&model={model}&valid={time_str}")

r = requests.get(url, timeout=20)
r.raise_for_status()
payload = r.json()

if "error" in payload:
    raise ValueError(f"IEM error: {payload['error']}")

levels   = []
sfc_hght = None

for lv in payload.get("profile", []):
    try:
        pres = float(lv["pres"]); hght = float(lv["hght"])
        temp = float(lv["tmpc"]); dwpt = float(lv["dwpc"])
        wdir = float(lv["drct"]); wspd = float(lv["sknt"])
    except (KeyError, TypeError, ValueError):
        continue
    if sfc_hght is None:
        sfc_hght = hght
    levels.append({"pres": pres, "hght": hght - sfc_hght,
                    "hght_msl": hght, "temp": temp, "dwpt": dwpt,
                    "wdir": wdir, "wspd": wspd})

if len(levels) < 5:
    raise ValueError("Insufficient sounding levels returned.")

levels.sort(key=lambda x: -x["pres"])
return levels
```

# ══════════════════════════════════════════════════════════════════════════════

# GUI

# ══════════════════════════════════════════════════════════════════════════════

class WeatherApp:
def **init**(self, root):
self.root = root
root.title(“Severe Weather Risk Checker V2.1”)
root.configure(bg=”#0d0d1f”)
root.geometry(“960x860”)
root.resizable(True, True)
self._build_ui()

```
def _build_ui(self):
    style = ttk.Style()
    style.theme_use("clam")
    for s in ("TFrame", "TLabel"):
        style.configure(s, background="#0d0d1f", foreground="#e0e0e0",
                        font=("Consolas", 10))
    style.configure("Header.TLabel", font=("Consolas", 17, "bold"),
                    foreground="#00d4ff", background="#0d0d1f")
    style.configure("Sub.TLabel", foreground="#7799bb",
                    background="#0d0d1f", font=("Consolas", 9))
    style.configure("TButton", font=("Consolas", 11, "bold"),
                    background="#0055cc", foreground="white")
    style.configure("TCombobox", fieldbackground="#151530",
                    foreground="white", background="#0d0d1f")
    style.configure("TEntry", fieldbackground="#151530",
                    foreground="white", insertcolor="white")

    # Header
    hdr = ttk.Frame(self.root)
    hdr.pack(fill="x", pady=(12, 2))
    ttk.Label(hdr, text="⛈  SEVERE WEATHER RISK CHECKER V2.1",
              style="Header.TLabel").pack()
    ttk.Label(hdr, text="GFS • NAM • HRRR • RAP  |  WFO Louisville Tip Sheet",
              style="Sub.TLabel").pack()

    # Input row
    inp = ttk.Frame(self.root)
    inp.pack(fill="x", padx=20, pady=6)

    ttk.Label(inp, text="City:").grid(row=0, column=0, sticky="w", padx=4)
    self.city_var = tk.StringVar(value="Oklahoma City, OK")
    ttk.Entry(inp, textvariable=self.city_var, width=26).grid(
        row=0, column=1, padx=4)

    ttk.Label(inp, text="OR Lat:").grid(row=0, column=2, padx=(10, 2))
    self.lat_var = tk.StringVar()
    ttk.Entry(inp, textvariable=self.lat_var, width=9).grid(row=0, column=3)
    ttk.Label(inp, text="Lon:").grid(row=0, column=4, padx=(6, 2))
    self.lon_var = tk.StringVar()
    ttk.Entry(inp, textvariable=self.lon_var, width=9).grid(row=0, column=5)

    ttk.Label(inp, text="Model:").grid(row=0, column=6, padx=(12, 4))
    self.model_var = tk.StringVar(value="HRRR")
    ttk.Combobox(inp, textvariable=self.model_var,
                 values=list(MODEL_MAP.keys()),
                 width=7, state="readonly").grid(row=0, column=7, padx=4)

    self.run_btn = ttk.Button(inp, text="▶  ANALYZE",
                              command=self._run_threaded)
    self.run_btn.grid(row=0, column=8, padx=(14, 4))

    self.status_var = tk.StringVar(value="Ready.")
    ttk.Label(self.root, textvariable=self.status_var,
              foreground="#556688").pack(anchor="w", padx=22)

    # Results
    frame = tk.Frame(self.root, bg="#080818", relief="sunken", bd=2)
    frame.pack(fill="both", expand=True, padx=20, pady=(2, 10))

    self.txt = scrolledtext.ScrolledText(
        frame, bg="#080818", fg="#d0d0e8",
        font=("Consolas", 10), state="disabled",
        relief="flat", padx=14, pady=12, wrap="word")
    self.txt.pack(fill="both", expand=True)

    # Colour tags
    T = self.txt.tag_config
    T("title",    foreground="#00d4ff", font=("Consolas", 13, "bold"))
    T("risk_HIGH",     foreground="#ff2020", font=("Consolas", 15, "bold"))
    T("risk_MODERATE", foreground="#ff8800", font=("Consolas", 14, "bold"))
    T("risk_MARGINAL", foreground="#ffee00", font=("Consolas", 13, "bold"))
    T("risk_LOW",      foreground="#44cc66", font=("Consolas", 13, "bold"))
    T("mode",     foreground="#ff9944", font=("Consolas", 11, "bold"))
    T("flag",     foreground="#ffcc55", font=("Consolas", 10, "bold"))
    T("issue",    foreground="#ff5555")
    T("section",  foreground="#66bb88", font=("Consolas", 10, "underline"))
    T("param",    foreground="#aaddff")
    T("hodo",     foreground="#cc88ff")
    T("muted",    foreground="#555577")
    T("good",     foreground="#44cc66")

# ── Threading ──────────────────────────────────────────────────────────────
def _run_threaded(self):
    self.run_btn.config(state="disabled")
    self.status_var.set("Fetching sounding…")
    threading.Thread(target=self._analyze, daemon=True).start()

def _analyze(self):
    try:
        lat_s = self.lat_var.get().strip()
        lon_s = self.lon_var.get().strip()
        city  = self.city_var.get().strip()

        if lat_s and lon_s:
            lat, lon  = float(lat_s), float(lon_s)
            loc_label = f"{lat:.3f}°N, {lon:.3f}°E"
        elif city:
            self.status_var.set("Geocoding…")
            lat, lon  = geocode_city(city)
            loc_label = city
        else:
            raise ValueError("Enter a city or lat/lon.")

        model_key = self.model_var.get()
        model     = MODEL_MAP[model_key]

        self.status_var.set(f"Downloading {model_key} sounding…")
        levels = fetch_sounding(lat, lon, model)

        self.status_var.set("Computing parameters…")

        # Mixed-layer parcel
        ml_t, ml_td, ml_p = mixed_layer_parcel(levels)
        mlcape, mlcin, mllcl_m, lfc_p, el_p = cape_cin(
            levels, ml_t, ml_td, ml_p)

        # Bunkers
        rm_u, rm_v, lm_u, lm_v, mw_u, mw_v = bunkers_storm_motion(levels)
        rm_spd = uv_to_spd(rm_u, rm_v); rm_dir = uv_to_dir(rm_u, rm_v)
        lm_spd = uv_to_spd(lm_u, lm_v); lm_dir = uv_to_dir(lm_u, lm_v)

        # SRH
        srh_01 = compute_srh(levels, rm_u, rm_v, 1000)
        srh_03 = compute_srh(levels, rm_u, rm_v, 3000)

        # Effective inflow layer
        base_idx, top_idx = effective_inflow_layer(levels)
        eff_shear_kt = effective_bulk_shear(levels, base_idx, top_idx)
        if base_idx is not None and top_idx is not None:
            eff_base_h = levels[base_idx]["hght"]
            eff_top_h  = levels[top_idx]["hght"]
            eff_srh    = compute_srh(levels, rm_u, rm_v, eff_top_h)
        else:
            eff_base_h = 0; eff_top_h = 1000; eff_srh = srh_01

        # Bulk shear layers
        sh_01 = bulk_shear(levels, 1000)
        sh_03 = bulk_shear(levels, 3000)
        sh_06 = bulk_shear(levels, 6000)
        sh_08 = bulk_shear(levels, 8000)

        # BRN
        brn, brn_shear = compute_brn(mlcape, levels)

        # Lapse rates
        lr_03     = lapse_rate_layer(levels, 0, 3000)
        lr_700_500 = lapse_rate_pressure(levels, 700, 500)

        # BL RH + PW
        bl_rh = boundary_layer_rh(levels)
        pw    = precipitable_water(levels)

        # Hodograph
        hodo = hodograph_shape(levels)

        # DCAPE
        dcape = compute_dcape(levels)

        # Composites
        scp = supercell_composite(mlcape, eff_srh, eff_shear_kt)
        stp = significant_tornado(mlcape, mllcl_m, mlcin,
                                  eff_srh, eff_shear_kt)

        # Storm mode + risk
        risk, modes, flags, issues = classify_storm_mode(
            mlcape, mlcin, mllcl_m,
            sh_01, sh_03, sh_06, sh_08,
            srh_01, srh_03, scp, stp,
            brn, brn_shear,
            lr_03, lr_700_500,
            bl_rh, pw, hodo, dcape)

        self.root.after(0, lambda: self._render(
            loc_label, model_key, levels,
            mlcape, mlcin, mllcl_m, lfc_p, el_p,
            rm_spd, rm_dir, lm_spd, lm_dir,
            srh_01, srh_03, eff_srh, eff_base_h, eff_top_h,
            sh_01, sh_03, sh_06, sh_08, eff_shear_kt,
            brn, brn_shear,
            lr_03, lr_700_500, bl_rh, pw,
            hodo, dcape, scp, stp,
            risk, modes, flags, issues))

        self.status_var.set("Analysis complete.")

    except Exception as e:
        self.root.after(0, lambda: self._show_error(str(e)))
    finally:
        self.root.after(0, lambda: self.run_btn.config(state="normal"))

# ── Render ────────────────────────────────────────────────────────────────
def _w(self, text, tag=None):
    self.txt.config(state="normal")
    self.txt.insert("end", text, tag) if tag else self.txt.insert("end", text)
    self.txt.config(state="disabled")

def _clear(self):
    self.txt.config(state="normal")
    self.txt.delete("1.0", "end")
    self.txt.config(state="disabled")

def _render(self, loc, model, levels,
            mlcape, mlcin, mllcl_m, lfc_p, el_p,
            rm_spd, rm_dir, lm_spd, lm_dir,
            srh_01, srh_03, eff_srh, eff_base_h, eff_top_h,
            sh_01, sh_03, sh_06, sh_08, eff_shear_kt,
            brn, brn_shear, lr_03, lr_700_500, bl_rh, pw,
            hodo, dcape, scp, stp,
            risk, modes, flags, issues):

    self._clear()
    ts = datetime.now().strftime("%Y-%m-%d %H:%M")

    self._w(f"\n  ⛈  SEVERE WEATHER ANALYSIS  —  {loc}  [{model}]  {ts}\n\n", "title")
    self._w(f"  OVERALL RISK: {risk}\n\n", f"risk_{risk}")

    # Storm modes
    if modes:
        self._w("  ── STORM MODE DIAGNOSIS ─────────────────────────────\n", "section")
        for m in modes:
            self._w(f"  ▸ {m}\n", "mode")
        self._w("\n")

    # Instability
    self._w("  ── INSTABILITY ──────────────────────────────────────\n", "section")
    lfc_s = f"{lfc_p:.0f} hPa" if lfc_p else "—"
    el_s  = f"{el_p:.0f} hPa"  if el_p  else "—"
    self._w(f"  MLCAPE : {mlcape:>7.0f} J/kg   MLCIN  : {mlcin:>6.0f} J/kg   MLLCL : {mllcl_m:>5.0f} m\n", "param")
    self._w(f"  LFC    : {lfc_s:>10}   EL     : {el_s:>10}\n", "param")

    # Lapse rates
    lr03_s     = f"{lr_03:.1f} °C/km"     if lr_03      is not None else "—"
    lr700_s    = f"{lr_700_500:.1f} °C/km" if lr_700_500 is not None else "—"
    self._w(f"  LR 0-3 km : {lr03_s:>12}   LR 700-500 mb : {lr700_s}\n", "param")

    # Moisture
    self._w("\n  ── MOISTURE ─────────────────────────────────────────\n", "section")
    self._w(f"  BL RH : {bl_rh:>5.0f}%   PW : {pw:.2f} in\n", "param")

    # Shear
    self._w("\n  ── WIND SHEAR ───────────────────────────────────────\n", "section")
    self._w(f"  0-1 km : {sh_01:>5.1f} kt   0-3 km : {sh_03:>5.1f} kt\n", "param")
    self._w(f"  0-6 km : {sh_06:>5.1f} kt   0-8 km : {sh_08:>5.1f} kt\n", "param")
    self._w(f"  Eff Bulk Shear : {eff_shear_kt:>5.1f} kt  (layer {eff_base_h:.0f}–{eff_top_h:.0f} m AGL)\n", "param")

    # SRH
    self._w("\n  ── HELICITY ─────────────────────────────────────────\n", "section")
    self._w(f"  0-1 km SRH : {srh_01:>6.0f} m²/s²   0-3 km SRH : {srh_03:>6.0f} m²/s²\n", "param")
    self._w(f"  Eff SRH    : {eff_srh:>6.0f} m²/s²\n", "param")

    # Bunkers
    self._w("\n  ── BUNKERS STORM MOTION ─────────────────────────────\n", "section")
    self._w(f"  Right Mover : {rm_dir:.0f}° @ {rm_spd:.0f} kt\n", "hodo")
    self._w(f"  Left Mover  : {lm_dir:.0f}° @ {lm_spd:.0f} kt\n", "hodo")

    # BRN
    self._w("\n  ── BULK RICHARDSON NUMBER ───────────────────────────\n", "section")
    brn_s = f"{brn:.0f}" if brn is not None else "—"
    self._w(f"  BRN : {brn_s:>6}   BRN Shear : {brn_shear:.0f} m²/s²\n", "param")

    # Hodograph
    self._w("\n  ── HODOGRAPH ────────────────────────────────────────\n", "section")
    if hodo:
        tag = "flag" if hodo["tornadic_sig"] else "hodo"
        self._w(f"  Shape : {hodo['shape']}\n", tag)
        for n in hodo["notes"]:
            self._w(f"  ↳ {n}\n", tag)

    # Composites
    self._w("\n  ── COMPOSITE INDICES ────────────────────────────────\n", "section")
    self._w(f"  SCP (Supercell Composite — Eff) : {scp:>6.1f}  {'[≥2]' if scp>=2 else ''}\n", "param")
    self._w(f"  STP (Sig Tornado — Eff Layer)   : {stp:>6.2f}  {'[≥1 ELEVATED]' if stp>=1 else ''}\n", "param")
    self._w(f"  DCAPE (400-500 hPa)              : {dcape:>6.0f} J/kg\n", "param")

    # Flags
    if flags:
        self._w("\n  ── ACTIVE FLAGS ─────────────────────────────────────\n", "section")
        for f in flags:
            self._w(f"  {f}\n", "flag")

    # Issues
    if issues:
        self._w("\n  ── LIMITING FACTORS / FAIL MODES ────────────────────\n", "section")
        for iss in issues:
            self._w(f"  ✗ {iss}\n", "issue")

    # Surface
    sfc = levels[0]
    self._w("\n  ── SURFACE ──────────────────────────────────────────\n", "section")
    self._w(
        f"  T: {sfc['temp']:.1f}°C  Td: {sfc['dwpt']:.1f}°C  "
        f"Wind: {sfc['wdir']:.0f}°@{sfc['wspd']:.0f}kt  "
        f"Pres: {sfc['pres']:.0f} hPa\n", "muted")

    self._w("\n")

def _show_error(self, msg):
    self._clear()
    self._w(f"\n  ERROR: {msg}\n\n"
            f"  Check location, model availability, or network.\n"
            f"  HRRR/RAP are CONUS-only. Try GFS if others fail.\n", "issue")
```

# ── Entry ──────────────────────────────────────────────────────────────────────

if **name** == “**main**”:
root = tk.Tk()
WeatherApp(root)
root.mainloop()
