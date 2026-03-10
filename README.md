# ⛈ Severe Weather Risk Checker V2

A desktop severe weather analysis tool that fetches model soundings and computes a full suite of meteorological parameters to identify tornadic supercell, bow echo, MCS, hail, and other severe weather setups.

Built on the **WFO Louisville Severe Weather Forecasting Tip Sheet** and standard composite index formulations (Thompson 2004).

-----

## Features

- **Model sounding analysis** — fetches vertical profiles from GFS, NAM, HRRR, or RAP via the Iowa Environmental Mesonet (IEM) API
- **Location input** — city name (geocoded via OpenStreetMap Nominatim), lat/lon, or both
- **Full thermodynamic computation** — all parameters calculated locally from the raw sounding profile using Bolton (1980) formulas
- **Storm mode classifier** — explicitly diagnoses likely storm mode based on environmental parameters

-----

## Parameters Computed

### Instability

|Parameter            |Description                                 |
|---------------------|--------------------------------------------|
|MLCAPE               |Mixed-layer CAPE (lowest 100 hPa average)   |
|MLCIN                |Mixed-layer CIN                             |
|MLLCL                |Mixed-layer LCL height (m AGL)              |
|LFC / EL             |Level of Free Convection / Equilibrium Level|
|Lapse Rate 0–3 km    |Environmental lapse rate, °C/km             |
|Lapse Rate 700–500 mb|Mid-level lapse rate, °C/km                 |

### Wind Shear

|Parameter           |Description                                        |
|--------------------|---------------------------------------------------|
|0–1 km bulk shear   |Low-level shear (kt)                               |
|0–3 km bulk shear   |Used for bow echo / bowing segment potential       |
|0–6 km bulk shear   |Primary supercell vs. multicell discriminator      |
|0–8 km bulk shear   |Long-lived supercell / long-track tornado indicator|
|Effective bulk shear|Shear through the effective inflow layer           |
|BRN / BRN Shear     |Bulk Richardson Number and shear term              |

### Helicity

|Parameter    |Description                                  |
|-------------|---------------------------------------------|
|0–1 km SRH   |Storm-relative helicity (Bunkers right-mover)|
|0–3 km SRH   |Storm-relative helicity (Bunkers right-mover)|
|Effective SRH|SRH within the effective inflow layer        |

### Composite Indices

|Parameter|Formulation                                                         |
|---------|--------------------------------------------------------------------|
|SCP      |Supercell Composite Parameter — effective layer (Thompson 2004)     |
|STP      |Significant Tornado Parameter — effective layer (Thompson 2004/2012)|
|DCAPE    |Downdraft CAPE from the 400–500 hPa layer                           |

### Storm Motion

- **Bunkers Right Mover** — direction and speed (kt)
- **Bunkers Left Mover** — direction and speed (kt)
- Internal dynamics method (Bunkers 2000)

### Moisture

- Boundary layer mean RH (lowest 100 hPa)
- Precipitable water (inches)

### Hodograph Shape

- Straight / Curved-CW / Curved-CCW / Backed
- Tornadic hodograph signature detection: straight 0–1 km spike with curvature above

-----

## Storm Mode Classifier

The tool diagnoses which storm modes the environment supports, using WFO Louisville tip sheet thresholds:

- **Tornadic Supercell** — SCP, STP, SRH, LCL, BL RH, hodograph signature
- **Supercell** — 0–6 km shear ≥ 40 kt, SCP ≥ 2; or marginal shear with CAPE ≥ 2500 J/kg
- **Bow Echo / QLCS** — 0–3 km shear ≥ 30–40 kt
- **Derecho** — high CAPE + fast unidirectional low/mid-level flow
- **Organized Multicell / MCS** — 0–6 km shear 20–35 kt
- **Disorganized Multicell / Pulse** — 0–6 km shear < 20 kt
- **Mini-Supercell** — low CAPE + sufficient low-level shear
- **Non-Mesocyclone Tornado** — steep 0–3 km lapse rates + weak CIN + high LCL

Each output also includes **limiting factors / fail modes** — explicit notes on what is suppressing or capping the severe weather potential (e.g. strong CIN, high LCL, insufficient SRH).

### Risk Thresholds

|Risk Level|Criteria                                      |
|----------|----------------------------------------------|
|HIGH      |Tornadic supercell mode + STP ≥ 2             |
|MODERATE  |Tornadic supercell mode OR SCP ≥ 4 + SRH ≥ 120|
|MARGINAL  |SCP ≥ 2 OR SRH ≥ 120 OR derecho environment   |
|LOW       |Below all thresholds                          |

Base watch thresholds: **SCP > 2, 0–1 km SRH > 120 m²/s²**

-----

## Data Sources

|Source                                                                |Use                                |
|----------------------------------------------------------------------|-----------------------------------|
|[Iowa Environmental Mesonet (IEM)](https://mesonet.agron.iastate.edu/)|GFS, NAM, HRRR, RAP model soundings|
|[OpenStreetMap Nominatim](https://nominatim.openstreetmap.org/)       |City name geocoding                |


> **Note:** HRRR and RAP are CONUS-only. If a sounding is unavailable for the current model run, try GFS or NAM as a fallback.

-----

## Installation

**Requirements:** Python 3.8+

```bash
pip install requests
```

`tkinter` is included in the Python standard library. On some Linux systems you may need:

```bash
sudo apt install python3-tk
```

**Run:**

```bash
python weather_checker_v2.1.py
```

-----

## Usage

1. Enter a **city name** (e.g. `Oklahoma City, OK`) or **lat/lon** coordinates
1. Select a **model** — HRRR is recommended for short-range (< 18h), GFS for longer range
1. Click **▶ ANALYZE**
1. Review the storm mode diagnosis, parameter table, active flags, and limiting factors

-----

## References

- **Bolton, D. (1980)** — The computation of equivalent potential temperature. *Monthly Weather Review*, 108, 1046–1053.
- **Bunkers, M.J. et al. (2000)** — Predicting supercell motion using a new hodograph technique. *Weather and Forecasting*, 15, 61–79.
- **Davies-Jones, R. et al. (1990)** — An overview of tornado observations. *AMS Monograph*.
- **Thompson, R.L. et al. (2004)** — Close proximity soundings within supercell environments. *Weather and Forecasting*, 19, 1–21.
- **WFO Louisville Severe Weather Forecasting Tip Sheet** — Vertical wind shear, SRH, storm mode, tornado, hail, bow echo, and MCS criteria.

-----

## Version History

|Version|Notes                                                                                                                                                   |
|-------|--------------------------------------------------------------------------------------------------------------------------------------------------------|
|V1     |Open-Meteo API, basic CAPE/shear/SRH, NiceGUI interface                                                                                                 |
|V2     |IEM model soundings, full sounding thermodynamics, Bunkers, effective layer STP/SCP, DCAPE, Tkinter GUI                                                 |
|V2.1   |WFO Louisville tip sheet integration: BRN, bulk shear layers, lapse rates, BL RH, PW, hodograph shape classifier, storm mode diagnosis, fail mode output|

-----

## Disclaimer

This tool is intended to assist forecasters familiar with severe weather meteorology. Output should always be interpreted alongside NWS products, current observations, and radar data. Do not use as a sole basis for severe weather decisions.