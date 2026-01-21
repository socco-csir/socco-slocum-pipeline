import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt 
from glob import glob
from pathlib import Path

data_dir='./data/Alg307-CTD-Data-Processed/down_casts/'
out_dir = './data/Alg307-CTD-Data-Processed/'

# merge ctd casts
files = sorted(glob(data_dir+"*.asc"))

dfs = []
for f in files:
    d = pd.read_csv(f,delim_whitespace=True)

    # optional: parse time if you have a time column
    if "time" in d.columns:
        d["time"] = pd.to_datetime(d["time"], errors="coerce", utc=True)

    d["cast_id"] = Path(f).stem   
    dfs.append(d)

df = pd.concat(dfs, ignore_index=True, sort=False)

# sort first columns to cast_id, timeS, Latitude, Longitude, Depth, then the rest
cols = df.columns.tolist()
sorted_cols = ['cast_id', 'TimeS', 'Latitude', 'Longitude', 'DepSM','PrDM'] + [c for c in cols if c not in ['cast_id', 'TimeS', 'Latitude', 'Longitude', 'DepSM','PrDM']]
df = df[sorted_cols]

# Rename columns for easier access
# TODO: update to follow standard naming conventions
df = df.rename(columns={'TimeS': 'time', 
                        'Latitude': 'lat',
                        'Longitude': 'lon',
                        'DepSM': 'depth_m',
                        'PrDM': 'pressure_dbar',
                        'COS/m':'conductivity_s/m',
                        'AvgsvWM':'sound velocity',
                        'FlECO-AFL': 'flr',
                        'Sbox0Mm/Kg': 'oxygen_concentration_mmol_kg',
                        'OxsatMm/Kg': 'oxygen_saturation_mmol_kg',
                        'Par':'PAR',
                        'T090C':'temperature_C',
                        'Density00':'density',
                        'Sal00':'practical_salinity_PSU',
                        'Gsw_saA0':'absolute_salinity',
                        'Gsw_ctA0':'conservative_temperature'
                        })

# save as csv
# TODO: save as netcdf with metadata

df.to_csv(out_dir+'Alg_dcasts.csv')


