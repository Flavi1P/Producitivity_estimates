from erddapy import ERDDAP
import pandas as pd

# Define server
e = ERDDAP(server='https://erddap.ifremer.fr/erddap',
           protocol='tabledap')

# Use the global Argo profile index dataset
e.dataset_id = 'ArgoFloats-index'

# Define constraints
constraints = {
    'longitude>=': -45,
    'longitude<=': -10,
    'latitude>=': 54,
    'latitude<=': 63,
    'parameter_has_do': '=1',  # Only floats that have DOXY
}

# Select variables of interest
e.variables = [
    'platform_number',
    'latitude',
    'longitude',
    'time'
]

# Fetch data as pandas DataFrame
df = e.to_pandas(constraints=constraints)

# Clean WMO list
wmo_list = sorted(df['platform_number'].dropna().unique().astype(int))

print(f"Found {len(wmo_list)} floats with DOXY in region:")
print(wmo_list)
