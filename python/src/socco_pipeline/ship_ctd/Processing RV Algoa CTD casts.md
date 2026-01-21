## Processing RV Algoa CTD casts

Processed by I.Giddy, 19 Jan 2026

Software SBE Data Processing
https://www.seabird.com/sbe-911plus-ctd/product-downloads?id=60761421595

Data processed from .hex, using alg307.xmlcon config file. 

1. Data Conversion
    - Include Header file (.hdr)
    - Select variables:
        - Latitude 
        - Longitude
        - C0 (S/m) (Conductivity)
        - Avgsv (WM) (sound velocity)
        - FlECO-AFL (*wetlabs, need to check this!)
        - Sbox0 (Mm/Kg) Oxygen
        - Oxsat (Mm/Kg) Oxygen Sat
        - Par
        - T090C (Temperature)
        - PrDM (Binned depths, m)
        - Time (time elapsed in seconds)
        - Density00 (density) 
        - DepSM  (depth, m)
        - Sal00 (Practical Salinity, PSU)
        - Gsw_saA0 (Gibbs seawater Absolute Salinity)
        - Gsw_ctA0 (Gibbs seawater Conservative Temperature)
    
2. Align CTD
    - Temperature: 0.000 s
    - Conductivity: +0.073 s
    - Oxygen (SBE 43): +2 s
    - Optics: +0.5 s
 
3. Cell Thermal Mass
    - use defaults
        - Alpha: 0.03
        - Tau: 7.0

4. Loop edit
    - minimum fall rate 0.25m/s
    
5. Derive additional variables

6. Filter 
    - use defaults

7. Bin average 
    - 1db

8. Split
    - For physics use downcasts only

9. Output ACSII




