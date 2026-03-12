# so_biological_pump

Python scripts used to produce diatom transfer functions submitted to Climate of the Past:
"Southern Ocean biological pump over the last glacial cycle from new diatom transfer functions" by Rembauville and Pichat. https://egusphere.copernicus.org/preprints/2026/egusphere-2025-6347/

**1_map.py** : Produces the Southern Ocean Map with fronts, sediment traps and sediment cores locations

**2_calibration.py** : Calibrates the transfer functions (MLR : multiple linear regression, PLSR : partial least square regression, GBR : gradient boosting regression) using the sediment trap data to reconstruct POC fluxes and PIC:POC export ratio.

**3_application.py** : Applies the transfer functions to the sediment cores diatom assemblages. 

**4_comparison.py** : Compares the reconstructed POC flux and PIC:POC export ratio with previously published data for cores PS1768-8, PS1786-1, PS2606-6 and PS75/072-4.
