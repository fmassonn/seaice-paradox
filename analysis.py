import numpy as np
from netCDF4 import Dataset

# functions to compute sea ice area/extent
from seaice_commondiags import * 

REarth = 6356e3 # Earth'radius in meters, to compute grid cell areas

# Analysis parameters
# -------------------

# Where the source data is located [string]
dataRootDir = "/cofast/fmasson/ECMWF_SEAICE/"

# The ensemble members to be selected [list of strings]
members = [str(j) for j in np.arange(50 + 1)]

# Months of initialization [list of strings]
monthsInit = [str(m).zfill(2) for m in [2, 5, 8, 11]]

# Years of initialization [list of strings]
yearsInit = [str(y) for y in np.arange(1993, 2022 + 1)]

# Definition of masks [list of lists each containing a name, and a list of coordinates defining the boundaries]
            # Name              lonW        lonE    latS    latN
masks = [[  "Arctic",           [-180.0,    180.0,  0.0,    90.0]], \
         [  "Antarctic",        [-180.0,    180.0,  -90.0,  90.0]], \
        ]

nDiags = 2 # Sea ice extent and sea ice area, others can be added depending on how many diagnostics are calculated in the loop below (if so, update this integer)
nMasks = len(masks)
nMonthsInit = len(monthsInit)
nYearsInit = len(yearsInit)
nMembers = len(members)
nTime = 7 # The lead time (in months)

# We start looping over months of initialization
# We set a boolean variable to "True" to create grid cell areas
createCellArea = True
createMasks    = True
createData     = True

for jmI, monthInit in enumerate(monthsInit):
    # Then we loop over the years of initialization for a given month of initializatio
    for jyI, yearInit in enumerate(yearsInit):
        # Then we loop over the ensemble members
        for jmb, member in enumerate(members):
            fileName = "seas5_ens_" + member + "_sic_" + yearInit + monthInit + "01" + ".nc"
            fileIn = dataRootDir + fileName

            try:

                print("Reading file " + fileIn)
                f = Dataset(fileIn, mode = "r")
                sic = f.variables["sic"][:] * 100.0 # Convert to %
                lat = f.variables["latitude"][:]
                lon = f.variables["longitude"][:]
                f.close()

                if createCellArea:
                    if len(set(np.diff(lon))) != 1 or len(set(np.diff(lat))) != 1:
                                # If the spacing is not regular
                        stop("Spacing not regular")
                    else:
                        dLon = np.min(np.diff(lon))
                        dLat = np.min(np.diff(lat)) # min or max does not
                                                    # since all equal

                        # Reset lon to [-180.0, 180.0]
                        lon[lon > 180.0] = lon[lon > 180.0] - 360.0

                        # mesh lon and lat from 1D to 2D
                        LON, LAT = np.meshgrid(lon, lat)

                        cellArea = dLon / 360.0 * 2 * np.pi * REarth \
                                * np.cos(LAT / 360.0 * 2 * np.pi) * \
                                REarth * dLat / 360.0 * 2 * np.pi

                        createCellArea = False
        
                if createMasks:
                    # Create masks of 0's and 1 corresponding to the masks
                    for jms, mask in enumerate(masks):
                        lonW, lonE, latS, latN = mask[1][0], mask[1][1], mask[1][2], mask[1][3]

                        # Two cases depending on whether the region straddles the 180° meridian
                        if lonW < lonE:
                            maskValue = (LON >= lonW) * (LON <= lonE) * (LAT >= latS) * (LAT <= latN) * 1
                        else:
                            maskValue = ((LON >= lonW) + (LON <= lonE)) * (LAT >= latS) * (LAT <= latN) * 1 

                        masks[jms].append(maskValue)
                    
                    createMasks = False


                # Loop over masks
                for jms, mask in enumerate(masks):
                    # Compute the sea ice area and sea ice extents
                    sia = compute_area(sic, cellArea, mask = mask[2])
                    sie = compute_extent(sic, cellArea, mask = mask[2])


                    # Save that in an array that has as dimensions:
                    # 1) diagnostic (extent then area) 
                    # 2) regions (masks)
                    # 3) month of initialization
                    # 4) year of initialization
                    # 5) ensemble member
                    # 6) time

                    if createData:
                        data = np.full((nDiags, nMasks, nMonthsInit, nYearsInit, nMembers, nTime), np.nan)
                        createData = False

                    data[0, jms, jmI, jyI, jmb, :] = sie
                    data[1, jms, jmI, jyI, jmb, :] = sie


                
                # Plot the raw data
                # One figure per region and per month of initialization

            except FileNotFoundError:
                print("File not found: " + fileIn)

