import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import pickle
import os

# functions to compute sea ice area/extent
from seaice_commondiags import * 



REarth = 6356e3 # Earth'radius in meters, to compute grid cell areas
monthsNames = ["January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"]
# Our goal in this script is first to create an array called "forecasts" that will record several diagnostics out of the native
# sea ice data provided by ECMWF. The array will have the following dimensions:

# 1st dimension: the *d*iagnostic (e.g., extent, area of sea ice)
# 2nd dimension: the *r*egion (e.g., Arctic, Ross Sea, ...)     
# 3rd dimension: the *m*onth of initialization (e.g., Jan, Feb, ...)
# 4th dimension: the *y*ear of initialization (e.g., 1993, 1994, ...)
# 5th dimension: the mem*b*er (e.g., 1, 2, ...)
# 6th dimension: the *l*ead time (e.g., 1, 2, ...)

# Each dimension has a corresponding running index called j${letter} where $letter is d, r, m, y, b, l following the dimensions above. Each index runs from 0 to the length of the dimension (not included, as per Python conventions)


# A switch to decide whether we re-compute the data or if we load them from a pkl file
createData = True

if createData:

    # Analysis parameters
    # -------------------
    
    # Where the source data is located [string]
    dataRootDir = "/cofast/fmasson/ECMWF_SEAICE/"
    
    # 1. Definition of diagnostics Sea ice extent and sea ice area, others can be added depending on how many diagnostics are calculated in the loop below (if so, update this integer)
    diagnostics = ["extent",]# "area"]
    unitsDiagnostics = ["10$^6$ km$^2$",]# "10$^6$ km$^2$"]
    nd = len(diagnostics)
    
    # 2. Definition of regions [list of lists each containing a name, and a list of coordinates defining the boundaries]
                  # Name              lonW        lonE    latS    latN
    regions = [[  "Arctic",           [-180.0,    180.0,  0.0,    90.0]], \
               [  "Antarctic",        [-180.0,    180.0,  -90.0,  90.0]], \
              ]
    nr = len(regions)
    
    # 3. Months of initialization [list of strings]
    monthsInit = [2, 5, 8, 11] #! NON-Pythonic: 1 = January.
    nm = len(monthsInit)
    
    # 4. Years of initialization [list of strings]
    yearsInit = np.arange(1993, 2022 + 1)
    #yearsInit = np.arange(1993, 1995 + 1)
    ny = len(yearsInit)
    
    # 5. The ensemble members to be selected [list of strings]
    members = [str(j) for j in np.arange(50 + 1)]
    nb = len(members)
    
    # 6. The lead times (in months from initialization)
    leadTimes = [0, 1, 2, 3, 4, 5, 6]
    nl = len(leadTimes)
    
    
    # Phase 1: producing the diagnostics
    
    # For each of the diagnostics, we start looping over months of initialization
    # We set a boolean variable to "True" to create grid cell areas only once
    
    createCellArea = True
    createMasks    = True
    createDataArray= True
    
    
    # Our goal here is to create an array called "forecasts" of dimensions nd, nr, nm, ny, nb, nl
    
    # Each dimension has a running index called j${letter} where $letter is d, r, m, y, b, l following the dimensions above
    # The order of the loops should be organized to avoid opening and closing the same file many times. The files are organized
    # as one file =  one member for one year+month initial time seas5_ens_26_sic_20050801.nc. The region, lead time should therefore
    # arrive as last loops
    
    for jm, monthInit in enumerate(monthsInit):
        # Then we loop over the years of initialization for a given month of initialization
        for jy, yearInit in enumerate(yearsInit):
            # Then we loop over the ensemble members
            for jb, member in enumerate(members):
                fileName = "seas5_ens_" + member + "_sic_" + str(yearInit) + str(monthInit).zfill(2) + "01" + ".nc"
                fileIn = dataRootDir + fileName
    
                try:
                    #print("Reading file " + fileIn)
                    f = Dataset(fileIn, mode = "r")
                    sic = f.variables["sic"][:] * 100.0 # Convert to %
                    lat = f.variables["latitude"][:]
                    lon = f.variables["longitude"][:]
                    f.close()
    
                    if createCellArea:
                        print("Creating cell area")
    
                        if len(set(np.diff(lon))) != 1 or len(set(np.diff(lat))) != 1:
                                # If the spacing is not regular
                            stop("Spacing not regular")
                        else:
                            dLon = np.min(np.diff(lon))
                            dLat = np.min(np.diff(lat)) # min or max does not
                                                    # since all equal
    
                            # Reset lon to [-180.0, 180.0] to follow the region mask conventions
                            lon[lon > 180.0] = lon[lon > 180.0] - 360.0
    
                            # mesh lon and lat from 1D to 2D
                            LON, LAT = np.meshgrid(lon, lat)
    
                            cellArea = dLon / 360.0 * 2 * np.pi * REarth \
                                    * np.cos(LAT / 360.0 * 2 * np.pi) * \
                                    REarth * dLat / 360.0 * 2 * np.pi
    
                            createCellArea = False
        
                    if createMasks:
                        # Create masks of 0's and 1 corresponding to the masks
                        for jr, region in enumerate(regions):
                            lonW, lonE, latS, latN = region[1][0], region[1][1], region[1][2], region[1][3]
    
                            # Two cases depending on whether the region straddles the 180° meridian
                            if lonW < lonE:
                                maskValue = (LON >= lonW) * (LON <= lonE) * (LAT >= latS) * (LAT <= latN) * 1
                            else:
                                maskValue = ((LON >= lonW) + (LON <= lonE)) * (LAT >= latS) * (LAT <= latN) * 1 
    
                            # Update the regions variable with the mask
                            regions[jr].append(maskValue)
                    
                        createMasks = False
    
    
                    # Loop over masks
                    for jr, region in enumerate(regions):
                        # Compute the sea ice area and sea ice extents
                        sia = compute_area(  sic, cellArea, mask = region[2])
                        sie = compute_extent(sic, cellArea, mask = region[2])
    
    
                        # Save that in an array that has as dimensions:
                        # 1) diagnostic (extent then area) 
                        # 2) regions (defined by their masks)
                        # 3) month of initialization
                        # 4) year of initialization
                        # 5) ensemble member
                        # 6) lead time
    
                        if createDataArray:
                            forecastData = np.full((nd, nr, nm, ny, nb, nl), np.nan)
                            createDataArray = False
    
                        # Run through the diagnostics
                        for jd in range(nd):
                            if diagnostics[jd] == "extent":
                                forecastData[jd, jr, jm, jy, jb, :] = sie[leadTimes]
                            elif diagnostics[jd] == "area":
                                forecastData[jd, jr, jm, jy, jb, :] = sia[leadTimes]
    
    
                
                    # Plot the raw data
                    # One figure per region and per month of initialization
    
                except FileNotFoundError:
                    print("File not found: " + fileIn)
    
    # Save forecastData to a file
    saved_data = {
            'forecastData'     : forecastData,
            'diagnostics'      : diagnostics,
            'unitsDiagnostics' : unitsDiagnostics,
            'nd'               : nd,
            'regions'          : regions,
            'nr'               : nr,
            'monthsInit'       : monthsInit,
            'nm'               : nm,
            'yearsInit'        : yearsInit,
            'ny'               : ny,
            'members'          : members,
            'nb'               : nb,
            'leadTimes'        : leadTimes,
            'nl'               : nl,
            }
    with open('saved_data.pkl', 'wb') as f:
        pickle.dump(saved_data, f)
        print("saved_data.pkl saved to disk")

else: # If we don't create data, we read it from saved_data.pkl
    if os.path.exists('saved_data.pkl'):
        # Load the data from the file
        with open('saved_data.pkl', 'rb') as f:
            saved_data = pickle.load(f)
            print("Loaded forecastData from file.")

            forecastData = saved_data['forecastData']
            diagnostics = saved_data['diagnostics']
            unitsDiagnostics = saved_data['unitsDiagnostics']
            nd = saved_data['nd']
            regions = saved_data['regions']
            nr = saved_data['nr']
            monthsInit = saved_data['monthsInit']
            nm = saved_data['nm']
            yearsInit = saved_data['yearsInit']
            ny = saved_data['ny']
            members = saved_data['members']
            nb = saved_data['nb']
            leadTimes = saved_data['leadTimes']
            nl = saved_data['nl']
    else:
        # If the file doesn't exist, perform the calculations
        print("ERROR: saved_data.pkl not found")
        stop()



# Load observations (currently done very sloppy)
print("WARNING WARNING WARNING: OBS SLOPPY")
f = Dataset("siextentn_r1i1p1_mon_197801-202312.nc", mode= "r")
obs = f.variables["siextentn"][(yearsInit[0]- 1978) * 12 + 8:(yearsInit[-1] - 1978) * 12 + 12: 12]
f.close()


# Make statistics
# This array has dimensions nd, nr, nm, ny, nl
forecastDataEnsembleMean = np.mean(forecastData, axis = 4)

# Compute temporal standard deviation of ensemble mean
sigEnsembleMean = np.std(forecastDataEnsembleMean, axis = 3)
# Compute sig2_tot as in Eade et al: time average of per-year ensemble variance
sigEnsembleMembers = np.sqrt(np.mean(np.var(forecastData, axis = 4), axis = 3))
sigObs = np.std(obs)
# Make figures to check that everything. One figure per region, per lead time, per initialization month, and per diagnostic

for jr in range(nr):
    regionNoSpace = regions[jr][0].replace(" ", "")
    for jm in range(nm):
        for jd in range(nd):
            for jl in range(nl):
                fig, ax = plt.subplots(1, 1, figsize = (8, 5))
                labelFlag = True
                for jb in range(nb):
                    if labelFlag:
                        label = "Members"
                        labelFlag = False
                    else:
                        label = None
                    ax.scatter(yearsInit, forecastData[jd, jr, jm, :, jb, jl], 1, color  = [0.3, 0.3, 0.3], label = label)
                # Ensemble mean
                ax.scatter(yearsInit + 0.2, forecastDataEnsembleMean[jd, jr, jm, :, jl], 30, color = [1, 0.5, 0], label = "Ensemble mean")


                # Plots obs
                ax.scatter(yearsInit + 0.4, obs, 30, color = [0.0, 0.5, 0.0], marker  = "x", lw = 2, label = "OBS")
                ax.grid()
                ax.legend()
                ax.set_ylabel(unitsDiagnostics[jd])
                ax.set_axisbelow(True)
                ax.set_title(regions[jr][0] + " sea ice " + diagnostics[jd] + "\n" + "seas 5 - " + monthsNames[monthsInit[jm] - 1] + " initialization - " + monthsNames[(monthsInit[jm] + leadTimes[jl]) % 12 - 1] + " target" + \
                        "\n$(\sigma_{sig})^2$ = " + "(" + str(np.round(sigEnsembleMean[jd, jr, jm, jl], 2)) + " " + unitsDiagnostics[jd] + ")$^2$"  +  \
                        "\n$(\sigma_{tot})^2$ = " + "(" + str(np.round(sigEnsembleMembers[jd, jr, jm, jl], 2))  + " " + unitsDiagnostics[jd]  + ")$^2$"  \
                        "\n$(\sigma_{obs})^2$ = " + "(" + str(np.round(sigObs, 2))  + " " + unitsDiagnostics[jd]   + ")$^2$" \
                        "\n r(obs, EnsMean) = PC_obs =" + str(np.round(np.corrcoef(obs, forecastDataEnsembleMean[jd, jr, jm, :, jl])[0, 1], 2)))
                figName = "seas5_" + str(yearsInit[0]) + "-" + str(yearsInit[-1]) + "_" + diagnostics[jd] + "_m" + str(monthsInit[jm]).zfill(2) + "_l" + str(leadTimes[jl]).zfill(2) + "_" + regions[jr][0] + ".png"
                fig.tight_layout()
                plt.savefig("./figs/" + figName, dpi = 300)
                print(figName + "  printed")
                plt.close(fig)

