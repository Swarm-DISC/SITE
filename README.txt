April 18, 2021.

This Swarm Ion Temperature Estimation (SITE) software package is designed to generate ion temperatures along the orbits of Swarm satellite. It was developed as a pat of ESA Swarm DISC project at the University of Calgary (PI: Johnathan Burchill, project scientist: Levan Lomidze).

The temperature estimation is based on ion heat balance equation and requires inputs of various neutral and plasma parameters.

The electron densities and temperatures are taken from Swarm Langmuir probe measurements, neutral composition and temperature from NRLMSISE-00 model, and neutral winds from HWM14 model. The model also requires ion drift velocities at high latitudes to calculate ion frictional heating. These velocities have two sources: Swarm Thermal Ion Imager measurements and Weimer 2005 empirical model, and consequently the SITE software generates two different versions of ion temperatures at high latitudes.

The software is written in MATLAB (tested to work on 2018 and later versions) and can be run on Windows (tested on Windows 10). To generate MSIS and HWM executable user needs a Fortran compiler (e.g. Gfortran). For run Weimer 2005 model it is also required to have IDL.

The following sourcecode was developed by UCalgary under contract to DTU and is released under an MIT license:

function_adjust_swarm_ne_te.m
function_geomagn_indexes.m
function_get_tii_data.m
function_hwm14_matlab_3hrAP.m
function_msis00_matlab.m
function_position_magn_fordate_B.m
function_prepare_geomagn_indexes_celest.m
function_Ti_SITE.m
main_generate_SITE_Ti_ver0101.m
hwm14_matlab_driver.f90
create_new_omi_nc_for_weimer05.m
nrlmsise00_driver_fileIO.for
E_weimer05.pro
tilt_gen.pro.

The SITE processor has the following external dependencies, which can be obtained as follows under their own respective software licenses:

The NRLMSISE-00 FORTRAN source is available at https://ccmc.gsfc.nasa.gov/pub/modelweb/atmospheric/msis/nrlmsise00/. 
The HWM14 FORTRAN source code is available at http://onlinelibrary.wiley.com/doi/10.1002/2014EA000089/full (supplemental information: ess224-sup-0002-supinfo.tgz). 
The Weimer 2005 IDL source code is available at https://zenodo.org/record/2530324.

SITE input data:

Swarm LP data are available at ftp://swarm-diss.eo.esa.int/Advanced/Plasma_Data/2_Hz_Langmuir_Probe_Extended_Dataset/.

Swarm TII data are available at ftp://swarm-diss.eo.esa.int/Advanced/Plasma_Data/2Hz_TII_Cross-track_Dataset/New_baseline/. 

The empirical models used in SITE also require various indexes, such as solar, geomagnetic, and solar wind parameters.  They are included in this package until certain date. User can obtain these indexes for newer dates at https://celestrak.com/SpaceData/sw19571001.txt and https://spdf.gsfc.nasa.gov/pub/data/omni/high_res_omni/, and update data files ‘celest.txt’ and ‘omni_1min_2013_2020_25trailavg_20min_lag.nc’, respectively (see below).

To calculate Altitude adjusted corrected geomagnetic (AACGM) coordinates, the model uses dynamic link library (.dll) file created by A. Kouznetsov (University of Calgary).

Questions about the SITE model can be sent to levan.lomidze@ucalgary.ca.

Below are detailed steps interested user needs to undertake to set up and run the SITE code and generate ion temperatures.

1.	Check and if necessary update solar/geomagnetic indexes in the file ‘celest.txt’ (see above) for the date of interest. In the provided file the date range between Jan/01/2013 - Jan/31/2021.

2.	Update ‘omni_1min_2013_2020_25trailavg_20min_lag.nc’ file if needed. The provided files covers data between 2013-2020. It contains 1 min resolution data of IMF and solar wind parameters, and is needed for Weimer 2005 model. To update, obtain yearly files from the link above, copy them in folder ‘IMF_and_SW_data_for_weimer05’ and run ‘create_new_omi_nc_for_weimer05.m’ in MATLAB.

3.	Download Swarm LP and TII data (from the link above) and copy data into corresponding folders ‘LP_data’ and ‘TII_data’. This software comes with Swarm-A data for 25-Jun-2014 to run the test case.

4.	Download NRLMSISE-00 main subroutine file (‘nrlmsise00_sub.for’) and copy it to ‘nrlmsise00_matlab’ directory. Run 
“gfortran nrlmsise00_sub.for nrlmsise00_driver_fileIO.for -w -o NRLMSIS00.exe” 
on windows terminal (or similar). Note, Gfortran needs to be installed first. One option is MINGW64 (GCC compiler on Windows systems).

5.	Download and extract HWM14 code to ‘hwm14_matlab’ directory. Run the following command to generate executable:
“gfortran hwm14.f90 hwm14_matlab_driver.f90 -o hwm14_3hrAP.exe” (see note above).

6.	To use Weimer 2005 model in SITE active IDL license is required. SITE allows to generate ion temperature only with Swarm TII data (see corresponding setting in ‘main_generate_SITE_Ti_ver0101.m’). In this case the W05-based estimates will not include ion frictional heating at high-latitudes. If not using Weimer 2005 skip steps 7-9.

7.	Obtain Weimer2005 IDL code and copy content to folder ‘weimer05’.

8.	Update file ‘tilt_2013-2023.txt’ if needed. It contains hourly values of geomagnetic dipole tilt angle, which is the angle of the north magnetic pole to the GSM z axis (positive when the north magnetic pole is tilted towards the Sun). The file provided contains data between Jan/1/2013-Dec/31/2023. To generate new file run ‘tilt_gen.pro’ in IDL (note, update path in the IDL file as required).

9.	Update path in ‘E_weimer05.pro’ file.

10.	In MATLAB run ‘main_generate_SITE_Ti_ver0101.m’ to generate daily CDF file containing estimated ion temperatures and other parameters defined in SITE product definition document (SW-TN-UoC-GS-001_2-3_SITE_Product_Definition). Lines 13-18 specify satellite (A, B, or C), start date (d1), end date (d2) of analysis, software version, and option to include/exclude Weimer 2005 model.

11.	Included LP and TII  data generates file “example_SW_OPER_EFIATIE_2__20140725T000000_20140725T235959_0101.cdf”, which can be used to test the SITE model.


