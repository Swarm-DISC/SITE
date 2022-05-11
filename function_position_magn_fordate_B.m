function [mlat, mlon, mlt, qdlat, qdlon, D, I, Bn, Be, Bc, errMSG]=...
    function_position_magn_fordate_B(doy2000, sec_of_day, alt, lat, lon, path)

i_Time_Points_Calculated = Orbit_TLE_Ephemeris.C_Orbit_TLE_Ephemeris.Calc_QD_AACGM_ZA_ML_by_Ephemeris...
    (doy2000, sec_of_day, alt, lat, lon, [path,'\apexsh-1955-2025.dat'], [path,'\IGRF13-1900-2020.COF']);

errMSG='No Errors';
mm=length(sec_of_day);

if i_Time_Points_Calculated~=mm
   errMSG=Orbit_TLE_Ephemeris.C_Orbit_TLE_Ephemeris.sErrMsg_Ephemris_and_IGRF;
   disp(errMSG)
   stop
end

mlat=zeros(mm,1)-999;
mlon=zeros(mm,1)-999;
mlt=zeros(mm,1)-999;
qdlat=zeros(mm,1)-999;
qdlon=zeros(mm,1)-999;
I=zeros(mm,1)-999;
D=zeros(mm,1)-999;
Bn=zeros(mm,1)-999;
Be=zeros(mm,1)-999;
Bc=zeros(mm,1)-999;

temp_mlat=Orbit_TLE_Ephemeris.C_Orbit_TLE_Ephemeris.ad_mGLat_Satellite_m90_90;
temp_mlon=Orbit_TLE_Ephemeris.C_Orbit_TLE_Ephemeris.ad_mGLong_Satellite_m180_180;
temp_mlt=Orbit_TLE_Ephemeris.C_Orbit_TLE_Ephemeris.ad_Magnetic_Local_Time_0_24_hr;

temp_qdlat=Orbit_TLE_Ephemeris.C_Orbit_TLE_Ephemeris.ad_QD_Lat_deg;
temp_qdlon=Orbit_TLE_Ephemeris.C_Orbit_TLE_Ephemeris.ad_QD_Long_deg;

temp_east=Orbit_TLE_Ephemeris.C_Orbit_TLE_Ephemeris.ad_IGRF_EAST_nT;
temp_north=Orbit_TLE_Ephemeris.C_Orbit_TLE_Ephemeris.ad_IGRF_NORTH_nT;
temp_down=Orbit_TLE_Ephemeris.C_Orbit_TLE_Ephemeris.ad_IGRF_DOWN_nT;

for i=1:mm
    
mlat(i,1)=temp_mlat.Get(i-1);
mlon(i,1)=temp_mlon.Get(i-1);
mlt(i,1)=temp_mlt.Get(i-1);

qdlat(i,1)=temp_qdlat.Get(i-1);
qdlon(i,1)=temp_qdlon.Get(i-1);

D(i,1)=atan2(temp_east.Get(i-1),temp_north.Get(i-1));
I(i,1)=atan2(  temp_down.Get(i-1),(temp_north.Get(i-1)^2+temp_east.Get(i-1)^2)^0.5  );

Bn(i,1)=temp_north.Get(i-1);
Be(i,1)=temp_east.Get(i-1);
Bc(i,1)=temp_down.Get(i-1);

end
