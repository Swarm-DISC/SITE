% This software calculates ion temperatue along the orbit of
% Swarm satellite

clear variables;
clear global;
fclose all;
clc;

%----user can make changes here---------
%---------------------------------------
%---------------------------------------
%---------------------------------------
main_path='M:\SITE_software'; % path to software package (change accordingly)
sat='A'; % specify Swarm satellite
d1=datenum(2014,07,25); % start date (year, month, day)
d2=datenum(2014,07,25); % end date
vers='0101'; % SITE Ti dataset version
include_weimer=true; % 'false' if Weimer 2005 is not used (i.e. no frictional heating in 'Ti_model_drift')
%---------------------------------------
%---------------------------------------
%---------------------------------------
%---------------------------------------
%---------------------------------------






%###########################################
%###########################################
%###########################################
%###########################################
cd (main_path)

dllPath = fullfile(main_path,'DB_Model_Signals1900.dll'); % dll for magnetic coordinate calculation
asm1=NET.addAssembly(dllPath);

% load solar/geomagnetic indexes
[ap_date,ap,Ap,f10y_array,f10yb_array]=function_prepare_geomagn_indexes_celest(pwd);

if include_weimer==true
cd ([main_path,'\weimer05'])
file_omni='omni_1min_2013_2021_25trailavg_20min_lag.nc';
swvel0=ncread(file_omni,'swvel');
swden0=ncread(file_omni,'swden');
by0_imf=ncread(file_omni,'by');
bz0_imf=ncread(file_omni,'bz');
datte0_omni=ncread(file_omni,'date');

% load pre-calculated tilt angle info (angle bettween GSM and dipole)
a=load('tilt_2013-2023.txt');
yr_t=a(:,1);
mo_t=a(:,2);
dd_t=a(:,3);
hr_t=a(:,4);
tilt=a(:,5);
dat_tilt=datenum(yr_t,mo_t,dd_t,hr_t, 0, 0);
clear a
end

cd (main_path)

 for ij=1:(d2-d1)+1
     clearvars -except d1 d2 sat vers main_path ij ap Ap ap_date asm1 bx0_imf by0_imf bz0_imf dat_tilt ...
         datte0_omni dd_t dllpath f10y_array f10yb_array file_omni hr_t mo_t path swden0 swvel0 tilt yr_t include_weimer
        ii=d1+ij-1;
        datt=datestr(ii,'yyyymmdd');


% - Read LP data
cd ([main_path,'\LP_data\Sat_',sat])
s3zp=['SW_EXTD_EFI',sat,'_LP_HM_',datt,'*.CDF.ZIP'];
    f=dir(s3zp);
    if isempty(f)~=1 % check if file exists
    fname_zip=f(end).name;    
    mkdir(fname_zip(1:end-4))
    unzip(fname_zip,fname_zip(1:end-4))
    cd(fname_zip(1:end-4))
    s3='*.CDF';
    fname_cdf=dir(s3);
    disp(['working on ',fname_cdf.name])
    
    % read specific variables----
    data=cdfread(fname_cdf.name,'Variable',{'Flagbits','Height','Latitude','Longitude','n','Te_hgn','Te_lgn','Timestamp','Radius'});
    Flagbits=data(:,1);
    Height=data(:,2);
    gLatitude=data(:,3);
    gLongitude=data(:,4);
    Ne=data(:,5);
    Te_hgn=data(:,6);
    Te_lgn=data(:,7);
    Timestamp=data(:,8);
    Radius = data(:,9);
    clear data
        
    nx=length(gLatitude);
    Ne_flag=ones(nx,1);
    Te_lgn_flag=ones(nx,1);
    Te_hgn_flag=ones(nx,1);
    
    yr=NaN(nx,1);
    mo=NaN(nx,1);
    d=NaN(nx,1);
    hr=NaN(nx,1);
    mn=NaN(nx,1);
    sc=NaN(nx,1);
    for j=1:nx
    [x1,x2,x3,x4,x5,x6]=datevec(todatenum(Timestamp{j}));
    yr(j,1)=x1;
    mo(j,1)=x2;
    d(j,1)=x3;
    hr(j,1)=x4;
    mn(j,1)=x5;
    sc(j,1)=x6;
    end
    
    msec_of_day=(hr*3600+mn*60+sc)*1000;
    sec_of_day=(hr*3600+mn*60+sc);
    
    f0=ones(nx,1)-99990;
    f1=ones(nx,1)-99990;
    f2=ones(nx,1)-99990;
    f3=ones(nx,1)-99990;
    f4=ones(nx,1)-99990;
    f5=ones(nx,1)-99990;
    f6=ones(nx,1)-99990;
    f7=ones(nx,1)-99990;
    f8=ones(nx,1)-99990;
    f13=ones(nx,1)-99990;
    f14=ones(nx,1)-99990;
    f15=ones(nx,1)-99990;
    f16=ones(nx,1)-99990;
    f21=ones(nx,1)-99990;
    
    flag=decimalToBinaryVector(cell2mat(Flagbits),24,'LSBFirst');
    
    f0(:,1)=flag(:,1);
    f1(:,1)=flag(:,2);
    f2(:,1)=flag(:,3);
    f3(:,1)=flag(:,4);
    f4(:,1)=flag(:,5);
    f5(:,1)=flag(:,6);
    f6(:,1)=flag(:,7);
    f7(:,1)=flag(:,8);
    f8(:,1)=flag(:,9);
    f13(:,1)=flag(:,14);
    f14(:,1)=flag(:,15);
    f15(:,1)=flag(:,16);
    f16(:,1)=flag(:,17);
    f21(:,1)=flag(:,22);
    
    Ne=cell2mat(Ne);
    Te_hgn=cell2mat(Te_hgn);
    Te_lgn=cell2mat(Te_lgn);
    
    gLatitude=cell2mat(gLatitude);
    gLongitude=cell2mat(gLongitude);
    Height=cell2mat(Height);
    Radius=cell2mat(Radius);

%---------------------------------------------------------        
    Ne_flag(f0==0 & f1==0)=4000; % Ne from low gain    
    Ne_flag(f0==1 & f1==1)=1; 
    Ne_flag(f0==0 & f1==1)=1; % normal configuration (pr1 high, pr2 low)
    Ne_flag(f0==1 & f1==0)=1; 
    Ne_flag(f8==1)=-99999; % don't use 
    Ne_flag(Ne>10^7)=-99998; % don't use 
    Ne_flag(Ne<0)=-99997; % don't use 
    
    Ne_flag(f21==1)=4000; % Ne from low gain
    
    
    Ne(strcmp(f8,'1')==1)   =-10^38;
    
    Ne(Ne_flag<0 | Ne<0)=NaN; % ignores all bad/unusable data
%---------------------------------------------------------

    Te_hgn_flag(f0==0 & f1==1)= 1; % normal configuration (pr1 high, pr2 low)
    Te_hgn_flag(f0==1 & f1==0)= 4000; % pr1 in low-gain, pr2 in high-gain mode
    Te_hgn_flag(f0==0 & f1==1 & f2==1)= 20; % normal configuration (pr1 high, pr2 low) with Overflow at linear bias
    Te_hgn_flag(f0==1 & f1==0 & f2==1)= 4020; % pr1 in low-gain, pr2 in high-gain mode, with Overflow at linear bias
    
    Te_hgn_flag(f0==0 & f1==0)=-99999; % don't use
    Te_hgn_flag(f0==1 & f1==1)=-99999; % don't use

    Te_hgn_flag(f4==1)  =-99999; % don't use
    Te_hgn_flag(f6==1)  =-99999; % don't use
    Te_hgn_flag(f8==1)  =-99999; % don't use
    Te_hgn_flag(f13==1) =-99999; % don't use
    Te_hgn_flag(f15==1) =-99999; % don't use
    Te_hgn_flag(Te_hgn>2*10^4)=-99998; % don't use
    Te_hgn_flag(Te_hgn<0)=-99997; % don't use
    
    
    Te_hgn(f0==0 & f1==0)=-10^38;
    Te_hgn(f0==1 & f1==1)=-10^38;
    Te_hgn(f4==1)  =-10^38;
    Te_hgn(f6==1)  =-10^38;
    Te_hgn(f8==1)  =-10^38;
    Te_hgn(f13==1) =-10^38;
    Te_hgn(f15==1) =-10^38;
    
    Te_hgn(Te_hgn_flag<0 | Te_hgn<0)=NaN; %ignores bad/unusable data
    Te_hgn(Te_hgn_flag==20 | Te_hgn_flag==4020)=NaN; %ignores overflows
%--------------------------------------------------------- 

    Te_lgn_flag(f0==0 & f1==0)=-99999; % don't use
    Te_lgn_flag(f0==1 & f1==1)=-99999; % don't use
    Te_lgn_flag(f0==0 & f1==1)= 1; % normal configuration (pr1 high, pr2 low)
    Te_lgn_flag(f0==1 & f1==0)= 4000; %  pr1 in low-gain, pr2 in high-gain mode
    Te_lgn_flag(f0==0 & f1==1 & f3==1)= 20; % normal configuration (pr1 high, pr2 low) with Overflow at linear bias
    Te_lgn_flag(f0==1 & f1==0 & f3==1)= 4020; %  pr1 in low-gain, pr2 in high-gain mode with Overflow at linear bias
    
    Te_lgn_flag(f5==1)  =-99999; % don't use
    Te_lgn_flag(f7==1)  =-99999; % don't use
    Te_lgn_flag(f8==1)  =-99999; % don't use
    Te_lgn_flag(f14==1) =-99999; % don't use
    Te_lgn_flag(f16==1) =-99999; % don't use
    Te_lgn_flag(Te_lgn>2*10^4)=-99998; % don't use
    Te_lgn_flag(Te_lgn<0)=-99997; % don't use
    
    Te_lgn(f0==0 & f1==0)=-10^38;
    Te_lgn(f0==1 & f1==1)=-10^38;
    Te_lgn(f5==1)  =-10^38;
    Te_lgn(f7==1)  =-10^38;
    Te_lgn(f8==1)  =-10^38;
    Te_lgn(f14==1) =-10^38;
    Te_lgn(f16==1) =-10^38;
    
    Te_lgn(Te_lgn_flag<0 | Te_lgn<0)=NaN; %ignores bad data
    Te_lgn(Te_lgn_flag==20 | Te_lgn_flag==4020)=NaN; %ignores overflows
    %----------------------------------------------------

    %---end Read LP data-----
    cd (main_path)
        
    doy=datenum(yr,mo,d)-datenum(yr-1,12,31);
    datte=yr*1000+doy+(hr+mn/60+sc/3600)/24;
    doy2000=datenum(yr, mo, d)-datenum(2000,01,01);
    if length(unique(doy2000))>2
    disp('error in doy2000')
    stop
    end
    doy2000=unique(doy2000);
   
%----get solar wind data------------- 
if include_weimer==true
k=find(datte0_omni>=min(datte) & datte0_omni<=max(datte));
swvel=swvel0(k);
swden=swden0(k);
by_imf=by0_imf(k);
bz_imf=bz0_imf(k);
datte_omni=datte0_omni(k);

run_weimer=0;
if length(k)>1
run_weimer=1;    
%interpolate
by_imf=interp1(datte_omni,by_imf,datte,'spline','extrap');
bz_imf=interp1(datte_omni,bz_imf,datte,'spline','extrap');
swvel=interp1(datte_omni,swvel,datte,'spline','extrap');
swden=interp1(datte_omni,swden,datte,'spline','extrap');
end
end
        
%--------------------------------------------------------- 
doy2000_file=unique(doy2000);
    if length(doy2000_file)~=1
        disp('ERROR IN DATES')
        stop
    end 
 %-------get solar/geomagn indexes-------

    [ap_3hr,ap_avg,Ap_daily,f107_obs,f107a_obs,f107_adj,f107a_adj]=....
        function_geomagn_indexes(yr,mo,d,hr,mn,sc,24,pwd,ap_date,ap,Ap,f10y_array,f10yb_array);

    pF107=(f107_obs+f107a_obs)/2;
    
    time_prev_day=datevec(datenum(yr,mo,d,hr,mn,sc)-1);
    
        [~,~,~,f107_obs_prev_day,~,~,~]=....
        function_geomagn_indexes(time_prev_day(:,1),time_prev_day(:,2),time_prev_day(:,3),...
        time_prev_day(:,4),time_prev_day(:,5),time_prev_day(:,6),...
        24,pwd,ap_date,ap,Ap,f10y_array,f10yb_array);
    
    %%get magn data %%%%%%%
    clear mlt qdlat qdlon

[mlat_aacgm, ~, mlt_aacgm, qdlat, qdlon, D, I, Bnn, Bee, Bcc, errMSG]=...
    function_position_magn_fordate_B(doy2000_file, sec_of_day, Height, gLatitude, gLongitude, main_path);

if length(qdlat)~=length(gLatitude)
    disp('lengths do not match')
    stop
end

if include_weimer==true
if run_weimer==1
%--generate input file for weimer05 model
cd ([main_path,'\weimer05'])
nn=length(gLatitude);
fspec=[repmat('%5.2f \t ',1,7),repmat('%d \t ',1,4),'%5.2f \n'];
fid=fopen('imf.txt','wt');
fprintf(fid,'%d\n',nn);
for i=1:nn
    tlt=tilt(yr(i)==yr_t & mo(i)==mo_t & d(i)==dd_t & hr(i)==hr_t);
    
fprintf(fid,fspec,mlat_aacgm(i),mlt_aacgm(i),Height(i),by_imf(i),bz_imf(i),swvel(i),swden(i),...
    yr(i), mo(i), d(i), hr(i), tlt);
    clear tlt
end
fclose(fid);
% run weimer05 (IDL)
disp('running W05')
[~,~]=system (['idl -e ".RUN ',main_path,'\weimer05\E_weimer05.pro" ']);
delete 'imf.txt'

%get E field (magn north and east components)
a=load('Efield_out_with_pot.txt');
Emn=a(:,1)*10^-3; % in V/m
Eme=a(:,2)*10^-3;
Elpot=a(:,3); % electric potential, kV
clear a
delete 'Efield_out_with_pot.txt'

Btot=((Bnn.^2+Bee.^2+Bcc.^2).^0.5)*10^-9; % in Tesla
Bn=Btot.*cos(I).*cos(D);
Be=Btot.*cos(I).*sin(D);
Bc=Btot.*sin(I);


vime= -( Emn.*sin(I) )./ Btot; % magn east
vimn=  ( Eme.*sin(I) )./ Btot; % magn north
vic= - ( Eme.*cos(I) )./ Btot; % vertical, towards the earth

vim_ge=vimn.*sin(D) + vime.*cos(D); % Geographic east, Weimer05 drift
vim_gn=vimn.*cos(D) - vime.*sin(D); % Geographic north, Weimer05 drift
vim_z = vic; % vertical, towards the earth
end
end

cd (main_path)

% get MSIS data
[H,O,N2,O2,HE,Tn,rho,N]=function_msis00_matlab(doy, sec_of_day, gLatitude, gLongitude, Height, Ap_daily, f107_obs_prev_day, f107a_obs, pwd);

%get HWM data
[u_merid,u_zon]=function_hwm14_matlab_3hrAP(yr, doy, sec_of_day, gLatitude, gLongitude, Height, Ap_daily, ap_3hr, f107_obs, f107a_obs, pwd);
disp('MSIS/HWM runs complete')

%--------------------
% read CT data
clear i_tii
cd ([main_path,'\TII_data\Sat_',sat])
s1=['SW_EXPT_EFI',sat,'_TCT02_',datt,'*.ZIP'];
    ftii=dir(s1);
    i_tii=0; % 0 if no tii data exists
    if isempty(ftii)~=1 % check if file exists
    i_tii=1; % 1 if tii data exists
    vi_tii_gn=[];
    vi_tii_ge=[];
    vi_tii_z=[];
    hr_ct=[];
    mn_ct=[];
    sc_ct=[];
    
    for it=1:length(ftii) % neede because multiple TII files on a given day are present
    [vi_tii_gn_temp, vi_tii_ge_temp, vi_tii_z_temp,hr_ct_temp,mn_ct_temp,sc_ct_temp]=function_get_tii_data(ftii(it),sat);
    vi_tii_gn=[vi_tii_gn;vi_tii_gn_temp];
    vi_tii_ge=[vi_tii_ge;vi_tii_ge_temp];
    vi_tii_z =[vi_tii_z;vi_tii_z_temp];
    hr_ct=[hr_ct; hr_ct_temp];
    mn_ct=[mn_ct; mn_ct_temp];
    sc_ct=[sc_ct; sc_ct_temp];
    cd ([main_path,'\TII_data\Sat_',sat])    
    end
   
    end % check if TII data exists


%--------------------------------------------
cd (main_path)
% adjust Ne/Te
[ne_adj,~,teh_adj,~,tel_adj,~]=function_adjust_swarm_ne_te(Ne,Te_hgn,Te_lgn,sat);

% if Te values after adjustement get <Tn, set them to NaN
teh_adj(teh_adj<Tn)=NaN;
tel_adj(tel_adj<Tn)=NaN;

% create Ti flags
Flag_ti_meas=zeros(nx,3);
Flag_ti_model=zeros(nx,3);

% find missing high-gain Te and replace with low-gain Te, record this in Ti flag
clear k k1 k2
k1=isnan(teh_adj);
k2=isnan(tel_adj);
k=find(k1==1 & k2~=1);
teh_adj(k)=tel_adj(k);
Flag_ti_meas(k,2)=1;
Flag_ti_model(k,2)=1;

vi_tii_GN=zeros(nx,1);
vi_tii_GE=zeros(nx,1);
vi_tii_Z=zeros(nx,1);

% make timing of TII and LP data consistent (there is ~0.03 sec shift)
if i_tii==1 % do when TII data exist, otherwise zeros 

temp = fix(10*(hr*3600 + mn * 60 + sc));
temp_ct = fix(10*(hr_ct*3600 + mn_ct * 60 + sc_ct - 0.03));
 [~,ix,ix_ct]=intersect(temp,temp_ct, 'stable'); % indexes of common elements
    vi_tii_GN(ix)=vi_tii_gn(ix_ct);
    vi_tii_GE(ix)=vi_tii_ge(ix_ct);
    vi_tii_Z(ix)=vi_tii_z(ix_ct);
       
    if length(ix)<2 || length(ix_ct)<2
        disp('no mathching LP and TII measurements found')
    end
    
end

% find when TII drift is missing and set wind to zero
% (sets frictional heating to zero)
k=find(abs(qdlat)>=44 & abs(qdlat)<=60 & vi_tii_GE==0);
u_zon(k)=0;
k=find(abs(qdlat)>=44 & abs(qdlat)<=60 & vi_tii_GN==0); 
u_merid(k)=0;

dvx=u_zon-vi_tii_GE;
dvy=u_merid-vi_tii_GN;
dvz=vi_tii_Z;

% no frictional heating at mid/low latitudes because no drifts are expected to be available
dvx(abs(qdlat)<44)=0;
dvy(abs(qdlat)<44)=0;
dvz(abs(qdlat)<44)=0;


% identify when high-lat drifts are zero (happens outside W05 lower boundary) and set neutral winds
% also sets winds to zero to eliminate frictional heating;
if include_weimer==true
if run_weimer==1
u_zon_mod=u_zon;
u_merid_mod=u_merid;
k=find(abs(qdlat)>=44 & abs(qdlat)<=60 & vim_ge==0); 
u_zon_mod(k)=0;
k=find(abs(qdlat)>=44 & abs(qdlat)<=60 & vim_gn==0); 
u_merid_mod(k)=0;

dvx_m=u_zon_mod-vim_ge;
dvy_m=u_merid_mod-vim_gn;
dvz_m=vim_z;

dvx_m(abs(qdlat)<44)=0;
dvy_m(abs(qdlat)<44)=0;
dvz_m(abs(qdlat)<44)=0;
end


if run_weimer==0
dvx_m=u_zon.*0;
dvy_m=u_merid.*0;
dvz_m=u_zon.*0;
vim_ge=u_zon.*0;
vim_gn=u_zon.*0;
vim_z=u_zon.*0;
end

else
dvx_m=u_zon.*0;
dvy_m=u_merid.*0;
dvz_m=u_zon.*0;
vim_ge=u_zon.*0;
vim_gn=u_zon.*0;
vim_z=u_zon.*0;    
end


% Generate ion temperature
% TII-based

[Ti]=function_Ti_SITE(H, O, N2, O2, HE, Tn, ne_adj, teh_adj, dvx, dvy, dvz);

Ti(Ti>teh_adj & abs(qdlat)<44)=NaN;
Ti(Ti<Tn)=NaN;

% Weimer based
[Ti_m]=function_Ti_SITE(H, O, N2, O2, HE, Tn, ne_adj, teh_adj, dvx_m, dvy_m, dvz_m);

Ti_m(Ti_m>teh_adj & abs(qdlat)<44)=NaN;
Ti_m(Ti_m<Tn)=NaN;

%--------------------------------------------


% identify when high-lat frict heating is included:
k=find(abs(qdlat)>=44 & (vi_tii_GE~=0 | vi_tii_GN~=0 | vi_tii_Z~=0)); 
Flag_ti_meas(k,1)=1;

k=find(abs(qdlat)>=44 & (vim_ge~=0 | vim_gn~=0 | vim_z~=0)); 
Flag_ti_model(k,1)=1;

% identify if ion temperature data is available
Flag_ti_meas(:,3)=1;
Flag_ti_model(:,3)=1;

clear k k1 k2
k=isnan(Ti);
Ti(k)=-9999;
Flag_ti_meas(k,3)=0;

clear k
k=isnan(Ti_m);
Ti_m(k)=-9999;
Flag_ti_model(k,3)=0;

% set NaN Te to -9999
clear k
k=isnan(teh_adj);
teh_adj(k)=-9999;


flag_ti_meas=bin2dec(num2str(Flag_ti_meas(:,:)));
flag_ti_model=bin2dec(num2str(Flag_ti_model(:,:)));


datatime = double(86400*1000*datenum(yr,mo,d,hr,mn,sc))-86400*1000; % in msec and 1 day back to be consistent with cdf epoch definition

% Create a CDF file
fname_out_cdf=['SW_OPER_EFI',sat,'TIE_2__',datt,'T000000_',datt,'T235959_',vers,'.cdf'];
cdfid = cdflib.create(fname_out_cdf);

% Create Variables 
time_id = cdflib.createVar(cdfid,'Timestamp','CDF_EPOCH',1,[],true,[]);
lat_id = cdflib.createVar(cdfid,'Latitude','CDF_REAL8',1,[],true,[]);
lon_id = cdflib.createVar(cdfid,'Longitude','CDF_REAL8',1,[],true,[]);
ht_id = cdflib.createVar(cdfid,'Height','CDF_REAL8',1,[],true,[]);
rad_id = cdflib.createVar(cdfid,'Radius','CDF_REAL8',1,[],true,[]);
qd_id = cdflib.createVar(cdfid,'QDLatitude','CDF_REAL8',1,[],true,[]);
mlt_id = cdflib.createVar(cdfid,'MLT','CDF_REAL8',1,[],true,[]);
Tn_id = cdflib.createVar(cdfid,'Tn_msis','CDF_REAL8',1,[],true,[]);
Te_id = cdflib.createVar(cdfid,'Te_adj_LP','CDF_REAL8',1,[],true,[]);
Ti_mes_id = cdflib.createVar(cdfid,'Ti_meas_drift','CDF_REAL8',1,[],true,[]);
Ti_mod_id = cdflib.createVar(cdfid,'Ti_model_drift','CDF_REAL8',1,[],true,[]);
flag_ti_ms_id = cdflib.createVar(cdfid,'Flag_ti_meas','CDF_UINT1',1,[],true,[]);
flag_ti_mod_id = cdflib.createVar(cdfid,'Flag_ti_model','CDF_UINT1',1,[],true,[]);

%set variable compression
for i_comp=1:13
cdflib.setVarCompression(cdfid,i_comp-1,'GZIP_COMPRESSION',4)
end


% Write data to other variables.
for i=1:nx
cdflib.putVarRecordData(cdfid,time_id,i-1,datatime(i));
cdflib.putVarRecordData(cdfid,lat_id,i-1, round(double(gLatitude(i)),13) );
cdflib.putVarRecordData(cdfid,lon_id,i-1, round(double(gLongitude(i)),13) );
cdflib.putVarRecordData(cdfid,rad_id,i-1, round(double(1000*Radius(i)),13) ); % in meters
cdflib.putVarRecordData(cdfid,ht_id,i-1, round(double(1000*Height(i)),13) ); % in meters
cdflib.putVarRecordData(cdfid,qd_id,i-1, round(double(qdlat(i)),13) );
cdflib.putVarRecordData(cdfid,mlt_id,i-1, round(double(mlt_aacgm(i)),13) );
cdflib.putVarRecordData(cdfid,Tn_id,i-1, round(double(Tn(i)),2) );
cdflib.putVarRecordData(cdfid,Te_id,i-1, round(double(teh_adj(i)),2) );
cdflib.putVarRecordData(cdfid,Ti_mes_id,i-1, round(double(Ti(i)),2) );
cdflib.putVarRecordData(cdfid,Ti_mod_id,i-1, round(double(Ti_m(i)),2) );
cdflib.putVarRecordData(cdfid,flag_ti_ms_id,i-1,uint8(flag_ti_meas(i)));
cdflib.putVarRecordData(cdfid,flag_ti_mod_id,i-1,uint8(flag_ti_model(i)));
end


% Write to Global Attribute
titleAttrNum = cdflib.createAttr(cdfid,'TITLE','global_scope');
% Write values to entries in the global attribute.
cdflib.putAttrEntry(cdfid,titleAttrNum,0,'CDF_CHAR',['Swarm ',sat,' Ion Temperature Model Estimate']);

% ----
glattr=cdflib.createAttr(cdfid,'File_Name','global_scope');
cdflib.putAttrEntry(cdfid,glattr,0,'CDF_CHAR',fname_out_cdf);

glattr=cdflib.createAttr(cdfid,'Creator','global_scope');
cdflib.putAttrEntry(cdfid,glattr,0,'CDF_CHAR','University of Calgary, Alberta, Canada');

glattr=cdflib.createAttr(cdfid,'Creator_Software','global_scope');
cdflib.putAttrEntry(cdfid,glattr,0,'CDF_CHAR','SITE');

glattr=cdflib.createAttr(cdfid,'Creator_Version','global_scope');
cdflib.putAttrEntry(cdfid,glattr,0,'CDF_CHAR','1.0');

glattr=cdflib.createAttr(cdfid,'Generation_date','global_scope');
gen_date=datetime(datetime,'TimeZone','local','Format','yyyy-MM-dd hh:mm:ss');
gen_date.TimeZone = 'UTC';
gen_date=datestr(gen_date,'yyyy-mm-ddTHH:MM:SS');
cdflib.putAttrEntry(cdfid,glattr,0,'CDF_CHAR',['UTC=',gen_date]);

glattr=cdflib.createAttr(cdfid,'Input_files','global_scope');
cdflib.putAttrEntry(cdfid,glattr,0,'CDF_CHAR',[f.name]);
for it=1:length(ftii)
cdflib.putAttrEntry(cdfid,glattr,it,'CDF_CHAR',ftii(it).name);
end

%-- metainfo for input file ----
fid=fopen([fname_out_cdf(1:end-4),'.metainfo'],'wt');
fprintf(fid,'%s\n','ProcessingCenter:UOC');
fprintf(fid,'%s\n','Processor:UOC_SITE');
fprintf(fid,'%s\n','ProcessorVersion:01.00');
fprintf(fid,'%s\n','ProductError:0');
fprintf(fid,'%s\n',['Input:',f.name]);
for it=1:length(ftii)
fprintf(fid,'%s\n',['Input:',ftii(it).name]);
end
fclose(fid);
%--------------------------

glattr=cdflib.createAttr(cdfid,'Models_Used','global_scope');
cdflib.putAttrEntry(cdfid,glattr,0,'CDF_CHAR','NRLMSISE-00');
cdflib.putAttrEntry(cdfid,glattr,1,'CDF_CHAR','HWM14');
cdflib.putAttrEntry(cdfid,glattr,2,'CDF_CHAR','Weimer 2005');
% -----

glattr=cdflib.createAttr(cdfid,'PI_name','global_scope');
cdflib.putAttrEntry(cdfid,glattr,0,'CDF_CHAR','Johnathan Burchill');

glattr=cdflib.createAttr(cdfid,'Project_scientist','global_scope');
cdflib.putAttrEntry(cdfid,glattr,0,'CDF_CHAR','Levan Lomidze');

glattr=cdflib.createAttr(cdfid,'PI_affiliation','global_scope');
cdflib.putAttrEntry(cdfid,glattr,0,'CDF_CHAR','UCalgary');

glattr=cdflib.createAttr(cdfid,'Mission_group','global_scope');
cdflib.putAttrEntry(cdfid,glattr,0,'CDF_CHAR','Swarm');

glattr=cdflib.createAttr(cdfid,'Acknowledgement','global_scope');
cdflib.putAttrEntry(cdfid,glattr,0,'CDF_CHAR','ESA Swarm EFI LP and TII data are available from https://swarm-diss.eo.esa.int');

glattr=cdflib.createAttr(cdfid,'Data_version','global_scope');
cdflib.putAttrEntry(cdfid,glattr,0,'CDF_CHAR',vers);

glattr=cdflib.createAttr(cdfid,'Discipline','global_scope');
cdflib.putAttrEntry(cdfid,glattr,0,'CDF_CHAR','Space Physics>Ionospheric Science');

glattr=cdflib.createAttr(cdfid,'Project','global_scope');
cdflib.putAttrEntry(cdfid,glattr,0,'CDF_CHAR','ESA Living Planet Programme');

glattr=cdflib.createAttr(cdfid,'TEXT','global_scope');
cdflib.putAttrEntry(cdfid,glattr,0,'CDF_CHAR',['Knudsen, D.J., Burchill, J.K., Buchert, S.C., Eriksson, A.I., Gill, R., Wahlund, J.E., Ahlen, L., Smith, M. and Moffat, B., 2017. Thermal ion imagers and Langmuir probes in the Swarm electric field instruments. Journal of Geophysical Research: Space Physics, 122(2), pp.2655-2673.']);

glattr=cdflib.createAttr(cdfid,'Time_resolution','global_scope');
cdflib.putAttrEntry(cdfid,glattr,0,'CDF_CHAR','0.5 seconds');

glattr=cdflib.createAttr(cdfid,'MODS','global_scope');
cdflib.putAttrEntry(cdfid,glattr,0,'CDF_CHAR','Draft release');

% Write to Attributes Associated with Variables
% Create attributes associated with variables in the CDF file.

fieldAttrNum = cdflib.createAttr(cdfid,'FIELDNAM','variable_scope');
unitsAttrNum = cdflib.createAttr(cdfid,'UNITS','variable_scope');
catdescAttrNum = cdflib.createAttr(cdfid,'CATDESC','variable_scope');
formatAttrNum = cdflib.createAttr(cdfid,'FORMAT','variable_scope');

% Write to attributes of the Time variable.
cdflib.putAttrEntry(cdfid,fieldAttrNum,time_id,...
    'CDF_CHAR','Time of observation');
cdflib.putAttrEntry(cdfid,unitsAttrNum,time_id,...
    'CDF_CHAR',' ');
cdflib.putAttrEntry(cdfid,catdescAttrNum,time_id,...
    'CDF_CHAR','UT');
cdflib.putAttrEntry(cdfid,formatAttrNum,time_id,...
    'CDF_CHAR','I14');

% Write to attributes of other variables.
% Lat
cdflib.putAttrEntry(cdfid,fieldAttrNum,lat_id,...
    'CDF_CHAR','Latitude');
cdflib.putAttrEntry(cdfid,unitsAttrNum,lat_id,...
    'CDF_CHAR','degrees');
cdflib.putAttrEntry(cdfid,catdescAttrNum,lat_id,...
    'CDF_CHAR','Geodetic latitude');
cdflib.putAttrEntry(cdfid,formatAttrNum,lat_id,...
    'CDF_CHAR','F17.13');

% Long
cdflib.putAttrEntry(cdfid,fieldAttrNum,lon_id,...
    'CDF_CHAR','Longitude');
cdflib.putAttrEntry(cdfid,unitsAttrNum,lon_id,...
    'CDF_CHAR','degrees');
cdflib.putAttrEntry(cdfid,catdescAttrNum,lon_id,...
    'CDF_CHAR','Geodetic longitude');
cdflib.putAttrEntry(cdfid,formatAttrNum,lon_id,...
    'CDF_CHAR','F18.13');

% height
cdflib.putAttrEntry(cdfid,fieldAttrNum,ht_id,...
    'CDF_CHAR','Height');
cdflib.putAttrEntry(cdfid,unitsAttrNum,ht_id,...
    'CDF_CHAR','m');
cdflib.putAttrEntry(cdfid,catdescAttrNum,ht_id,...
    'CDF_CHAR','Height above WGS84 reference ellipsoid');
cdflib.putAttrEntry(cdfid,formatAttrNum,ht_id,...
    'CDF_CHAR','F20.13');

% Radius
cdflib.putAttrEntry(cdfid,fieldAttrNum,rad_id,...
    'CDF_CHAR','Radius');
cdflib.putAttrEntry(cdfid,unitsAttrNum,rad_id,...
    'CDF_CHAR','m');
cdflib.putAttrEntry(cdfid,catdescAttrNum,rad_id,...
    'CDF_CHAR','Geocentric radius');
cdflib.putAttrEntry(cdfid,formatAttrNum,rad_id,...
    'CDF_CHAR','F21.13');

% qd lat
cdflib.putAttrEntry(cdfid,fieldAttrNum,qd_id,...
    'CDF_CHAR','QDLatitude');
cdflib.putAttrEntry(cdfid,unitsAttrNum,qd_id,...
    'CDF_CHAR','degrees');
cdflib.putAttrEntry(cdfid,catdescAttrNum,qd_id,...
    'CDF_CHAR','Quasi-dipole magnetic latitude');
cdflib.putAttrEntry(cdfid,formatAttrNum,qd_id,...
    'CDF_CHAR','F17.13');

% MLT
cdflib.putAttrEntry(cdfid,fieldAttrNum,mlt_id,...
    'CDF_CHAR','MLT');
cdflib.putAttrEntry(cdfid,unitsAttrNum,mlt_id,...
    'CDF_CHAR','hour');
cdflib.putAttrEntry(cdfid,catdescAttrNum,mlt_id,...
    'CDF_CHAR','Magnetic local time');
cdflib.putAttrEntry(cdfid,formatAttrNum,mlt_id,...
    'CDF_CHAR','F16.13');

% Tn_msis
cdflib.putAttrEntry(cdfid,fieldAttrNum,Tn_id,...
    'CDF_CHAR','Tn_msis');
cdflib.putAttrEntry(cdfid,unitsAttrNum,Tn_id,...
    'CDF_CHAR','K');
cdflib.putAttrEntry(cdfid,catdescAttrNum,Tn_id,...
    'CDF_CHAR','Neutral temperature from NRLMSISE00 model');
cdflib.putAttrEntry(cdfid,formatAttrNum,Tn_id,...
    'CDF_CHAR','F7.2');

% Te_adj_LP
cdflib.putAttrEntry(cdfid,fieldAttrNum,Te_id,...
    'CDF_CHAR',['Te_adj_LP']);
cdflib.putAttrEntry(cdfid,unitsAttrNum,Te_id,...
    'CDF_CHAR','K');
cdflib.putAttrEntry(cdfid,catdescAttrNum,Te_id,...
    'CDF_CHAR',['Corrected Swarm ',sat,' LP electron temperature']);
cdflib.putAttrEntry(cdfid,formatAttrNum,Te_id,...
    'CDF_CHAR','F8.2');

% Ti_meas_drift
cdflib.putAttrEntry(cdfid,fieldAttrNum,Ti_mes_id,...
    'CDF_CHAR',['Ti_meas_drift']);
cdflib.putAttrEntry(cdfid,unitsAttrNum,Ti_mes_id,...
    'CDF_CHAR','K');
cdflib.putAttrEntry(cdfid,catdescAttrNum,Ti_mes_id,...
    'CDF_CHAR','Ion temperature estimated using Swarm TII drift at high latitudes');
cdflib.putAttrEntry(cdfid,formatAttrNum,Ti_mes_id,...
    'CDF_CHAR','F8.2');

% Ti_model_drift
cdflib.putAttrEntry(cdfid,fieldAttrNum,Ti_mod_id,...
    'CDF_CHAR',['Ti_model_drift']);
cdflib.putAttrEntry(cdfid,unitsAttrNum,Ti_mod_id,...
    'CDF_CHAR','K');
cdflib.putAttrEntry(cdfid,catdescAttrNum,Ti_mod_id,...
    'CDF_CHAR','Ion temperature estimated using Weimer 2005 model drifts at high latitudes');
cdflib.putAttrEntry(cdfid,formatAttrNum,Ti_mod_id,...
    'CDF_CHAR','F8.2');

% Flag_ti_meas
cdflib.putAttrEntry(cdfid,fieldAttrNum,flag_ti_ms_id,...
    'CDF_CHAR',['Flags characterising TII-based Ti']);
cdflib.putAttrEntry(cdfid,unitsAttrNum,flag_ti_ms_id,...
    'CDF_CHAR',' ');
cdflib.putAttrEntry(cdfid,catdescAttrNum,flag_ti_ms_id,...
    'CDF_CHAR',['Bit0 = 1/0: high-latitude frictional heating included/omitted, Bit1=1/0: electron temperature from low/high-gain probe, Bit2=1/0: ion temperature data available/data value is set to -9999']);
cdflib.putAttrEntry(cdfid,formatAttrNum,flag_ti_ms_id,...
    'CDF_CHAR','I2');

% Flag_ti_mod
cdflib.putAttrEntry(cdfid,fieldAttrNum,flag_ti_mod_id,...
    'CDF_CHAR',['Flags characterising model-based Ti']); %Bitwise flag for Ti_model_drift
cdflib.putAttrEntry(cdfid,unitsAttrNum,flag_ti_mod_id,...
    'CDF_CHAR',' ');
cdflib.putAttrEntry(cdfid,catdescAttrNum,flag_ti_mod_id,...
    'CDF_CHAR',['Bit0 = 1/0: high-latitude frictional heating included/omitted, Bit1=1/0: electron temperature from low/high-gain probe, Bit2=1/0: ion temperature data available/data value is set to -9999']);
cdflib.putAttrEntry(cdfid,formatAttrNum,flag_ti_mod_id,...
    'CDF_CHAR','I2');

cdflib.close(cdfid);
disp([datt,' - CDF file generated'])
    end % file date check
 end

