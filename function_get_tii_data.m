function[vi_tii_gn,vi_tii_ge,viz,hr_ct,mn_ct,sc_ct]=function_get_tii_data(ftii,sat)

narginchk(2,2)
nargoutchk(6,6);
    
    fnametii_zip=ftii.name;    
    mkdir(fnametii_zip(1:end-4))
    unzip(fnametii_zip,fnametii_zip(1:end-4))

    cd(fnametii_zip(1:end-4))
    s3tii='*.CDF';
    fnametii_cdf=dir(s3tii);
    
    if isempty(fnametii_cdf)==1
      cd(['databases\TCT\',fnametii_zip(end-7:end-4),'\TCT02']) % some TII data need this
          s3tii='*.CDF';
    fnametii_cdf=dir(s3tii);
    end
    
    disp(['working on ',fnametii_cdf.name])
    
    data=cdfread(fnametii_cdf.name,'Variable',{'Timestamp','Quality_flags','Vixh','Vixv','Viy','Viz','VsatN','VsatE'});
    Timestamp=data(:,1);
    Quality_flags=data(:,2);
    Vixh=data(:,3);
    Vixv=data(:,4);
    Viy=data(:,5);
    Viz=data(:,6);
    VsatN=data(:,7);
    VsatE=data(:,8);
    clear data
%----read specific variables (fast)---

 
    ndata=length(Viy);
     
    y_ct=NaN(ndata,1);
    mo_ct=NaN(ndata,1);
    d_ct=NaN(ndata,1);
    hr_ct=NaN(ndata,1);
    mn_ct=NaN(ndata,1);
    sc_ct=NaN(ndata,1);
 for j=1:ndata
    [x1,x2,x3,x4,x5,x6]=datevec(todatenum(Timestamp{j}));
    y_ct(j,1)=x1;
    mo_ct(j,1)=x2;
    d_ct(j,1)=x3;
    hr_ct(j,1)=x4;
    mn_ct(j,1)=x5;
    sc_ct(j,1)=x6;
 end
 
%   qxh=NaN(ndata,1);
%   qxv=NaN(ndata,1);
  qy=NaN(ndata,1);
%   qz=NaN(ndata,1);
  
  flag=decimalToBinaryVector(cell2mat(Quality_flags),4,'LSBFirst');
  
%      qxh(:,1)=flag(:,1);
%      qxv(:,1)=flag(:,2);
     qy(:,1)=flag(:,3);
%      qz(:,1)=flag(:,4);
  

 clear x* v

%  doy_ct=datenum(y_ct,mo_ct,d_ct)-datenum(y_ct-1,12,31);
%  datte_ct=y_ct*1000+doy_ct+(hr_ct+mn_ct/60+sc_ct/3600)/24;
%  doy2000_ct=datenum(y_ct,mo_ct,d_ct)-datenum(2000,01,01);
%  dat_fortilt_ct=datenum(y_ct,mo_ct,d_ct,hr_ct,0,0);
 
vixh=cell2mat(Vixh);
vixv=cell2mat(Vixv);
viy=cell2mat(Viy);
viz=cell2mat(Viz);

%  latitude_ct=cell2mat(Latitude);
%  longitude_ct=cell2mat(Longitude);
%  radius_ct=cell2mat(Radius);
%  QDLatitude_ct=cell2mat(QDLatitude);

 vsatnorth=cell2mat(VsatN);
 vsateast=cell2mat(VsatE);
%  vsatcentre=cell2mat(VsatC);
 
 alfa=atan2(vsateast,vsatnorth);

 
if strcmp(sat,'C')==1 % for Sw-C, all qy is allowed
viy(qy~=1 & qy~=0)=NaN;    
viy(abs(viy)>4000)=NaN;
% set vix and viz to zero untill they are calibrated
vix=viy.*0;
viz=viy.*0;
vix(isnan(viy))=NaN;
viz(isnan(viy))=NaN;
else
    
% take flags into account
% currently only CT flag(qy) matters

%  vixh(qxh~=1)=NaN;
%  vixv(qxv~=1)=NaN;
viy(qy~=1)=NaN;
%  viz(qz~=1)=NaN;

% set AT and Vert to NaN when CT flag=0
vixh(qy~=1)=NaN;
vixv(qy~=1)=NaN;
viz(qy~=1)=NaN;

vix=(vixh+vixv)/2; % average of two sensors, the along-track component

% remove |values| > 4000 m/s
vix(abs(vix)>4000)=NaN;
viy(abs(viy)>4000)=NaN;
viz(abs(viz)>4000)=NaN;

% set vix and viz to zero untill they are calibrated
vix(:)=0;
viz(:)=0;
vix(isnan(viy))=NaN;
viz(isnan(viy))=NaN;
end

vi_tii_gn = vix.*cos(alfa) - viy.*sin(alfa); % geographic N component of TII drift
vi_tii_ge = vix.*sin(alfa) + viy.*cos(alfa); % geographic E component of TII drift

