clear variables;
clear global;
file_new_name='omni_1min_2013_2020_25trailavg_20min_lag.nc'; % specify name here

 file_source='schema.nc';
 
 files = dir('*.asc');
    y=[];
    doy=[];
    d=[];
    hr=[];
    mn=[];
    bx=[];
    by=[];
    bz=[];
    swvel=[];
    swden=[];
    
for i=1:length(files)
    
ax=load(files(i).name);

    yx=ax(:,1);
    doyx=ax(:,2);
    hrx=ax(:,3);
    mnx=ax(:,4);
    bxx=ax(:,15); % GSM
    byx=ax(:,18); % GSM
    bzx=ax(:,19); % GSM
    swvelx=ax(:,22); % GSE
    swdenx=ax(:,26);

y=[y;yx];
doy=[doy;doyx];
hr=[hr;hrx];
mn=[mn;mnx];
bx=[bx;bxx];
by=[by;byx];
bz=[bz;bzx];
swvel=[swvel;swvelx];
swden=[swden;swdenx];

disp(i)

end

 datte=y*1000+doy+(hr+mn/60)/24;

 %----needed to avoid NaNs in the interpolation----
k1=find(bz~=(9999.99));
k2=find(swden~=(999.99));

k=min(k1(end),k2(end));

k=1:k;
datte=datte(k);
bx=bx(k);
by=by(k);
bz=bz(k);
swvel=swvel(k);
swden=swden(k);
clear k
%------------------------------------------------

 ndata=length(swvel);

 k=find(bx~=(9999.99));
 temp_1=bx(k);
 datte_1=datte(k);
 temp_i=interp1(datte_1,temp_1,datte,'linear','extrap');
 bxMask=zeros(ndata,1);
 bxMask(k)=1;
 bxMask=int8(bxMask);
 bx=double(temp_i);
 clear k datte_1 temp_i temp_1
 
 k=find(by~=(9999.99));
 temp_1=by(k);
 datte_1=datte(k);
 temp_i=interp1(datte_1,temp_1,datte,'linear','extrap');
 byMask=zeros(ndata,1);
 byMask(k)=1;
 byMask=int8(byMask);
 by=double(temp_i);
 clear k datte_1 temp_i temp_1
 
 k=find(bz~=(9999.99));
 temp_1=bz(k);
 datte_1=datte(k);
 temp_i=interp1(datte_1,temp_1,datte,'linear','extrap');
 bzMask=zeros(ndata,1);
 bzMask(k)=1;
 bzMask=int8(bzMask);
 bz=double(temp_i);
 clear k datte_1 temp_i temp_1
 
 k=find(swden~=(999.99));
 temp_1=swden(k);
 datte_1=datte(k);
 temp_i=interp1(datte_1,temp_1,datte,'linear','extrap');
 denMask=zeros(ndata,1);
 denMask(k)=1;
 denMask=int8(denMask);
 swden=double(temp_i);
 clear k datte_1 temp_i temp_1
 
 k=find(swvel~=(99999.9));
 temp_1=swvel(k);
 datte_1=datte(k);
 temp_i=interp1(datte_1,temp_1,datte,'linear','extrap');
 velMask=zeros(ndata,1);
 velMask(k)=1;
 velMask=int8(velMask);
 swvel=double(temp_i);
 clear k datte_1 temp_i temp_1
 

t_delay=20; % 20 min time delay 
mov_avg=25; %25 min trailing average 

 bx=circshift(bx,t_delay); 
 bx=movmean(bx,[mov_avg-1,0]); 
 
 bxMask=circshift(bxMask,t_delay); 

 by=circshift(by,t_delay); 
 by=movmean(by,[mov_avg-1,0]); 

 byMask=circshift(byMask,t_delay); 
 
 bz=circshift(bz,t_delay); 
 bz=movmean(bz,[mov_avg-1,0]); 

 bzMask=circshift(bzMask,t_delay); 

 swvel=circshift(swvel,t_delay); 
 swvel=movmean(swvel,[mov_avg-1,0]); 

 velMask=circshift(velMask,t_delay);  

 swden=circshift(swden,t_delay); 
 swden=movmean(swden,[mov_avg-1,0]); 

 denMask=circshift(denMask,t_delay); 

 
figure
subplot(3,2,1)
plot(datte,bx,'.')
subplot(3,2,2)
plot(datte,by,'.')
subplot(3,2,3)
plot(datte,bz,'.')
subplot(3,2,4)
plot(datte,swvel,'.')
subplot(3,2,5)
plot(datte,swden,'.')

S = ncinfo(file_source);
S.Dimensions.Length=ndata;
S.Dimensions.Unlimited=true;
file_new=file_new_name;
ncwriteschema(file_new,S);
ncwrite(file_new,'bx',bx);
ncwrite(file_new,'by',by);
ncwrite(file_new,'bz',bz);
ncwrite(file_new,'swvel',swvel);
ncwrite(file_new,'swden',swden);
ncwrite(file_new,'date',datte);

ncwrite(file_new,'bxMask',bxMask);
ncwrite(file_new,'byMask',byMask);
ncwrite(file_new,'bzMask',bzMask);
ncwrite(file_new,'velMask',velMask);
ncwrite(file_new,'denMask',denMask);


%check data
swvel_n=ncread(file_new,'swvel');
swden_n=ncread(file_new,'swden');
bx_n=ncread(file_new,'bx');
by_n=ncread(file_new,'by');
bz_n=ncread(file_new,'bz');
datte_n=ncread(file_new,'date');

if sum(bx-bx_n)+sum(by-by_n)+sum(bz-bz_n)+sum(swden-swden_n)+ sum(swvel-swvel_n)==0
    disp('porting success')
else
    disp('FAIL')
    stop
end

