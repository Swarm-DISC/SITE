function[app,ap_avg,AP_daily,f107_obs,f107a_obs,f107_adj,f107a_adj]=function_geomagn_indexes(y,mo,d,hr,mn,sec,ap_hr_before,w_path,...
    ap_date,ap,Ap,f10y_array,f10yb_array)

apy=ap_date(:,1);
apm=ap_date(:,2);
apd=ap_date(:,3);

f10y=f10y_array(:,1);
f10m=f10y_array(:,2);
f10d=f10y_array(:,3);
f10=f10y_array(:,4);
f10a=f10y_array(:,5);

f10yb=f10yb_array(:,1);
f10mb=f10yb_array(:,2);
f10db=f10yb_array(:,3);
f10b =f10yb_array(:,4);
f10ba=f10yb_array(:,5);

if nargout~=7
    stop
end


if ap_hr_before==24 || ap_hr_before==12
    
% apdoy=datenum(apy,apm,apd);


nn=length(y);

f107_obs=NaN(nn,1);
f107a_obs=NaN(nn,1);
f107_adj=NaN(nn,1);
f107a_adj=NaN(nn,1);
app=NaN(nn,1);
ap_avg=NaN(nn,1);
AP_daily=NaN(nn,1);

for k=1:length(y)

% Geophys index part    
    k1=find(y(k)==f10y & mo(k)==f10m & d(k)==f10d);
    f107_obs(k)=f10(k1); % daily observed F10.7 flux
    f107a_obs(k)=f10a(k1); % 81 day centered avg
%     f107a_obs(k)=nanmean(f10(k1-40:k1+40));
       
    kb1=find(y(k)==f10yb & mo(k)==f10mb & d(k)==f10db);
    f107_adj(k)=f10b(kb1); % daily adjusted F10.7 flux
    f107a_adj(k)=f10ba(kb1); % 81 day centered avg
%     f107a_adj(k)=nanmean(f10b(kb1-40:kb1+40));
    
    
    
    k2=find(y(k)==apy & mo(k)==apm & d(k)==apd);
    AP_daily(k)=Ap(k2);
    UT=hr(k)+mn(k)/60+sec(k)/3600;
    UT=mod(UT,24); % to exclude UT=24 h
    apti=fix(fix(UT)/3)+1; %index that determines the part of the day (1-8) for ap-index purposes
    in_ap=apti; % hourly ap index (column) for a given measurement
    if ap_hr_before==12
    sn_ap=sign(in_ap-(mod((in_ap+4:in_ap+7)-1,8)+1)); %- 12h before
    end
    if ap_hr_before==24
    sn_ap=sign(in_ap-(mod((in_ap:in_ap+7)-1,8)+1)); %- 24h before
    end
    sn_ap(sn_ap==0)=-1;
    sn_ap(sn_ap==1)=0; %if same day (or 1) then 0.
    if ap_hr_before==12
    aph_ind=mod((in_ap+4:in_ap+7)-1,8)+1; %12h before
    ap_avg=mean([ap(k2+sn_ap(1),aph_ind(1)),ap(k2+sn_ap(2),aph_ind(2)),ap(k2+sn_ap(3),aph_ind(3)),...
            ap(k2+sn_ap(4),aph_ind(4))]); % 12h average ap
    end
    if ap_hr_before==24
    aph_ind=mod((in_ap:in_ap+7)-1,8)+1; %24h before
    ap_avg(k)=mean([ap(k2+sn_ap(1),aph_ind(1)),ap(k2+sn_ap(2),aph_ind(2)),ap(k2+sn_ap(3),aph_ind(3)),...
            ap(k2+sn_ap(4),aph_ind(4)),ap(k2+sn_ap(5),aph_ind(5)),ap(k2+sn_ap(6),aph_ind(6)),...
            ap(k2+sn_ap(7),aph_ind(7)),ap(k2+sn_ap(8),aph_ind(8))]); % 24h average ap
    end       
   app(k)=ap(k2,in_ap); %3-hour ap  
% end Geophys index part
end

else
    disp('ap time average error: should be 12 or 24')
    stop
end

cd(w_path)
