function[ne_adj,err_ne,teh_adj,err_teh,tel_adj,err_tel]=function_adjust_swarm_ne_te(ne,teh,tel,satellite)

% Coefficient update: Jan-2021
% correct Ne and Te based on ISR measurements

narginchk(4,4)
nargoutchk(6,6);

f = ( ne./(1.24*10^4) ).^0.5; % convert to frequency

if strcmp(satellite, 'A')==1 
    
a_f=1.1069; da_f=0.0060;

a_teh=1.2844; b_teh=-1083; c_teh=0;
da_teh=0.0464; db_teh=127; dc_teh=0;

a_tel=1.0; b_tel=-723; c_tel=0.0;
da_tel=0.0; db_tel=35; dc_tel=0.0;

end

if strcmp(satellite, 'B')==1 
    
a_f=1.0952; da_f=0.0071;
    
a_teh=1.1626; b_teh=-827; c_teh=0;
da_teh=0.0622; db_teh=176; dc_teh=0;

a_tel=1.0; b_tel=-698; c_tel=0.0;
da_tel=0.0; db_tel=41; dc_tel=0.0;

end

if strcmp(satellite, 'C')==1 
    
a_f=1.1087; da_f=0.0060;
    
a_teh=1.2153; b_teh=-916; c_teh=0;
da_teh=0.0410; db_teh=111; dc_teh=0;

a_tel=1.0; b_tel=-682; c_tel=0.0;
da_tel=0.0; db_tel=32; dc_tel=0.0;

end


    f1=(a_f-da_f) * f;
    f2=(a_f+da_f) * f;
    ne1=1.24*10^4 * f1.^2; % convert back to density
    ne2=1.24*10^4 * f2.^2;
    ne_adj=(ne1+ne2)/2;
    err_ne=abs(ne2-ne1)/2;

ah1=a_teh-da_teh;
ah2=a_teh+da_teh;

bh1=b_teh-db_teh;
bh2=b_teh+db_teh;

ch1=c_teh-dc_teh;
ch2=c_teh+dc_teh;

a_hgn=combvec([ah1,ah2],[bh1,bh2],[ch1,ch2]);

al1=a_tel-da_tel;
al2=a_tel+da_tel;

bl1=b_tel-db_tel;
bl2=b_tel+db_tel;

cl1=c_tel-dc_tel;
cl2=c_tel+dc_tel;

a_lgn=combvec([al1,al2],[bl1,bl2],[cl1,cl2]);

tehx_adj=NaN(length(teh),length(a_hgn));
telx_adj=NaN(length(teh),length(a_hgn));

for i=1:length(a_hgn)
    tehx_adj(:,i) = a_hgn(1,i) * teh + a_hgn(2,i) + a_hgn(3,i) * ne;
    telx_adj(:,i) = a_lgn(1,i) * tel + a_lgn(2,i) + a_lgn(3,i) * ne;
end

teh_adj=a_teh * teh + b_teh + c_teh * ne;
err_teh=( max(tehx_adj,[],2) - min(tehx_adj,[],2) )/2;

tel_adj=a_tel * tel + b_tel + c_tel * ne;
err_tel=( max(telx_adj,[],2) - min(telx_adj,[],2) )/2;

err_tel=err_tel';
err_teh=err_teh';



    
    
    
    
    
