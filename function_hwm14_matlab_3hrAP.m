function[u_merid,u_zon]=function_hwm14_matlab_3hrAP(yy, doy, sec_of_day, glat, glon, alt, Ap_daily, ap_3hr, f107_obs, f107a_obs, w_path)
% (m/sec + Northward)
% (m/sec + Eastward)
narginchk(11,11)
nargoutchk(2,2);
k=length(doy);

yyddd=yy*1000+doy;
yyddd=yyddd-fix(yyddd/10^5)*10^5;

LT=mod(sec_of_day/3600+glon/15,24);

%---run msis-------
formatSpec = 'hwm14_3hrAP.exe%10d';
cd([w_path,'\hwm14_matlab'])
delete 'input_hwm14.dat' 'output_hwm14.dat'

        fid = fopen('input_hwm14.dat','wt');
       
        for j=1:k           
        fprintf(fid,'%s%6.0f%9.2f%9.2f%9.2f%9.2f%9.2f%9.2f%9.2f%9.2f\n',...
        num2str(yyddd(j),'%05.f'), sec_of_day(j), LT(j), glat(j), glon(j), alt(j),...
        Ap_daily(j), ap_3hr(j), f107_obs(j), f107a_obs(j) );
        end
        
        fclose (fid);
        
      clear inputdata   
      inputdata = sprintf(formatSpec,k);
      [~,~] = system(inputdata);
      
        clear a
        a=load('output_hwm14.dat');
        if k ~=length(a(:,1))
            disp('Wrong I/O')
            stop
        end
        
        u_merid=a(:,1);
        u_zon=a(:,2);

cd(w_path)