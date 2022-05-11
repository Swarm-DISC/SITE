function[H,O,N2,O2,HE,Tn,rho,N]=function_msis00_matlab(doy, sec_of_day, glat, glon, alt, Ap_daily, f107_obs, f107a_obs, w_path)
narginchk(9,9)
nargoutchk(8,8);
k=length(doy);
%---run msis-------
formatSpec = 'NRLMSIS00.exe%10d';
cd([w_path,'\nrlmsise00_matlab'])
delete 'input_msis.dat' 'output_msis.dat'

        fid = fopen('input_msis.dat','wt');
       
        for j=1:k           
        fprintf(fid,'%4.0f%6.0f%9.2f%9.2f%9.2f%9.2f%9.2f%9.2f\n',...
            doy(j), sec_of_day(j), glat(j), glon(j), alt(j), Ap_daily(j), f107_obs(j), f107a_obs(j) );
        end
        
        fclose (fid);

      clear inputdata   
      inputdata = sprintf(formatSpec,k);
      [~,~] = system(inputdata);
      
             
        clear a
        a=load('output_msis.dat');
        if k ~=length(a(:,1))
            disp('Wrong I/O')
            stop
        end
        
        HE  =a(:,1);
        H   =a(:,7);
        O   =a(:,2);
        N2  =a(:,3);
        O2  =a(:,4);
        N   =a(:,8);
        Tn  =a(:,11);
        rho =a(:,6);

cd(w_path)