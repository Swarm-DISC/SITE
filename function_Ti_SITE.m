function[Ti_th3]=function_Ti_SITE(Hn,On,N2n,O2n,Hen,tn,dene,te,dvx,dvy,dvz)

narginchk(11,11)
nargoutchk(1,1);

[n1,n2,n3,n4]=size(dene);
Hn=reshape(Hn,1,[]);
On=reshape(On,1,[]);
N2n=reshape(N2n,1,[]);
O2n=reshape(O2n,1,[]);
Hen=reshape(Hen,1,[]);
tn=reshape(tn,1,[]);
dene=reshape(dene,1,[]);
te=reshape(te,1,[]);

dvx=reshape(dvx,1,[]);
dvy=reshape(dvy,1,[]);
dvz=reshape(dvz,1,[]);

iz=length(On);

tiOp=tn; % assume Ti=Tn initially

%---------------ion temperature model-------------------------
amu_m=1.6605*10^-27; % kg
kb   = 8.6*10^-5; % Boltzman constant

rmsd_ti=100; % some random large number
 % Coll. frequencies - Bank and Kockarts 1973, Huba 2000, Schunk 2009

      n=5; % number of neutrals
      An(1)  = 1.; %H
      An(2)  = 16.;%O
      An(3) = 28.; %N2
      An(4) = 32.; %O2
      An(5) = 4.;  %HE
      
      io=2; %( O+ index)
      Ai(io)  = 16.;% O+

      
      % Polarizabilities (Banks and Kockarts, Part A, pg 219)
      alpha0(1)  = 0.667;
      alpha0(2)  = 0.79;
      alpha0(3) = 1.76;
      alpha0(4) = 1.59;
      alpha0(5) = 0.21;
      
      denn(1,:)=Hn;
      denn(2,:)=On;
      denn(3,:)=N2n;
      denn(4,:)=O2n;
      denn(5,:)=Hen;
      
      
Qei=7.7*10^-6/Ai(io)*dene.*(te.^(-1.5)); % electron - ion heating rate

Ti_th3=NaN(1,iz);
      nuin=NaN(n,1);
      Q1=NaN(n,1);
      gamma1=NaN(n,1);  

for i=1:iz
while rmsd_ti>1 % defines # of iterations, untill RMSD <= 1K
ti=tiOp(i);
Tr    = 0.5 * ( ti + tn(i) );

   
      
     %nuin is O+  to neutral coll. freq.
 

         for j=1:n
             
            if j==1 % O+ - H resonant collision
            
            nuin(j)  = 4.63*10^-12 * denn(j,i).* ( tn(i) + ti/16 ).^0.5; % Schunk 2009, pg 107
             
            elseif j==2 % O+ - O resonant collision

%             nuin(j)  = 3.67*10^-11 * denn(j,i).* sqrt(Tr).* ( 1.0 - 0.064 * log10(Tr) ).^2; % Schunk 2009, pg 107
%             nuin(j)  = 4.45*10^-11 * denn(j,i).* sqrt(Tr).* ( 1.04 - .067  * log10(Tr) ).^2; % Bailey and Balan 1996, Huba 2000 (equivalent to Raitt 1975 multiplied by 1.3)
              nuin(j)  = 3.0*10^-11 * denn(j,i).* sqrt(Tr).* ( 1.0 - 0.135 * log10(Tr/1000) ).^2; % Pesnell 1993
%             nuin(j)  = 5.9*10^-11 * denn(j,i).* sqrt(Tr).* ( 1.0 - 0.096 * log10(Tr) ).^2; %Alterenative version of % Pesnell 1993 given in Salah 1993
            
            else % non-resonant interactions
            amu_redmass    = Ai(io) * An(j) / ( Ai(io) + An(j) );
            nuin_tilde = 2.69*10^-9 * denn(j,i) * sqrt( alpha0(j)/amu_redmass ) ;
            nuin(j) = nuin_tilde * An(j) / ( Ai(io) + An(j) ); 
            end
            
          end

       
      
    for j = 1:n
          Q1(j) = 3. * kb * nuin(j) * Ai(io) / ( Ai(io) + An(j) ) ;
    end

    
    Qin=sum(Q1); % ion neutral cooling rate.

% calculate factor for frictional heating term
    for j = 1:n
          gamma1(j) = nuin(j) * amu_m * Ai(io) * An(j) / ( Ai(io) + An(j) );
    end
   gamma_in=sum(gamma1) * 0.624 * 10^19; % this factor comes from SI to eV convertion
    
    
    ti_prev=ti; % save previous step
    ti_new=(Qei(i)*te(i) + Qin*tn(i) + gamma_in*(dvx(i)^2+dvy(i)^2+dvz(i)^2))/ ( Qei(i) + Qin ); % calculate new
    dtemp=ti_new - ti_prev;
    rmsd_ti=rms( dtemp(~isnan(dtemp)) ); 

tiOp(i)=ti_new; % update ion temperature and iterate

end % while ends here
rmsd_ti=100; % reset rmsd
Ti_th3(1,i)=ti_new;
end
    
   
Ti_th3=reshape(Ti_th3,[n1,n2,n3,n4]); % restore array
   
