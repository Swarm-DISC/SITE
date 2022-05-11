function[ap_date,ap,Ap,f10y_array,f10yb_array]=function_prepare_geomagn_indexes_celest(w_path)

b=load('celest.txt');

ap_date(:,1)=b(:,1); % year
ap_date(:,2)=b(:,2); %month
ap_date(:,3)=b(:,3); %day

ap(:,1)=b(:,15); 
ap(:,2)=b(:,16); 
ap(:,3)=b(:,17); 
ap(:,4)=b(:,18); 
ap(:,5)=b(:,19);
ap(:,6)=b(:,20); 
ap(:,7)=b(:,21); 
ap(:,8)=b(:,22);
Ap=b(:,23);

%-----read F10.7-----------------------


f10y_array(:,1)=b(:,1); % year
f10y_array(:,2)=b(:,2); %month
f10y_array(:,3)=b(:,3); %day
f10y_array(:,4)=b(:,31);%F10.7-obs
f10y_array(:,5)=b(:,32); % Centered 81-day arithmetic average of F10.7 (observed).

f10yb_array(:,1)=b(:,1); % year
f10yb_array(:,2)=b(:,2); %month
f10yb_array(:,3)=b(:,3); %day
f10yb_array(:,4)=b(:,27);%F10.7 -adj
f10yb_array(:,5)=b(:,29);%Centered 81-day arithmetic average of F10.7 (adjusted).

cd(w_path)
