;pro E_weimer05

Forward_Function EpotVal

CD, 'M:\SITE_software\weimer05'
RESTORE, 'W05plot1.sav'
file = 'imf.txt'
openr,lu1,file,/get_lun
readf,lu1,ndata
data = FLTARR(12, ndata)
Enorth = FLTARR(1,ndata)
Eeast = FLTARR(1,ndata)
fi = FLTARR(1,ndata)
readf,lu1,data
mlat=data(0,*)
mlt=data(1,*)
alt=data(2,*)
by=data(3,*)
bz=data(4,*)
swvel=data(5,*)
swden=data(6,*)
yr=data(7,*)
mo=data(8,*)
dd=data(9,*)
hr=data(10,*)
tlt=data(11,*)
;AL=data(12,*)


openw, lu2, 'Efield_out_with_pot.txt',/get_lun

for i=0,ndata-1,1 do begin

if (mlat(i) le 0) then begin
  tlt(i)=-tlt(i)
  by(i)=-by(i)
endif
;SetModel, by(i), bz(i), tlt(i), swvel(i), swden(i),[AL(i)], /YZ
SetModel, by(i), bz(i), tlt(i), swvel(i), swden(i),[], /YZ
fi(i)=EpotVal (abs(mlat(i)), mlt(i))
if (mlat(i) gt 0) then efield, abs(mlat(i)), mlt(i), Enorth_temp, Eeast_temp, ALT=alt(i)
if (mlat(i) le 0) then efield, abs(mlat(i)), mlt(i), Enorth_temp, Eeast_temp, ALT=alt(i),/south
Enorth(i)= Enorth_temp 
Eeast(i)=Eeast_temp 
printf, lu2, Enorth(i), Eeast(i), fi(i), FORMAT='(F11.5,1X,F11.5,1X,F11.5)'
endfor


close, lu1
close, lu2

end