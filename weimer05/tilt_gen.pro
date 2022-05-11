;pro tilt_gen

Forward_Function Tilt_Angle

CD, 'M:\SITE_software\weimer05'
RESTORE, 'W05plot1.sav' ; may need to manually execute this line
file = 'tilt_2013-2023.txt'
openw, lu2, file,/get_lun
for yr=2013,2023 do begin
  for mo=1,12 do begin
    for dd=1,31 do begin
      for hr=0,23 do begin
      EARTHTRANS, yr,mo,dd,hr
      tlt=Tilt_Angle()
      printf, lu2, yr, mo, dd, hr, tlt, FORMAT='(F9.0,1X,F8.0,1X,F8.0,1X,F8.0,1X,F8.4)' ;FORMAT='(5(F8.4,2X))'
      endfor
    endfor
  endfor
endfor

close, lu2

end
