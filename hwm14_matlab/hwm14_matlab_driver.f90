!
!  Matlab driver for HWM93
!  L.L. Jan/2020
!
!*******************************************************************

      program matlab_hwm14

      implicit none
      INTEGER      :: IYD, ip, ixx
      REAL(4)      :: SEC,ALT,GLAT,GLON,STL,F107A,F107,AP(2),zAP,zap3hr
      REAL(4)      :: W(2)
      REAL		   :: SW(25)
      CHARACTER*10 :: IxxVAL

      open(9005,file='input_hwm14.dat')
      open(9002,file='output_hwm14.dat')

      CALL GETARG(1, IxxVAL)
      READ (IxxVAL,'(I10)') ixx

      do 9001 ip=1,ixx
      read(9005,*) iyd,sec,stl,glat,glon,alt,zAP,zap3hr,f107,f107a


      ! iyd  - YEAR AND DAY AS YYDDD



      ap(1) = zAP ! not used in HWM14
      ap(2) = zap3hr !CURRENT 3HR ap INDEX

    
      call hwm14(iyd,sec,alt,glat,glon,stl,f107a,f107,ap,w)
 
      write(9002,101) w(1), w(2) 

 9001 continue
  101 FORMAT (2E13.5)
!STOP
!END

 
      end program matlab_hwm14