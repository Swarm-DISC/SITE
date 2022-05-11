C      Generate exexutable

      Integer zIDAY, ip, ixx
      REAL zUT,zALT,zXLAT,zXLONG,zF107A,zF107,zAP,zXLST


      DIMENSION D(9,1),T(2,1),SW(25),APH(1)
         
      DIMENSION IDAY(1),UT(1),ALT(1),XLAT(1),XLONG(1),XLST(1),
     & F107A(1),F107(1),AP(1)    
      COMMON/GTS3C/DL(16)
      COMMON/DATIME/ISDATE(3),ISTIME(2),NAME(2)

      CHARACTER*10 IxxVAL

      open(9005,file='input_msis.dat')
      open(9002,file='output_msis.dat')

      !read (*,*) ixx
      
      CALL GETARG(1, IxxVAL)
      READ (IxxVAL,'(I10)') ixx

      do 9001 ip=1,ixx
      read(9005,*) zIDAY, zUT, zXLAT, zXLONG,
     $ zALT, zAP, zF107, zF107A
	  
      IDAY=zIDAY
      UT=zUT
      ALT=zALT
      XLAT=zXLAT
      XLONG=zXLONG
      F107A=zF107A
      F107=zF107
      AP=zAP


      XLST=UT/3600+XLONG/15

      DATA SW/25*1./
      CALL TSELEC(SW)

      DO I=1,1
         CALL GTD7(IDAY(I),UT(I),ALT(I),XLAT(I),XLONG(I),XLST(I),
     &             F107A(I),F107(I),AP(I),48,D(1,I),T(1,I))


	  write(9002,101) D(1,I),D(2,I),D(3,I),D(4,I),D(5,I),
     $    D(6,I),D(7,I),D(8,I),D(9,I),T(1,I),T(2,I)
     
C        D(1) - HE NUMBER DENSITY(CM-3)
C        D(2) - O NUMBER DENSITY(CM-3)
C        D(3) - N2 NUMBER DENSITY(CM-3)
C        D(4) - O2 NUMBER DENSITY(CM-3)
C        D(5) - AR NUMBER DENSITY(CM-3)                       
C        D(6) - TOTAL MASS DENSITY(GM/CM3)
C        D(7) - H NUMBER DENSITY(CM-3)
C        D(8) - N NUMBER DENSITY(CM-3)
C        D(9) - Anomalous oxygen NUMBER DENSITY(CM-3)
C        T(1) - EXOSPHERIC TEMPERATURE
C        T(2) - TEMPERATURE AT ALT
      ENDDO

 9001 continue
      
  101 FORMAT(9E13.5,2F11.2)
      

  200 FORMAT(I10,8F10.3)

 39100 FORMAT(F6.1,F7.1,F7.3,F9.2,F7.2,F8.2,F9.2,1P5E12.3,
     $1P5E12.3,1P5E12.3,1P5E12.3,1P5E12.3)
  

      STOP
      END
