*                                                                               
*=== evevap ===========================================================*        
*                                                                               
      SUBROUTINE EVEVAP ( WEE )                                                 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      RETURN                                                                    
*=== End of subroutine Evevap =========================================*        
      END                                                                       
*                                                                      *        
*=== eberttp ===========================================================*        
*                                                                      *        
C     SUBROUTINE BERTTP                                                         
C     RETURN                                                                    
C     END                                                                       
*                                                                      *        
*=== sincini ===========================================================*        
*                                                                      *        
C     SUBROUTINE INCINI                                                         
C     RETURN                                                                    
*=== End of subroutine incini =========================================*        
C     END                                                                       
*                                                                               

C     eva2he ===========================================================*        
*                                                                      *        
C     SUBROUTINE EVA2HE(MO,IREJ)
C     RETURN
C     END

C     energy===========================================================
C      DOUBLE PRECISION FUNCTION ENERGY(AI,AIZ)
C      ENERGY=0.D0
C      RETURN
C      END
*=== raco =============================================================*
*
      SUBROUTINE RACO (WX,WY,WZ)
 
C     INCLUDE '(DBLPRC)'
*$ CREATE DBLPRC.ADD
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( KALGNM = 2 )
      PARAMETER ( ANGLGB = 5.0D-16 )
      PARAMETER ( ANGLSQ = 2.5D-31 )
      PARAMETER ( AXCSSV = 0.2D+16 )
      PARAMETER ( ANDRFL = 1.0D-38 )
      PARAMETER ( AVRFLW = 1.0D+38 )
      PARAMETER ( AINFNT = 1.0D+30 )
      PARAMETER ( AZRZRZ = 1.0D-30 )
      PARAMETER ( EINFNT = +69.07755278982137 D+00 )
      PARAMETER ( EZRZRZ = -69.07755278982137 D+00 )
      PARAMETER ( ONEMNS = 0.999999999999999  D+00 )
      PARAMETER ( ONEPLS = 1.000000000000001  D+00 )
      PARAMETER ( CSNNRM = 2.0D-15 )
      PARAMETER ( DMXTRN = 1.0D+08 )
      PARAMETER ( ZERZER = 0.D+00 )
      PARAMETER ( ONEONE = 1.D+00 )
      PARAMETER ( TWOTWO = 2.D+00 )
      PARAMETER ( THRTHR = 3.D+00 )
      PARAMETER ( FOUFOU = 4.D+00 )
      PARAMETER ( FIVFIV = 5.D+00 )
      PARAMETER ( SIXSIX = 6.D+00 )
      PARAMETER ( SEVSEV = 7.D+00 )
      PARAMETER ( EIGEIG = 8.D+00 )
      PARAMETER ( ANINEN = 9.D+00 )
      PARAMETER ( TENTEN = 10.D+00 )
      PARAMETER ( HLFHLF = 0.5D+00 )
      PARAMETER ( ONETHI = ONEONE / THRTHR )
      PARAMETER ( TWOTHI = TWOTWO / THRTHR )
      PARAMETER ( ONEFOU = ONEONE / FOUFOU )
      PARAMETER ( THRTWO = THRTHR / TWOTWO )
      PARAMETER ( PIPIPI = 3.141592653589793238462643383279D+00 )
      PARAMETER ( TWOPIP = 6.283185307179586476925286766559D+00 )
      PARAMETER ( PIP5O2 = 7.853981633974483096156608458199D+00 )
      PARAMETER ( PIPISQ = 9.869604401089358618834490999876D+00 )
      PARAMETER ( PIHALF = 1.570796326794896619231321691640D+00 )
      PARAMETER ( ERFA00 = 0.886226925452758013649083741671D+00 )
      PARAMETER ( ENEPER = 2.718281828459045235360287471353D+00 )
      PARAMETER ( SQRENT = 1.648721270700128146848650787814D+00 )
      PARAMETER ( SQRSIX = 2.449489742783178098197284074706D+00 )
      PARAMETER ( SQRSEV = 2.645751311064590590501615753639D+00 )
      PARAMETER ( SQRT12 = 3.464101615137754587054892683012D+00 )
      PARAMETER ( CLIGHT = 2.99792458         D+10 )
      PARAMETER ( AVOGAD = 6.0221367          D+23 )
      PARAMETER ( BOLTZM = 1.380658           D-23 )
      PARAMETER ( AMELGR = 9.1093897          D-28 )
      PARAMETER ( PLCKBR = 1.05457266         D-27 )
      PARAMETER ( ELCCGS = 4.8032068          D-10 )
      PARAMETER ( ELCMKS = 1.60217733         D-19 )
      PARAMETER ( AMUGRM = 1.6605402          D-24 )
      PARAMETER ( AMMUMU = 0.113428913        D+00 )
      PARAMETER ( AMPRMU = 1.007276470        D+00 )
      PARAMETER ( AMNEMU = 1.008664904        D+00 )
      PARAMETER ( ALPFSC = 7.2973530791728595 D-03 )
      PARAMETER ( FSCTO2 = 5.3251361962113614 D-05 )
      PARAMETER ( FSCTO3 = 3.8859399018437826 D-07 )
      PARAMETER ( FSCTO4 = 2.8357075508200407 D-09 )
      PARAMETER ( PLABRC = 0.197327053        D+00 )
      PARAMETER ( AMELCT = 0.51099906         D-03 )
      PARAMETER ( AMUGEV = 0.93149432         D+00 )
      PARAMETER ( AMMUON = 0.105658389        D+00 )
      PARAMETER ( AMPRTN = 0.93827231         D+00 )
      PARAMETER ( AMNTRN = 0.93956563         D+00 )
      PARAMETER ( AMDEUT = 1.87561339         D+00 )
      PARAMETER ( COUGFM = ELCCGS * ELCCGS / ELCMKS * 1.D-07 * 1.D+13
     &                   * 1.D-09 )
      PARAMETER ( RCLSEL = 2.8179409183694872 D-13 )
      PARAMETER ( BLTZMN = 8.617385           D-14 )
      PARAMETER ( GEVMEV = 1.0                D+03 )
      PARAMETER ( EMVGEV = 1.0                D-03 )
      PARAMETER ( ALGVMV = 6.90775527898214   D+00 )
      PARAMETER ( RADDEG = 180.D+00 / PIPIPI )
      PARAMETER ( DEGRAD = PIPIPI / 180.D+00 )
      LOGICAL LGBIAS, LGBANA
      COMMON / GLOBAL / LGBIAS, LGBANA
C     INCLUDE '(DIMPAR)'
*$ CREATE DIMPAR.ADD
      PARAMETER ( MXXRGN = 5000 )
      PARAMETER ( MXXMDF = 56   )
      PARAMETER ( MXXMDE = 50   )
      PARAMETER ( MFSTCK = 1000 )
      PARAMETER ( MESTCK = 100  )
      PARAMETER ( NALLWP = 39   )
      PARAMETER ( MPDPDX = 8    )
      PARAMETER ( ICOMAX = 180  )
      PARAMETER ( NSTBIS = 304  )
      PARAMETER ( IDMAXP = 210  )
      PARAMETER ( IDMXDC = 620  )
      PARAMETER ( MKBMX1 = 1    )
      PARAMETER ( MKBMX2 = 1    )
C     INCLUDE '(IOUNIT)'
*$ CREATE IOUNIT.ADD
      PARAMETER ( LUNIN  = 5  )
      PARAMETER ( LUNOUT = 6  )
      PARAMETER ( LUNERR = 15 )
      PARAMETER ( LUNBER = 14 )
      PARAMETER ( LUNECH = 8  )
      PARAMETER ( LUNFLU = 13 )
      PARAMETER ( LUNGEO = 16 )
      PARAMETER ( LUNPGS = 12 )
      PARAMETER ( LUNRAN = 2  )
      PARAMETER ( LUNXSC = 9  )
      PARAMETER ( LUNDET = 17 )
      PARAMETER ( LUNRAY = 10 )
      PARAMETER ( LUNRDB = 1  )
C********************************************************************
C     VERSION JUNE 81 BY             PERTTI AARNIO
C     LAST CHANGE 22. JUNE 81 BY     PERTTI AARNIO
C                                    HELSINKI UNIVERSITY OF
C                                    TECHNOLOGY, FINLAND
C
C
C     SUBROUTINE OF FLUKA TO GIVE THE DIRECTION COSINES OF RANDOM
C     UNIFORM (ISOTROPIC) DIRECTION IN THREE DIMENSIONAL SPACE.
C********************************************************************
C
  10  CONTINUE
         X=TWOTWO*RNDM(X)-ONEONE
         Y=RNDM(Y)
         X2=X*X
         Y2=Y*Y
      IF ( X2+Y2 .GT. ONEONE ) GO TO 10
      CFE=(X2-Y2)/(X2+Y2)
      SFE=TWOTWO*X*Y/(X2+Y2)
* z = 1/2 [ 1 + cos (theta) ]
      Z =RNDM(Z)
* 1/2 sin (theta)
      WZ=SQRT(Z*(ONEONE-Z))
      WX=TWOTWO*WZ*CFE
      WY=TWOTWO*WZ*SFE
      WZ=TWOTWO*Z-ONEONE
      RETURN
      END
 
C      SUBROUTINE HISTOG(I)
C      RETURN
C      END
       SUBROUTINE STALIN
       RETURN
       END
       SUBROUTINE FRBKIN(L,LP)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       LOGICAL L,LP
       RETURN
       END
