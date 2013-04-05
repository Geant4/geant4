
#ifdef batch

      PROGRAM main

      PARAMETER (NG=200000,NH=10000)
      COMMON/GCBANK/Q(NG)
      COMMON/PAWC/H(NH)
C     ALLOCATE MEMORY FOR ZEBRA/HBOOK
C
      CALL GZEBRA(NG)
      CALL HLIMIT(-NH)
C
      CALL HPLINT(0)
C
C     INITIALISATION
C
      CALL UGINIT
C
C     PROCESSING
C
      CALL GRUN
C
C     TERMINATION
C
      CALL UGLAST
C
      STOP
      END
      
#else
      
      PROGRAM main
*
*     GEANT main program. To link with the MOTIF user interface
*     the routine GPAWPP(NWGEAN,NWPAW) should be called, whereas
*     the routine GPAW(NWGEAN,NWPAW) gives access to the basic
*     graphics version.
*
      PARAMETER (NWGEAN=3000000, NWPAW=1000000)
      COMMON/GCBANK/GEANT(NWGEAN)
      COMMON/PAWC  /PAW  (NWPAW)
*
*
      CALL GPAW (NWGEAN,NWPAW)
*
      END
*
      SUBROUTINE qnext
      END
*
      SUBROUTINE czopen
      END
*
      SUBROUTINE cztcp
      END
*
      SUBROUTINE czclos
      END
*
      SUBROUTINE czputa
      END
      
#endif           
