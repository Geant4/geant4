*    
*     ********************************************************************
*     * DISCLAIMER                                                       *
*     *                                                                  *
*     * The following disclaimer summarizes all the specific disclaimers *
*     * of contributors to this software. The specific disclaimers,which *
*     * govern, are listed with their locations in:                      *
*     *   http://cern.ch/geant4/license                                  *
*     *                                                                  *
*     * Neither the authors of this software system, nor their employing *
*     * institutes,nor the agencies providing financial support for this *
*     * work  make  any representation or  warranty, express or implied, *
*     * regarding  this  software system or assume any liability for its *
*     * use.                                                             *
*     *                                                                  *
*     * This  code  implementation is the  intellectual property  of the *
*     * GEANT4 collaboration.                                            *
*     * By copying,  distributing  or modifying the Program (or any work *
*     * based  on  the Program)  you indicate  your  acceptance of  this *
*     * statement, and all its terms.                                    *
*     ********************************************************************
*    
#ifdef batch  

      PROGRAM main
*
*
      PARAMETER (NGBANK=100000, NHBOOK=10000)
      COMMON/GCBANK/Q(NGBANK)
      COMMON/PAWC  /H(NHBOOK)
*
      CALL GZEBRA( NGBANK)
      CALL HLIMIT(-NHBOOK)
*
* *** initialize HIGZ
      CALL HPLINT(0)
*
* *** GEANT initialisation
      CALL UGINIT
*
* *** End of RUN
      CALL UGLAST
*
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
