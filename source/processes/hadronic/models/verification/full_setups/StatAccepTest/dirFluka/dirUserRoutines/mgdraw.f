C----------------------------------------------------------------------
C
C As far as the analysis of calorimeter observables is concerned, 
C the only two assumptions made on this code are the following:
C
C   1) NO magnetic field;
C
C   2) active layers have EVEN region numbers
C      (i.e., the first region of the calorimeter, with region
C       number = 3, is the passive layer).
C 
C Look at  ***LOOKHERE***  :
C  -  to decide if you want the full counting of all particles, 
C     which is quite CPU time consuming. By default you get the
C     counting of all particles, except for electrons, positrons,
C     and gammas, with negligible CPU time overhead.
C  -  to decide whether you want Birks quenching switched on.
C     (The Birks coefficients should be specified in the card
C      USERDUMP, with "UDQUENCH" as SDUM .)
C
C----------------------------------------------------------------------

*$ CREATE MGDRAW.FOR
*COPY MGDRAW
*                                                                      *
*=== mgdraw ===========================================================*
*                                                                      *
      SUBROUTINE MGDRAW ( ICODE, MREG )

      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'
*
*----------------------------------------------------------------------*
*                                                                      *
*     Copyright (C) 1990-2005      by        Alfredo Ferrari           *
*     All Rights Reserved.                                             *
*                                                                      *
*                                                                      *
*     MaGnetic field trajectory DRAWing: actually this entry manages   *
*                                        all trajectory dumping for    *
*                                        drawing                       *
*                                                                      *
*     Created by                     Alfredo Ferrari                   *
*                                    INFN - Milan                      *
*     last change  09-jul-05  by     Alfredo Ferrari                   *
*                                    INFN - Milan                      *
*                                                                      *
*----------------------------------------------------------------------*
*
      INCLUDE '(CASLIM)'
      INCLUDE '(COMPUT)'
      INCLUDE '(SOURCM)'
      INCLUDE '(FHEAVY)'
      INCLUDE '(FLKSTK)'
      INCLUDE '(GENSTK)'
      INCLUDE '(MGDDCM)'
      INCLUDE '(PAPROP)'
      INCLUDE '(SUMCOU)'
      INCLUDE '(TRACKR)'

      INCLUDE '(FLKMAT)'
      INCLUDE '(RESNUC)'
      INCLUDE '(NUCDAT)'

      INCLUDE '(QUEMGD)'
*
      INCLUDE 'myinclude.inc'
*
      DIMENSION DTQUEN( MXTRCK, MAXQMG )
*
      CHARACTER*20 FILNAM
      LOGICAL LFCOPE
      SAVE LFCOPE
      DATA LFCOPE / .FALSE. /

      INTEGER iReadoutLayer, iBinRadius
      REAL radius, currentRadius
      LOGICAL found
      REAL Ekin

      REAL stepEdep

      LOGICAL isFirstCall
      SAVE isFirstCall
      DATA isFirstCall/.TRUE./

C  ***LOOKHERE*** : decide whether to switch on/off Birks quenching. 
      LOGICAL isBirksOn
      PARAMETER ( isBirksOn = .FALSE. )

C  ***LOOKHERE*** : if you set true the variable  isFullParticleCountOn
C                   then the counting of the particles includes all
C                   the particles, but the CPU time overhead can be
C                   more than a factor of 10.
C                   If you leave instead that variable to its default
C                   value (false), then you get the counting of all
C                   particles excluding only electrons, positrons, and
C                   gammas, with negligible CPU time overhead.
      LOGICAL isFullParticleCountOn
      PARAMETER ( isFullParticleCountOn = .FALSE. )  

*
*----------------------------------------------------------------------*
*                                                                      *
*     Icode = 1: call from Kaskad                                      *
*     Icode = 2: call from Emfsco                                      *
*     Icode = 3: call from Kasneu                                      *
*     Icode = 4: call from Kashea                                      *
*     Icode = 5: call from Kasoph                                      *
*                                                                      *
*----------------------------------------------------------------------*
*            
      IF ( isFirstCall .AND.  LTRACK.EQ.1 ) THEN
         isFirstCall = .FALSE.
         print *, '----------------------------------'
         IF ( isBirksOn ) THEN
            print *, '*** Birks ON ***'             
         ENDIF
         print *, 'Beam Particle = ', JTRACK, '   ',
     &            PRNAME(JTRACK)
         print *, 'Beam Particle Kinetic Energy =', 
     &            ETRACK - AM(JTRACK), ' GeV'
         print *, '----------------------------------'
      ENDIF

C---- If the continuous energy deposition happens inside the 
C     calorimeter then its value is added to the variable
C     edepAllCalo, which keeps the sum of all energy depositions
C     in the current event.
C     For edepAllCalo, we do not apply the Birks qunching even 
C     if it is active and the energy deposition happens in an 
C     active layer.
C     Notice that the active layers correspond to regions with
C     even number, excluding region number 2 (experimental hall).

      IF ( ( Ntrack.EQ.1 ).AND.( Mtrack.EQ.1 ) ) THEN
         IF ( MREG.GT.2 ) THEN

CCC            IF ( LTRACK.eq.1 ) then   ! Only for primary
CCC               print *, ' MGDRAW : mreg=', MREG, ' (', 
CCC     &              XTRACK(1), ',', YTRACK(1), ',', ZTRACK(1), ')'
CCC            ENDIF

            stepEdep = DTRACK(1)

            edepAllCalo = edepAllCalo + stepEdep
            IF ( JTRACK.EQ.3 .OR. JTRACK.EQ.4 ) THEN
               edepAllCalo_electron = edepAllCalo_electron + stepEdep
            ELSE IF ( JTRACK.EQ.10 .OR. JTRACK.EQ.11 ) THEN
               edepAllCalo_muon = edepAllCalo_muon + stepEdep
            ELSE IF ( JTRACK.EQ.15 .OR. JTRACK.EQ.16 ) THEN
               edepAllCalo_kaon = edepAllCalo_kaon + stepEdep
            ELSE IF ( JTRACK.EQ.1 .OR. JTRACK.EQ.2 ) THEN
               edepAllCalo_proton = edepAllCalo_proton + stepEdep
            ELSE IF ( JTRACK.EQ.13 .OR. JTRACK.EQ.14 ) THEN
               edepAllCalo_pion = edepAllCalo_pion + stepEdep
            ELSE IF ( JTRACK.EQ.8 .OR. JTRACK.EQ.9 .OR. ! neutrons/antineutrons
     &                JTRACK.LE.-2 .OR. ! ions
     &                JTRACK.EQ.20 .OR. JTRACK.EQ.21 .OR.   ! Sigma-/+
     &                JTRACK.EQ.31 .OR. JTRACK.EQ.33 .OR.   ! AntiSigma-/+
     &                JTRACK.EQ.36 .OR. JTRACK.EQ.37 .OR.   ! Xsi-, AntiXsi-
     &                JTRACK.EQ.38 .OR. JTRACK.EQ.39 ) THEN ! Omega-, AntiOmega-
               edepAllCalo_ion = edepAllCalo_ion + stepEdep
            ENDIF

           IF ( MOD(MREG,2).EQ.0 ) THEN  

              IF ( isBirksOn ) THEN
                 RULLL = ZERZER
                 CALL QUENMG( ICODE, MREG, RULLL, DTQUEN ) 
                 stepEdep = DTQUEN(1,1)
              ELSE
                 stepEdep = DTRACK(1)
              ENDIF   

CCC               print *, ' continuous: stepEdep=', stepEdep

               edepVisCalo = edepVisCalo + stepEdep
               IF ( JTRACK.EQ.3 .OR. JTRACK.EQ.4 ) THEN
                  edepVisCalo_electron = edepVisCalo_electron + 
     &                 stepEdep
               ELSE IF ( JTRACK.EQ.10 .OR. JTRACK.EQ.11 ) THEN
                  edepVisCalo_muon = edepVisCalo_muon + 
     &                 stepEdep
               ELSE IF ( JTRACK.EQ.15 .OR. JTRACK.EQ.16 ) THEN
                  edepVisCalo_kaon = edepVisCalo_kaon + 
     &                 stepEdep
               ELSE IF ( JTRACK.EQ.1 .OR. JTRACK.EQ.2 ) THEN
                  edepVisCalo_proton = edepVisCalo_proton + 
     &                 stepEdep
               ELSE IF ( JTRACK.EQ.13 .OR. JTRACK.EQ.14 ) THEN
                  edepVisCalo_pion = edepVisCalo_pion + 
     &                 stepEdep
               ELSE IF ( JTRACK.EQ.8 .OR. JTRACK.EQ.9 .OR. ! neutrons/antineutrons
     &                   JTRACK.LE.-2 .OR. ! ions
     &                   JTRACK.EQ.20 .OR. JTRACK.EQ.21 .OR.   ! Sigma-/+
     &                   JTRACK.EQ.31 .OR. JTRACK.EQ.33 .OR.   ! AntiSigma-/+
     &                   JTRACK.EQ.36 .OR. JTRACK.EQ.37 .OR.   ! Xsi-, AntiXsi-
     &                   JTRACK.EQ.38 .OR. JTRACK.EQ.39 ) THEN ! Omega-, AntiOmega-
                  edepVisCalo_ion = edepVisCalo_ion + 
     &                 stepEdep
               ENDIF

               iReadoutLayer = 1 + 
     &              (MREG - 4) / ( 2 * numberOfLayersPerReadoutLayer )
CCC               print *, ' MGDRAW : MREG=', MREG, 
CCC     &              ' iReadouLayer=', iReadoutLayer

               longitudinalProfile( iReadoutLayer ) =
     &              longitudinalProfile( iReadoutLayer ) + stepEdep
               IF ( JTRACK.EQ.3 .OR. JTRACK.EQ.4 ) THEN
                  longitudinalProfile_electron( iReadoutLayer ) =
     &                 longitudinalProfile_electron( iReadoutLayer ) + 
     &                 stepEdep                  
               ELSE IF ( JTRACK.EQ.10 .OR. JTRACK.EQ.11 ) THEN
                  longitudinalProfile_muon( iReadoutLayer ) =
     &                 longitudinalProfile_muon( iReadoutLayer ) + 
     &                 stepEdep                  
               ELSE IF ( JTRACK.EQ.15 .OR. JTRACK.EQ.16 ) THEN
                  longitudinalProfile_kaon( iReadoutLayer ) =
     &                 longitudinalProfile_kaon( iReadoutLayer ) + 
     &                 stepEdep                  
               ELSE IF ( JTRACK.EQ.1 .OR. JTRACK.EQ.2 ) THEN
                  longitudinalProfile_proton( iReadoutLayer ) =
     &                 longitudinalProfile_proton( iReadoutLayer ) + 
     &                 stepEdep                  
               ELSE IF ( JTRACK.EQ.13 .OR. JTRACK.EQ.14 ) THEN
                  longitudinalProfile_pion( iReadoutLayer ) =
     &                 longitudinalProfile_pion( iReadoutLayer ) + 
     &                 stepEdep                  
               ELSE IF ( JTRACK.EQ.8 .OR. JTRACK.EQ.9 .OR.     ! neutrons/antineutrons
     &                   JTRACK.LE.-2 .OR.                     ! ions
     &                   JTRACK.EQ.20 .OR. JTRACK.EQ.21 .OR.   ! Sigma-/+
     &                   JTRACK.EQ.31 .OR. JTRACK.EQ.33 .OR.   ! AntiSigma-/+
     &                   JTRACK.EQ.36 .OR. JTRACK.EQ.37 .OR.   ! Xsi-, AntiXsi-
     &                   JTRACK.EQ.38 .OR. JTRACK.EQ.39 ) THEN ! Omega-, AntiOmega-
                  longitudinalProfile_ion( iReadoutLayer ) =
     &                 longitudinalProfile_ion( iReadoutLayer ) + 
     &                 stepEdep                  
               ENDIF

               radius = SQRT( XTRACK(1)**2 + YTRACK(1)**2 )
               iBinRadius = 1
               currentRadius = radiusBin;
               DO WHILE ( ( radius.GT.currentRadius ).AND.
     &              ( iBinRadius.LT.numberOfRadiusBins ) )
                  currentRadius = currentRadius + iBinRadius*radiusBin
                  iBinRadius = iBinRadius + 1
               ENDDO
CCC               print *, ' MGDRAW : r=', radius, 
CCC     &              ' iBinRadius=', iBinRadius

               transverseProfile( iBinRadius ) = 
     &              transverseProfile( iBinRadius ) + stepEdep
               IF ( JTRACK.EQ.3 .OR. JTRACK.EQ.4 ) THEN
                  transverseProfile_electron( iBinRadius ) =
     &                 transverseProfile_electron( iBinRadius ) + 
     &                 stepEdep                  
               ELSE IF ( JTRACK.EQ.10 .OR. JTRACK.EQ.11 ) THEN
                  transverseProfile_muon( iBinRadius ) =
     &                 transverseProfile_muon( iBinRadius ) + 
     &                 stepEdep                  
               ELSE IF ( JTRACK.EQ.15 .OR. JTRACK.EQ.16 ) THEN
                  transverseProfile_kaon( iBinRadius ) =
     &                 transverseProfile_kaon( iBinRadius ) + 
     &                 stepEdep                  
               ELSE IF ( JTRACK.EQ.1 .OR. JTRACK.EQ.2 ) THEN
                  transverseProfile_proton( iBinRadius ) =
     &                 transverseProfile_proton( iBinRadius ) + 
     &                 stepEdep                  
               ELSE IF ( JTRACK.EQ.13 .OR. JTRACK.EQ.14 ) THEN
                  transverseProfile_pion( iBinRadius ) =
     &                 transverseProfile_pion( iBinRadius ) + 
     &                 stepEdep                  
               ELSE IF ( JTRACK.EQ.8 .OR. JTRACK.EQ.9 .OR.     ! neutrons/antineutrons
     &                   JTRACK.LE.-2 .OR.                     ! ions
     &                   JTRACK.EQ.20 .OR. JTRACK.EQ.21 .OR.   ! Sigma-/+
     &                   JTRACK.EQ.31 .OR. JTRACK.EQ.33 .OR.   ! AntiSigma-/+
     &                   JTRACK.EQ.36 .OR. JTRACK.EQ.37 .OR.   ! Xsi-, AntiXsi-
     &                   JTRACK.EQ.38 .OR. JTRACK.EQ.39 ) THEN ! Omega-, AntiOmega-
                  transverseProfile_ion( iBinRadius ) =
     &                 transverseProfile_ion( iBinRadius ) + 
     &                 stepEdep                  
               ENDIF

            ENDIF
         ENDIF
      ENDIF

C---- Count the tracks created inside the calorimeter.
      IF ( MREG.GT.2 ) THEN
C------- Because the proper counting of all particles is quite 
C        CPU consuming (the overhead can be more than a factor of 10),
C        by default we skip from the counting the following
C        electromagnetic particles:  e- , e+ , gamma .
C        Excluding these particles, the overhead for counting all
C        other particles is negligible.
C        We only do the full counting when  isFullParticleCountOn  
C        is true.
         found = .TRUE.
         IF ( isFullParticleCountOn .OR.
     &        ( JTRACK.NE.3 .AND. JTRACK.NE.4 .AND. JTRACK.NE.7 ) 
     &      ) THEN
C---------- Check whether the track has been already considered, 
C           by looking at its track number in the vector vecTrkNum, 
C           whose effective length is indVecTracks. 
            found = .FALSE.
            DO I = 1, indVecTracks
               IF ( vecTrkNum( I ).EQ.ISPUSR(MKBMX2) ) THEN
                  found = .TRUE.
               ENDIF
            ENDDO
         ENDIF
         IF ( .NOT.found ) THEN
C---------- The track has not yet been considered: store its
C           track number in the vector, and then update the counting
C           of particle types in the vector vecParticles.
            indVecTracks = indVecTracks + 1
            IF ( indVecTracks.LE.maxNumTracksPerEvent ) THEN
               vecTrkNum( indVecTracks ) = ISPUSR(MKBMX2)
C------------- We use the Fluka particle code as index of the vector.
C              There are 64 particles with positive index, from 1 to 64
C              (including some reserved codes), but there are some with
C              non positive index: for values from 0 to -6, we use the
C              last 7 elements of the  vecParticles  (that has 100 
C              elements). We use the element  90  to collect all the
C              particles with Fluka codes < -6 .
               indVecParticles = JTRACK
               IF ( JTRACK.LE.0 .AND. JTRACK.GE.-6 ) THEN
                  indVecParticles = 100 + JTRACK
               ELSE IF ( JTRACK.LT.-6 ) THEN
                  indVecParticles = 90
               ENDIF
               vecParticles( indVecParticles ) = 
     &              vecParticles( indVecParticles ) + 1.0
            ELSE
               indVecTracks = maxNumTracksPerEvent
               print *, 
     &          '***IGNORED track above maxNumTracksPerEvent=',
     &          maxNumTracksPerEvent
            ENDIF
         ENDIF
      ENDIF

C---- Debugging for the primary
CCC      IF ( LTRACK.eq.1 ) THEN
CCC         print *, ' MGDRAW : mreg=', MREG, ' (', 
CCC     &        XTRACK(1), ',', YTRACK(1), ',', ZTRACK(1), ')',
CCC     &        ' trk_id=', ISPUSR(MKBMX2), ' L=', TTRACK(1),
CCC     &        ' edep=', DTRACK(1)           
CCC      ENDIF



      IF ( .NOT. LFCOPE ) THEN
         LFCOPE = .TRUE.
         IF ( KOMPUT .EQ. 2 ) THEN
            FILNAM = '/'//CFDRAW(1:8)//' DUMP A'
         ELSE
            FILNAM = CFDRAW
         END IF
CCCALB         OPEN ( UNIT = IODRAW, FILE = FILNAM, STATUS = 'NEW', FORM =
CCCALB     &          'UNFORMATTED' )
      END IF
CCCALB      WRITE (IODRAW) NTRACK, MTRACK, JTRACK, SNGL (ETRACK),
CCCALB     &               SNGL (WTRACK)
CCCALB      WRITE (IODRAW) ( SNGL (XTRACK (I)), SNGL (YTRACK (I)),
CCCALB     &                 SNGL (ZTRACK (I)), I = 0, NTRACK ),
CCCALB     &               ( SNGL (DTRACK (I)), I = 1,MTRACK ),
CCCALB     &                 SNGL (CTRACK)

C
      RETURN
*
*======================================================================*
*                                                                      *
*     Boundary-(X)crossing DRAWing:                                    *
*                                                                      *
*     Icode = 1x: call from Kaskad                                     *
*             19: boundary crossing                                    *
*     Icode = 2x: call from Emfsco                                     *
*             29: boundary crossing                                    *
*     Icode = 3x: call from Kasneu                                     *
*             39: boundary crossing                                    *
*     Icode = 4x: call from Kashea                                     *
*             49: boundary crossing                                    *
*     Icode = 5x: call from Kasoph                                     *
*             59: boundary crossing                                    *
*                                                                      *
*======================================================================*
*                                                                      *
      ENTRY BXDRAW ( ICODE, MREG, NEWREG, XSCO, YSCO, ZSCO )
 
      IF ( MREG.EQ.2 .AND. NEWREG.EQ.1 ) THEN

         Ekin = -999.9
         IF ( JTRACK.GE.-6 ) THEN
            Ekin = ETRACK-AM(JTRACK)
         ENDIF

C------- Also in this case, exactly as we did in MGDRAW for counting
C        the particles produced, we use the Fluka particle code as 
C        index of the vector.
         indVecParticles = JTRACK
         IF ( JTRACK.LE.0 .AND. JTRACK.GE.-6 ) THEN
            indVecParticles = 100 + JTRACK
         ELSE IF ( JTRACK.LT.-6 ) THEN
            indVecParticles = 90
         ENDIF

         vecCountExitingParticles( indVecParticles ) = 
     &        vecCountExitingParticles( indVecParticles ) + 1.0

         vecEkinExitingParticles( indVecParticles ) = 
     &        vecEkinExitingParticles( indVecParticles ) + Ekin
         
CCC         print *, 'Exit track=', ISPUSR(MKBMX2),
CCC     &        ' ', PRNAME( JTRACK ), 
CCC     &        '(', XSCO, ',', YSCO, ',', ZSCO, ')',
CCC     &        ' Ekin=', Ekin, ' GeV'

      ENDIF

      RETURN
*
*======================================================================*
*                                                                      *
*     Event End DRAWing:                                               *
*                                                                      *
*======================================================================*
*                                                                      *
      ENTRY EEDRAW ( ICODE )
      RETURN
*
*======================================================================*
*                                                                      *
*     ENergy deposition DRAWing:                                       *
*                                                                      *
*     Icode = 1x: call from Kaskad                                     *
*             10: elastic interaction recoil                           *
*             11: inelastic interaction recoil                         *
*             12: stopping particle                                    *
*             13: pseudo-neutron deposition                            *
*             14: escape                                               *
*             15: time kill                                            *
*     Icode = 2x: call from Emfsco                                     *
*             20: local energy deposition (i.e. photoelectric)         *
*             21: below threshold, iarg=1                              *
*             22: below threshold, iarg=2                              *
*             23: escape                                               *
*             24: time kill                                            *
*     Icode = 3x: call from Kasneu                                     *
*             30: target recoil                                        *
*             31: below threshold                                      *
*             32: escape                                               *
*             33: time kill                                            *
*     Icode = 4x: call from Kashea                                     *
*             40: escape                                               *
*             41: time kill                                            *
*             42: delta ray stack overflow                             *
*     Icode = 5x: call from Kasoph                                     *
*             50: optical photon absorption                            *
*             51: escape                                               *
*             52: time kill                                            *
*                                                                      *
*======================================================================*
*                                                                      *
      ENTRY ENDRAW ( ICODE, MREG, RULL, XSCO, YSCO, ZSCO )


C---- If the discrete energy deposition happens inside the 
C     calorimeter then its value is added to the variable
C     edepAllCalo, which keeps the sum of all energy depositions
C     in the current event.
C     For edepAllCalo, we do not apply the Birks quenching even
C     if it is active and the energy deposition happens in an 
C     active layer.
C     Notice that the active layers correspond to regions with
C     even number, excluding region number 2 (experimental hall).

      IF ( MREG.GT.2 ) THEN

CCC         IF ( LTRACK.eq.1 ) then ! Only for primary
CCC            print *, ' ENDRAW : mreg=', MREG, ' (', 
CCC     &           XSCO, ',', YSCO, ',', ZSCO, ')'
CCC         ENDIF
     
         stepEdep = RULL

         edepAllCalo = edepAllCalo + stepEdep
         IF ( JTRACK.EQ.3 .OR. JTRACK.EQ.4 ) THEN
            edepAllCalo_electron = edepAllCalo_electron + stepEdep
         ELSE IF ( JTRACK.EQ.10 .OR. JTRACK.EQ.11 ) THEN
            edepAllCalo_muon = edepAllCalo_muon + stepEdep
         ELSE IF ( JTRACK.EQ.15 .OR. JTRACK.EQ.16 ) THEN
            edepAllCalo_kaon = edepAllCalo_kaon + stepEdep
         ELSE IF ( JTRACK.EQ.1 .OR. JTRACK.EQ.2 ) THEN
            edepAllCalo_proton = edepAllCalo_proton + stepEdep
         ELSE IF ( JTRACK.EQ.13 .OR. JTRACK.EQ.14 ) THEN
            edepAllCalo_pion = edepAllCalo_pion + stepEdep
         ELSE IF ( JTRACK.EQ.8 .OR. JTRACK.EQ.9 .OR. ! neutrons/antineutrons
     &             JTRACK.LE.-2 .OR. ! ions
     &             JTRACK.EQ.20 .OR. JTRACK.EQ.21 .OR. ! Sigma-/+
     &             JTRACK.EQ.31 .OR. JTRACK.EQ.33 .OR. ! AntiSigma-/+
     &             JTRACK.EQ.36 .OR. JTRACK.EQ.37 .OR. ! Xsi-, AntiXsi-
     &             JTRACK.EQ.38 .OR. JTRACK.EQ.39 ) THEN ! Omega-, AntiOmega-
            edepAllCalo_ion = edepAllCalo_ion + stepEdep
         ENDIF

         IF ( MOD(MREG,2).EQ.0 ) THEN 

            IF ( isBirksOn ) THEN
               RULLL = RULL
               CALL QUENMG( ICODE, MREG, RULLL, DTQUEN )
               stepEdep = DTQUEN(1,1)
            ELSE 
               stepEdep = RULL
            ENDIF

CCC            print *, ' discrete: stepEdep=', stepEdep

            edepVisCalo = edepVisCalo + stepEdep
            IF ( JTRACK.EQ.3 .OR. JTRACK.EQ.4 ) THEN
               edepVisCalo_electron = edepVisCalo_electron + stepEdep
            ELSE IF ( JTRACK.EQ.10 .OR. JTRACK.EQ.11 ) THEN
               edepVisCalo_muon = edepVisCalo_muon + stepEdep
            ELSE IF ( JTRACK.EQ.15 .OR. JTRACK.EQ.16 ) THEN
               edepVisCalo_kaon = edepVisCalo_kaon + stepEdep
            ELSE IF ( JTRACK.EQ.1 .OR. JTRACK.EQ.2 ) THEN
               edepVisCalo_proton = edepVisCalo_proton + stepEdep
            ELSE IF ( JTRACK.EQ.13 .OR. JTRACK.EQ.14 ) THEN
               edepVisCalo_pion = edepVisCalo_pion + stepEdep
            ELSE IF ( JTRACK.EQ.8 .OR. JTRACK.EQ.9 .OR. ! neutrons/antineutrons
     &                JTRACK.LE.-2 .OR. ! ions
     &                JTRACK.EQ.20 .OR. JTRACK.EQ.21 .OR. ! Sigma-/+
     &                JTRACK.EQ.31 .OR. JTRACK.EQ.33 .OR. ! AntiSigma-/+
     &                JTRACK.EQ.36 .OR. JTRACK.EQ.37 .OR. ! Xsi-, AntiXsi-
     &                JTRACK.EQ.38 .OR. JTRACK.EQ.39 ) THEN ! Omega-, AntiOmega-
               edepVisCalo_ion = edepVisCalo_ion + stepEdep
            ENDIF

            iReadoutLayer = 1 + 
     &           (MREG - 4) / ( 2 * numberOfLayersPerReadoutLayer )
CCC            print *, ' ENDRAW : MREG=', MREG, 
CCC    &           ' iReadouLayer=', iReadoutLayer

            longitudinalProfile( iReadoutLayer ) =
     &           longitudinalProfile( iReadoutLayer ) + stepEdep
            IF ( JTRACK.EQ.3 .OR. JTRACK.EQ.4 ) THEN
               longitudinalProfile_electron( iReadoutLayer ) =
     &              longitudinalProfile_electron( iReadoutLayer ) + 
     &              stepEdep                  
            ELSE IF ( JTRACK.EQ.10 .OR. JTRACK.EQ.11 ) THEN
               longitudinalProfile_muon( iReadoutLayer ) =
     &              longitudinalProfile_muon( iReadoutLayer ) + 
     &              stepEdep                  
            ELSE IF ( JTRACK.EQ.15 .OR. JTRACK.EQ.16 ) THEN
               longitudinalProfile_kaon( iReadoutLayer ) =
     &              longitudinalProfile_kaon( iReadoutLayer ) + 
     &              stepEdep                  
            ELSE IF ( JTRACK.EQ.1 .OR. JTRACK.EQ.2 ) THEN
               longitudinalProfile_proton( iReadoutLayer ) =
     &              longitudinalProfile_proton( iReadoutLayer ) + 
     &              stepEdep                  
            ELSE IF ( JTRACK.EQ.13 .OR. JTRACK.EQ.14 ) THEN
               longitudinalProfile_pion( iReadoutLayer ) =
     &              longitudinalProfile_pion( iReadoutLayer ) + 
     &              stepEdep                  
            ELSE IF ( JTRACK.EQ.8 .OR. JTRACK.EQ.9 .OR. ! neutrons/antineutrons
     &                JTRACK.LE.-2 .OR. ! ions
     &                JTRACK.EQ.20 .OR. JTRACK.EQ.21 .OR. ! Sigma-/+
     &                JTRACK.EQ.31 .OR. JTRACK.EQ.33 .OR. ! AntiSigma-/+
     &                JTRACK.EQ.36 .OR. JTRACK.EQ.37 .OR. ! Xsi-, AntiXsi-
     &                JTRACK.EQ.38 .OR. JTRACK.EQ.39 ) THEN ! Omega-, AntiOmega-
               longitudinalProfile_ion( iReadoutLayer ) =
     &              longitudinalProfile_ion( iReadoutLayer ) + 
     &              stepEdep                  
            ENDIF

            radius = SQRT( XSCO**2 + YSCO**2 )
            iBinRadius = 1
            currentRadius = radiusBin;
            DO WHILE ( ( radius.GT.currentRadius ).AND.
     &           ( iBinRadius.LT.numberOfRadiusBins ) )
               currentRadius = currentRadius + iBinRadius*radiusBin
               iBinRadius = iBinRadius + 1
            ENDDO
CCC            print *, ' ENDRAW : r=', radius, 
CCC     &           ' iBinRadius=', iBinRadius

            transverseProfile( iBinRadius ) = 
     &           transverseProfile( iBinRadius ) + stepEdep
            IF ( JTRACK.EQ.3 .OR. JTRACK.EQ.4 ) THEN
               transverseProfile_electron( iBinRadius ) =
     &              transverseProfile_electron( iBinRadius ) + 
     &              stepEdep                  
            ELSE IF ( JTRACK.EQ.10 .OR. JTRACK.EQ.11 ) THEN
               transverseProfile_muon( iBinRadius ) =
     &              transverseProfile_muon( iBinRadius ) + 
     &              stepEdep                  
            ELSE IF ( JTRACK.EQ.15 .OR. JTRACK.EQ.16 ) THEN
               transverseProfile_kaon( iBinRadius ) =
     &              transverseProfile_kaon( iBinRadius ) + 
     &              stepEdep                  
            ELSE IF ( JTRACK.EQ.1 .OR. JTRACK.EQ.2 ) THEN
               transverseProfile_proton( iBinRadius ) =
     &              transverseProfile_proton( iBinRadius ) + 
     &              stepEdep                  
            ELSE IF ( JTRACK.EQ.13 .OR. JTRACK.EQ.14 ) THEN
               transverseProfile_pion( iBinRadius ) =
     &              transverseProfile_pion( iBinRadius ) + 
     &              stepEdep                  
            ELSE IF ( JTRACK.EQ.8 .OR. JTRACK.EQ.9 .OR. ! neutrons/antineutrons
     &                JTRACK.LE.-2 .OR. ! ions
     &                JTRACK.EQ.20 .OR. JTRACK.EQ.21 .OR. ! Sigma-/+
     &                JTRACK.EQ.31 .OR. JTRACK.EQ.33 .OR. ! AntiSigma-/+
     &                JTRACK.EQ.36 .OR. JTRACK.EQ.37 .OR. ! Xsi-, AntiXsi-
     &                JTRACK.EQ.38 .OR. JTRACK.EQ.39 ) THEN ! Omega-, AntiOmega-
               transverseProfile_ion( iBinRadius ) =
     &              transverseProfile_ion( iBinRadius ) + 
     &              stepEdep                  
            ENDIF   
         ENDIF
      ENDIF

C---- Debugging for the primary
CCC      IF ( LTRACK.eq.1 ) THEN
CCC         print *, ' ENDRAW : mreg=', MREG, ' (', 
CCC     &        XSCO, ',', YSCO, ',', ZSCO, ')',
CCC     &        ' trk_id=', ISPUSR(MKBMX2), ' L=', TTRACK(1),
CCC     &        ' edep=', RULL           
CCC      ENDIF

C---- Count the tracks created inside the calorimeter: this part
C     matters only for electrons, positrons, and gammas, so it is
C     done only if  isFullParticleCountOn  is true.
      IF ( MREG.GT.2 .AND. isFullParticleCountOn ) THEN
         found = .FALSE.
         DO I = 1, indVecTracks
            IF ( vecTrkNum( I ).EQ.ISPUSR(MKBMX2) ) THEN
               found = .TRUE.
            ENDIF
         ENDDO
         IF ( .NOT.found ) THEN 
            indVecTracks = indVecTracks + 1
            IF ( indVecTracks.LE.maxNumTracksPerEvent ) THEN
               vecTrkNum( indVecTracks ) = ISPUSR(MKBMX2)
               indVecParticles = JTRACK
               IF ( JTRACK.LE.0 .AND. JTRACK.GE.-6 ) THEN
                  indVecParticles = 100 + JTRACK
               ELSE IF ( JTRACK.LT.-6 ) THEN
                  indVecParticles = 90
               ENDIF
               vecParticles( indVecParticles ) = 
     &              vecParticles( indVecParticles ) + 1.0
            ELSE
               indVecTracks = maxNumTracksPerEvent
               print *, 
     &          '***IGNORED track above maxNumTracksPerEvent=',
     &          maxNumTracksPerEvent
            ENDIF
         ENDIF
      ENDIF

      IF ( .NOT. LFCOPE ) THEN
         LFCOPE = .TRUE.
         IF ( KOMPUT .EQ. 2 ) THEN
            FILNAM = '/'//CFDRAW(1:8)//' DUMP A'
         ELSE
            FILNAM = CFDRAW
         END IF
CCCALB         OPEN ( UNIT = IODRAW, FILE = FILNAM, STATUS = 'NEW', FORM =
CCCALB     &          'UNFORMATTED' )
      END IF
CCCALB      WRITE (IODRAW)  0, ICODE, JTRACK, SNGL (ETRACK), SNGL (WTRACK)
CCCALB      WRITE (IODRAW)  SNGL (XSCO), SNGL (YSCO), SNGL (ZSCO), SNGL (RULL)
      RETURN
*
*======================================================================*
*                                                                      *
*     SOurce particle DRAWing:                                         *
*                                                                      *
*======================================================================*
*
      ENTRY SODRAW
      IF ( .NOT. LFCOPE ) THEN
         LFCOPE = .TRUE.
         IF ( KOMPUT .EQ. 2 ) THEN
            FILNAM = '/'//CFDRAW(1:8)//' DUMP A'
         ELSE
            FILNAM = CFDRAW
         END IF
CCCALB         OPEN ( UNIT = IODRAW, FILE = FILNAM, STATUS = 'NEW', FORM =
CCCALB     &          'UNFORMATTED' )
      END IF
CCCALB      WRITE (IODRAW) -NCASE, NPFLKA, NSTMAX, SNGL (TKESUM),
CCCALB     &                SNGL (WEIPRI)
*  +-------------------------------------------------------------------*
*  |  Patch for heavy ions: it works only for 1 source particle on
*  |  the stack for the time being
      IF ( ABS (ILOFLK (NPFLKA)) .GE. 10000 ) THEN
         IONID = ILOFLK (NPFLKA)
         CALL DCDION ( IONID )
CCCALB         WRITE (IODRAW) ( IONID,SNGL(TKEFLK(I)+AMNHEA(-IONID)),
CCCALB     &                    SNGL (WTFLK(I)), SNGL (XFLK (I)),
CCCALB     &                    SNGL (YFLK (I)), SNGL (ZFLK (I)),
CCCALB     &                    SNGL (TXFLK(I)), SNGL (TYFLK(I)),
CCCALB     &                    SNGL (TZFLK(I)), I = 1, NPFLKA )
*  |
*  +-------------------------------------------------------------------*
*  |
      ELSE
CCCALB         WRITE (IODRAW) ( ILOFLK(I), SNGL (TKEFLK(I)+AM(ILOFLK(I))),
CCCALB     &                    SNGL (WTFLK(I)), SNGL (XFLK (I)),
CCCALB     &                    SNGL (YFLK (I)), SNGL (ZFLK (I)),
CCCALB     &                    SNGL (TXFLK(I)), SNGL (TYFLK(I)),
CCCALB     &                    SNGL (TZFLK(I)), I = 1, NPFLKA )
      END IF
*  |
*  +-------------------------------------------------------------------*
      RETURN
*
*======================================================================*
*                                                                      *
*     USer dependent DRAWing:                                          *
*                                                                      *
*     Icode = 10x: call from Kaskad                                    *
*             100: elastic   interaction secondaries                   *
*             101: inelastic interaction secondaries                   *
*             102: particle decay  secondaries                         *
*             103: delta ray  generation secondaries                   *
*             104: pair production secondaries                         *
*             105: bremsstrahlung  secondaries                         *
*     Icode = 20x: call from Emfsco                                    *
*             208: bremsstrahlung secondaries                          *
*             210: Moller secondaries                                  *
*             212: Bhabha secondaries                                  *
*             214: in-flight annihilation secondaries                  *
*             215: annihilation at rest   secondaries                  *
*             217: pair production        secondaries                  *
*             219: Compton scattering     secondaries                  *
*             221: photoelectric          secondaries                  *
*             225: Rayleigh scattering    secondaries                  *
*     Icode = 30x: call from Kasneu                                    *
*             300: interaction secondaries                             *
*     Icode = 40x: call from Kashea                                    *
*             400: delta ray  generation secondaries                   *
*  For all interactions secondaries are put on GENSTK common (kp=1,np) *
*  but for KASHEA delta ray generation where only the secondary elec-  *
*  tron is present and stacked on FLKSTK common for kp=lstack          *
*                                                                      *
*======================================================================*
*
      ENTRY USDRAW ( ICODE, MREG, XSCO, YSCO, ZSCO )
      IF ( .NOT. LFCOPE ) THEN
         LFCOPE = .TRUE.
         IF ( KOMPUT .EQ. 2 ) THEN
            FILNAM = '/'//CFDRAW(1:8)//' DUMP A'
         ELSE
            FILNAM = CFDRAW
         END IF
CCCALB         OPEN ( UNIT = IODRAW, FILE = FILNAM, STATUS = 'NEW', FORM =
CCCALB     &          'UNFORMATTED' )
      END IF

* No output by default:
      RETURN
*=== End of subrutine Mgdraw ==========================================*
      END
