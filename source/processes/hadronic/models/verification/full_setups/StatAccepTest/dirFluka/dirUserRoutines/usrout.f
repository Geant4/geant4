*$ CREATE USROUT.FOR
*COPY USROUT
*
*=== usrout ===========================================================*
*
      SUBROUTINE USROUT ( WHAT, SDUM )

      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'
      INCLUDE '(PAPROP)'
*
      INCLUDE 'myinclude.inc'
*
*----------------------------------------------------------------------*
*                                                                      *
*     Copyright (C) 1991-2005      by    Alfredo Ferrari & Paola Sala  *
*     All Rights Reserved.                                             *
*                                                                      *
*                                                                      *
*     USeR OUTput: this routine is called every time the USROCALL card *
*                  is found in the input stream                        *
*                                                                      *
*                                                                      *
*     Created on 01 january 1991   by    Alfredo Ferrari & Paola Sala  *
*                                                   Infn - Milan       *
*                                                                      *
*     Last change on 20-mar-05     by    Alfredo Ferrari               *
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*
      DIMENSION WHAT (6)
      CHARACTER SDUM*8

      INTEGER numEntriesForEnergyResolution
      REAL mu, sigma, mu_sigma
      REAL mu_Evis, mu_Evis_sigma
      REAL sum, sum2
 
      REAL sumAll, sumAll2, sumVis, sumVis2
      REAL vecL( numberOfReadoutLayers )
      REAL vecL2( numberOfReadoutLayers )
      REAL vecR( numberOfRadiusBins )
      REAL vecR2( numberOfRadiusBins )
      REAL totalVisibleEnergy, sigma_totalVisibleEnergy
      REAL totalDepositedEnergy
      CHARACTER*8 flukaParticleName
      CHARACTER*70 caseName

      REAL fractionLongitudinal1stQuarter
      REAL fractionLongitudinal2ndQuarter
      REAL fractionLongitudinal3rdQuarter
      REAL fractionLongitudinal4thQuarter
      REAL fractionLongitudinal1stQuarter_sigma
      REAL fractionLongitudinal2ndQuarter_sigma
      REAL fractionLongitudinal3rdQuarter_sigma
      REAL fractionLongitudinal4thQuarter_sigma

      REAL fractionTransverse1stThird 
      REAL fractionTransverse2ndThird
      REAL fractionTransverse3rdThird
      REAL fractionTransverse1stThird_sigma
      REAL fractionTransverse2ndThird_sigma
      REAL fractionTransverse3rdThird_sigma

      REAL width_Evis, width_Evis_sigma
      REAL energyResolution, energyResolution_sigma

      INTEGER flukaParticleCode

      INTEGER totNumExitingParticles
      REAL totEkinExitingParticles, averageEkin

C---- Print out the average total energy deposit per event
C     and its statistical error.

      print *, 'USROUT--- '

      if ( numEvent.LE.1 ) then
         numEvent = 2
      endif

      totalVisibleEnergy = 0.0
      sigma_totalVisibleEnergy = 0.0
      totalDepositedEnergy = 0.0

C---- Consider the following 7 cases:
C       1) total (no particle type distinction);
C       2) electrons & positrons
C       3) muon- & muon+
C       4) kaon- & kaon+
C       5) protons & antiprotons
C       6) pion- & pion+
C       7) ions (including neutrons and hyperones).
      DO iCase = 1, 7  

         IF ( iCase.EQ.1 ) THEN
            caseName = " "
            sumVis  = sumEdep_vis
            sumVis2 = sumEdep_vis2
            sumAll  = sumEdep_all
            sumAll2 = sumEdep_all2
            DO I = 1, numberOfReadoutLayers
               vecL( I )  = sumL( I )
               vecL2( I ) = sumL2( I )
            ENDDO            
            DO I = 1, numberOfRadiusBins
               vecR( I )  = sumR( I )  
               vecR2( I ) = sumR2( I )
            ENDDO
CCC            print *, ' ------------ vecEvis() in MeV -------------- '
CCC            DO I = 1, maxNumberOfEvents
CCC               if ( vecEvis( I ).GT.1.0E-09 ) then
CCC                  print *, ' evt=', I, ' Evis=', 1000.0*vecEvis( I )
CCC               endif
CCC            ENDDO
C         
         ELSE IF ( iCase.EQ.2 ) THEN
            caseName = " --- Electrons / Positrons --- "
            sumVis  = sumEdep_vis_electron
            sumVis2 = sumEdep_vis_electron2
            sumAll  = sumEdep_all_electron
            sumAll2 = sumEdep_all_electron2
            DO I = 1, numberOfReadoutLayers
               vecL( I )  = sumL_electron( I )
               vecL2( I ) = sumL_electron2( I )
            ENDDO            
            DO I = 1, numberOfRadiusBins
               vecR( I )  = sumR_electron( I )  
               vecR2( I ) = sumR_electron2( I )
            ENDDO
C
         ELSE IF ( iCase.EQ.3 ) THEN
            caseName = " --- Muon- / Muon+ --- "
            sumVis  = sumEdep_vis_muon
            sumVis2 = sumEdep_vis_muon2
            sumAll  = sumEdep_all_muon
            sumAll2 = sumEdep_all_muon2
            DO I = 1, numberOfReadoutLayers
               vecL( I )  = sumL_muon( I )
               vecL2( I ) = sumL_muon2( I )
            ENDDO            
            DO I = 1, numberOfRadiusBins
               vecR( I )  = sumR_muon( I )  
               vecR2( I ) = sumR_muon2( I )
            ENDDO
C
         ELSE IF ( iCase.EQ.4 ) THEN
            caseName = " --- Kaon+ / Kaon- --- "
            sumVis  = sumEdep_vis_kaon
            sumVis2 = sumEdep_vis_kaon2
            sumAll  = sumEdep_all_kaon
            sumAll2 = sumEdep_all_kaon2
            DO I = 1, numberOfReadoutLayers
               vecL( I )  = sumL_kaon( I )
               vecL2( I ) = sumL_kaon2( I )
            ENDDO            
            DO I = 1, numberOfRadiusBins
               vecR( I )  = sumR_kaon( I )  
               vecR2( I ) = sumR_kaon2( I )
            ENDDO
C
         ELSE IF ( iCase.EQ.5 ) THEN
            caseName = " --- Protons / Antiprotons --- "
            sumVis  = sumEdep_vis_proton
            sumVis2 = sumEdep_vis_proton2
            sumAll  = sumEdep_all_proton
            sumAll2 = sumEdep_all_proton2
            DO I = 1, numberOfReadoutLayers
               vecL( I )  = sumL_proton( I )
               vecL2( I ) = sumL_proton2( I )
            ENDDO            
            DO I = 1, numberOfRadiusBins
               vecR( I )  = sumR_proton( I )  
               vecR2( I ) = sumR_proton2( I )
            ENDDO
C
         ELSE IF ( iCase.EQ.6 ) THEN
            caseName = " --- Pion+ / Pion- --- "
            sumVis  = sumEdep_vis_pion
            sumVis2 = sumEdep_vis_pion2
            sumAll  = sumEdep_all_pion
            sumAll2 = sumEdep_all_pion2
            DO I = 1, numberOfReadoutLayers
               vecL( I )  = sumL_pion( I )
               vecL2( I ) = sumL_pion2( I )
            ENDDO            
            DO I = 1, numberOfRadiusBins
               vecR( I )  = sumR_pion( I )  
               vecR2( I ) = sumR_pion2( I )
            ENDDO
C
         ELSE IF ( iCase.EQ.7 ) THEN
            caseName = " --- Ions (and neutrons & hyperons ) --- "
            sumVis  = sumEdep_vis_ion
            sumVis2 = sumEdep_vis_ion2
            sumAll  = sumEdep_all_ion
            sumAll2 = sumEdep_all_ion2
            DO I = 1, numberOfReadoutLayers
               vecL( I )  = sumL_ion( I )
               vecL2( I ) = sumL_ion2( I )
            ENDDO            
            DO I = 1, numberOfRadiusBins
               vecR( I )  = sumR_ion( I )  
               vecR2( I ) = sumR_ion2( I )
            ENDDO
C
         ENDIF

         print *, caseName

         sum  = sumVis
         sum2 = sumVis2
         mu = sum / (1.0*numEvent)
         sigma = sqrt( abs( sum2 - sum*sum/(1.0*numEvent) ) 
     &        / ( 1.0*numEvent - 1.0 ) )
         mu_sigma = sigma / sqrt( 1.0*numEvent )
         IF ( iCase.EQ.1 ) THEN
            totalVisibleEnergy = mu
            sigma_totalVisibleEnergy = mu_sigma
            IF ( totalVisibleEnergy.LE.0.0 ) THEN
               totalVisibleEnergy = 1.0
            ENDIF
            print *, ' <E_vis> = ', 1000.0*mu, ' +/- ', 
     &           1000.0*mu_sigma, ' MeV' 
         ELSE
            print *, ' <E_vis> = ', 1000.0*mu, ' +/- ', 
     &           1000.0*mu_sigma, ' MeV  (', 
     &           100.0*mu/totalVisibleEnergy, ' %)' 
         ENDIF
         mu_Evis = mu;          ! For later usage.
         mu_Evis_sigma = mu_sigma; !  "    "     "
       
         sum  = sumAll
         sum2 = sumAll2 
         mu = sum / (1.0*numEvent)
         sigma = sqrt( abs( sum2 - sum*sum/(1.0*numEvent) ) 
     &        / ( 1.0*numEvent - 1.0 ) )
         mu_sigma = sigma / sqrt( 1.0*numEvent )
         IF ( iCase.EQ.1 ) THEN
            totalDepositedEnergy = mu
            IF ( totalDepositedEnergy.LE.0.0 ) THEN
               totalDepositedEnergy = 1.0
            ENDIF
            print *, ' <E_all> = ', 1000.0*mu, ' +/- ', 
     &           1000.0*mu_sigma, ' MeV' 
         ELSE
            print *, ' <E_all> = ', 1000.0*mu, ' +/- ', 
     &           1000.0*mu_sigma, ' MeV (', 
     &           100.0*mu/totalDepositedEnergy, ' %)' 
         ENDIF

CCC         print *, '  '
CCC         print *, ' ----------- longitudinalInfo() in MeV ---------- '
CCC         DO I = 1, numberOfReadoutLayers
CCC            print *, ' iL=', I, ' vecL=', 1000.0*vecL( I ), 
CCC     &           '  vecL2=', 1000.0*1000.0*vecL2( I ) 
CCC         ENDDO
CCC         
CCC         print *, '  '
CCC         print *, ' ----------- transverseInfo() in MeV ---------- '
CCC         DO I = 1, numberOfRadiusBins 
CCC            print *, ' iR=', I, ' vecR=', 1000.0*vecR( I ), 
CCC     &           '  vecR2=', 1000.0*1000.0*vecR2( I ) 
CCC         ENDDO
      
         fractionLongitudinal1stQuarter = 0.0
         fractionLongitudinal2ndQuarter = 0.0
         fractionLongitudinal3rdQuarter = 0.0
         fractionLongitudinal4thQuarter = 0.0
         fractionLongitudinal1stQuarter_sigma = 0.0
         fractionLongitudinal2ndQuarter_sigma = 0.0
         fractionLongitudinal3rdQuarter_sigma = 0.0
         fractionLongitudinal4thQuarter_sigma = 0.0
         print *, ' Average <E> [MeV] in each Layer ' 
         DO I = 1, numberOfReadoutLayers

            sum  = vecL( I )
            sum2 = vecL2( I )
            mu = sum / (1.0*numEvent)
            sigma = sqrt( abs( sum2 - sum*sum/(1.0*numEvent) ) 
     &           / ( 1.0*numEvent - 1.0 ) )
            mu_sigma = sigma / sqrt( 1.0*numEvent )
            print *, '   layer = ', I, '  E =', 1000.0*mu, ' +/- ', 
     &           1000.0*mu_sigma
            
            IF ( I.LE.(numberOfReadoutLayers/4) ) THEN
               fractionLongitudinal1stQuarter = 
     &              fractionLongitudinal1stQuarter + mu
               fractionLongitudinal1stQuarter_sigma = 
     &              fractionLongitudinal1stQuarter_sigma + 
     &              mu_sigma * mu_sigma
            ELSE IF ( I.LE.(2*numberOfReadoutLayers/4) ) THEN
               fractionLongitudinal2ndQuarter = 
     &              fractionLongitudinal2ndQuarter + mu
               fractionLongitudinal2ndQuarter_sigma = 
     &              fractionLongitudinal2ndQuarter_sigma + 
     &              mu_sigma * mu_sigma
            ELSE IF ( I.LE.(3*numberOfReadoutLayers/4) ) THEN
               fractionLongitudinal3rdQuarter = 
     &              fractionLongitudinal3rdQuarter + mu
               fractionLongitudinal3rdQuarter_sigma = 
     &              fractionLongitudinal3rdQuarter_sigma + 
     &              mu_sigma * mu_sigma
            ELSE
               fractionLongitudinal4thQuarter = 
     &              fractionLongitudinal4thQuarter + mu
               fractionLongitudinal4thQuarter_sigma = 
     &              fractionLongitudinal4thQuarter_sigma +
     &              mu_sigma * mu_sigma
            ENDIF

         ENDDO

         IF ( mu_Evis.GT.1.0E-09 ) THEN

            print *, '  sumL_1 = ', 
     &           1000.0*fractionLongitudinal1stQuarter, 
     &           ' +/- ', 
     &           1000.0*sqrt( fractionLongitudinal1stQuarter_sigma ),
     &           ' MeV'
            print *, '  sumL_2 = ', 
     &           1000.0*fractionLongitudinal2ndQuarter, 
     &           ' +/- ', 
     &           1000.0*sqrt( fractionLongitudinal2ndQuarter_sigma ),
     &           ' MeV'
            print *, '  sumL_3 = ', 
     &           1000.0*fractionLongitudinal3rdQuarter, 
     &           ' +/- ', 
     &           1000.0*sqrt( fractionLongitudinal3rdQuarter_sigma ),
     &           ' MeV'
            print *, '  sumL_4 = ', 
     &           1000.0*fractionLongitudinal4thQuarter, 
     &           ' +/- ', 
     &           1000.0*sqrt( fractionLongitudinal4thQuarter_sigma ),
     &           ' MeV'
            
            IF ( fractionLongitudinal1stQuarter.GT.1.0E-09 ) THEN 
               fractionLongitudinal1stQuarter_sigma =
     &              fractionLongitudinal1stQuarter_sigma /
     &              ( fractionLongitudinal1stQuarter * 
     &              fractionLongitudinal1stQuarter )
            ENDIF
            fractionLongitudinal1stQuarter = 
     &           fractionLongitudinal1stQuarter / mu_Evis
            fractionLongitudinal1stQuarter_sigma = 
     &           fractionLongitudinal1stQuarter *
     &           sqrt( fractionLongitudinal1stQuarter_sigma +
     &                 ( mu_Evis_sigma / mu_Evis ) * 
     &                 ( mu_Evis_sigma / mu_Evis ) )
            
            IF ( fractionLongitudinal2ndQuarter.GT.1.0E-09 ) THEN 
               fractionLongitudinal2ndQuarter_sigma =
     &              fractionLongitudinal2ndQuarter_sigma /
     &              ( fractionLongitudinal2ndQuarter * 
     &              fractionLongitudinal2ndQuarter )
            ENDIF
            fractionLongitudinal2ndQuarter = 
     &           fractionLongitudinal2ndQuarter / mu_Evis
            fractionLongitudinal2ndQuarter_sigma = 
     &           fractionLongitudinal2ndQuarter *
     &           sqrt( fractionLongitudinal2ndQuarter_sigma +
     &                 ( mu_Evis_sigma / mu_Evis ) * 
     &                 ( mu_Evis_sigma / mu_Evis ) )
            
            IF ( fractionLongitudinal3rdQuarter.GT.1.0E-09 ) THEN 
               fractionLongitudinal3rdQuarter_sigma =
     &              fractionLongitudinal3rdQuarter_sigma /
     &              ( fractionLongitudinal3rdQuarter * 
     &              fractionLongitudinal3rdQuarter )
            ENDIF
            fractionLongitudinal3rdQuarter = 
     &           fractionLongitudinal3rdQuarter / mu_Evis
            fractionLongitudinal3rdQuarter_sigma = 
     &           fractionLongitudinal3rdQuarter *
     &           sqrt( fractionLongitudinal3rdQuarter_sigma +
     &                 ( mu_Evis_sigma / mu_Evis ) * 
     &                 ( mu_Evis_sigma / mu_Evis ) )
            
            IF ( fractionLongitudinal4thQuarter.GT.1.0E-09 ) THEN 
               fractionLongitudinal4thQuarter_sigma =
     &              fractionLongitudinal4thQuarter_sigma /
     &              ( fractionLongitudinal4thQuarter * 
     &              fractionLongitudinal4thQuarter )
            ENDIF
            fractionLongitudinal4thQuarter = 
     &           fractionLongitudinal4thQuarter / mu_Evis
            fractionLongitudinal4thQuarter_sigma = 
     &           fractionLongitudinal4thQuarter *
     &           sqrt( fractionLongitudinal4thQuarter_sigma +
     &                 ( mu_Evis_sigma / mu_Evis ) * 
     &                 ( mu_Evis_sigma / mu_Evis ) )

         ENDIF

         print *, ' frac L 1st quarter = ',
     &        fractionLongitudinal1stQuarter*100.0, ' +/- ',
     &        fractionLongitudinal1stQuarter_sigma*100.0, ' % '
         print *, '        2nd quarter = ',
     &        fractionLongitudinal2ndQuarter*100.0, ' +/- ',
     &        fractionLongitudinal2ndQuarter_sigma*100.0, ' % '
         print *, '        3rd quarter = ',
     &        fractionLongitudinal3rdQuarter*100.0, ' +/- ',
     &        fractionLongitudinal3rdQuarter_sigma*100.0, ' % '
         print *, '        4th quarter = ',
     &        fractionLongitudinal4thQuarter*100.0, ' +/- ',
     &        fractionLongitudinal4thQuarter_sigma*100.0, ' % '
         
         fractionTransverse1stThird = 0.0
         fractionTransverse2ndThird = 0.0
         fractionTransverse3rdThird = 0.0
         fractionTransverse1stThird_sigma = 0.0
         fractionTransverse2ndThird_sigma = 0.0
         fractionTransverse3rdThird_sigma = 0.0
         
         print *, ' Average <E> [MeV] in each Radius bin '
         DO I = 1, numberOfRadiusBins

            sum = vecR( I ) 
            sum2 = vecR2( I )

            mu = sum / (1.0*numEvent)
            sigma = sqrt( abs( sum2 - sum*sum/(1.0*numEvent) ) 
     &           / ( 1.0*numEvent - 1.0 ) )
            mu_sigma = sigma / sqrt( 1.0*numEvent )
            print *, '   ring = ', I, '  E =', 1000.0*mu, ' +/- ', 
     &           1000.0*mu_sigma
            
            IF ( I.LE.(numberOfRadiusBins/3) ) THEN
               fractionTransverse1stThird = 
     &              fractionTransverse1stThird + mu
               fractionTransverse1stThird_sigma = 
     &              fractionTransverse1stThird_sigma + 
     &              mu_sigma * mu_sigma
            ELSE IF ( I.LE.(2*numberOfRadiusBins/3) ) THEN
               fractionTransverse2ndThird = 
     &              fractionTransverse2ndThird + mu
               fractionTransverse2ndThird_sigma = 
     &              fractionTransverse2ndThird_sigma + 
     &              mu_sigma * mu_sigma
            ELSE
               fractionTransverse3rdThird = 
     &              fractionTransverse3rdThird + mu
               fractionTransverse3rdThird_sigma = 
     &              fractionTransverse3rdThird_sigma + 
     &              mu_sigma * mu_sigma
            ENDIF

         ENDDO

         IF ( mu_Evis.GT.1.0E-09 ) THEN

            print *, '  sumR_1 = ', 1000.0*fractionTransverse1stThird, 
     &           ' +/- ',
     &           1000.0*sqrt( fractionTransverse1stThird_sigma ),
     &           ' MeV'
            print *, '  sumR_2 = ', 1000.0*fractionTransverse2ndThird, 
     &           ' +/- ',
     &           1000.0*sqrt( fractionTransverse2ndThird_sigma ),
     &           ' MeV'
            print *, '  sumR_3 = ', 1000.0*fractionTransverse3rdThird, 
     &           ' +/- ',
     &           1000.0*sqrt( fractionTransverse3rdThird_sigma ),
     &           ' MeV'

            IF ( fractionTransverse1stThird.GT.1.0E-09 ) THEN
               fractionTransverse1stThird_sigma =
     &              fractionTransverse1stThird_sigma /
     &              ( fractionTransverse1stThird * 
     &              fractionTransverse1stThird )
            ENDIF
            fractionTransverse1stThird = 
     &           fractionTransverse1stThird / mu_Evis
            fractionTransverse1stThird_sigma = 
     &           fractionTransverse1stThird *
     &           sqrt( fractionTransverse1stThird_sigma +  
     &                 ( mu_Evis_sigma / mu_Evis ) * 
     &                 ( mu_Evis_sigma / mu_Evis ) )
            
            IF ( fractionTransverse2ndThird.GT.1.0E-09 ) THEN
               fractionTransverse2ndThird_sigma =
     &              fractionTransverse2ndThird_sigma /
     &              ( fractionTransverse2ndThird * 
     &              fractionTransverse2ndThird )
            ENDIF
            fractionTransverse2ndThird = 
     &           fractionTransverse2ndThird / mu_Evis
            fractionTransverse2ndThird_sigma = 
     &           fractionTransverse2ndThird *
     &           sqrt( fractionTransverse2ndThird_sigma +  
     &                 ( mu_Evis_sigma / mu_Evis ) * 
     &                 ( mu_Evis_sigma / mu_Evis ) )
            
            IF ( fractionTransverse3rdThird.GT.1.0E-09 ) THEN
               fractionTransverse3rdThird_sigma =
     &              fractionTransverse3rdThird_sigma /
     &              ( fractionTransverse3rdThird * 
     &              fractionTransverse3rdThird )
            ENDIF
            fractionTransverse3rdThird = 
     &           fractionTransverse3rdThird / mu_Evis
            fractionTransverse3rdThird_sigma = 
     &           fractionTransverse3rdThird *
     &           sqrt( fractionTransverse3rdThird_sigma +  
     &                 ( mu_Evis_sigma / mu_Evis ) * 
     &                 ( mu_Evis_sigma / mu_Evis ) )
            
         ENDIF

         print *, ' frac R 1st third = ',
     &        fractionTransverse1stThird*100.0, ' +/- ',
     &        fractionTransverse1stThird_sigma*100.0, ' % '
         print *, '        2nd third = ',
     &        fractionTransverse2ndThird*100.0, ' +/- ',
     &        fractionTransverse2ndThird_sigma*100.0, ' % '
         print *, '        3rd third = ',
     &        fractionTransverse3rdThird*100.0, ' +/- ',
     &        fractionTransverse3rdThird_sigma*100.0, ' % '
         
         print *, '  '

      ENDDO ! End loop over iCase

      print *, ' ----------------------------------- '
      print *, '  '

C---- Print information on the energy resolution.
C     Because non-gaussian distributions of the visible energy 
C     are expected for non-compensating calorimeters, it is not 
C     possible to use the sigma of the gaussian fit as estimator 
C     of the width of the distribution. On the other hand, the  rms  
C     is not appropriate as estimator of the energy resolution, 
C     because it weighs too much the tails of the distribution 
C     (because of the square of the deviation of a value from the 
C      average). A more appropriate estimator is the 
C     "normalized deviation", i.e. the average of the absolute value
C     of the deviation from the mean (so each value has the same weight,
C     regardless whether it is in the central region or in the tail),
C     normalized in such a way to coincide with  sigma  in the case of
C     a perfect gaussian distribution.

      width_Evis = 0.0
      numEntriesForEnergyResolution = numEvent
      IF ( numEvent.GT.maxNumberOfEvents ) THEN
        numEntriesForEnergyResolution = maxNumberOfEvents 
      ENDIF
      DO I = 1, numEntriesForEnergyResolution   
         IF ( vecEvis(I).GT.1.0E-09 ) THEN
            width_Evis = width_Evis + 
     &           abs( vecEvis( I ) - totalVisibleEnergy )
         ENDIF
      ENDDO

      width_Evis = width_Evis * 
     &     sqrt( 3.141592654/2.0 ) / 
     &     (1.0*numEntriesForEnergyResolution)  

      width_Evis_sigma = width_Evis / 
     &     sqrt( 2.0*(numEntriesForEnergyResolution - 1) )

      energyResolution = 0.0
      energyResolution_sigma = 0.0
      IF ( totalVisibleEnergy.GT.1.0E-09 ) THEN
         energyResolution = width_Evis / totalVisibleEnergy
         IF ( width_Evis.GT.1.0E-09 ) THEN
            energyResolution_sigma = energyResolution *
     &           sqrt( ( width_Evis_sigma * width_Evis_sigma ) / 
     &                 ( width_Evis * width_Evis ) 
     &                +
     &                 ( sigma_totalVisibleEnergy * 
     &                   sigma_totalVisibleEnergy ) / 
     &                 ( totalVisibleEnergy + totalVisibleEnergy ) )
         ENDIF
      ENDIF
 
      print *, ' Visible energy information [MeV] '
      print *, '  mu_Evis    = ', 
     &     1000.0*totalVisibleEnergy, ' +/- ', 
     &     1000.0*sigma_totalVisibleEnergy
      print *, '  sigma_Evis = ', 1000.0*width_Evis, ' +/- ', 
     &      1000.0*width_Evis_sigma
      print *, '  energy resolution = ', energyResolution, ' +/- ',  
     &      energyResolution_sigma

      print *, '  '
      print *, ' Average number of particles produced per event: '
      DO I = 1, 100
         IF ( vecParticles( I ).GT.0.0 ) THEN
            flukaParticleCode = I
            IF ( I.GE.90 ) THEN
               flukaParticleCode = I - 100
            ENDIF
            flukaParticleName = '< -6'
            IF ( flukaParticleCode.GE.-6 ) THEN
               flukaParticleName = PRNAME( flukaParticleCode )
            ENDIF
            print *, '  code=', flukaParticleCode, 
     &               '  ', flukaParticleName,
     &               '  # ', vecParticles( I )/(1.0*numEvent)
         ENDIF
      ENDDO
      print *, '  '

      print *, ' Average number of particles Exiting,',
     &     ' and their Ekin, per event: '
      totNumExitingParticles = 0.0
      totEkinExitingParticles = 0.0
      DO I = 1, 100
         totNumExitingParticles = totNumExitingParticles +
     &        vecCountExitingParticles( I )  
         IF ( vecEkinExitingParticles( I ).GT.0.0 ) THEN
            totEkinExitingParticles = totEkinExitingParticles +
     &           vecEkinExitingParticles( I )   
         ENDIF
      ENDDO

      print *, ' # tot = ', (1.0*totNumExitingParticles) / 
     &     (1.0*numEvent)
      print *, ' Ekin tot = ', totEkinExitingParticles / 
     &     (1.0*numEvent), ' GeV'
      DO I = 1, 100
         IF ( vecCountExitingParticles( I ).GT.0.0 ) THEN
            flukaParticleCode = I
            IF ( I.GE.90 ) THEN
               flukaParticleCode = I - 100
            ENDIF
            flukaParticleName = '< -6'
            IF ( flukaParticleCode.GE.-6 ) THEN
               flukaParticleName = PRNAME( flukaParticleCode )
            ENDIF
            averageEkin = -999.9
            IF ( vecEkinExitingParticles( I ).GT.0.0 ) THEN
               averageEkin = vecEkinExitingParticles( I ) / 
     &              (1.0*numEvent)
            ENDIF
            print *, '  code=', flukaParticleCode, 
     &               ' ', flukaParticleName,
     &               ' # ', vecCountExitingParticles( I )
     &           /(1.0*numEvent), 
     &           ' (', INT( 100.0*vecCountExitingParticles( I ) / 
     &           totNumExitingParticles ), ' %)', 
     &           ' Ekin=', averageEkin,
     &           ' (', INT( 100.0*averageEkin * numEvent / 
     &           totEkinExitingParticles ), ' %)'
         ENDIF
      ENDDO
      print *, '  '

      RETURN
*=== End of subroutine Usrout =========================================*
      END

