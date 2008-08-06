*$ CREATE USREOU.FOR
*COPY USREOU
*
*=== Usreou ===========================================================*
*
      SUBROUTINE USREOU

      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'
*
      INCLUDE 'myinclude.inc'
*
*----------------------------------------------------------------------*
*                                                                      *
*     Copyright (C) 1991-2005      by    Alfredo Ferrari & Paola Sala  *
*     All Rights Reserved.                                             *
*                                                                      *
*                                                                      *
*     USeR Event OUtput: this routine is called at the end of each     *
*     event                                                            *
*                                                                      *
*     Created on 01 january 1991   by    Alfredo Ferrari & Paola Sala  *
*                                                   Infn - Milan       *
*                                                                      *
*     Last change on 09-apr-99     by    Alfredo Ferrari               *
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*
CCC      PRINT *, ' USREOU--  edepAllCalo=', edepAllCalo,
CCC     &     '  edepVisCalo=', edepVisCalo

C---- At the end of the event, add the total energy deposition of
C     the current event to the sum of all the energy depositions
C     for all the events so far. Update also the square of it
C     (it will be useful at the end of the Run to calculate the
C      statistical error of the average total energy deposition
C      per event).
      sumEdep_all  = sumEdep_all  + edepAllCalo
      sumEdep_all2 = sumEdep_all2 + edepAllCalo*edepAllCalo
      sumEdep_vis  = sumEdep_vis  + edepVisCalo
      sumEdep_vis2 = sumEdep_vis2 + edepVisCalo*edepVisCalo

      sumEdep_all_electron  = sumEdep_all_electron  + 
     &     edepAllCalo_electron
      sumEdep_all_electron2 = sumEdep_all_electron2 + 
     &     edepAllCalo_electron*edepAllCalo_electron
      sumEdep_vis_electron  = sumEdep_vis_electron  + 
     &     edepVisCalo_electron
      sumEdep_vis_electron2 = sumEdep_vis_electron2 + 
     &     edepVisCalo_electron*edepVisCalo_electron
CCC      print *, ' USREOU-- edepAllCalo_electron=', 
CCC     &     edepAllCalo_electron, 
CCC     &     ' edepVisCalo_electron=', edepVisCalo_electron

      sumEdep_all_muon  = sumEdep_all_muon  + 
     &     edepAllCalo_muon
      sumEdep_all_muon2 = sumEdep_all_muon2 + 
     &     edepAllCalo_muon*edepAllCalo_muon
      sumEdep_vis_muon  = sumEdep_vis_muon  + 
     &     edepVisCalo_muon
      sumEdep_vis_muon2 = sumEdep_vis_muon2 + 
     &     edepVisCalo_muon*edepVisCalo_muon
CCC      print *, ' USREOU-- edepAllCalo_muon=', 
CCC     &     edepAllCalo_muon, 
CCC     &     ' edepVisCalo_muon=', edepVisCalo_muon

      sumEdep_all_kaon  = sumEdep_all_kaon  + 
     &     edepAllCalo_kaon
      sumEdep_all_kaon2 = sumEdep_all_kaon2 + 
     &     edepAllCalo_kaon*edepAllCalo_kaon
      sumEdep_vis_kaon  = sumEdep_vis_kaon  + 
     &     edepVisCalo_kaon
      sumEdep_vis_kaon2 = sumEdep_vis_kaon2 + 
     &     edepVisCalo_kaon*edepVisCalo_kaon
CCC      print *, ' USREOU-- edepAllCalo_kaon=', 
CCC     &     edepAllCalo_kaon, 
CCC     &     ' edepVisCalo_kaon=', edepVisCalo_kaon

      sumEdep_all_proton  = sumEdep_all_proton  + 
     &     edepAllCalo_proton
      sumEdep_all_proton2 = sumEdep_all_proton2 + 
     &     edepAllCalo_proton*edepAllCalo_proton
      sumEdep_vis_proton  = sumEdep_vis_proton  + 
     &     edepVisCalo_proton
      sumEdep_vis_proton2 = sumEdep_vis_proton2 + 
     &     edepVisCalo_proton*edepVisCalo_proton
CCC      print *, ' USREOU-- edepAllCalo_proton=', 
CCC     &     edepAllCalo_proton, 
CCC     &     ' edepVisCalo_proton=', edepVisCalo_proton

      sumEdep_all_pion  = sumEdep_all_pion  + 
     &     edepAllCalo_pion
      sumEdep_all_pion2 = sumEdep_all_pion2 + 
     &     edepAllCalo_pion*edepAllCalo_pion
      sumEdep_vis_pion  = sumEdep_vis_pion  + 
     &     edepVisCalo_pion
      sumEdep_vis_pion2 = sumEdep_vis_pion2 + 
     &     edepVisCalo_pion*edepVisCalo_pion
CCC      print *, ' USREOU-- edepAllCalo_pion=', 
CCC     &     edepAllCalo_pion, 
CCC     &     ' edepVisCalo_pion=', edepVisCalo_pion

      sumEdep_all_ion  = sumEdep_all_ion  + 
     &     edepAllCalo_ion
      sumEdep_all_ion2 = sumEdep_all_ion2 + 
     &     edepAllCalo_ion*edepAllCalo_ion
      sumEdep_vis_ion  = sumEdep_vis_ion  + 
     &     edepVisCalo_ion
      sumEdep_vis_ion2 = sumEdep_vis_ion2 + 
     &     edepVisCalo_ion*edepVisCalo_ion
CCC      print *, ' USREOU-- edepAllCalo_ion=', 
CCC     &     edepAllCalo_ion, 
CCC     &     ' edepVisCalo_ion=', edepVisCalo_ion

C---- Update the information for the longitudinal and transverse
C     shower profiles.
      DO I = 1, numberOfReadoutLayers
         sumL( I )  = sumL( I )  + longitudinalProfile( I )
         sumL2( I ) = sumL2( I ) +
     &        longitudinalProfile( I ) * longitudinalProfile( I )
CCC         print *, ' USREOUT: Longitudinal profile: layer=', I,
CCC     &        ' energy=', longitudinalProfile( I )

         sumL_electron( I )  = sumL_electron( I )  + 
     &        longitudinalProfile_electron( I )
         sumL_electron2( I ) = sumL_electron2( I ) +
     &        longitudinalProfile_electron( I ) * 
     &        longitudinalProfile_electron( I )

         sumL_muon( I )  = sumL_muon( I )  + 
     &        longitudinalProfile_muon( I )
         sumL_muon2( I ) = sumL_muon2( I ) +
     &        longitudinalProfile_muon( I ) * 
     &        longitudinalProfile_muon( I )

         sumL_kaon( I )  = sumL_kaon( I )  + 
     &        longitudinalProfile_kaon( I )
         sumL_kaon2( I ) = sumL_kaon2( I ) +
     &        longitudinalProfile_kaon( I ) * 
     &        longitudinalProfile_kaon( I )

         sumL_proton( I )  = sumL_proton( I )  + 
     &        longitudinalProfile_proton( I )
         sumL_proton2( I ) = sumL_proton2( I ) +
     &        longitudinalProfile_proton( I ) * 
     &        longitudinalProfile_proton( I )

         sumL_pion( I )  = sumL_pion( I )  + 
     &        longitudinalProfile_pion( I )
         sumL_pion2( I ) = sumL_pion2( I ) +
     &        longitudinalProfile_pion( I ) * 
     &        longitudinalProfile_pion( I )

         sumL_ion( I )  = sumL_ion( I )  + 
     &        longitudinalProfile_ion( I )
         sumL_ion2( I ) = sumL_ion2( I ) +
     &        longitudinalProfile_ion( I ) * 
     &        longitudinalProfile_ion( I )
      ENDDO

      DO I = 1, numberOfRadiusBins 
         sumR( I )  = sumR( I )  + transverseProfile( I )
         sumR2( I ) = sumR2( I ) + 
     &        transverseProfile( I ) * transverseProfile( I )
CCC         print *, ' USREOUT: Transverse profile: ring=', I,
CCC     &        ' energy=', transverseProfile( I )

         sumR_electron( I )  = sumR_electron( I )  + 
     &        transverseProfile_electron( I )
         sumR_electron2( I ) = sumR_electron2( I ) + 
     &        transverseProfile_electron( I ) * 
     &        transverseProfile_electron( I )

         sumR_muon( I )  = sumR_muon( I )  + 
     &        transverseProfile_muon( I )
         sumR_muon2( I ) = sumR_muon2( I ) + 
     &        transverseProfile_muon( I ) * 
     &        transverseProfile_muon( I )

         sumR_kaon( I )  = sumR_kaon( I )  + 
     &        transverseProfile_kaon( I )
         sumR_kaon2( I ) = sumR_kaon2( I ) + 
     &        transverseProfile_kaon( I ) * 
     &        transverseProfile_kaon( I )

         sumR_proton( I )  = sumR_proton( I )  + 
     &        transverseProfile_proton( I )
         sumR_proton2( I ) = sumR_proton2( I ) + 
     &        transverseProfile_proton( I ) * 
     &        transverseProfile_proton( I )

         sumR_pion( I )  = sumR_pion( I )  + 
     &        transverseProfile_pion( I )
         sumR_pion2( I ) = sumR_pion2( I ) + 
     &        transverseProfile_pion( I ) * 
     &        transverseProfile_pion( I )

         sumR_ion( I )  = sumR_ion( I )  + 
     &        transverseProfile_ion( I )
         sumR_ion2( I ) = sumR_ion2( I ) + 
     &        transverseProfile_ion( I ) * 
     &        transverseProfile_ion( I )
      ENDDO

C---- Store information of the visible energy, for later computing
C     of the energy resolution.
      IF ( numEvent.LE.maxNumberOfEvents ) THEN
         vecEvis( numEvent ) = edepVisCalo
      ELSE
         print *, 
     & '***IGNORED event (for energy res.) above maxNumberOfEvents=',
     & maxNumberOfEvents 
      ENDIF
*
      RETURN
*=== End of subroutine Usreou =========================================*
      END

