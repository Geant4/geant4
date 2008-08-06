*$ CREATE USRINI.FOR
*COPY USRINI
*
*=== usrini ===========================================================*
*
      SUBROUTINE USRINI ( WHAT, SDUM )

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
*     USeR INItialization: this routine is called every time the       *
*                          USRICALL card is found in the input stream  *
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

*  Don't change the following line:
      LUSRIN = .TRUE.
* *** Write from here on *** *
*
      PRINT *, 'USRINI---' 

C---- Initialize the variables.
      numEvent = 0
      sumEdep_all = 0.0
      sumEdep_vis = 0.0
      sumEdep_all2 = 0.0
      sumEdep_vis2 = 0.0

      sumEdep_all_electron = 0.0
      sumEdep_vis_electron = 0.0
      sumEdep_all_electron2 = 0.0
      sumEdep_vis_electron2 = 0.0

      sumEdep_all_muon = 0.0
      sumEdep_vis_muon = 0.0
      sumEdep_all_muon2 = 0.0
      sumEdep_vis_muon2 = 0.0

      sumEdep_all_kaon = 0.0
      sumEdep_vis_kaon = 0.0
      sumEdep_all_kaon2 = 0.0
      sumEdep_vis_kaon2 = 0.0

      sumEdep_all_proton = 0.0
      sumEdep_vis_proton = 0.0
      sumEdep_all_proton2 = 0.0
      sumEdep_vis_proton2 = 0.0

      sumEdep_all_pion = 0.0
      sumEdep_vis_pion = 0.0
      sumEdep_all_pion2 = 0.0
      sumEdep_vis_pion2 = 0.0

      sumEdep_all_ion = 0.0
      sumEdep_vis_ion = 0.0
      sumEdep_all_ion2 = 0.0
      sumEdep_vis_ion2 = 0.0

      DO I = 1, numberOfReadoutLayers
         sumL( I )  = 0.0
         sumL2( I ) = 0.0

         sumL_electron( I )  = 0.0
         sumL_electron2( I ) = 0.0

         sumL_muon( I )  = 0.0
         sumL_muon2( I ) = 0.0

         sumL_kaon( I )  = 0.0
         sumL_kaon2( I ) = 0.0

         sumL_proton( I )  = 0.0
         sumL_proton2( I ) = 0.0

         sumL_pion( I )  = 0.0
         sumL_pion2( I ) = 0.0

         sumL_ion( I )  = 0.0
         sumL_ion2( I ) = 0.0
      ENDDO

      DO I = 1, numberOfRadiusBins 
         sumR( I )  = 0.0
         sumR2( I ) = 0.0

         sumR_electron( I )  = 0.0
         sumR_electron2( I ) = 0.0

         sumR_muon( I )  = 0.0
         sumR_muon2( I ) = 0.0

         sumR_kaon( I )  = 0.0
         sumR_kaon2( I ) = 0.0

         sumR_proton( I )  = 0.0
         sumR_proton2( I ) = 0.0

         sumR_pion( I )  = 0.0
         sumR_pion2( I ) = 0.0

         sumR_ion( I )  = 0.0
         sumR_ion2( I ) = 0.0
      ENDDO
   
      DO I = 1, maxNumberOfEvents
         vecEvis( I ) = 0.0
      ENDDO

      DO I = 1, 100
         vecParticles( I ) = 0.0
      ENDDO

      DO I = 1, 100
         vecCountExitingParticles( I ) = 0.0
         vecEkinExitingParticles( I ) = 0.0
      ENDDO

      RETURN
*=== End of subroutine Usrini =========================================*
      END

