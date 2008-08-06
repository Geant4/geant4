*$ CREATE USREIN.FOR
*COPY USREIN
*
*=== Usrein ===========================================================*
*
      SUBROUTINE USREIN

      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'
*
*----------------------------------------------------------------------*
*                                                                      *
*     Copyright (C) 1991-2005      by    Alfredo Ferrari & Paola Sala  *
*     All Rights Reserved.                                             *
*                                                                      *
*                                                                      *
*     USeR Event INitialization: this routine is called before the     *
*     showering of an event is started, but after the source particles *
*     of that event have been already loaded on the stack              *
*                                                                      *
*     Created on 01 january 1991   by    Alfredo Ferrari & Paola Sala  *
*                                                   Infn - Milan       *
*                                                                      *
*     Last change on 09-apr-99     by    Alfredo Ferrari               *
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*
      INCLUDE '(CASLIM)'

      INCLUDE 'myinclude.inc'

 
      PRINT *, ' USREIN-- event=', NCASE

C---- Reset the variable that keeps the total deposit energy
C     per event, also for different particle types, and increase 
C     the counter of the events.
      edepAllCalo = 0.0
      edepVisCalo = 0.0
      numEvent = numEvent + 1

      edepAllCalo_electron = 0.0
      edepVisCalo_electron = 0.0

      edepAllCalo_muon = 0.0
      edepVisCalo_muon = 0.0

      edepAllCalo_kaon = 0.0
      edepVisCalo_kaon = 0.0

      edepAllCalo_proton = 0.0
      edepVisCalo_proton = 0.0

      edepAllCalo_pion = 0.0
      edepVisCalo_pion = 0.0

      edepAllCalo_ion = 0.0
      edepVisCalo_ion = 0.0

C---- Reset the vectors for the event information on the visible
C     energy in each layer (longitudinal shower profile) and for 
C     each ring (transverse shower profile).
      DO I = 1, numberOfReadoutLayers
         longitudinalProfile( I ) = 0.0
         longitudinalProfile_electron( I ) = 0.0
         longitudinalProfile_muon( I ) = 0.0
         longitudinalProfile_kaon( I ) = 0.0
         longitudinalProfile_proton( I ) = 0.0
         longitudinalProfile_pion( I ) = 0.0
         longitudinalProfile_ion( I ) = 0.0
      ENDDO

      DO I = 1, numberOfRadiusBins 
         transverseProfile( I ) = 0.0
         transverseProfile_electron( I ) = 0.0
         transverseProfile_muon( I ) = 0.0
         transverseProfile_kaon( I ) = 0.0
         transverseProfile_proton( I ) = 0.0
         transverseProfile_pion( I ) = 0.0
         transverseProfile_ion( I ) = 0.0
      ENDDO

C---- Reset the index of tracks at the beginning of each event.
      indVecTracks = 0      

      RETURN
*=== End of subroutine Usrein =========================================*
      END

