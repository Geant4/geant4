C-----------------------------------------------------------------
C Keep in a include file all the varibles that are in my COMMONs.
C-----------------------------------------------------------------

C---- The common block MY_COM_1 is made of the variables which are 
C     needed in order to evaluate, at the end of the Run, the
C     average energy deposit per event in the calorimeter, and
C     its statistical error.

      INTEGER numEvent
      REAL sumEdep_all, sumEdep_vis
      REAL sumEdep_all2, sumEdep_vis2
      REAL edepAllCalo, edepVisCalo
      COMMON/MY_COM_1/numEvent, sumEdep_all, sumEdep_vis,
     &     sumEdep_all2, sumEdep_vis2, edepAllCalo, edepVisCalo


C---- The common block MY_COM_2 is made of the variables which are
C     needed for the longitudinal & transverse shower profile, and
C     for the energy resolution.    

      INCLUDE 'mysmallinclude.inc'

      REAL vecEvis( maxNumberOfEvents ) 
      REAL longitudinalProfile( numberOfReadoutLayers )  
      REAL sumL( numberOfReadoutLayers )
      REAL sumL2( numberOfReadoutLayers )
      REAL transverseProfile( numberOfRadiusBins )  
      REAL sumR( numberOfRadiusBins )
      REAL sumR2( numberOfRadiusBins )
      COMMON/MY_COM_2/vecEvis,
     &     longitudinalProfile, sumL, sumL2,
     &     transverseProfile, sumR, sumR2


C---- The common block MY_COM_3 is made of the variables which are
C     needed for counting the number of particles produced.    
C     We use the Fluka track ID as index of the vector vecTrkNum.
C     The index of this vector is indVecTracks, which is reset
C     to 0 at the beginning of each event.
C     We use the Fluka particle code as index of the vector vecParticles.
C     There are 64 particles with positive index, from 1 to 64
C     (including some reserved codes), but there are some with non 
C     positive index: for values from 0 to -6, we use the last 7 elements
C     of vecParticles; We use the element  90  to collect all the
C     particles with Fluka codes < -6 .
C     The elements of vecParticles corresponds at the total number
C     of tracks (as a real, instead of an integer, to avoid overflow
C     problems in runs with high beam energies and many events), 
C     created during the whole Run, of type given by
C     the index number (according to the above simple rule).

      INTEGER maxNumTracksPerEvent
      PARAMETER ( maxNumTracksPerEvent = 10000000 )
      INTEGER vecTrkNum( maxNumTracksPerEvent ), indVecTracks
      DOUBLE PRECISION vecParticles(100)
      COMMON/MY_COM_3/vecTrkNum, vecParticles, indVecTracks


C---- The common block MY_COM_4 is made of the variables which are
C     needed for counting the number of particles (as a real, 
C     instead of an integer, to avoid overflow problems in runs 
C     with high beam energies and many events) exiting the
C     calorimeter, and their kinetic energy.

      REAL vecCountExitingParticles(100)
      REAL vecEkinExitingParticles(100)
      COMMON/MY_COM_4/vecCountExitingParticles, 
     &	   vecEkinExitingParticles


C---- The common block MY_COM_5 is made of the variables which are
C     needed for the longitudinal & transverse shower profile for
C     different types of particles.    

      REAL longitudinalProfile_electron( numberOfReadoutLayers )  
      REAL sumL_electron( numberOfReadoutLayers )
      REAL sumL_electron2( numberOfReadoutLayers )
      REAL transverseProfile_electron( numberOfRadiusBins )  
      REAL sumR_electron( numberOfRadiusBins )
      REAL sumR_electron2( numberOfRadiusBins )
      REAL sumEdep_all_electron, sumEdep_all_electron2
      REAL sumEdep_vis_electron, sumEdep_vis_electron2
      REAL edepAllCalo_electron, edepVisCalo_electron
C
      REAL longitudinalProfile_muon( numberOfReadoutLayers )  
      REAL sumL_muon( numberOfReadoutLayers )
      REAL sumL_muon2( numberOfReadoutLayers )
      REAL transverseProfile_muon( numberOfRadiusBins )  
      REAL sumR_muon( numberOfRadiusBins )
      REAL sumR_muon2( numberOfRadiusBins )
      REAL sumEdep_all_muon, sumEdep_all_muon2
      REAL sumEdep_vis_muon, sumEdep_vis_muon2
      REAL edepAllCalo_muon, edepVisCalo_muon
C
      REAL longitudinalProfile_kaon( numberOfReadoutLayers )  
      REAL sumL_kaon( numberOfReadoutLayers )
      REAL sumL_kaon2( numberOfReadoutLayers )
      REAL transverseProfile_kaon( numberOfRadiusBins )  
      REAL sumR_kaon( numberOfRadiusBins )
      REAL sumR_kaon2( numberOfRadiusBins )
      REAL sumEdep_all_kaon, sumEdep_all_kaon2
      REAL sumEdep_vis_kaon, sumEdep_vis_kaon2
      REAL edepAllCalo_kaon, edepVisCalo_kaon
C
      REAL longitudinalProfile_proton( numberOfReadoutLayers )  
      REAL sumL_proton( numberOfReadoutLayers )
      REAL sumL_proton2( numberOfReadoutLayers )
      REAL transverseProfile_proton( numberOfRadiusBins )  
      REAL sumR_proton( numberOfRadiusBins )
      REAL sumR_proton2( numberOfRadiusBins )
      REAL sumEdep_all_proton, sumEdep_all_proton2
      REAL sumEdep_vis_proton, sumEdep_vis_proton2
      REAL edepAllCalo_proton, edepVisCalo_proton
C
      REAL longitudinalProfile_pion( numberOfReadoutLayers )  
      REAL sumL_pion( numberOfReadoutLayers )
      REAL sumL_pion2( numberOfReadoutLayers )
      REAL transverseProfile_pion( numberOfRadiusBins )  
      REAL sumR_pion( numberOfRadiusBins )
      REAL sumR_pion2( numberOfRadiusBins )
      REAL sumEdep_all_pion, sumEdep_all_pion2
      REAL sumEdep_vis_pion, sumEdep_vis_pion2
      REAL edepAllCalo_pion, edepVisCalo_pion
C
      REAL longitudinalProfile_ion( numberOfReadoutLayers )  
      REAL sumL_ion( numberOfReadoutLayers )
      REAL sumL_ion2( numberOfReadoutLayers )
      REAL transverseProfile_ion( numberOfRadiusBins )  
      REAL sumR_ion( numberOfRadiusBins )
      REAL sumR_ion2( numberOfRadiusBins )
      REAL sumEdep_all_ion, sumEdep_all_ion2
      REAL sumEdep_vis_ion, sumEdep_vis_ion2
      REAL edepAllCalo_ion, edepVisCalo_ion
C
      COMMON/MY_COM_5/
     &     longitudinalProfile_electron, 
     &	   sumL_electron, sumL_electron2,
     &     longitudinalProfile_muon, 
     &	   sumL_muon, sumL_muon2,
     &     longitudinalProfile_kaon, 
     &	   sumL_kaon, sumL_kaon2,
     &     longitudinalProfile_proton, 
     &	   sumL_proton, sumL_proton2,
     &     longitudinalProfile_pion, 
     &	   sumL_pion, sumL_pion2,
     &     longitudinalProfile_ion, 
     &	   sumL_ion, sumL_ion2,
C
     &     transverseProfile_electron, 
     &     sumR_electron, sumR_electron2,
     &     transverseProfile_muon, 
     &     sumR_muon, sumR_muon2,
     &     transverseProfile_kaon, 
     &     sumR_kaon, sumR_kaon2,
     &     transverseProfile_proton, 
     &     sumR_proton, sumR_proton2,
     &     transverseProfile_pion, 
     &     sumR_pion, sumR_pion2,
     &     transverseProfile_ion, 
     &     sumR_ion, sumR_ion2,
C
     & 	   sumEdep_all_electron, sumEdep_all_electron2,
     &     sumEdep_vis_electron, sumEdep_vis_electron2,
     & 	   edepAllCalo_electron, edepVisCalo_electron,
     & 	   sumEdep_all_muon, sumEdep_all_muon2,
     &     sumEdep_vis_muon, sumEdep_vis_muon2,
     & 	   edepAllCalo_muon, edepVisCalo_muon,
     & 	   sumEdep_all_kaon, sumEdep_all_kaon2,
     &     sumEdep_vis_kaon, sumEdep_vis_kaon2,
     & 	   edepAllCalo_kaon, edepVisCalo_kaon,
     & 	   sumEdep_all_proton, sumEdep_all_proton2,
     &     sumEdep_vis_proton, sumEdep_vis_proton2,
     & 	   edepAllCalo_proton, edepVisCalo_proton,
     & 	   sumEdep_all_pion, sumEdep_all_pion2,
     &     sumEdep_vis_pion, sumEdep_vis_pion2,
     & 	   edepAllCalo_pion, edepVisCalo_pion,
     & 	   sumEdep_all_ion, sumEdep_all_ion2,
     &     sumEdep_vis_ion, sumEdep_vis_ion2,
     & 	   edepAllCalo_ion, edepVisCalo_ion


