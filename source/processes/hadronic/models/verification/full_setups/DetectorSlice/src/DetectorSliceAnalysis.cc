#include "DetectorSliceAnalysis.hh"
#include <string>
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"

#ifdef G4ANALYSIS_USE
#include <AIDA/AIDA.h>
#endif


DetectorSliceAnalysis* DetectorSliceAnalysis::instance = 0;


DetectorSliceAnalysis::DetectorSliceAnalysis() : 
  numberOfEvents( 0 )
  // The following integers need to be constant, otherwise you cannot
  // use them in switch statements.
  , electronId( G4Electron::ElectronDefinition()->GetPDGEncoding() )
  , positronId( G4Positron::PositronDefinition()->GetPDGEncoding() )
  , gammaId( G4Gamma::GammaDefinition()->GetPDGEncoding() )
  , muonMinusId( G4MuonMinus::MuonMinusDefinition()->GetPDGEncoding() )
  , muonPlusId( G4MuonPlus::MuonPlusDefinition()->GetPDGEncoding() )
  , tauMinusId( G4TauMinus::TauMinusDefinition()->GetPDGEncoding() )
  , tauPlusId( G4TauPlus::TauPlusDefinition()->GetPDGEncoding() )
  , eNeutrinoId( G4NeutrinoE::NeutrinoEDefinition()->GetPDGEncoding() )
  , antiENeutrinoId( G4AntiNeutrinoE::AntiNeutrinoEDefinition()->GetPDGEncoding() ) 
  , muNeutrinoId( G4NeutrinoMu::NeutrinoMuDefinition()->GetPDGEncoding() ) 
  , antiMuNeutrinoId( G4AntiNeutrinoMu::AntiNeutrinoMuDefinition()->GetPDGEncoding() )
  , tauNeutrinoId( G4NeutrinoTau::NeutrinoTauDefinition()->GetPDGEncoding() )
  , antiTauNeutrinoId( G4AntiNeutrinoTau::AntiNeutrinoTauDefinition()->GetPDGEncoding() )
  , pionPlusId( G4PionPlus::PionPlusDefinition()->GetPDGEncoding() )
  , pionMinusId( G4PionMinus::PionMinusDefinition()->GetPDGEncoding() )
  , pionZeroId( G4PionZero::PionZeroDefinition()->GetPDGEncoding() ) 
  , kaonMinusId( G4KaonMinus::KaonMinusDefinition()->GetPDGEncoding() )
  , kaonPlusId( G4KaonPlus::KaonPlusDefinition()->GetPDGEncoding() )
  , kaonZeroId( G4KaonZero::KaonZeroDefinition()->GetPDGEncoding() )
  , antiKaonZeroId( G4AntiKaonZero::AntiKaonZeroDefinition()->GetPDGEncoding() )
  , kaonShortId( G4KaonZeroShort::KaonZeroShortDefinition()->GetPDGEncoding() )
  , kaonLongId( G4KaonZeroLong::KaonZeroLongDefinition()->GetPDGEncoding() )
  , protonId( G4Proton::ProtonDefinition()->GetPDGEncoding() )
  , antiProtonId( G4AntiProton::AntiProtonDefinition()->GetPDGEncoding() )
  , neutronId( G4Neutron::NeutronDefinition()->GetPDGEncoding() )
  , antiNeutronId( G4AntiNeutron::AntiNeutronDefinition()->GetPDGEncoding() )
  //
#ifdef G4ANALYSIS_USE
  , analysisFactory( 0 ), tree( 0 ), tuple( 0 ), histoFactory( 0 )
#endif

{
#ifdef G4ANALYSIS_USE
  analysisFactory = AIDA_createAnalysisFactory();
  if ( analysisFactory ) {    
    AIDA::ITreeFactory* treeFactory = analysisFactory->createTreeFactory();
    if ( treeFactory ) {
      // Trees in memory: create a "tree" associated to an hbook file,
      // which will be filled with an ntuple, and several other trees
      // associated to an XML file, for data sets.
      bool readOnly = false;  // we want to write.
      bool createNew = true ; // create file if it doesn't exist.
      tree = treeFactory->create( "ntuple.hbook", "hbook", readOnly, createNew );
      if ( tree ) {
	// Get a tuple factory :
	AIDA::ITupleFactory* tupleFactory = analysisFactory->createTupleFactory( *tree );
	if ( tupleFactory ) {
	  // Create a tuple :
	  std::string description = "float ID, E, ED_TRA, ED_EM_S, ED_EM,";
          description += " ED_HAD_S, ED_HAD, ED_MU, R_MU;"; 
	  tuple = tupleFactory->create( "1", "Event info", description );
	  assert( tuple );
	  delete tupleFactory;
	}
        // Create a factory for histograms :
        histoFactory = analysisFactory->createHistogramFactory( *tree );  
      }
      delete treeFactory; // It will not delete the ITree.
    }
  }
#endif
}


DetectorSliceAnalysis::~DetectorSliceAnalysis() {}


void DetectorSliceAnalysis::close() {
#ifdef G4ANALYSIS_USE
  if ( tree ) {
    tree->commit();
    tree->close();
  }
  if ( tree ) {
    delete tree;
    tree = 0;
  }
  if ( analysisFactory ) {
    delete analysisFactory;
    analysisFactory = 0;
  }
  if ( histoFactory ) {
    delete histoFactory;
    histoFactory = 0;
  }
#endif
}


DetectorSliceAnalysis* DetectorSliceAnalysis::getInstance() {
  if ( instance == 0 ) instance = new DetectorSliceAnalysis();
  return instance;
}


void DetectorSliceAnalysis::init() {
  //G4cout << " DetectorSliceAnalysis::init() : Cleaning up..." << G4endl;

  // We need to reset the content of the tuple and the
  // profile containers at the beginning of a new Run,
  // because they can be incompatible between two different
  // Runs. Notice that we never execute jobs with more than
  // one Run to really get data, so in practice it is not
  // a problem to loose the data of the previous runs.
  // In practice, jobs with more runs are used only sometimes
  // to check that there is no crashes and the parameters are
  // properly set, but not when you are really interested in
  // the result of the simulation, i.e. the content of the
  // ntuple in the file ntuple.hbook.

#ifdef G4ANALYSIS_USE  
  if ( tuple ) tuple->reset();
#endif

  primaryParticleId = 0;
  beamEnergy      = 0.0;
  sumEdepTracker  = 0.0;
  sumEdepTracker2 = 0.0;
  sumEdepEmAct    = 0.0;
  sumEdepEmAct2   = 0.0;
  sumEdepEmTot   = 0.0;
  sumEdepEmTot2  = 0.0;
  sumEdepHadAct   = 0.0;
  sumEdepHadAct2  = 0.0;
  sumEdepHadTot   = 0.0;
  sumEdepHadTot2  = 0.0;
  sumEdepMuon     = 0.0;
  sumEdepMuon2    = 0.0;

  numberOfEvents  = 0;

  vecEvisEm.clear();
  vecEvisHad.clear();

  sumEdepEmAct_electron  = 0.0;
  sumEdepEmAct_electron2 = 0.0;
  sumEdepEmAct_muon  = 0.0;
  sumEdepEmAct_muon2 = 0.0;
  sumEdepEmAct_pion  = 0.0;
  sumEdepEmAct_pion2 = 0.0;
  sumEdepEmAct_kaon  = 0.0;
  sumEdepEmAct_kaon2 = 0.0;
  sumEdepEmAct_proton  = 0.0;
  sumEdepEmAct_proton2 = 0.0;
  sumEdepEmAct_nuclei  = 0.0;
  sumEdepEmAct_nuclei2 = 0.0;

  sumEdepHadAct_electron  = 0.0;
  sumEdepHadAct_electron2 = 0.0;
  sumEdepHadAct_muon  = 0.0;
  sumEdepHadAct_muon2 = 0.0;
  sumEdepHadAct_pion  = 0.0;
  sumEdepHadAct_pion2 = 0.0;
  sumEdepHadAct_kaon  = 0.0;
  sumEdepHadAct_kaon2 = 0.0;
  sumEdepHadAct_proton  = 0.0;
  sumEdepHadAct_proton2 = 0.0;
  sumEdepHadAct_nuclei  = 0.0;
  sumEdepHadAct_nuclei2 = 0.0;

  eventTimeSet.clear();

#ifdef G4ANALYSIS_USE
  assert( histoFactory );
  if ( histoFactory ) {
    // Kinetic spectra of some particles when they are created.
    // The x-axis is LOG10(energy/MeV), from inf=-4.0 to sup=+6.0
    // (i.e. from 0.1 keV up to 1000 GeV), with 1000 (logarithmic) 
    // bins. 
    gammaSpectrumInEmCal = histoFactory->
      createHistogram1D( "71", "log10(energy/MeV) for gamma in EM Cal", 
			 1000, -4.0, 6.0 );
    neutronSpectrumInEmCal = histoFactory->
      createHistogram1D( "72", "log10(energy/MeV) for neutron in EM Cal", 
			 1000, -4.0, 6.0 );
    protonSpectrumInEmCal = histoFactory->
      createHistogram1D( "73", "log10(energy/MeV) for proton in EM Cal",
			 1000, -4.0, 6.0 );
    pionZeroSpectrumInEmCal = histoFactory->
      createHistogram1D( "74", "log10(energy/MeV) for pionZero in EM Cal", 
			 1000, -4.0, 6.0 );
    pionPlusSpectrumInEmCal = histoFactory->
      createHistogram1D( "75", "log10(energy/MeV) for pionPlus in EM Cal",
			 1000, -4.0, 6.0 );
    pionMinusSpectrumInEmCal = histoFactory->
      createHistogram1D( "76", "log10(energy/MeV) for pionMinus in EM Cal", 
			 1000, -4.0, 6.0 );

    gammaSpectrumInHadCal = histoFactory->
      createHistogram1D( "81", "log10(energy/MeV) for gamma in HAD Cal", 
			 1000, -4.0, 6.0 );
    neutronSpectrumInHadCal = histoFactory->
      createHistogram1D( "82", "log10(energy/MeV) for neutron in HAD Cal", 
			 1000, -4.0, 6.0 );
    protonSpectrumInHadCal = histoFactory->
      createHistogram1D( "83", "log10(energy/MeV) for proton in HAD Cal",
			 1000, -4.0, 6.0 );
    pionZeroSpectrumInHadCal = histoFactory->
      createHistogram1D( "84", "log10(energy/MeV) for pionZero in HAD Cal", 
			 1000, -4.0, 6.0 );
    pionPlusSpectrumInHadCal = histoFactory->
      createHistogram1D( "85", "log10(energy/MeV) for pionPlus in HAD Cal",
			 1000, -4.0, 6.0 );
    pionMinusSpectrumInHadCal = histoFactory->
      createHistogram1D( "86", "log10(energy/MeV) for pionMinus in HAD Cal", 
			 1000, -4.0, 6.0 );

  }
#endif
}                       


void DetectorSliceAnalysis::fillNtuple( float incidentParticleId, 
					float incidentParticleEnergy, 
					float totalEnergyDepositedInTracker,
					float totalEnergyDepositedInActiveEmCal,
					float totalEnergyDepositedInEmCal, 
					float totalEnergyDepositedInActiveHadCal,
					float totalEnergyDepositedInHadCal, 
					float totalEnergyDepositedInMuonDet,
					float exitingRadiusPrimaryMuon ) {

  primaryParticleId = static_cast< G4int >( incidentParticleId );
  beamEnergy = incidentParticleEnergy;

  //G4cout << " DetectorSliceAnalysis::fillNtuple : DEBUG Info " << G4endl
  //       << "\t incidentParticleId = " << incidentParticleId << G4endl
  //       << "\t incidentParticleEnergy = " << incidentParticleEnergy << G4endl
  //	   << "\t totalEnergyDepositedInTracker = " 
  //	   << totalEnergyDepositedInTracker << G4endl
  //	   << "\t totalEnergyDepositedInActiveEmCal = " 
  //	   << totalEnergyDepositedInActiveEmCal << G4endl 
  //       << "\t totalEnergyDepositedInEmCal = "
  //       << totalEnergyDepositedInEmCal << G4endl
  //	   << "\t totalEnergyDepositedInActiveHadCal = "
  //       << totalEnergyDepositedInActiveHadCal << G4endl
  //	   << "\t totalEnergyDepositedInHadCal = " 
  //       << totalEnergyDepositedInHadCal << G4endl
  //	   << "\t totalEnergyDepositedInMuonDet = "
  //       << totalEnergyDepositedInMuonDet << G4endl
  //	   << "\t exitingRadiusPrimaryMuon = "
  //       << exitingRadiusPrimaryMuon << G4endl;       // ***DEBUG***

#ifdef G4ANALYSIS_USE
  if ( tuple ) {
    tuple->fill( tuple->findColumn( "ID" ), incidentParticleId );
    tuple->fill( tuple->findColumn( "E" ), incidentParticleEnergy );
    tuple->fill( tuple->findColumn( "ED_TRA" ), totalEnergyDepositedInTracker );
    tuple->fill( tuple->findColumn( "ED_EM_S" ), totalEnergyDepositedInActiveEmCal );
    tuple->fill( tuple->findColumn( "ED_EM" ), totalEnergyDepositedInEmCal );
    tuple->fill( tuple->findColumn( "ED_HAD_S" ), totalEnergyDepositedInActiveHadCal );
    tuple->fill( tuple->findColumn( "ED_HAD" ), totalEnergyDepositedInHadCal );
    tuple->fill( tuple->findColumn( "ED_MU" ), totalEnergyDepositedInMuonDet );
    tuple->fill( tuple->findColumn( "R_MU" ), exitingRadiusPrimaryMuon );
    tuple->addRow();
  }
#endif

  sumEdepTracker += totalEnergyDepositedInTracker;
  sumEdepTracker2 += totalEnergyDepositedInTracker * totalEnergyDepositedInTracker;

  sumEdepEmAct  += totalEnergyDepositedInActiveEmCal;
  sumEdepEmAct2 += totalEnergyDepositedInActiveEmCal * totalEnergyDepositedInActiveEmCal;
  sumEdepEmTot  += totalEnergyDepositedInEmCal;
  sumEdepEmTot2 += totalEnergyDepositedInEmCal * totalEnergyDepositedInEmCal;

  sumEdepHadAct  += totalEnergyDepositedInActiveHadCal;
  sumEdepHadAct2 += 
    totalEnergyDepositedInActiveHadCal * totalEnergyDepositedInActiveHadCal;
  sumEdepHadTot  += totalEnergyDepositedInHadCal;
  sumEdepHadTot2 += totalEnergyDepositedInHadCal * totalEnergyDepositedInHadCal;

  sumEdepMuon += totalEnergyDepositedInMuonDet;
  sumEdepMuon2 += totalEnergyDepositedInMuonDet * totalEnergyDepositedInMuonDet;

  numberOfEvents++;

#ifdef G4ANALYSIS_USE
  // This method is called at each event, so it is useful to commit
  // the tree from time to time, for instance every 10 events, in
  // such a way that it is possible to see the ntuple  ntuple.hbook
  // while the job is running, for instance for a quick check that
  // it makes sense. Or, if the job crashes after many events, we
  // have anyhow some data already stored in the ntuple to be checked.
  // Remember that when you look the ntuple with PAW, while is running,
  // you need to close the session and open a new one to see the updates.
  if ( numberOfEvents % 1000 == 0 ) {
    if ( tree ) {
      tree->commit();
      //G4cout << " tree commit ,  at event=" << numberOfEvents-1 
      //      << G4endl; //***DEBUG***
    }
  }
#endif

  // Store information of the visible energy, for later computing
  // of the energy resolution.
  vecEvisEm.push_back( totalEnergyDepositedInActiveEmCal );
  vecEvisHad.push_back( totalEnergyDepositedInActiveHadCal );

}


void DetectorSliceAnalysis::
particleDepositingEnergy( const G4String caloType,
			  const G4double edep, const G4int particlePDG ) {

  // Consider now the separate contribution due to the following particles:
  //    -  electron (e-  and  e+    together)
  //    -  muons    (mu- and  mu+   together)
  //    -  pions    (pi- and  pi+   together)
  //    -  kaons    (k-  and  k+    together)
  //    -  protons  (p   and  pbar  together)
  //    -  nuclei   (all particles with PDG code = 0 and neutrons together)

  if ( particlePDG == electronId  ||  particlePDG == positronId ) {
    if ( caloType == "EM" ) {
      sumEdepEmAct_electron += edep;
    } else if ( caloType == "HAD" ) {
      sumEdepHadAct_electron += edep;
    }
  } else if ( particlePDG == muonMinusId  ||  particlePDG == muonPlusId ) {
    if ( caloType == "EM" ) {
      sumEdepEmAct_muon += edep;
    } else if ( caloType == "HAD" ) {
      sumEdepHadAct_muon += edep;
    }
  } else if ( particlePDG == pionPlusId  ||  particlePDG == pionMinusId ) {
    if ( caloType == "EM" ) {
      sumEdepEmAct_pion += edep;
    } else if ( caloType == "HAD" ) {
      sumEdepHadAct_pion += edep;
    }
  } else if ( particlePDG == kaonMinusId  ||  particlePDG == kaonPlusId ) {
    if ( caloType == "EM" ) {
      sumEdepEmAct_kaon += edep;
    } else if ( caloType == "HAD" ) {
      sumEdepHadAct_kaon += edep;
    }
  } else if ( particlePDG == protonId  ||  particlePDG == antiProtonId ) {
    if ( caloType == "EM" ) {
      sumEdepEmAct_proton += edep;
    } else if ( caloType == "HAD" ) {
      sumEdepHadAct_proton += edep;
    }
  } else if ( particlePDG == 0  || 
	      particlePDG / 1000000000 >= 1  ||
	      particlePDG == neutronId ) {
    // Before 8.1.ref04, the PDG code for nuclei was 0; 
    // from 8.1.ref04, the new PDG code convention for nuclei has been
    // supported: in this schema, the PDG code of nuclei is characterized
    // by 10 digits (all other particles have lower codes).
    // Starting with 8.1, neutrons can also deposit energy, as an 
    // effective way to take into account the recoil of nuclei 
    // below a certain threshold.
    if ( caloType == "EM" ) {
      sumEdepEmAct_nuclei += edep;
    } else if ( caloType == "HAD" ) {
      sumEdepHadAct_nuclei += edep;
    }
  }

}


void DetectorSliceAnalysis::infoTrack( const G4Track* aTrack ) {
#ifdef G4ANALYSIS_USE
  G4int particleId = aTrack->GetDefinition()->GetPDGEncoding();
  G4String volumeName = aTrack->GetVolume()->GetName();
  G4double trackEkin = aTrack->GetKineticEnergy();
  if ( trackEkin < 1.0*eV ) {
    // To avoid problem with the logarithm below.
    trackEkin = 1.0*eV;  
  }
  if ( volumeName.find( "Em" ) != std::string::npos ) {
    if ( particleId == gammaId ) {
      gammaSpectrumInEmCal->fill( std::log10( trackEkin/MeV ) );
    } else if ( particleId == neutronId ) {
      neutronSpectrumInEmCal->fill( std::log10( trackEkin/MeV ) );      
    } else if ( particleId == protonId ) {
      protonSpectrumInEmCal->fill( std::log10( trackEkin/MeV ) );
    } else if ( particleId == pionZeroId ) {
      pionZeroSpectrumInEmCal->fill( std::log10( trackEkin/MeV ) );
    } else if ( particleId == pionPlusId ) {
      pionPlusSpectrumInEmCal->fill( std::log10( trackEkin/MeV ) );
    } else if ( particleId == pionMinusId ) {
      pionMinusSpectrumInEmCal->fill( std::log10( trackEkin/MeV ) );
    }
  } else if ( volumeName.find( "Had" ) != std::string::npos ) {
    if ( particleId == gammaId ) {
      gammaSpectrumInHadCal->fill( std::log10( trackEkin/MeV ) );
    } else if ( particleId == neutronId ) {
      neutronSpectrumInHadCal->fill( std::log10( trackEkin/MeV ) );      
    } else if ( particleId == protonId ) {
      protonSpectrumInHadCal->fill( std::log10( trackEkin/MeV ) );
    } else if ( particleId == pionZeroId ) {
      pionZeroSpectrumInHadCal->fill( std::log10( trackEkin/MeV ) );
    } else if ( particleId == pionPlusId ) {
      pionPlusSpectrumInHadCal->fill( std::log10( trackEkin/MeV ) );
    } else if ( particleId == pionMinusId ) {
      pionMinusSpectrumInHadCal->fill( std::log10( trackEkin/MeV ) );
    }
  }
#endif
}


void DetectorSliceAnalysis::endOfEvent( const G4double timeEventInSec ) {
  // This method is useful to update the "squared" event variables
  // which are used at the end of the Run to compute the statistical
  // uncertainties of various quantities.
  // Notice that only the quantities that are meaningful on an event
  // by event basis, like the number of steps, the number of tracks,
  // the exiting kinetic energy, the number of exiting particles;
  // this does not apply to the track length.

  // Apply the same trick for the energy deposits and shower profiles
  // of the following groups of particles:
  // e-  and  e+  together
   static G4double sumEdepEmAct_electron_previous = 0.0;
   sumEdepEmAct_electron2 += 
     ( sumEdepEmAct_electron - sumEdepEmAct_electron_previous ) *
     ( sumEdepEmAct_electron - sumEdepEmAct_electron_previous );
   sumEdepEmAct_electron_previous = sumEdepEmAct_electron;

   static G4double sumEdepHadAct_electron_previous = 0.0;
   sumEdepHadAct_electron2 += 
     ( sumEdepHadAct_electron - sumEdepHadAct_electron_previous ) *
     ( sumEdepHadAct_electron - sumEdepHadAct_electron_previous );
   sumEdepHadAct_electron_previous = sumEdepHadAct_electron;
 
   // mu-  and  mu+  together
   static G4double sumEdepEmAct_muon_previous = 0.0;
   sumEdepEmAct_muon2 += ( sumEdepEmAct_muon - sumEdepEmAct_muon_previous ) *
                         ( sumEdepEmAct_muon - sumEdepEmAct_muon_previous );
   sumEdepEmAct_muon_previous = sumEdepEmAct_muon;

   static G4double sumEdepHadAct_muon_previous = 0.0;
   sumEdepHadAct_muon2 += ( sumEdepHadAct_muon - sumEdepHadAct_muon_previous ) *
                          ( sumEdepHadAct_muon - sumEdepHadAct_muon_previous );
   sumEdepHadAct_muon_previous = sumEdepHadAct_muon;
 
   // pi+  and  pi-  together
   static G4double sumEdepEmAct_pion_previous = 0.0;
   sumEdepEmAct_pion2 += ( sumEdepEmAct_pion - sumEdepEmAct_pion_previous ) *
                         ( sumEdepEmAct_pion - sumEdepEmAct_pion_previous );
   sumEdepEmAct_pion_previous = sumEdepEmAct_pion;

   static G4double sumEdepHadAct_pion_previous = 0.0;
   sumEdepHadAct_pion2 += ( sumEdepHadAct_pion - sumEdepHadAct_pion_previous ) *
                          ( sumEdepHadAct_pion - sumEdepHadAct_pion_previous );
   sumEdepHadAct_pion_previous = sumEdepHadAct_pion;
 
   // k+  and  k-  together
   static G4double sumEdepEmAct_kaon_previous = 0.0;
   sumEdepEmAct_kaon2 += ( sumEdepEmAct_kaon - sumEdepEmAct_kaon_previous ) *
                         ( sumEdepEmAct_kaon - sumEdepEmAct_kaon_previous );
   sumEdepEmAct_kaon_previous = sumEdepEmAct_kaon;

   static G4double sumEdepHadAct_kaon_previous = 0.0;
   sumEdepHadAct_kaon2 += ( sumEdepHadAct_kaon - sumEdepHadAct_kaon_previous ) *
                          ( sumEdepHadAct_kaon - sumEdepHadAct_kaon_previous );
   sumEdepHadAct_kaon_previous = sumEdepHadAct_kaon;
 
   // p  and  pbar  together
   static G4double sumEdepEmAct_proton_previous = 0.0;
   sumEdepEmAct_proton2 += ( sumEdepEmAct_proton - sumEdepEmAct_proton_previous ) *
                           ( sumEdepEmAct_proton - sumEdepEmAct_proton_previous );
   sumEdepEmAct_proton_previous = sumEdepEmAct_proton;

   static G4double sumEdepHadAct_proton_previous = 0.0;
   sumEdepHadAct_proton2 += ( sumEdepHadAct_proton - sumEdepHadAct_proton_previous ) *
                            ( sumEdepHadAct_proton - sumEdepHadAct_proton_previous );
   sumEdepHadAct_proton_previous = sumEdepHadAct_proton;
 
   // Neutrons and nuclei
   static G4double sumEdepEmAct_nuclei_previous = 0.0;
   sumEdepEmAct_nuclei2 += ( sumEdepEmAct_nuclei - sumEdepEmAct_nuclei_previous ) *
                           ( sumEdepEmAct_nuclei - sumEdepEmAct_nuclei_previous );
   sumEdepEmAct_nuclei_previous = sumEdepEmAct_nuclei;

   static G4double sumEdepHadAct_nuclei_previous = 0.0;
   sumEdepHadAct_nuclei2 += ( sumEdepHadAct_nuclei - sumEdepHadAct_nuclei_previous ) *
                            ( sumEdepHadAct_nuclei - sumEdepHadAct_nuclei_previous );
   sumEdepHadAct_nuclei_previous = sumEdepHadAct_nuclei;
 
   // Keep the time per event in a multiset.
   eventTimeSet.insert( timeEventInSec );

}


void DetectorSliceAnalysis::finish() {

   // Print results. 
   G4cout << " Primary particle PDG Id = " << primaryParticleId << G4endl;
   G4cout << " Beam energy [MeV] = " << beamEnergy / MeV << G4endl << G4endl;
 
   G4double n = static_cast< G4double >( numberOfEvents );
   if ( n <= 1.0 ) {
     n = 2.0;         // To avoid division by zero in  sigma
   }
   //G4cout << " n=" << n << G4endl;                             //***DEBUG***

   G4double sum, sum2, mu, sigma, mu_sigma;

   G4cout << " Average <E> [MeV] deposited in the subdetctors " << G4endl;

   sum  = sumEdepTracker;
   sum2 = sumEdepTracker2;
   mu       = sum / n;
   sigma    = std::sqrt( std::abs( ( sum2 - sum*sum/n ) ) / (n - 1.0) );
   mu_sigma = sigma / std::sqrt( n );
   G4cout << "\t Tracker                = " << mu << " +/- " << mu_sigma << G4endl;
 
   sum  = sumEdepEmAct;
   sum2 = sumEdepEmAct2;
   mu       = sum / n;
   sigma    = std::sqrt( std::abs( ( sum2 - sum*sum/n ) ) / (n - 1.0) );
   mu_sigma = sigma / std::sqrt( n );
   G4cout << "\t EM Calo active layers  = " << mu << " +/- " << mu_sigma << G4endl;
   G4double mu_EvisEm = mu;             // For later usage.
   G4double mu_EvisEm_sigma = mu_sigma; //  "    "     "
 
   sum  = sumEdepEmTot;
   sum2 = sumEdepEmTot2;
   mu       = sum / n;
   sigma    = std::sqrt( std::abs( ( sum2 - sum*sum/n ) ) / (n - 1.0) );
   mu_sigma = sigma / std::sqrt( n );
   G4cout << "\t EM Calo                = " << mu << " +/- " << mu_sigma << G4endl;

   sum  = sumEdepHadAct;
   sum2 = sumEdepHadAct2;
   mu       = sum / n;
   sigma    = std::sqrt( std::abs( ( sum2 - sum*sum/n ) ) / (n - 1.0) );
   mu_sigma = sigma / std::sqrt( n );
   G4cout << "\t HAD Calo active layers = " << mu << " +/- " << mu_sigma << G4endl;
   G4double mu_EvisHad = mu;             // For later usage.
   G4double mu_EvisHad_sigma = mu_sigma; //  "    "     "
 
   sum  = sumEdepHadTot;
   sum2 = sumEdepHadTot2;
   mu       = sum / n;
   sigma    = std::sqrt( std::abs( ( sum2 - sum*sum/n ) ) / (n - 1.0) );
   mu_sigma = sigma / std::sqrt( n );
   G4cout << "\t HAD Calo               = " << mu << " +/- " << mu_sigma << G4endl;

   sum  = sumEdepMuon;
   sum2 = sumEdepMuon2;
   mu       = sum / n;
   sigma    = std::sqrt( std::abs( ( sum2 - sum*sum/n ) ) / (n - 1.0) );
   mu_sigma = sigma / std::sqrt( n );
   G4cout << "\t Muon detector          = " << mu << " +/- " << mu_sigma << G4endl;
  
   // Print information useful for the linearity of the energy response,
   // energy resolutio, e/pi .
   // Because non-gaussian distribution of the visible energy (i.e. the
   // sum of the deposited energies in all active layers) are expected
   // for a non-compensating calorimeter, it is not possible to use the
   //  sigma  of the gaussian fit as estimator of the width of the
   // distribution. On the other hand, the  rms  is not appropriate as
   // estimator of the energy resolution, because it weighs too much the
   // tails of the distribution (because of the square of the deviation
   // of a value from the average). A more appropriate estimator is the
   // "normalized deviation", i.e. the average of the absolute value of
   // the deviation from the mean (so each value has the same weight,
   // regardless whether it is in the central region or in the tail),
   // normalized in such a way to coincide with  sigma  in the case of
   // a perfect gaussian distribution.

   G4double width_EvisEm = 0.0;
   for ( std::vector< G4double >::const_iterator cit = vecEvisEm.begin();
 	cit != vecEvisEm.end() ; ++cit ) {
     width_EvisEm += std::abs( *cit - mu_EvisEm );
   }
   width_EvisEm *= std::sqrt( 3.141592654/2.0 ) / n ;
   G4double width_EvisEm_sigma = width_EvisEm / std::sqrt( 2.0*(n - 1) );
   G4double energyResolutionEm = 0.0;
   G4double energyResolutionEm_sigma = 0.0;
   G4double samplingFractionEm = 0.0;
   G4double samplingFractionEm_sigma = 0.0;
   if ( mu_EvisEm > 1.0E-06) {
     energyResolutionEm = width_EvisEm / mu_EvisEm;
     if ( width_EvisEm > 1.0E-06 ) {
       energyResolutionEm_sigma = energyResolutionEm *
 	std::sqrt( ( width_EvisEm_sigma * width_EvisEm_sigma ) / 
		   ( width_EvisEm * width_EvisEm ) 
 		   +
 		   ( mu_EvisEm_sigma * mu_EvisEm_sigma ) / 
		   ( mu_EvisEm + mu_EvisEm ) );
     }
     samplingFractionEm = mu_EvisEm / beamEnergy;
     samplingFractionEm_sigma = samplingFractionEm * mu_EvisEm_sigma / mu_EvisEm;
   }
 
   G4cout << G4endl
          << " Visible energy information [MeV] for EM Calo " << G4endl
	  << "\t mu_EvisEm     = " << mu_EvisEm 
          << " +/- " << mu_EvisEm_sigma << G4endl
          << "\t sigma_EvisEm = " << width_EvisEm 
          << " +/- " << width_EvisEm_sigma << G4endl
          << "\t energy resolution = " << energyResolutionEm 
          << " +/- "  << energyResolutionEm_sigma << G4endl
          << "\t sampling fraction = " << samplingFractionEm 
          << " +/- " << samplingFractionEm_sigma << G4endl;
 
   G4double width_EvisHad = 0.0;
   for ( std::vector< G4double >::const_iterator cit = vecEvisHad.begin();
 	cit != vecEvisHad.end() ; ++cit ) {
     width_EvisHad += std::abs( *cit - mu_EvisHad );
   }
   width_EvisHad *= std::sqrt( 3.141592654/2.0 ) / n ;
   G4double width_EvisHad_sigma = width_EvisHad / std::sqrt( 2.0*(n - 1) );
   G4double energyResolutionHad = 0.0;
   G4double energyResolutionHad_sigma = 0.0;
   G4double samplingFractionHad = 0.0;
   G4double samplingFractionHad_sigma = 0.0;
   if ( mu_EvisHad > 1.0E-06) {
     energyResolutionHad = width_EvisHad / mu_EvisHad;
     if ( width_EvisHad > 1.0E-06 ) {
       energyResolutionHad_sigma = energyResolutionHad *
 	std::sqrt( ( width_EvisHad_sigma * width_EvisHad_sigma ) / 
		   ( width_EvisHad * width_EvisHad ) 
 		   +
 		   ( mu_EvisHad_sigma * mu_EvisHad_sigma ) / 
		   ( mu_EvisHad + mu_EvisHad ) );
     }
     samplingFractionHad = mu_EvisHad / beamEnergy;
     samplingFractionHad_sigma = samplingFractionHad * mu_EvisHad_sigma / mu_EvisHad;
   }
 
   G4cout << G4endl
          << " Visible energy information [MeV] for HAD Calo " << G4endl
	  << "\t mu_EvisHad     = " << mu_EvisHad 
          << " +/- " << mu_EvisHad_sigma << G4endl
          << "\t sigma_EvisHad = " << width_EvisHad 
          << " +/- " << width_EvisHad_sigma << G4endl
          << "\t energy resolution = " << energyResolutionHad 
          << " +/- "  << energyResolutionHad_sigma << G4endl
          << "\t sampling fraction = " << samplingFractionHad 
          << " +/- " << samplingFractionHad_sigma << G4endl;

 #ifdef G4ANALYSIS_USE
   if ( tree ) tree->commit();
 #endif
 
   // Print information on the different particle contributions to
   // the visible energy in EM and HAD calorimeters.
 
   G4cout << G4endl << " Contributions of the main particle types [MeV] in EM Calo" 
	  << G4endl;
   for ( int iCase = 0; iCase < 6; iCase++ ) {
 
     std::string caseName = "";
     G4double sumVis = 0.0;
     G4double sumVis2 = 0.0;
 
     switch ( iCase ) {
     case 0 : {
       caseName = "electron";
       sumVis = sumEdepEmAct_electron;
       sumVis2 = sumEdepEmAct_electron2;
       break;
     }
     case 1 : {
       caseName = "muon";
       sumVis = sumEdepEmAct_muon;
       sumVis2 = sumEdepEmAct_muon2;
       break;
     }
     case 2 : {
       caseName = "pion";
       sumVis = sumEdepEmAct_pion;
       sumVis2 = sumEdepEmAct_pion2;
       break;
     }
     case 3 : {
       caseName = "kaon";
       sumVis = sumEdepEmAct_kaon;
       sumVis2 = sumEdepEmAct_kaon2;
       break;
     }
     case 4 : {
       caseName = "proton";
       sumVis = sumEdepEmAct_proton;
       sumVis2 = sumEdepEmAct_proton2;
       break;
     }
     case 5 : {
       caseName = "nuclei";
       sumVis = sumEdepEmAct_nuclei;
       sumVis2 = sumEdepEmAct_nuclei2;
       break;
     }
     } // End of switch
 
     G4cout << "\t Particle type: " << caseName << G4endl;
  
     G4double f, f_sigma;
 
     sum  = sumVis;
     sum2 = sumVis2;
     mu       = sum / n;
     sigma    = std::sqrt( std::abs( sum2 - sum*sum/n ) / (n - 1.0) );
     mu_sigma = sigma / std::sqrt( n );
     f = 0.0;
     f_sigma = 0.0;
     if ( mu_EvisEm > 1.0E-06 ) {
       f = mu / mu_EvisEm;
       if ( mu > 1.0E-06 ) {
 	f_sigma = f * std::sqrt( ( (mu_sigma*mu_sigma) / (mu*mu) ) +
 				 ( (mu_EvisEm_sigma*mu_EvisEm_sigma) / 
				   (mu_EvisEm*mu_EvisEm) ) );
       }
     }
     G4cout << "\t \t <E_vis> = " << mu << " +/- " << mu_sigma 
 	   << "  ( " << 100.0*f << " +/- " << 100.0*f_sigma << " % )" << G4endl;
 
   } // End of the loop over the particle types.

   G4cout << G4endl << " Contributions of the main particle types [MeV] in HAD Calo" 
	  << G4endl;
   for ( int iCase = 0; iCase < 6; iCase++ ) {
 
     std::string caseName = "";
     G4double sumVis = 0.0;
     G4double sumVis2 = 0.0;
 
     switch ( iCase ) {
     case 0 : {
       caseName = "electron";
       sumVis = sumEdepHadAct_electron;
       sumVis2 = sumEdepHadAct_electron2;
       break;
     }
     case 1 : {
       caseName = "muon";
       sumVis = sumEdepHadAct_muon;
       sumVis2 = sumEdepHadAct_muon2;
       break;
     }
     case 2 : {
       caseName = "pion";
       sumVis = sumEdepHadAct_pion;
       sumVis2 = sumEdepHadAct_pion2;
       break;
     }
     case 3 : {
       caseName = "kaon";
       sumVis = sumEdepHadAct_kaon;
       sumVis2 = sumEdepHadAct_kaon2;
       break;
     }
     case 4 : {
       caseName = "proton";
       sumVis = sumEdepHadAct_proton;
       sumVis2 = sumEdepHadAct_proton2;
       break;
     }
     case 5 : {
       caseName = "nuclei";
       sumVis = sumEdepHadAct_nuclei;
       sumVis2 = sumEdepHadAct_nuclei2;
       break;
     }
     } // End of switch
 
     G4cout << "\t Particle type: " << caseName << G4endl;
  
     G4double f, f_sigma;
 
     sum  = sumVis;
     sum2 = sumVis2;
     mu       = sum / n;
     sigma    = std::sqrt( std::abs( sum2 - sum*sum/n ) / (n - 1.0) );
     mu_sigma = sigma / std::sqrt( n );
     f = 0.0;
     f_sigma = 0.0;
     if ( mu_EvisHad > 1.0E-06 ) {
       f = mu / mu_EvisHad;
       if ( mu > 1.0E-06 ) {
 	f_sigma = f * std::sqrt( ( (mu_sigma*mu_sigma) / (mu*mu) ) +
 				 ( (mu_EvisHad_sigma*mu_EvisHad_sigma) / 
				   (mu_EvisHad*mu_EvisHad) ) );
       }
     }
     G4cout << "\t \t <E_vis> = " << mu << " +/- " << mu_sigma 
 	   << "  ( " << 100.0*f << " +/- " << 100.0*f_sigma << " % )" << G4endl;
 
   } // End of the loop over the particle types.
 
   // Print information about the CPU time per event.
   // We are interested in the mininum and maximum time, and on the
   // average, its error, and the rms of the distribution. However,
   // to avoid that an eventual long tail in the distribution would
   // bias the mean to high values, we also calculate the reduced
   // mean and rms taking into consideration only those values which
   // are within 5% and 95% of the ordered distribution of values
   // (kept automatically in a multiset object).
   G4cout << G4endl << " Information on Time per Event " << G4endl;
   G4int numberOfTimes = eventTimeSet.size();
   G4double tmin = 0.0, tmax = 0.0;
   sum = 0.0; sum2 = 0.0;
   G4double sumReduced = 0.0, sumReduced2 = 0.0;
   G4double t5perCent = 0.0, t95perCent = 0.0;
   const G4int num5perCent = numberOfTimes * 5 / 100 ;
   mu = sigma = mu_sigma = 0.0; 
   G4double muReduced = 0.0, sigmaReduced = 0.0, mu_sigmaReduced = 0.0;
   if ( numberOfTimes > 0 ) {
     G4int count = 0, countReduced = 0;
     for ( std::multiset<G4double>::const_iterator cit = eventTimeSet.begin();
 	  cit != eventTimeSet.end(); ++cit ) {
       count++;
       G4double t = *cit;
       //G4cout << "\t" << count - 1 << "\t" << t << "s" << G4endl;
       if ( count == 1 ) tmin = t;
       if ( count == num5perCent ) t5perCent = t;
       if ( count == ( numberOfTimes - num5perCent + 1 ) ) t95perCent = t;
       if ( count == numberOfTimes ) tmax = t;
       sum += t;
       sum2 += t*t;
       if ( count > num5perCent  &&  count <= ( numberOfTimes - num5perCent ) ) {
 	countReduced++;
 	sumReduced += t;
 	sumReduced2 += t*t;
       }
     }
     //G4cout << "\t countReduced = " << countReduced << G4endl; 
     n = static_cast< G4double >( numberOfTimes );
     mu = sum / n;
     if ( n > 1.0 ) {
       sigma = std::sqrt( std::abs( ( sum2 - sum*sum/n ) ) / (n - 1.0) );
     }
     mu_sigma = sigma / std::sqrt( n );
 
     n = static_cast< G4double >( countReduced );
     muReduced = sumReduced / n;
     if ( n > 1.0 ) {
       sigmaReduced = std::sqrt( std::abs( ( sumReduced2 - sumReduced*sumReduced/n ) )
 				/ (n - 1.0) );
     }
     mu_sigmaReduced = sigmaReduced / std::sqrt( n );
   }
   G4cout << "\t eventTimeSet.size() = " << numberOfTimes << G4endl
          << "\t num5perCent = " << num5perCent << G4endl
          << "\t tmin = " << tmin << "s" << G4endl
          << "\t t5perCent = " << t5perCent << "s" << G4endl
          << "\t t95perCent = " << t95perCent << "s" << G4endl
          << "\t tmax = " << tmax << "s" << G4endl
          //   << "\t sum = " << sum << " sec   sum2 = " << sum2 << " sec^2" << G4endl 
          << "\t mean = " << mu << " +/- " << mu_sigma << " sec" << G4endl
          << "\t rms = " << sigma << " sec" << G4endl
          //   << "\t sumReduced = " << sumReduced << " sec   sumReduced2 = " 
 	 //   << sumReduced2 << " sec^2" << G4endl
          << "\t meanReduced = " << muReduced << " +/- " << mu_sigmaReduced 
	  << " sec" << G4endl 
          << "\t rmsReduced = " << sigmaReduced << " sec" << G4endl;

}

