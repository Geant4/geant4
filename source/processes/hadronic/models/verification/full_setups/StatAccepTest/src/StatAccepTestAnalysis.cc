#include "StatAccepTestAnalysis.hh"
#include <string>
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"

#include <AIDA/AIDA.h>

//***TEMPORARY WORK-AROUND*** : need  AIDA_Dev/  subdirectory at the
//                              same level as  AIDA/  in PI/PI_1_3_0/include .
//                              But it does not work! 
// // //#include <AIDA_Dev/IDevHistogram1D.h>


StatAccepTestAnalysis* StatAccepTestAnalysis::instance = 0;

 
StatAccepTestAnalysis::StatAccepTestAnalysis() : 
  analysisFactory( 0 ), tree( 0 ), tuple( 0 ), histoFactory( 0 ), 
  numberOfEvents( 0 ), numberOfReplicas( 0 ), 
  numberOfReadoutLayers( 0 ), numberOfActiveLayersPerReadoutLayer( 1 ), 
  numberOfRadiusBins( 0 ), radiusBin( 0.0 ), 
  longitudinalProfileHisto( 0 ), transverseProfileHisto( 0 ) {

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
	  std::string description = "float ID, E, EDEP_ACT, EDEP_CAL; "; 
          description += "int nLayers, nBinR; ";
          description += "Tuple{ float L} L; ";
          description += "Tuple{ float R} R"; 
	  std::string option = "nLayers[0,100] L(nLayers) nBinR[0,30] R(nBinR)";
	  tuple = tupleFactory->create("1","Event info", description, option);
	  assert( tuple );
	  delete tupleFactory;
	}
        // Create a factory for histograms :
        histoFactory = analysisFactory->createHistogramFactory( *tree );  
      }
      delete treeFactory; // It will not delete the ITree.
    }
  }

}


StatAccepTestAnalysis::~StatAccepTestAnalysis() {}


void StatAccepTestAnalysis::close() {
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
}


StatAccepTestAnalysis* StatAccepTestAnalysis::getInstance() {
  if ( instance == 0 ) instance = new StatAccepTestAnalysis();
  return instance;
}


void StatAccepTestAnalysis::init() {
  // G4cout << " StatAccepTestAnalysis::init() : Cleaning up..." << G4endl;

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

  if ( tuple ) tuple->reset();
  longitudinalProfile.clear();
  beamEnergy  = 0.0;
  sumEdepAct  = 0.0;
  sumEdepAct2 = 0.0;
  sumEdepTot  = 0.0;
  sumEdepTot2 = 0.0;
  maxEdepTot  = 0.0;
  sumL.clear();
  sumL2.clear();
  for ( int layer = 0; layer < numberOfReadoutLayers; layer++ ) {
    longitudinalProfile.push_back( 0.0 );
    sumL.push_back( 0.0 );
    sumL2.push_back( 0.0 );
  }
  transverseProfile.clear();
  sumR.clear();
  sumR2.clear();
  for ( int ir = 0; ir < numberOfRadiusBins; ir++ ) {
    transverseProfile.push_back( 0.0 );
    sumR.push_back( 0.0 );
    sumR2.push_back( 0.0 );
  } 
  numberOfEvents = 0;
  vecEvis.clear();

  numStep = 0;
  numStepPositive = numStepNeutral = numStepNegative = 0;
  numStepPDGCodeZero = numStepPDGCodeUnrecognized = 0;
  numStepEM = 0;
  numStepEWK = 0;
  numStepHAD = 0; 
  numStepMeson = numStepBaryon = 0;     
  numStepMesonLight = numStepBaryonLight = 0;     
  numStepMesonStrange = numStepBaryonStrange = 0;     
  numStepMesonHeavy = numStepBaryonHeavy = 0;
  numStepElectron = numStepGamma = numStepPositron = 0;
  numStepMuMinus = numStepMuPlus = 0;
  numStepTauMinus = numStepTauPlus = 0;
  numStepNeutrino = 0;
  numStepPiPlus = numStepPi0 = numStepPiMinus = 0;
  numStepKPlus = numStepKNeutral = numStepKMinus = 0; 
  numStepProton = numStepAntiProton = 0;
  numStepNeutron = numStepAntiNeutron = 0;

  numTrack = 0;
  numTrackPositive = numTrackNeutral = numTrackNegative = 0;
  numTrackPDGCodeZero = numTrackPDGCodeUnrecognized = 0;
  numTrackEM = 0;
  numTrackEWK = 0;
  numTrackHAD = 0; 
  numTrackMeson = numTrackBaryon = 0;     
  numTrackMesonLight = numTrackBaryonLight = 0;     
  numTrackMesonStrange = numTrackBaryonStrange = 0;     
  numTrackMesonHeavy = numTrackBaryonHeavy = 0;
  numTrackElectron = numTrackGamma = numTrackPositron = 0;
  numTrackMuMinus = numTrackMuPlus = 0;
  numTrackTauMinus = numTrackTauPlus = 0;
  numTrackNeutrino = 0;
  numTrackPiPlus = numTrackPi0 = numTrackPiMinus = 0;
  numTrackKPlus = numTrackKNeutral = numTrackKMinus = 0; 
  numTrackProton = numTrackAntiProton = 0;
  numTrackNeutron = numTrackAntiNeutron = 0;

}                       


void StatAccepTestAnalysis::init( const G4int numberOfReplicasIn, 
				  const G4int numberOfReadoutLayersIn,
				  const G4int numberOfRadiusBinsIn,
				  const G4double radiusBinIn ) {
  if ( numberOfReplicasIn > 0 ) {
    numberOfReplicas = numberOfReplicasIn;
  }
  if ( numberOfReadoutLayersIn > 0 ) {
    if ( numberOfReplicasIn % numberOfReadoutLayersIn != 0 ) {
      G4cout << " ***WARNING*** StatAccepTestAnalysis::init(...) " << G4endl
             << " \t numberOfReplicas = " << numberOfReplicasIn
	     << "  is NOT compatible with  numberOfReadoutLayers = " 
	     << numberOfReadoutLayersIn << G4endl;
    } else {
      numberOfReadoutLayers = numberOfReadoutLayersIn;
      numberOfActiveLayersPerReadoutLayer = numberOfReplicas / numberOfReadoutLayers;
    }
  }
  if ( numberOfRadiusBinsIn > 0 ) {
    numberOfRadiusBins = numberOfRadiusBinsIn;
  }
  if ( radiusBinIn > 0 ) {
    radiusBin = radiusBinIn;
  }
  //G4cout << " StatAccepTestAnalysis::init( , , , ) : DEBUG Info " << G4endl
  //       << "\t numberOfReplicas      = " << numberOfReplicas << G4endl
  //       << "\t numberOfReadoutLayers = " << numberOfReadoutLayers << G4endl
  //       << "\t numberOfRadiusBins    = " << numberOfRadiusBins << G4endl
  //       << "\t radiusBin             = " << radiusBin/mm << " mm" 
  //       << G4endl;  //***DEBUG***

  // Create two histograms: one for the longitudinal shower profile,
  // and one for the transverse shower profile.
  assert( histoFactory );
  if ( histoFactory ) {
      // // // if ( ! tree->find( "50" ) ) {
      longitudinalProfileHisto = 
	histoFactory->createHistogram1D("50", "Longitudinal shower profile", 
					numberOfReadoutLayers, 0.0, 
					1.0*numberOfReadoutLayers );
      assert( longitudinalProfileHisto );
      // G4cout << " Created longitudinalProfileHisto " << G4endl;
      // // // if ( ! tree->find( "60" ) ) {
      transverseProfileHisto = 
	histoFactory->createHistogram1D("60", "Transverse shower profile", 
					numberOfRadiusBins, 0.0, 1.0*numberOfRadiusBins );
      assert( transverseProfileHisto );
      // G4cout << " Created transverseProfileHisto " << G4endl;
  }
  
}                       


void StatAccepTestAnalysis::fillNtuple( float incidentParticleId, 
					float incidentParticleEnergy, 
					float totalEnergyDepositedInActiveLayers,
					float totalEnergyDepositedInCalorimeter ) {
  beamEnergy = incidentParticleEnergy;
  if ( totalEnergyDepositedInCalorimeter > maxEdepTot ) {
    maxEdepTot = totalEnergyDepositedInCalorimeter;
  }
  if (tuple) {
    // G4cout << " StatAccepTestAnalysis::fillNtuple : DEBUG Info " << G4endl
    //        << "\t incidentParticleId = " << incidentParticleId << G4endl
    //        << "\t incidentParticleEnergy = " << incidentParticleEnergy << G4endl
    //        << "\t totalEnergyDepositedInActiveLayers = " 
    //        << totalEnergyDepositedInActiveLayers << G4endl
    //        << "\t totalEnergyDepositedInCalorimeter = " 
    //        << totalEnergyDepositedInCalorimeter << G4endl;       // ***DEBUG***

    tuple->fill( tuple->findColumn( "ID" ), incidentParticleId );
    tuple->fill( tuple->findColumn( "E" ), incidentParticleEnergy );
    tuple->fill( tuple->findColumn( "EDEP_ACT" ), totalEnergyDepositedInActiveLayers );
    tuple->fill( tuple->findColumn( "EDEP_CAL" ), totalEnergyDepositedInCalorimeter );
    tuple->fill( tuple->findColumn( "nLayers" ), numberOfReadoutLayers );
    tuple->fill( tuple->findColumn( "nBinR" ), numberOfRadiusBins );
    sumEdepAct  += totalEnergyDepositedInActiveLayers;
    sumEdepAct2 += 
      totalEnergyDepositedInActiveLayers * totalEnergyDepositedInActiveLayers;
    sumEdepTot  += totalEnergyDepositedInCalorimeter;
    sumEdepTot2 += 
      totalEnergyDepositedInCalorimeter * totalEnergyDepositedInCalorimeter;
    AIDA::ITuple* tpL = tuple->getTuple( tuple->findColumn( "L" ) );
    AIDA::ITuple* tpR = tuple->getTuple( tuple->findColumn( "R" ) );
    for ( int iLayer = 0; iLayer < numberOfReadoutLayers; iLayer++ ) {
      tpL->fill( 0, longitudinalProfile[ iLayer ] );
      tpL->addRow();
      sumL[ iLayer ]  += longitudinalProfile[ iLayer ];
      sumL2[ iLayer ] += longitudinalProfile[ iLayer ] * longitudinalProfile[ iLayer ];  
    }
    for ( int iBinR = 0; iBinR < numberOfRadiusBins; iBinR++ ) {
      tpR->fill( 0, transverseProfile[ iBinR ] );
      tpR->addRow();
      sumR[ iBinR ]  += transverseProfile[ iBinR ];
      sumR2[ iBinR ] += transverseProfile[ iBinR ] * transverseProfile[ iBinR ];  
    }
    tuple->addRow();
  }

  // Reset the longitudinal and transverse profiles, for the next event.
  for ( int layer = 0; layer < numberOfReadoutLayers; layer++ ) {
    // G4cout << " StatAccepTestAnalysis::fillNtuple : DEBUG Info " << G4endl
    //        << "\t Longitudinal profile: layer = " << layer
    //        << "   energy = " << longitudinalProfile[ layer ] / MeV 
    //        << " MeV " << G4endl;                                 //***DEBUG***
    longitudinalProfile[ layer ] = 0.0;
  }
  for ( int ir = 0; ir < numberOfRadiusBins; ir++ ) {
    // G4cout << " StatAccepTestAnalysis::fillNtuple : DEBUG Info " << G4endl
    //        << "\t Transverse profile: iBinRadius = " << ir / mm
    //        << " mm   energy = " << transverseProfile[ ir ] / MeV 
    //        << " MeV " << G4endl;                                 //***DEBUG***
    transverseProfile[ ir ] = 0.0;
  }

  // This method is called at each event, so it is useful to commit
  // the tree from time to time, for instance every 10 events, in
  // such a way that it is possible to see the ntuple  ntuple.hbook
  // while the job is running, for instance for a quick check that
  // it makes sense. Or, if the job crashes after many events, we
  // have anyhow some data already stored in the ntuple to be checked.
  // Remember that when you look the ntuple with PAW, while is running,
  // you need to close the session and open a new one to see the updates.
  numberOfEvents++;
  if ( numberOfEvents % 1000 == 0 ) {
    if ( tree ) {
      tree->commit();
      // G4cout << " tree commit ,  at event=" << numberOfEvents-1 
      //       << G4endl; //***DEBUG***
    }
  }

  // Store information of the visible energy, for later computing
  // of the energy resolution.
  vecEvis.push_back( totalEnergyDepositedInActiveLayers );

}


void StatAccepTestAnalysis::
fillShowerProfile( G4int replica, G4double radius, G4double edep ) {

  if ( replica >= numberOfReplicas ) {
    G4cout << " StatAccepTestAnalysis::fillShowerProfile : ***ERROR*** " << G4endl
           << "\t replica=" << replica 
	   << "  >=  numberOfReplicas=" << numberOfReplicas << G4endl; 
    replica = numberOfReplicas - 1;  // Just to avoid a crash
  }

  G4int readoutLayer = replica / numberOfActiveLayersPerReadoutLayer;

  longitudinalProfile[ readoutLayer ] += edep;

  // The last bin of the transverse profile includes all the hits with 
  // remaining radius. 
  // The bins are not constants: the specified radius refers to the first
  // (narrow) bin, then the second one has a width which is double the first
  // one, then the third has a width which is three time the first time,
  // and so on. The reason for this is to keep a reasonable statistics
  // in each bin, given the fast decrease as the radius increases.
  int iBinRadius = 0;
  G4double currentRadius = radiusBin;
  while ( radius > currentRadius  &&  iBinRadius < numberOfRadiusBins-1 ) {
    iBinRadius++;
    currentRadius += iBinRadius*radiusBin;
  }

  transverseProfile[ iBinRadius ] += edep;

  //G4cout << " StatAccepTestAnalysis::fillShowerProfile : DEBUG Info " << G4endl
  //       << " \t replica = " << replica 
  //       << "  readoutLayer = " << readoutLayer << G4endl
  //       << " \t radius = " << radius / mm 
  //       << " mm   iBinRadius = " << iBinRadius << G4endl
  //       << " \t edep = " << edep << " MeV "  << G4endl;  //***DEBUG***

}


void StatAccepTestAnalysis::infoStep( const G4Step* aStep ) {
  classifyParticle( false , aStep->GetTrack()->GetDefinition() );
}


void StatAccepTestAnalysis::infoTrack( const G4Track* aTrack ) {
  classifyParticle( true , aTrack->GetDefinition() );
}


void StatAccepTestAnalysis::
classifyParticle( const bool isTrack, const G4ParticleDefinition* particleDef ) {

  if ( isTrack ) {
    numTrack++;
  } else {
    numStep++;
  }

  G4int id = particleDef->GetPDGEncoding();

  G4int aid = std::abs( id );
  
  if ( id == 0 ) {      // Currently, nuclear fragments have zero PDG code.
    if ( isTrack ) {
      numTrackPDGCodeZero++;
    } else {
      numStepPDGCodeZero++;
    }
  } else if ( particleDef == G4Gamma::GammaDefinition() ) {
    if ( isTrack ) {
      numTrackNeutral++;
      numTrackEM++;
      numTrackGamma++;
    } else {
      numStepNeutral++;
      numStepEM++;
      numStepGamma++;
    }
  } else if ( particleDef == G4Electron::ElectronDefinition() ) {
    if ( isTrack ) {
      numTrackNegative++;
      numTrackEM++;
      numTrackElectron++;
    } else {
      numStepNegative++;
      numStepEM++;
      numStepElectron++;
    }
  } else if ( particleDef == G4Positron::PositronDefinition() ) {
    if ( isTrack ) {
      numTrackPositive++;
      numTrackEM++;
      numTrackPositron++;
    } else {
      numStepPositive++;
      numStepEM++;
      numStepPositron++;
    }
  } else if ( particleDef == G4MuonMinus::MuonMinusDefinition() ) {
    if ( isTrack ) {
      numTrackNegative++;
      numTrackEWK++;
      numTrackMuMinus++;
    } else {
      numStepNegative++;
      numStepEWK++;
      numStepMuMinus++;
    } 
  } else if ( particleDef == G4MuonPlus::MuonPlusDefinition() ) {
    if ( isTrack ) {
      numTrackPositive++;
      numTrackEWK++;
      numTrackMuPlus++;
    } else {
      numStepPositive++;
      numStepEWK++;
      numStepMuPlus++;
    }
  } else if ( particleDef == G4TauMinus::TauMinusDefinition() ) {
    if ( isTrack ) {
      numTrackNegative++;
      numTrackEWK++;
      numTrackTauMinus++;
    } else {
      numStepNegative++;
      numStepEWK++; 
      numStepTauMinus++;
    } 
  } else if ( particleDef == G4TauPlus::TauPlusDefinition() ) {
    if ( isTrack ) {
      numTrackPositive++;
      numTrackEWK++;
      numTrackTauPlus++;
    } else {
      numStepPositive++;
      numStepEWK++; 
      numStepTauPlus++;
    } 
  } else if ( particleDef == G4NeutrinoE::NeutrinoEDefinition()  ||
	      particleDef == G4NeutrinoMu::NeutrinoMuDefinition()  ||
	      particleDef == G4NeutrinoTau::NeutrinoTauDefinition()  ||
	      particleDef == G4AntiNeutrinoE::AntiNeutrinoEDefinition()  ||
	      particleDef == G4AntiNeutrinoMu::AntiNeutrinoMuDefinition()  ||
	      particleDef == G4AntiNeutrinoTau::AntiNeutrinoTauDefinition() ) {
    if ( isTrack ) {
      numTrackNeutral++;
      numTrackEWK++;
      numTrackNeutrino++;
    } else {
      numStepNeutral++;
      numStepEWK++;
      numStepNeutrino++;
    }
  } else if ( particleDef == G4PionMinus::PionMinusDefinition() ) {
    if ( isTrack ) {
      numTrackNegative++;
      numTrackHAD++;
      numTrackMeson++;
      numTrackMesonLight++;
      numTrackPiMinus++;
    } else {
      numStepNegative++;
      numStepHAD++;
      numStepMeson++;
      numStepMesonLight++;
      numStepPiMinus++;
    }
  } else if ( particleDef == G4PionZero::PionZeroDefinition() ) {
    if ( isTrack ) {
      numTrackNeutral++;
      numTrackHAD++;
      numTrackMeson++;
      numTrackMesonLight++;
      numTrackPi0++;
    } else {
      numStepNeutral++;
      numStepHAD++;
      numStepMeson++;
      numStepMesonLight++;
      numStepPi0++;
    }
  } else if ( particleDef == G4PionPlus::PionPlusDefinition() ) {
    if ( isTrack ) {
      numTrackPositive++;
      numTrackHAD++;
      numTrackMeson++;
      numTrackMesonLight++;
      numTrackPiPlus++;
    } else {
      numStepPositive++;
      numStepHAD++;
      numStepMeson++;
      numStepMesonLight++;
      numStepPiPlus++;
    }
  } else if ( particleDef == G4KaonMinus::KaonMinusDefinition() ) {
    if ( isTrack ) {
      numTrackNegative++;
      numTrackHAD++;
      numTrackMeson++;
      numTrackMesonStrange++;
      numTrackKMinus++;
    } else {
      numStepNegative++;
      numStepHAD++;
      numStepMeson++;
      numStepMesonStrange++;
      numStepKMinus++;
    }
  } else if ( particleDef == G4KaonZero::KaonZeroDefinition() ||
	      particleDef == G4AntiKaonZero::AntiKaonZeroDefinition() || 
	      particleDef == G4KaonZeroLong::KaonZeroLongDefinition() ||
	      particleDef == G4KaonZeroShort::KaonZeroShortDefinition() ) {
    if ( isTrack ) {
      numTrackNeutral++;
      numTrackHAD++;
      numTrackMeson++;
      numTrackMesonStrange++;
      numTrackKNeutral++;
    } else {
      numStepNeutral++;
      numStepHAD++;
      numStepMeson++;
      numStepMesonStrange++;
      numStepKNeutral++;
    }
  } else if ( particleDef == G4KaonPlus::KaonPlusDefinition() ) {
    if ( isTrack ) {
      numTrackPositive++;
      numTrackHAD++;
      numTrackMeson++;
      numTrackMesonStrange++;
      numTrackKPlus++;
    } else {
      numStepPositive++;
      numStepHAD++;
      numStepMeson++;
      numStepMesonStrange++;
      numStepKPlus++;
    }
  } else if ( particleDef == G4Neutron::NeutronDefinition() ) {
    if ( isTrack ) {
      numTrackNeutral++;
      numTrackHAD++;
      numTrackBaryon++;
      numTrackBaryonLight++;
      numTrackNeutron++;
    } else {
      numStepNeutral++;
      numStepHAD++;
      numStepBaryon++;
      numStepBaryonLight++;
      numStepNeutron++;
    }
  } else if ( particleDef == G4AntiNeutron::AntiNeutronDefinition() ) {
    if ( isTrack ) {
      numTrackNeutral++;
      numTrackHAD++;
      numTrackBaryon++;
      numTrackBaryonLight++;
      numTrackAntiNeutron++;
    } else {
      numStepNeutral++;
      numStepHAD++;
      numStepBaryon++;
      numStepBaryonLight++;
      numStepAntiNeutron++;
    }
  } else if ( particleDef == G4Proton::ProtonDefinition() ) {
    if ( isTrack ) {
      numTrackPositive++;
      numTrackHAD++;
      numTrackBaryon++;
      numTrackBaryonLight++;
      numTrackProton++;
    } else {
      numStepPositive++;
      numStepHAD++;
      numStepBaryon++;
      numStepBaryonLight++;
      numStepProton++;
    }
  } else if ( particleDef == G4AntiProton::AntiProtonDefinition() ) {
    if ( isTrack ) {
      numTrackNegative++;
      numTrackHAD++;
      numTrackBaryon++;
      numTrackBaryonLight++;
      numTrackAntiProton++;
    } else {
      numStepNegative++;
      numStepHAD++;
      numStepBaryon++;
      numStepBaryonLight++;
      numStepAntiProton++;
    }
    // Other (light) mesons
  } else if ( aid == 9000211  ||  aid == 100211  ||  // Light I = 1 mesons
	      aid == 10211  ||  aid == 200211  ||
	      aid == 213  ||  aid == 10213  ||
	      aid == 20213  ||  aid == 9000213  ||
	      aid == 100213  ||  aid == 9010213  ||
	      aid == 9020213  ||  aid == 30213  ||
              aid == 9030213  ||  
	      aid == 215  ||  aid == 9000215  ||
	      aid == 10215  ||  aid == 9010215  ||
	      aid == 9020215 ||  
	      aid == 217  ||  aid == 9000217  ||
	      aid == 219 ) {
    if ( isTrack ) {
      if ( id > 0 ) {
	numTrackPositive++;
      } else {
	numTrackNegative++;
      }
      numTrackHAD++;
      numTrackMeson++;
      numTrackMesonLight++;
    } else {
      if ( id > 0 ) {
	numStepPositive++;
      } else {
	numStepNegative++;
      }
      numStepHAD++;
      numStepMeson++;
      numStepMesonLight++;
    }
  } else if ( aid == 9000111  ||  aid == 100111  ||  // Light I = 1 mesons
	      aid == 10111  ||  aid == 200111  ||
	      aid == 113  ||  aid == 10113  ||
	      aid == 20113  ||  aid == 9000113  ||
	      aid == 100113  ||  aid == 9010113  ||
	      aid == 9020113  ||  aid == 30113  ||
              aid == 9030113  ||  
	      aid == 115  ||  aid == 9000115  ||
	      aid == 10115  ||  aid == 9010115  ||
	      aid == 9020115 ||  
	      aid == 117  ||  aid == 9000117  ||
	      aid == 119  ||
	      aid == 221  ||  aid == 331  ||         // Light I = 0 mesons
	      aid == 9000221  ||  aid == 9010221  ||
	      aid == 100221  ||  aid == 10221  ||
	      aid == 100331  ||  aid == 9020221  ||
	      aid == 10331  ||  aid == 200221  ||
	      aid == 9030221  ||  aid == 9040221  ||
	      aid == 9050221  ||  aid == 9060221  ||
	      aid == 223  ||  aid == 333  ||
	      aid == 10223  ||  aid == 20223  ||
	      aid == 10333  ||  aid == 20333  ||
	      aid == 100223  ||  aid == 9000223  ||
	      aid == 30223  ||  aid == 100333  ||
	      aid == 225  ||  aid == 9000225  ||
	      aid == 335  ||  aid == 9010225  ||
	      aid == 9020225  ||  aid == 10225  ||
	      aid == 100225  ||  aid == 10335  ||
	      aid == 9030225  ||  aid == 100335  ||
	      aid == 9040225  ||  aid == 9050225  ||
	      aid == 9060225  ||  aid == 227  ||
	      aid == 337  ||  aid == 229  ||
	      aid == 9000339  ||  aid == 9000229 ) {
    if ( isTrack ) {
      numTrackNeutral++;
      numTrackHAD++;
      numTrackMeson++;
      numTrackMesonLight++;
    } else {
      numStepNeutral++;
      numStepHAD++;
      numStepMeson++;
      numStepMesonLight++;
    }
    // Other light baryons.
  } else if ( aid == 2224  ||  aid == 2214 ) {   // Delta++ , Delta+
    if ( isTrack ) {
      if ( id > 0 ) {
	numTrackPositive++;
      } else {
	numTrackNegative++;
      }
      numTrackHAD++;
      numTrackBaryon++;
      numTrackBaryonLight++;
    } else {
      if ( id > 0 ) {
	numStepPositive++;
      } else {
	numStepNegative++;
      }
      numStepHAD++;
      numStepBaryon++;
      numStepBaryonLight++;
    }
  } else if ( aid == 2114 ) {   // Delta0
    if ( isTrack ) {
      numTrackNeutral++;
      numTrackHAD++;
      numTrackBaryon++;
      numTrackBaryonLight++;
    } else {
      numStepNeutral++;
      numStepHAD++;
      numStepBaryon++;
      numStepBaryonLight++;
    }
  } else if ( aid == 1114 ) {   // Delta-
    if ( isTrack ) {
      if ( id > 0 ) {
	numTrackNegative++;
      } else {
	numTrackPositive++;
      }
      numTrackHAD++;
      numTrackBaryon++;
      numTrackBaryonLight++;
    } else {
      if ( id > 0 ) {
	numStepNegative++;
      } else {
	numStepPositive++;
      }
      numStepHAD++;
      numStepBaryon++;
      numStepBaryonLight++;
    }
    // Other Strange mesons
  } else if ( aid == 10321  ||  aid == 100321  ||
	      aid == 200321  ||  aid == 9000321  ||
	      aid == 323  ||  aid == 10323  ||
	      aid == 20323  ||  aid == 100323  ||
	      aid == 9000323  ||  aid == 30323  ||
	      aid == 325  ||  aid == 9000325  ||
	      aid == 10325  ||  aid == 20325  ||
	      aid == 100325  ||  aid == 9010325  ||
	      aid == 327  ||  aid == 9010327  ||
	      aid == 329  ||  aid == 9000329 ) {
    if ( isTrack ) {
      if ( id > 0 ) {
	numTrackPositive++;
      } else {
	numTrackNegative++;
      }
      numTrackHAD++;
      numTrackMeson++;
      numTrackMesonStrange++;
    } else {
      if ( id > 0 ) {
	numStepPositive++;
      } else {
	numStepNegative++;
      }
      numStepHAD++;
      numStepMeson++;
      numStepMesonStrange++;
    }
  } else if ( aid == 10311  ||  aid == 100311  ||
	      aid == 200311  ||  aid == 9000311  ||
	      aid == 313  ||  aid == 10313  ||
	      aid == 20313  ||  aid == 100313  ||
	      aid == 9000313  ||  aid == 30313  ||
	      aid == 315  ||  aid == 9000315  ||
	      aid == 10315  ||  aid == 20315  ||
	      aid == 100315  ||  aid == 9010315  ||
	      aid == 317  ||  aid == 9010317  ||
	      aid == 319  ||  aid == 9000319 ) {
    if ( isTrack ) {
      numTrackNeutral++;
      numTrackHAD++;
      numTrackMeson++;
      numTrackMesonStrange++;
    } else {
      numStepNeutral++;
      numStepHAD++;
      numStepMeson++;
      numStepMesonStrange++;
    } 
    // Strange baryons
  } else if ( aid == 3222  ||  aid == 3224 ) {   // Sigma+ , Sigma*+
    if ( isTrack ) {
      if ( id > 0 ) {
	numTrackPositive++;
      } else {
	numTrackNegative++;
      }
      numTrackHAD++;
      numTrackBaryon++;
      numTrackBaryonStrange++;
    } else {
      if ( id > 0 ) {
	numStepPositive++;
      } else {
	numStepNegative++;
      }
      numStepHAD++;
      numStepBaryon++;
      numStepBaryonStrange++;
    }
  } else if ( aid == 3122  ||                     // Lambda0
              aid == 3212  ||  aid == 3214  ||    // Sigma0 , Sigma*0
	      aid == 3322  ||  aid == 3324 ) {    // Csi0 , Csi*0
    if ( isTrack ) {
      numTrackNeutral++;
      numTrackHAD++;
      numTrackBaryon++;
      numTrackBaryonStrange++;
    } else {
      numStepNeutral++;
      numStepHAD++;
      numStepBaryon++;
      numStepBaryonStrange++;
    }
  } else if ( aid == 3334  ||                     // Omega-
              aid == 3112  ||  aid == 3114  ||    // Sigma- , Sigma*-
	      aid == 3312  ||  aid == 3314 ) {    // Csi- , Csi*-
    if ( isTrack ) {
      if ( id > 0 ) {
	numTrackNegative++;
      } else {
	numTrackPositive++;
      }
      numTrackHAD++;
      numTrackBaryon++;
      numTrackBaryonStrange++;
    } else {
      if ( id > 0 ) {
	numStepNegative++;
      } else {
	numStepPositive++;
      }
      numStepHAD++;
      numStepBaryon++;
      numStepBaryonStrange++;
    }
    // Heavy mesons
  } else if ( aid == 411  ||  aid == 10411  ||   // D+ , D0*+
	      aid == 413  ||  aid == 10413  ||   // D*(2010)+ , D1(2420)+ 
	      aid == 20413  ||                   // D1(H)+
	      aid == 415  ||                     // D2*(2460)+
	      aid == 431  ||  aid == 10431  ||   // Ds+ , Ds0*+
	      aid == 433  ||  aid == 10433  ||   // Ds*+, Ds1(2536)+
	      aid == 20433  ||                   // Ds1(H)^+
	      aid == 435  ||                     // Ds2*+
	      aid == 521  ||  aid == 10521  ||   // B+ , B0*+
	      aid == 523  ||  aid == 10523  ||   // B*+ , B1(L)+
	      aid == 20523  ||                   // B1(H)+
	      aid == 525  ||                     // B2*+
	      aid == 541  ||  aid == 10541  ||   // Bc+, B0c*+
	      aid == 543  ||  aid == 10543  ||   // Bc*+, Bc1(L)+
	      aid == 20543  ||                   // Bc1(H)+
	      aid == 545 ) {                     // Bc2*+
    if ( isTrack ) {
      if ( id > 0 ) {
	numTrackPositive++;
      } else {
	numTrackNegative++;
      }
      numTrackHAD++;
      numTrackMeson++;
      numTrackMesonHeavy++;
    } else {
      if ( id > 0 ) {
	numStepPositive++;
      } else {
	numStepNegative++;
      }
      numStepHAD++;
      numStepMeson++;
      numStepMesonHeavy++;
    }
  } else if ( aid == 421  ||  aid == 10421  ||   // D0 , D0*0
	      aid == 423  ||  aid == 10423  ||   // D*(2007), D1(2420)0
	      aid == 20423  ||                   // D1(H)0
	      aid == 425  ||                     // D2*(2460)0
	      aid == 511  ||  aid == 10511  ||   // B0 , B0*0
	      aid == 513  ||  aid == 10513  ||   // B*0, B1(L)0
	      aid == 20513  ||                   // B1(H)0
	      aid == 515  ||                     // B2*0
	      aid == 531  ||  aid == 10531  ||   // Bs0 , Bs0*0
	      aid == 533  ||  aid == 10533  ||   // Bs*0 , Bs1(L)0
	      aid == 20533  ||                   // Bs1(H)0
	      aid == 535  ||                     // Bs2*0
	      aid == 441  ||  aid == 10441  ||            // ccbar mesons
	      aid == 100441  ||                 
	      aid == 443  ||  aid == 10443  ||  
	      aid == 20443  ||  aid == 100443  ||
	      aid == 30443  || aid == 9000443  ||
	      aid == 9010443  ||  aid == 9020443  ||
	      aid == 445  ||  aid == 9000445  ||
	      aid == 551  ||  aid == 10551  ||            // bbar mesons
	      aid == 100551  ||  aid == 110551  ||
	      aid == 200551  ||  aid == 210551  ||
	      aid == 553  ||  aid == 10553  ||
	      aid == 20553  ||  aid == 30553  ||
	      aid == 100553  ||  aid == 110553  ||
	      aid == 120553  ||  aid == 130553  ||
	      aid == 200553  ||  aid == 210553  ||
	      aid == 220553  ||  aid == 300553  ||
	      aid == 9000553  ||  aid == 9010553  ||
	      aid == 555  ||  aid == 10555  ||
	      aid == 20555  ||  aid == 100555  ||
	      aid == 110555  ||  aid == 120555  ||
	      aid == 200555  ||
	      aid == 557  ||  aid == 100557 ) {
    if ( isTrack ) {
      numTrackNeutral++;
      numTrackHAD++;
      numTrackMeson++;
      numTrackMesonHeavy++;
    } else {
      numStepNeutral++;
      numStepHAD++;
      numStepMeson++;
      numStepMesonHeavy++;
    }
    // Heavy baryons
  } else if ( aid == 4122  ||                    // Lambda_c+
	      aid == 4222  ||  aid == 4224  ||   // Sigma_c++ , Sigma_c*++
	      aid == 4212  ||  aid == 4214  ||   // Sigma_c+ , Sigma_c*+
              aid == 4232  ||  aid == 4324  ||   // Csi_c+ , Csi_c*+
	      aid == 4322  ||                    // Csi_c'+
              aid == 4412  ||  aid == 4414  ||   // Csi_cc+ , Csi_cc*+
	      aid == 4422  ||  aid == 4424  ||   // Csi_cc++ , Csi_cc*++
	      aid == 4432  ||  aid == 4434  ||   // Omega_cc+ , Omega_cc*+
	      aid == 4444  ||                    // Omega_ccc++
	      aid == 5222  ||  aid == 5224  ||   // Sigma_b+ , Sigma_b*+
	      aid == 5242  ||  aid == 5424  ||   // Csi_bc+ , Csi_bc*+
              aid == 5422  ||                    // Csi_bc'+
              aid == 5442  ||  aid == 5444 ) {   // Omega_bcc*+
    if ( isTrack ) {
      if ( id > 0 ) {
	numTrackPositive++;
      } else {
	numTrackNegative++;
      }
      numTrackHAD++;
      numTrackBaryon++;
      numTrackBaryonHeavy++;
    } else {
      if ( id > 0 ) {
	numStepPositive++;
      } else {
	numStepNegative++;
      }
      numStepHAD++;
      numStepBaryon++;
      numStepBaryonHeavy++;
    }
  } else if ( aid == 4112  ||  aid == 4114  ||   // Sigma_c0 , Sigma_c*0
	      aid == 4132  ||  aid == 4314  ||   // Csi_c0 , Csi_c*0
	      aid == 4312  ||                    // Csi_c'0
	      aid == 4332  ||  aid == 4334  ||   // Omega_c0 , Omega_c*0
	      aid == 5122  ||                    // Lambda_b0
	      aid == 5212  ||  aid == 5214  ||   // Sigma_b0 , Sigma_b*0
	      aid == 5232  ||  aid == 5324  ||   // Csi_b0 , Csi_b*0
	      aid == 5322  ||                    // Csi_b'0
	      aid == 5142  ||  aid == 5414  ||   // Csi_bc0 , Csi_bc*0
	      aid == 5412  ||                    // Csi_bc'0
	      aid == 5342  ||  aid == 5434  ||   // Omega_bc0 , Omega_bc*0
	      aid == 5432  ||                    // Omega_bc'0
	      aid == 5522  ||  aid == 5524  ||   // Csi_bb0, Csi_bb*0
	      aid == 5542  ||  aid == 5544 ) {   // Omega_bbc0, Omega_bbc*0
    if ( isTrack ) {
      numTrackNeutral++;
      numTrackHAD++;
      numTrackBaryon++;
      numTrackBaryonHeavy++;
    } else {
      numStepNeutral++;
      numStepHAD++;
      numStepBaryon++;
      numStepBaryonHeavy++;
    }
  } else if ( aid == 5112  ||  aid == 5114  ||   // Sigma_b- , Sigma_b*-
              aid == 5132  ||  aid == 5314  ||   // Csi_b- , Csi_b*-
	      aid == 5312  ||                    // Csi_b'-
	      aid == 5332  ||  aid == 5334  ||   // Omega_b- , Omega_b*-
	      aid == 5512  ||  aid == 5514  ||   // Csi_bb- , Csi_bb*-
	      aid == 5532  ||  aid == 5534  ||   // Omega_bb-, Omega_bb*-
	      aid == 5554 ) {                    // Omega_bbb-
    if ( isTrack ) {
      if ( id > 0 ) {
	numTrackNegative++;
      } else {
	numTrackPositive++;
      }
      numTrackHAD++;
      numTrackBaryon++;
      numTrackBaryonHeavy++;
    } else {
      if ( id > 0 ) {
	numStepNegative++;
      } else {
	numStepPositive++;
      }
      numStepHAD++;
      numStepBaryon++;
      numStepBaryonHeavy++;
    }
  } else {
    if ( isTrack ) {
      G4cout << " UNRECOGNIZED PARTICLE: id = " << id << G4endl;
      numTrackPDGCodeUnrecognized++;
    } else {
      numStepPDGCodeUnrecognized++;
    }
  }
}


void StatAccepTestAnalysis::finish() {

  // Notice that the errors that are calculated here are not the
  // precisely correct ones, which would need also the pure weights
  // besides the values  x * weight . Only the latter are used 
  // below to calculate approximate errors, using the standard
  // formulas, as if the values  x * weight  were pure values
  // without a weight.

  // For the summary histograms for the longitudinal and transverse 
  // shower shapes (longitudinalProfileHisto, transverseProfileHisto
  // respectively), for the time being, in AIDA is not yet possible
  // to set the errors, so the errors are automatically set to the 
  // square root of the y-value of that bin (which is what is called 
  // "mu" below, rather than "mu_sigma"), which does not make any
  // sense! So, ignore eventual failures in the statistical tests
  // only for these two variables.

  // Print results. 
  G4cout << " Beam energy [MeV] = " << beamEnergy << G4endl;

  // Check for the energy conservation by comparing the beam energy
  // with the maximum total energy deposit.
  // We write a message of energy-non-conservation in the case that
  // the maximum total energy deposit is greater than the beam energy.
  // However, in this case, one cannot always conclude that there is
  // some violation of the energy conservation in some of the Geant4
  // processes, because different physics effects could explain such
  // "anomaly", for example:
  //    - biasing techniques could produce big energy non conservation;
  //    - decays of particles can transform part of their masses in
  //      kinetic energy of the decay products: for example 
  //      K+ -> mu+ nu_mu about half of the mass of the kaon (i.e.
  //      about 250 MeV) goes to the muon, whose mass is only a part
  //      of it (about 100 MeV), and so the remaining goes in kinetic
  //      energy of the muon that can be deposited in a large calorimeter.
  //    - fission of heavy elements, like Uranium, W, Pb, etc...
  // Notice that, a part the first effect, the other two are usually
  // small effects: you cannot obtain large energy violation!
  // So, it is at the level of post-processing, that the user should
  // decide if an eventual message of energy-non-conservation makes sense
  // or not.
  G4cout << " maxEdepTot  [MeV] = " << maxEdepTot;
  if ( maxEdepTot - beamEnergy > 0.001*MeV ) {
    G4cout << "\t ***ENERGY-NON-CONSERVATION*** ";
  } 
  G4cout << G4endl;

  G4double n = static_cast< G4double >( numberOfEvents );
  if ( n <= 1.0 ) {
    n = 2.0;         // To avoid division by zero in  sigma
  }
  // G4cout << " n=" << n << G4endl;                             //***DEBUG***
  G4double sum, sum2, mu, sigma, mu_sigma;
  sum  = sumEdepAct;
  sum2 = sumEdepAct2;
  mu       = sum / n;
  sigma    = sqrt( ( sum2 - sum*sum/n ) / (n - 1.0) );
  mu_sigma = sigma / sqrt( n );
  G4cout << " Average <E> [MeV] deposited in all active layers = " 
         << mu << " +/- " << mu_sigma << G4endl;
  G4double mu_Evis = mu;  // For later usage.
  sum  = sumEdepTot;
  sum2 = sumEdepTot2;
  mu       = sum / n;
  sigma    = sqrt( ( sum2 - sum*sum/n ) / (n - 1.0) );
  mu_sigma = sigma / sqrt( n );
  G4cout << " Average <E> [MeV] deposited in the whole calorimeter = " 
         << mu << " +/- " << mu_sigma << G4endl;
  G4cout << " Average <E> [MeV] in each Layer " << G4endl; 
  for ( int iLayer = 0; iLayer < numberOfReadoutLayers; iLayer++ ) {
    sum  = sumL[ iLayer ];
    sum2 = sumL2[ iLayer ];
    mu       = sum / n;
    sigma    = sqrt( ( sum2 - sum*sum/n ) / (n - 1.0) );
    mu_sigma = sigma / sqrt( n );
    //if ( mu > 1.0E-06 ) {
    G4cout << "\t layer = " << iLayer << "\t <E> = " 
	   << mu << " +/- " << mu_sigma << G4endl;
    //}
    if ( longitudinalProfileHisto ) {
      //***LOOKHERE*** :  mu_sigma  should be set as error.
      longitudinalProfileHisto->fill( 1.0*iLayer, mu );
    }
  }
  // // // std::vector< G4double > rmsTransverseProfile;   //***TEMPORARY WORK-AROUND***
  G4cout << " Average <E> [MeV] in each Radius bin " << G4endl; 
  for ( int iBinR = 0; iBinR < numberOfRadiusBins; iBinR++ ) {
    sum  = sumR[ iBinR ];
    sum2 = sumR2[ iBinR ];
    mu       = sum / n;
    sigma    = sqrt( ( sum2 - sum*sum/n ) / (n - 1.0) );
    mu_sigma = sigma / sqrt( n );
    //if ( mu > 1.0E-06 ) {
    G4cout << "\t iBinR = " << iBinR << "\t <E> = " 
	   << mu << " +/- " << mu_sigma << G4endl;
    //}
    if ( transverseProfileHisto ) {
      //***LOOKHERE*** :  mu_sigma  should be set as error.
      transverseProfileHisto->fill( 1.0*iBinR, mu );
    }
    // // // rmsTransverseProfile.push_back( mu_sigma );   //***TEMPORARY WORK-AROUND***
  }  

  // ***TEMPORARY WORK-AROUND*** in order to set properly the error bars
  // of a weighted histogram, until this possibility will be added 
  // in AIDA::IHistogram1D . But it does not work!
  // // // AIDA::Dev::IDevHistogram1D * devTransverseProfileHisto = 
  // // //   dynamic_cast< AIDA::Dev::IDevHistogram1D * > ( transverseProfileHisto ); 
  // // // assert( devTransverseProfileHisto ); 
  // // // for (int iBinR = 0; iBinR < devTransverseProfileHisto->axis().bins() ; iBinR++ ) { 
  // // //   devTransverseProfileHisto->setBinContents( iBinR, 
  // // // 					       transverseProfileHisto->binEntries( iBinR ),
  // // // 					       transverseProfileHisto->binHeight( iBinR ),
  // // // 					       rmsTransverseProfile[ iBinR ], 
  // // // 					       transverseProfileHisto->binMean( iBinR ) );
  // // // }

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
  G4double width_Evis = 0.0;
  for ( std::vector< G4double >::const_iterator cit = vecEvis.begin();
	cit != vecEvis.end() ; ++cit ) {
    width_Evis += std::abs( *cit - mu_Evis );
  }
  width_Evis *= sqrt( 3.141592654/2.0 ) / n ;
  G4cout << " Visible energy information [MeV] " << G4endl; 
  G4cout << "\t mu_Evis    = " << mu_Evis << G4endl
         << "\t sigma_Evis = " << width_Evis << G4endl
         << "\t energy resolution = " << width_Evis/mu_Evis << G4endl
         << "\t sampling fraction = " << mu_Evis/beamEnergy << G4endl;

  if ( tree ) tree->commit();

  // Print information regarding the average number of steps and tracks
  // per event.
  for ( int i = 0; i < 2; i++ ) {  // 0 : steps ; 1 : tracks. 
    if ( i == 0 ) {
      G4cout << G4endl << " Average number of STEPS per Event " << G4endl;
    } else {
      G4cout << G4endl << " Average number of TRACKS per Event " << G4endl;
    }
    for ( int j = 0; j < 35; j++ ) {  // Loop over the 35 cases considered.
      std::string caseName = "";
      sum = 0.0, sum2 = 0.0;
      switch ( j ) {
      case 0 : {
	caseName += "# total";
	if ( i == 0 ) {
	  sum = static_cast< G4double >( numStep );
	} else {
	  sum = static_cast< G4double >( numTrack );
	}
	break;
      }
      case 1 : {
	caseName += "# positives";
	if ( i == 0 ) {
	  sum = static_cast< G4double >( numStepPositive );
	} else {
	  sum = static_cast< G4double >( numTrackPositive );
	}
	break;
      }
      case 2 : {
	caseName += "# neutrals";
	if ( i == 0 ) {
	  sum = static_cast< G4double >( numStepNeutral );
	} else {
	  sum = static_cast< G4double >( numTrackNeutral );
	}
	break;
      }
      case 3 : {
	caseName += "# negatives";
	if ( i == 0 ) {
	  sum = static_cast< G4double >( numStepNegative );
	} else {
	  sum = static_cast< G4double >( numTrackNegative );
	}
	break;
      }
      case 4 : {
	caseName += "# particles with 0 PDG code";
	if ( i == 0 ) {
	  sum = static_cast< G4double >( numStepPDGCodeZero );
	} else {
	  sum = static_cast< G4double >( numTrackPDGCodeZero );
	}
	break;
      }
      case 5 : {
	caseName += "# particles with Unrecognized PDG code";
	if ( i == 0 ) {
	  sum = static_cast< G4double >( numStepPDGCodeUnrecognized );
	} else {
	  sum = static_cast< G4double >( numTrackPDGCodeUnrecognized );
	}
	break;
      }
      case 6 : {
	caseName += "# electromagnetic (e+ , e- , gammas)";
	if ( i == 0 ) {
	  sum = static_cast< G4double >( numStepEM );
	} else {
	  sum = static_cast< G4double >( numTrackEM );
	}
	break;
      }
      case 7 : {
	caseName += "# electroweak (mu+, mu-, tau+, tau-, neutrinos)";
	if ( i == 0 ) {
	  sum = static_cast< G4double >( numStepEWK );
	} else {
	  sum = static_cast< G4double >( numTrackEWK );
	}
	break;
      }
      case 8 : {
	caseName += "# hadrons";
	if ( i == 0 ) {
	  sum = static_cast< G4double >( numStepHAD );
	} else {
	  sum = static_cast< G4double >( numTrackHAD );
	}
	break;
      }
      case 9 : {
	caseName += "# mesons";
	if ( i == 0 ) {
	  sum = static_cast< G4double >( numStepMeson );
	} else {
	  sum = static_cast< G4double >( numTrackMeson );
	}
	break;
      }
      case 10 : {
	caseName += "# baryons";
	if ( i == 0 ) {
	  sum = static_cast< G4double >( numStepBaryon );
	} else {
	  sum = static_cast< G4double >( numTrackBaryon );
	}
	break;
      }
      case 11 : {
	caseName += "# light mesons (u/ubar/d/dbar)";
	if ( i == 0 ) {
	  sum = static_cast< G4double >( numStepMesonLight );
	} else {
	  sum = static_cast< G4double >( numTrackMesonLight );
	}
	break;
      }
      case 12 : {
	caseName += "# light baryons (u/ubar/d/dbar)";
	if ( i == 0 ) {
	  sum = static_cast< G4double >( numStepBaryonLight );
	} else {
	  sum = static_cast< G4double >( numTrackBaryonLight );
	}
	break;
      }
      case 13 : {
	caseName += "# strange (s/sbar) mesons";
	if ( i == 0 ) {
	  sum = static_cast< G4double >( numStepMesonStrange );
	} else {
	  sum = static_cast< G4double >( numTrackMesonStrange );
	}
	break;
      }
      case 14 : {
	caseName += "# strange (s/sbar) baryons";
	if ( i == 0 ) {
	  sum = static_cast< G4double >( numStepBaryonStrange );
	} else {
	  sum = static_cast< G4double >( numTrackBaryonStrange );
	}
	break;
      }
      case 15 : {
	caseName += "# heavy (c/cbar or b/bbar) mesons";
	if ( i == 0 ) {
	  sum = static_cast< G4double >( numStepMesonHeavy );
	} else {
	  sum = static_cast< G4double >( numTrackMesonHeavy );
	}
	break;
      }
      case 16 : {
	caseName += "# heavy (c/cbar or b/bbar) baryons";
	if ( i == 0 ) {
	  sum = static_cast< G4double >( numStepBaryonHeavy );
	} else {
	  sum = static_cast< G4double >( numTrackBaryonHeavy );
	}
	break;
      }
      case 17 : {
	caseName += "# electrons";
	if ( i == 0 ) {
	  sum = static_cast< G4double >( numStepElectron );
	} else {
	  sum = static_cast< G4double >( numTrackElectron );
	}
	break;
      }
      case 18 : {
	caseName += "# gammas";
	if ( i == 0 ) {
	  sum = static_cast< G4double >( numStepGamma );
	} else {
	  sum = static_cast< G4double >( numTrackGamma );
	}
	break;
      }
      case 19 : {
	caseName += "# positrons";
	if ( i == 0 ) {
	  sum = static_cast< G4double >( numStepPositron );
	} else {
	  sum = static_cast< G4double >( numTrackPositron );
	}
	break;
      }
      case 20 : {
	caseName += "# mu-";
	if ( i == 0 ) {
	  sum = static_cast< G4double >( numStepMuMinus );
	} else {
	  sum = static_cast< G4double >( numTrackMuMinus );
	}
	break;
      }
      case 21 : {
	caseName += "# mu+";
	if ( i == 0 ) {
	  sum = static_cast< G4double >( numStepMuPlus );
	} else {
	  sum = static_cast< G4double >( numTrackMuPlus );
	}
	break;
      }
      case 22 : {
	caseName += "# tau-";
	if ( i == 0 ) {
	  sum = static_cast< G4double >( numStepTauMinus );
	} else {
	  sum = static_cast< G4double >( numTrackTauMinus );
	}
	break;
      }
      case 23 : {
	caseName += "# tau+";
	if ( i == 0 ) {
	  sum = static_cast< G4double >( numStepTauPlus );
	} else {
	  sum = static_cast< G4double >( numTrackTauPlus );
	}
	break;
      }
      case 24 : {
	caseName += "# neutrinos";
	if ( i == 0 ) {
	  sum = static_cast< G4double >( numStepNeutrino );
	} else {
	  sum = static_cast< G4double >( numTrackNeutrino );
	}
	break;
      }
      case 25 : {
	caseName += "# pi+";
	if ( i == 0 ) {
	  sum = static_cast< G4double >( numStepPiPlus );
	} else {
	  sum = static_cast< G4double >( numTrackPiPlus );
	}
	break;
      }
      case 26 : {
	caseName += "# pi0";
	if ( i == 0 ) {
	  sum = static_cast< G4double >( numStepPi0 );
	} else {
	  sum = static_cast< G4double >( numTrackPi0 );
	}
	break;
      }
      case 27 : {
	caseName += "# pi-";
	if ( i == 0 ) {
	  sum = static_cast< G4double >( numStepPiMinus );
	} else {
	  sum = static_cast< G4double >( numTrackPiMinus );
	}
	break;
      }
      case 28 : {
	caseName += "# K+";
	if ( i == 0 ) {
	  sum = static_cast< G4double >( numStepKPlus );
	} else {
	  sum = static_cast< G4double >( numTrackKPlus );
	}
	break;
      }
      case 29 : {
	caseName += "# K-neutral (K0/K0bar or K0_S/K0_L)";
	if ( i == 0 ) {
	  sum = static_cast< G4double >( numStepKNeutral );
	} else {
	  sum = static_cast< G4double >( numTrackKNeutral );
	}
	break;
      }
      case 30 : {
	caseName += "# K-";
	if ( i == 0 ) {
	  sum = static_cast< G4double >( numStepKMinus );
	} else {
	  sum = static_cast< G4double >( numTrackKMinus );
	}
	break;
      }
      case 31 : {
	caseName += "# protons";
	if ( i == 0 ) {
	  sum = static_cast< G4double >( numStepProton );
	} else {
	  sum = static_cast< G4double >( numTrackProton );
	}
	break;
      }
      case 32 : {
	caseName += "# anti-protons";
	if ( i == 0 ) {
	  sum = static_cast< G4double >( numStepAntiProton );
	} else {
	  sum = static_cast< G4double >( numTrackAntiProton );
	}
	break;
      }
      case 33 : {
	caseName += "# neutrons";
	if ( i == 0 ) {
	  sum = static_cast< G4double >( numStepNeutron );
	} else {
	  sum = static_cast< G4double >( numTrackNeutron );
	}
	break;
      }
      case 34 : {
	caseName += "# anti-neutrons";
	if ( i == 0 ) {
	  sum = static_cast< G4double >( numStepAntiNeutron );
	} else {
	  sum = static_cast< G4double >( numTrackAntiNeutron );
	}
	break;
      }
      default : 
	{
	  G4cout << "\t ***WRONG*** I should NOT be inside the  default: case" 
		 << G4endl; 
	  break;
	}
      }
      mu       = sum / n;
      mu_sigma = sqrt( sum ) / n;
      G4cout << "\t" << caseName << " = " << mu  
	     << " +/- " << mu_sigma << G4endl;
    }
  }
}

