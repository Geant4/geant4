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


//***LOOKHERE***
bool StatAccepTestAnalysis::isHistogramOn = true;
bool StatAccepTestAnalysis::is2DHistogramStepLvsEOn = false;


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
  //G4cout << " StatAccepTestAnalysis::init() : Cleaning up..." << G4endl;

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
  primaryParticleId = 0;
  beamEnergy  = 0.0;
  sumEdepAct  = 0.0;
  sumEdepAct2 = 0.0;
  sumEdepTot  = 0.0;
  sumEdepTot2 = 0.0;
  maxEdepTot  = 0.0;
  countEnergyNonConservation = 0;
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

  sumEdepAct_electron  = 0.0;
  sumEdepAct_electron2 = 0.0;
  sumEdepTot_electron  = 0.0;
  sumEdepTot_electron2 = 0.0;
  sumL_electron.clear();
  sumL_electron2.clear();
  longitudinalProfile_electron.clear();
  for ( int layer = 0; layer < numberOfReadoutLayers; layer++ ) {
    longitudinalProfile_electron.push_back( 0.0 );
    sumL_electron.push_back( 0.0 );
    sumL_electron2.push_back( 0.0 );
  }
  sumR_electron.clear();
  sumR_electron2.clear();
  transverseProfile_electron.clear();
  for ( int ir = 0; ir < numberOfRadiusBins; ir++ ) {
    transverseProfile_electron.push_back( 0.0 );
    sumR_electron.push_back( 0.0 );
    sumR_electron2.push_back( 0.0 );
  } 

  sumEdepAct_muon  = 0.0;
  sumEdepAct_muon2 = 0.0;
  sumEdepTot_muon  = 0.0;
  sumEdepTot_muon2 = 0.0;
  sumL_muon.clear();
  sumL_muon2.clear();
  longitudinalProfile_muon.clear();
  for ( int layer = 0; layer < numberOfReadoutLayers; layer++ ) {
    longitudinalProfile_muon.push_back( 0.0 );
    sumL_muon.push_back( 0.0 );
    sumL_muon2.push_back( 0.0 );
  }
  sumR_muon.clear();
  sumR_muon2.clear();
  transverseProfile_muon.clear();
  for ( int ir = 0; ir < numberOfRadiusBins; ir++ ) {
    transverseProfile_muon.push_back( 0.0 );
    sumR_muon.push_back( 0.0 );
    sumR_muon2.push_back( 0.0 );
  } 

  sumEdepAct_pion  = 0.0;
  sumEdepAct_pion2 = 0.0;
  sumEdepTot_pion  = 0.0;
  sumEdepTot_pion2 = 0.0;
  sumL_pion.clear();
  sumL_pion2.clear();
  longitudinalProfile_pion.clear();
  for ( int layer = 0; layer < numberOfReadoutLayers; layer++ ) {
    longitudinalProfile_pion.push_back( 0.0 );
    sumL_pion.push_back( 0.0 );
    sumL_pion2.push_back( 0.0 );
  }
  sumR_pion.clear();
  sumR_pion2.clear();
  transverseProfile_pion.clear();
  for ( int ir = 0; ir < numberOfRadiusBins; ir++ ) {
    transverseProfile_pion.push_back( 0.0 );
    sumR_pion.push_back( 0.0 );
    sumR_pion2.push_back( 0.0 );
  } 

  sumEdepAct_kaon  = 0.0;
  sumEdepAct_kaon2 = 0.0;
  sumEdepTot_kaon  = 0.0;
  sumEdepTot_kaon2 = 0.0;
  sumL_kaon.clear();
  sumL_kaon2.clear();
  longitudinalProfile_kaon.clear();
  for ( int layer = 0; layer < numberOfReadoutLayers; layer++ ) {
    longitudinalProfile_kaon.push_back( 0.0 );
    sumL_kaon.push_back( 0.0 );
    sumL_kaon2.push_back( 0.0 );
  }
  sumR_kaon.clear();
  sumR_kaon2.clear();
  transverseProfile_kaon.clear();
  for ( int ir = 0; ir < numberOfRadiusBins; ir++ ) {
    transverseProfile_kaon.push_back( 0.0 );
    sumR_kaon.push_back( 0.0 );
    sumR_kaon2.push_back( 0.0 );
  } 

  sumEdepAct_proton  = 0.0;
  sumEdepAct_proton2 = 0.0;
  sumEdepTot_proton  = 0.0;
  sumEdepTot_proton2 = 0.0;
  sumL_proton.clear();
  sumL_proton2.clear();
  longitudinalProfile_proton.clear();
  for ( int layer = 0; layer < numberOfReadoutLayers; layer++ ) {
    longitudinalProfile_proton.push_back( 0.0 );
    sumL_proton.push_back( 0.0 );
    sumL_proton2.push_back( 0.0 );
  }
  sumR_proton.clear();
  sumR_proton2.clear();
  transverseProfile_proton.clear();
  for ( int ir = 0; ir < numberOfRadiusBins; ir++ ) {
    transverseProfile_proton.push_back( 0.0 );
    sumR_proton.push_back( 0.0 );
    sumR_proton2.push_back( 0.0 );
  } 

  sumEdepAct_pdg0  = 0.0;
  sumEdepAct_pdg02 = 0.0;
  sumEdepTot_pdg0  = 0.0;
  sumEdepTot_pdg02 = 0.0;
  sumL_pdg0.clear();
  sumL_pdg02.clear();
  longitudinalProfile_pdg0.clear();
  for ( int layer = 0; layer < numberOfReadoutLayers; layer++ ) {
    longitudinalProfile_pdg0.push_back( 0.0 );
    sumL_pdg0.push_back( 0.0 );
    sumL_pdg02.push_back( 0.0 );
  }
  sumR_pdg0.clear();
  sumR_pdg02.clear();
  transverseProfile_pdg0.clear();
  for ( int ir = 0; ir < numberOfRadiusBins; ir++ ) {
    transverseProfile_pdg0.push_back( 0.0 );
    sumR_pdg0.push_back( 0.0 );
    sumR_pdg02.push_back( 0.0 );
  } 

  numStep = 0.0;
  numStepPositive = numStepNeutral = numStepNegative = 0.0;
  numStepPDGCodeZero = numStepPDGCodeUnrecognized = 0.0;
  numStepEM = 0.0;
  numStepEWK = 0.0;
  numStepHAD = 0.0; 
  numStepMeson = numStepBaryon = 0.0;     
  numStepMesonLight = numStepBaryonLight = 0.0;     
  numStepMesonStrange = numStepBaryonStrange = 0.0;     
  numStepMesonHeavy = numStepBaryonHeavy = 0.0;
  numStepElectron = numStepGamma = numStepPositron = 0.0;
  numStepMuMinus = numStepMuPlus = 0.0;
  numStepTauMinus = numStepTauPlus = 0.0;
  numStepNeutrino = 0.0;
  numStepPiPlus = numStepPi0 = numStepPiMinus = 0.0;
  numStepKPlus = numStepKNeutral = numStepKMinus = 0.0; 
  numStepProton = numStepAntiProton = 0.0;
  numStepNeutron = numStepAntiNeutron = 0.0;

  numStep2 = 0.0;
  numStepPositive2 = numStepNeutral2 = numStepNegative2 = 0.0;
  numStepPDGCodeZero2 = numStepPDGCodeUnrecognized2 = 0.0;
  numStepEM2 = 0.0;
  numStepEWK2 = 0.0;
  numStepHAD2 = 0.0; 
  numStepMeson2 = numStepBaryon2 = 0.0;     
  numStepMesonLight2 = numStepBaryonLight2 = 0.0;     
  numStepMesonStrange2 = numStepBaryonStrange2 = 0.0;     
  numStepMesonHeavy2 = numStepBaryonHeavy2 = 0.0;
  numStepElectron2 = numStepGamma2 = numStepPositron2 = 0.0;
  numStepMuMinus2 = numStepMuPlus2 = 0.0;
  numStepTauMinus2 = numStepTauPlus2 = 0.0;
  numStepNeutrino2 = 0.0;
  numStepPiPlus2 = numStepPi02 = numStepPiMinus2 = 0.0;
  numStepKPlus2 = numStepKNeutral2 = numStepKMinus2 = 0.0; 
  numStepProton2 = numStepAntiProton2 = 0.0;
  numStepNeutron2 = numStepAntiNeutron2 = 0.0;

  numTrack = 0.0;
  numTrackPositive = numTrackNeutral = numTrackNegative = 0.0;
  numTrackPDGCodeZero = numTrackPDGCodeUnrecognized = 0.0;
  numTrackEM = 0.0;
  numTrackEWK = 0.0;
  numTrackHAD = 0.0; 
  numTrackMeson = numTrackBaryon = 0.0;     
  numTrackMesonLight = numTrackBaryonLight = 0.0;     
  numTrackMesonStrange = numTrackBaryonStrange = 0.0;     
  numTrackMesonHeavy = numTrackBaryonHeavy = 0.0;
  numTrackElectron = numTrackGamma = numTrackPositron = 0.0;
  numTrackMuMinus = numTrackMuPlus = 0.0;
  numTrackTauMinus = numTrackTauPlus = 0.0;
  numTrackNeutrino = 0.0;
  numTrackPiPlus = numTrackPi0 = numTrackPiMinus = 0.0;
  numTrackKPlus = numTrackKNeutral = numTrackKMinus = 0.0; 
  numTrackProton = numTrackAntiProton = 0.0;
  numTrackNeutron = numTrackAntiNeutron = 0.0;

  numTrack2 = 0.0;
  numTrackPositive2 = numTrackNeutral2 = numTrackNegative2 = 0.0;
  numTrackPDGCodeZero2 = numTrackPDGCodeUnrecognized2 = 0.0;
  numTrackEM2 = 0.0;
  numTrackEWK2 = 0.0;
  numTrackHAD2 = 0.0; 
  numTrackMeson2 = numTrackBaryon2 = 0.0;     
  numTrackMesonLight2 = numTrackBaryonLight2 = 0.0;     
  numTrackMesonStrange2 = numTrackBaryonStrange2 = 0.0;     
  numTrackMesonHeavy2 = numTrackBaryonHeavy2 = 0.0;
  numTrackElectron2 = numTrackGamma2 = numTrackPositron2 = 0.0;
  numTrackMuMinus2 = numTrackMuPlus2 = 0.0;
  numTrackTauMinus2 = numTrackTauPlus2 = 0.0;
  numTrackNeutrino2 = 0.0;
  numTrackPiPlus2 = numTrackPi02 = numTrackPiMinus2 = 0.0;
  numTrackKPlus2 = numTrackKNeutral2 = numTrackKMinus2 = 0.0; 
  numTrackProton2 = numTrackAntiProton2 = 0.0;
  numTrackNeutron2 = numTrackAntiNeutron2 = 0.0;

  electronTrackLength = 0.0;
  muonTrackLength = 0.0;
  pionChargedTrackLength = 0.0;
  protonTrackLength = 0.0;
  gammaTrackLength = 0.0;
  pion0TrackLength = 0.0;

  electronTrackLength2 = 0.0;
  muonTrackLength2 = 0.0;
  pionChargedTrackLength2 = 0.0;
  protonTrackLength2 = 0.0;
  gammaTrackLength2 = 0.0;
  pion0TrackLength2 = 0.0;
  neutronTrackLength2 = 0.0;

  kinEnergyExiting = 0.0;
  kinEnergyExitingGammas = 0.0; 
  kinEnergyExitingNeutrons = 0.0;
  kinEnergyExitingNeutrinos = 0.0;
  kinEnergyExitingMuons = 0.0;
  kinEnergyExitingElectrons = 0.0;
  kinEnergyExitingOthers = 0.0;
  numExiting = 0.0;
  numExitingGammas = 0.0;
  numExitingNeutrons = 0.0;
  numExitingNeutrinos = 0.0;
  numExitingMuons = 0.0;
  numExitingElectrons = 0.0;
  numExitingOthers = 0.0;

  kinEnergyExiting2 = 0.0;
  kinEnergyExitingGammas2 = 0.0; 
  kinEnergyExitingNeutrons2 = 0.0;
  kinEnergyExitingNeutrinos2 = 0.0;
  kinEnergyExitingMuons2 = 0.0;
  kinEnergyExitingElectrons2 = 0.0;
  kinEnergyExitingOthers2 = 0.0;
  numExiting2 = 0.0;
  numExitingGammas2 = 0.0;
  numExitingNeutrons2 = 0.0;
  numExitingNeutrinos2 = 0.0;
  numExitingMuons2 = 0.0;
  numExitingElectrons2 = 0.0;
  numExitingOthers2 = 0.0;

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
      //G4cout << " Created longitudinalProfileHisto " << G4endl;
      // // // if ( ! tree->find( "60" ) ) {
      transverseProfileHisto = 
	histoFactory->createHistogram1D("60", "Transverse shower profile", 
					numberOfRadiusBins, 0.0, 1.0*numberOfRadiusBins );
      //G4cout << " Created transverseProfileHisto " << G4endl;

  }

  // Create all the other histograms.
  if ( isHistogramOn && histoFactory ) {

      // Step Energy versus step Length.
      if ( is2DHistogramStepLvsEOn ) {
	h2stepEvsL_active = 
	  histoFactory->createHistogram2D( "80", "Step: Energy(MeV) vs. Length(mm), active", 100, 0.0, 10.0, 100, 0.0, 10.0 );
	h2stepEvsL_electron_active = 
	  histoFactory->createHistogram2D( "81", "Step: Energy(MeV) vs. Length(mm), active, ELECTRON", 100, 0.0, 10.0, 100, 0.0, 10.0 );
	h2stepEvsL_muon_active = 
	  histoFactory->createHistogram2D( "82", "Step: Energy(MeV) vs. Length(mm), active,  MUON", 100, 0.0, 10.0, 100, 0.0, 10.0 );
	h2stepEvsL_pionCharged_active = 
	  histoFactory->createHistogram2D( "83", "Step: Energy(MeV) vs. Length(mm), active, PION charged", 100, 0.0, 10.0, 100, 0.0, 10.0 );
	h2stepEvsL_proton_active = 
	  histoFactory->createHistogram2D( "84", "Step: Energy(MeV) vs. Length(mm), active, PROTON", 100, 0.0, 10.0, 100, 0.0, 10.0 );
	h2stepEvsL_gamma_active = 
	  histoFactory->createHistogram2D( "85", "Step: Energy(MeV) vs. Length(mm), active, GAMMA", 100, 0.0, 10.0, 100, 0.0, 10.0 );
	h2stepEvsL_neutron_active = 
	  histoFactory->createHistogram2D( "86", "Step: Energy(MeV) vs. Length(mm), active, NEUTRON", 100, 0.0, 10.0, 100, 0.0, 10.0 );
	h2stepEvsL_absorber = 
	  histoFactory->createHistogram2D( "90", "Step: Energy(MeV) vs. Length(mm), absorber", 100, 0.0, 10.0, 100, 0.0, 10.0 );
	h2stepEvsL_electron_absorber = 
	  histoFactory->createHistogram2D( "91", "Step: Energy(MeV) vs. Length(mm), absorber, ELECTRON", 100, 0.0, 10.0, 100, 0.0, 10.0 );
	h2stepEvsL_muon_absorber = 
	  histoFactory->createHistogram2D( "92", "Step: Energy(MeV) vs. Length(mm), absorber, MUON", 100, 0.0, 10.0, 100, 0.0, 10.0 );
	h2stepEvsL_pionCharged_absorber = 
	  histoFactory->createHistogram2D( "93", "Step: Energy(MeV) vs. Length(mm), absorber, PION charged", 100, 0.0, 10.0, 100, 0.0, 10.0 );
	h2stepEvsL_proton_absorber = 
	  histoFactory->createHistogram2D( "94", "Step: Energy(MeV) vs. Length(mm), absorber, PROTON", 100, 0.0, 10.0, 100, 0.0, 10.0 );
	h2stepEvsL_gamma_absorber = 
	  histoFactory->createHistogram2D( "95", "Step: Energy(MeV) vs. Length(mm), absorber, GAMMA", 100, 0.0, 10.0, 100, 0.0, 10.0 );
	h2stepEvsL_neutron_absorber = 
	  histoFactory->createHistogram2D( "96", "Step: Energy(MeV) vs. Length(mm), absorber, NEUTRON", 100, 0.0, 10.0, 100, 0.0, 10.0 );
      }

      // Kinetic spectra of some particles.
      // We evaluate the kinetic spectra in 10 active layers,
      // uniformely spread along the calorimeter.
      // Here is the decoding of the histograms:
      //    [particle type]*1000 + [energy type]*100 + [layer number]
      // where:
      //     particle type:  1 : electrons and positrons combined
      //                     2 : pi+ pi- pi0 combined
      //                     3 : protons
      //                     4 : neutrons
      //                     5 : gammas
      //                     6 : pi+
      //                     7 : pi-
      //     energy type: 
      //        1 : Particle flux up to 100 GeV (backward going have negative energy)
      //            -100. to +100. GeV  (2000 bins) - energy in GeV
      //        2 : Particle flux up to 10 GeV (backward going have negative energy)
      //            -10. to +10. GeV  (2000 bins) - energy in GeV
      //        3 : Particle flux up to 1 GeV (backward going have negative energy)
      //            -1. to +1. GeV  (2000 bins) - energy in GeV
      //        4 : Particle flux up to 0.1 GeV (backward going have negative energy)
      //            -0.1 to +0.1 GeV  (2000 bins) - energy in GeV
      //        5 : LOG10(energy/MeV) Particle flux from 0.1 keV up to 1000 GeV 
      //            (backward going are REMOVED)
      //            -4.0 to +6.0  (1000 logarithmic bins) 
      //            energy in MeV from  0.1 keV to 1000 GeV
      // Example:
      //   5201 is the histogram of energy of gammas passing into the 
      //   second layer (index #1) binned linearly between -10 and 10 GeV.
      // We put a "9" in front of the above ID, for the corresponding
      // weighted histograms, for example:
      //   95201 is the histogram of energy of gammas passing into the
      //   second layer (index #1) binned linearly between -10 and 10 GeV,
      //   and weighted with 1/momentum of the gamma.
      // Due to space problems, we keep only the histograms with "9"
      // in front of them.
      char id[5], histotag[60];
      for ( int iLayer = 0; iLayer < 10; iLayer++ ) {
	sprintf( id, "%d", iLayer+1100 );
	sprintf( histotag, "Elec./Pos. Lin. Energy Spectrum in Active Layer %d in GeV",
		 iLayer );
	emSpectrum1[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -100.0, 100.0 );
	sprintf( id, "%d", iLayer+91100 );
	emSpectrumWeighted1[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -100.0, 100.0 );

	sprintf( id, "%d", iLayer+1200 );
	sprintf( histotag, "Elec./Pos. Lin. Energy Spectrum in Active Layer %d in GeV",
		 iLayer );
	emSpectrum2[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -10.0, 10.0 );
	sprintf( id, "%d", iLayer+91200 );
	emSpectrumWeighted2[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -10.0, 10.0 );

	sprintf( id, "%d", iLayer+1300 );
	sprintf( histotag, "Elec./Pos. Lin. Energy Spectrum in Active Layer %d in GeV",
		 iLayer );
	emSpectrum3[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -1.0, 1.0 );
	sprintf( id, "%d", iLayer+91300 );
	emSpectrumWeighted3[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -1.0, 1.0 );

	sprintf( id, "%d", iLayer+1400 );
	sprintf( histotag, "Elec./Pos. Lin. Energy Spectrum in Active Layer %d in GeV",
		 iLayer );
	emSpectrum4[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -0.1, 0.1 );
	sprintf( id, "%d", iLayer+91400 );
	emSpectrumWeighted4[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -0.1, 0.1 );

	sprintf( id, "%d", iLayer+1500 );
	sprintf( histotag, "Elec./Pos. Log10 Energy Spectrum in Active Layer %d in MeV",
		 iLayer );
	emSpectrum5[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 1000, -4.0, 6.0 );
	sprintf( id, "%d", iLayer+91500 );
	emSpectrumWeighted5[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 1000, -4.0, 6.0 );

	sprintf( id, "%d", iLayer+2100 );
	sprintf( histotag, "Pion Lin. Energy Spectrum in Active Layer %d in GeV", 
		 iLayer) ;
	pionSpectrum1[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -100.0, 100.0 );
	sprintf( id, "%d", iLayer+92100 );
	pionSpectrumWeighted1[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -100.0, 100.0 );

	sprintf( id, "%d", iLayer+2200 );
	sprintf( histotag, "Pion Lin. Energy Spectrum in Active Layer %d in GeV", 
		 iLayer );
	pionSpectrum2[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -10.0, 10.0 );
	sprintf( id, "%d", iLayer+92200 );
	pionSpectrumWeighted2[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -10.0, 10.0 );

	sprintf( id, "%d", iLayer+2300 );
	sprintf( histotag, "Pion Lin. Energy Spectrum in Active Layer %d in GeV", 
		 iLayer );
	pionSpectrum3[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -1.0, 1.0 );
	sprintf( id, "%d", iLayer+92300 );
	pionSpectrumWeighted3[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -1.0, 1.0 );

	sprintf( id, "%d", iLayer+2400 );
	sprintf( histotag, "Pion Lin. Energy Spectrum in Active Layer %d in GeV", 
		 iLayer );
	pionSpectrum4[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -0.1, 0.1 );
	sprintf( id, "%d", iLayer+92400 );
	pionSpectrumWeighted4[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -0.1, 0.1 );

	sprintf( id, "%d", iLayer+2500 );
	sprintf( histotag, "Pion Log10 Energy Spectrum in Active Layer %d in MeV", 
		 iLayer );
	pionSpectrum5[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 1000, -4.0, 6.0 );
	sprintf( id, "%d", iLayer+92500 );
	pionSpectrumWeighted5[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 1000, -4.0, 6.0 );

	sprintf( id, "%d", iLayer+3100 );
	sprintf( histotag, "Proton Lin. Energy Spectrum in Active Layer %d in GeV",
		 iLayer );
	protonSpectrum1[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -100.0, 100.0 );
	sprintf( id, "%d", iLayer+93100 );
	protonSpectrumWeighted1[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -100.0, 100.0 );

	sprintf( id, "%d", iLayer+3200 );
	sprintf( histotag, "Proton Lin. Energy Spectrum in Active Layer %d in GeV", 
		 iLayer );
	protonSpectrum2[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -10.0, 10.0 );
	sprintf( id, "%d", iLayer+93200 );
	protonSpectrumWeighted2[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -10.0, 10.0 );

	sprintf( id, "%d", iLayer+3300 );
	sprintf( histotag, "Proton Lin. Energy Spectrum in Active Layer %d in GeV", 
		 iLayer );
	protonSpectrum3[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -1.0, 1.0 );
	sprintf( id, "%d", iLayer+93300 );
	protonSpectrumWeighted3[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -1.0, 1.0 );

	sprintf( id, "%d", iLayer+3400 );
	sprintf( histotag, "Proton Lin. Energy Spectrum in Active Layer %d in GeV", 
		 iLayer );
	protonSpectrum4[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -0.1, 0.1 );
	sprintf( id, "%d", iLayer+93400 );
	protonSpectrumWeighted4[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -0.1, 0.1 );

	sprintf( id, "%d", iLayer+3500 );
	sprintf( histotag, "Proton Log10 Energy Spectrum in Active Layer %d in MeV",
		 iLayer );
	protonSpectrum5[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 1000, -4.0, 6.0 );
	sprintf( id, "%d", iLayer+93500 );
	protonSpectrumWeighted5[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 1000, -4.0, 6.0 );

	sprintf( id, "%d", iLayer+4100 );
	sprintf( histotag, "Neutron Lin. Energy Spectrum in Active Layer %d in GeV", 
		 iLayer );
	neutronSpectrum1[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -100.0, 100.0 );
	sprintf( id, "%d", iLayer+94100 );
	neutronSpectrumWeighted1[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -100.0, 100.0 );

	sprintf( id, "%d", iLayer+4200 );
	sprintf( histotag, "Neutron Lin. Energy Spectrum in Active Layer%d   in GeV",
		 iLayer );
	neutronSpectrum2[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -10.0, 10.0 );
	sprintf( id, "%d", iLayer+94200 );
	neutronSpectrumWeighted2[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -10.0, 10.0 );

	sprintf( id, "%d", iLayer+4300 );
	sprintf( histotag, "Neutron Lin. Energy Spectrum in Active Layer%d   in GeV",
		 iLayer );
	neutronSpectrum3[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -1.0, 1.0 );
	sprintf( id, "%d", iLayer+94300 );
	neutronSpectrumWeighted3[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -1.0, 1.0 );

	sprintf( id, "%d", iLayer+4400 );
	sprintf( histotag, "Neutron Lin. Energy Spectrum in Active Layer%d   in GeV",
		 iLayer );
	neutronSpectrum4[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -0.1, 0.1 );
	sprintf( id, "%d", iLayer+94400 );
	neutronSpectrumWeighted4[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -0.1, 0.1 );

	sprintf( id, "%d", iLayer+4500 );
	sprintf( histotag, "Neutron Log10 Energy Spectrum in Active Layer %d in MeV",
		 iLayer );
	neutronSpectrum5[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 1000, -4.0, 6.0 );
	sprintf( id, "%d", iLayer+94500 );
	neutronSpectrumWeighted5[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 1000, -4.0, 6.0 );

	sprintf( id, "%d", iLayer+5100 );
	sprintf( histotag, "Gamma Lin. Energy Spectrum in Active Layer%d   in GeV", 
		 iLayer );
	gammaSpectrum1[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -100.0, 100.0 );
	sprintf( id, "%d", iLayer+95100 );
	gammaSpectrumWeighted1[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -100.0, 100.0 );

	sprintf( id, "%d", iLayer+5200 );
	sprintf( histotag, "Gamma Lin. Energy Spectrum in Active Layer %d in GeV", 
		 iLayer );
	gammaSpectrum2[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -10.0, 10.0 );
	sprintf( id, "%d", iLayer+95200 );
	gammaSpectrumWeighted2[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -10.0, 10.0 );

	sprintf( id, "%d", iLayer+5300 );
	sprintf( histotag, "Gamma Lin. Energy Spectrum in Active Layer %d in GeV", 
		 iLayer );
	gammaSpectrum3[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -1.0, 1.0 );
	sprintf( id, "%d", iLayer+95300 );
	gammaSpectrumWeighted3[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -1.0, 1.0 );

	sprintf( id, "%d", iLayer+5400 );
	sprintf( histotag, "Gamma Lin. Energy Spectrum in Active Layer %d in GeV", 
		 iLayer );
	gammaSpectrum4[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -0.1, 0.1 );
	sprintf( id, "%d", iLayer+95400 );
	gammaSpectrumWeighted4[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -0.1, 0.1 );

	sprintf( id, "%d", iLayer+5500 );
	sprintf( histotag, "Gamma Log10 Energy Spectrum in Active Layer %d in MeV",
		 iLayer );
	gammaSpectrum5[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 1000, -4.0, 6.0 );
	sprintf( id, "%d", iLayer+95500 );
	gammaSpectrumWeighted5[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 1000, -4.0, 6.0 );

	sprintf( id, "%d", iLayer+6100 );
	sprintf( histotag, "PionPlus Lin. Energy Spectrum in Active Layer %d in GeV", 
		 iLayer) ;
	pionPlusSpectrum1[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -100.0, 100.0 );
	sprintf( id, "%d", iLayer+96100 );
	pionPlusSpectrumWeighted1[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -100.0, 100.0 );

	sprintf( id, "%d", iLayer+6200 );
	sprintf( histotag, "Pion Plus Lin. Energy Spectrum in Active Layer %d in GeV", 
		 iLayer );
	pionPlusSpectrum2[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -10.0, 10.0 );
	sprintf( id, "%d", iLayer+96200 );
	pionPlusSpectrumWeighted2[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -10.0, 10.0 );

	sprintf( id, "%d", iLayer+6300 );
	sprintf( histotag, "PionPlus Lin. Energy Spectrum in Active Layer %d in GeV", 
		 iLayer );
	pionPlusSpectrum3[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -1.0, 1.0 );
	sprintf( id, "%d", iLayer+96300 );
	pionPlusSpectrumWeighted3[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -1.0, 1.0 );

	sprintf( id, "%d", iLayer+6400 );
	sprintf( histotag, "PionPlus Lin. Energy Spectrum in Active Layer %d in GeV", 
		 iLayer );
	pionPlusSpectrum4[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -0.1, 0.1 );
	sprintf( id, "%d", iLayer+96400 );
	pionPlusSpectrumWeighted4[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -0.1, 0.1 );

	sprintf( id, "%d", iLayer+6500 );
	sprintf( histotag, "PionPlus Log10 Energy Spectrum in Active Layer %d in MeV", 
		 iLayer );
	pionPlusSpectrum5[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 1000, -4.0, 6.0 );
	sprintf( id, "%d", iLayer+96500 );
	pionPlusSpectrumWeighted5[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 1000, -4.0, 6.0 );

	sprintf( id, "%d", iLayer+7100 );
	sprintf( histotag, "PionMinus Lin. Energy Spectrum in Active Layer %d in GeV", 
		 iLayer) ;
	pionMinusSpectrum1[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -100.0, 100.0 );
	sprintf( id, "%d", iLayer+97100 );
	pionMinusSpectrumWeighted1[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -100.0, 100.0 );

	sprintf( id, "%d", iLayer+7200 );
	sprintf( histotag, "Pion Minus Lin. Energy Spectrum in Active Layer %d in GeV", 
		 iLayer );
	pionMinusSpectrum2[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -10.0, 10.0 );
	sprintf( id, "%d", iLayer+97200 );
	pionMinusSpectrumWeighted2[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -10.0, 10.0 );

	sprintf( id, "%d", iLayer+7300 );
	sprintf( histotag, "PionMinus Lin. Energy Spectrum in Active Layer %d in GeV", 
		 iLayer );
	pionMinusSpectrum3[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -1.0, 1.0 );
	sprintf( id, "%d", iLayer+97300 );
	pionMinusSpectrumWeighted3[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -1.0, 1.0 );

	sprintf( id, "%d", iLayer+7400 );
	sprintf( histotag, "PionMinus Lin. Energy Spectrum in Active Layer %d in GeV", 
		 iLayer );
	pionMinusSpectrum4[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -0.1, 0.1 );
	sprintf( id, "%d", iLayer+97400 );
	pionMinusSpectrumWeighted4[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 2000, -0.1, 0.1 );

	sprintf( id, "%d", iLayer+7500 );
	sprintf( histotag, "PionMinus Log10 Energy Spectrum in Active Layer %d in MeV", 
		 iLayer );
	pionMinusSpectrum5[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 1000, -4.0, 6.0 );
	sprintf( id, "%d", iLayer+97500 );
	pionMinusSpectrumWeighted5[ iLayer ] = 
	  histoFactory->createHistogram1D( id, histotag, 1000, -4.0, 6.0 );
      }
  }
}                       


void StatAccepTestAnalysis::fillSpectrum( const G4ParticleDefinition* particleDef,
					  const G4int layerNumber, 
					  const G4double kinEnergy ) {
  // We evaluate the kinetic spectra in 10 active layers, whose
  // layer number is:  0, M, 2*M, ... , 9*M-1 ,
  // where M is an integer defined as: M = number of active layers / 10 .
  // For example, if the number of layers is 25, we consider 
  // the layers: 0, 2, 4, 6, 8, 10, 12, 14, 16, 18 .
  // To get Lorentz invariant spectra, we weigh the histograms with 
  //  1/|momentum|, exept for those in logarithmic scale in the x-axis
  // (i.e. log_10(kinEnergy/MeV) ), where the weight must be
  //  1/(|momentum*Ekin|) .

  if ( ! isHistogramOn ) return;

  if ( numberOfReplicas < 10  ||
       layerNumber % (numberOfReplicas/10) != 0  ||
       layerNumber / (numberOfReplicas/10) >= 10 ) return;
  G4int sampleLayer = layerNumber / (numberOfReplicas/10); 
  G4double p = std::abs( kinEnergy ); // Momentum for the weighting...
  if ( particleDef->GetPDGMass()/MeV > 1.0e-06 ) {
    p = std::sqrt( std::abs( kinEnergy ) * 
		   ( std::abs( kinEnergy ) + 2.0 * particleDef->GetPDGMass()/MeV ) );
    //G4cout << "  m=" << particleDef->GetPDGMass()/MeV
    //       << "  Ekin=" << std::abs( kinEnergy ) 
    //	     <<  "  p=" << p << " MeV" << G4endl;       //***DEBUG***
  }
  if ( particleDef == G4Electron::ElectronDefinition() || 
       particleDef == G4Positron::PositronDefinition() ) {
    emSpectrum1[ sampleLayer ]->fill( kinEnergy/GeV );
    emSpectrum2[ sampleLayer ]->fill( kinEnergy/GeV );
    emSpectrum3[ sampleLayer ]->fill( kinEnergy/GeV );
    emSpectrum4[ sampleLayer ]->fill( kinEnergy/GeV );
    if ( kinEnergy > 0.0 ) {
      emSpectrum5[ sampleLayer ]->fill( std::log10( kinEnergy/MeV ) );
    }
    if ( p > 0.0 ) {
      emSpectrumWeighted1[ sampleLayer ]->fill( kinEnergy/GeV, 1.0/p );
      emSpectrumWeighted2[ sampleLayer ]->fill( kinEnergy/GeV, 1.0/p );
      emSpectrumWeighted3[ sampleLayer ]->fill( kinEnergy/GeV, 1.0/p );
      emSpectrumWeighted4[ sampleLayer ]->fill( kinEnergy/GeV, 1.0/p );
      if ( kinEnergy > 0.0 ) {
	emSpectrumWeighted5[ sampleLayer ]->fill( std::log10( kinEnergy/MeV ), 
						  1.0/(p*std::abs(kinEnergy)) );
      }
    }
  } else if ( particleDef == G4Gamma::GammaDefinition() ) {
    gammaSpectrum1[ sampleLayer ]->fill( kinEnergy/GeV );
    gammaSpectrum2[ sampleLayer ]->fill( kinEnergy/GeV );
    gammaSpectrum3[ sampleLayer ]->fill( kinEnergy/GeV );
    gammaSpectrum4[ sampleLayer ]->fill( kinEnergy/GeV );
    if ( kinEnergy > 0.0 ) {
      gammaSpectrum5[ sampleLayer ]->fill( std::log10( kinEnergy/MeV ) );  
    }
    if ( p > 0.0 ) {
      gammaSpectrumWeighted1[ sampleLayer ]->fill( kinEnergy/GeV, 1.0/p );
      gammaSpectrumWeighted2[ sampleLayer ]->fill( kinEnergy/GeV, 1.0/p );
      gammaSpectrumWeighted3[ sampleLayer ]->fill( kinEnergy/GeV, 1.0/p );
      gammaSpectrumWeighted4[ sampleLayer ]->fill( kinEnergy/GeV, 1.0/p );
      if ( kinEnergy > 0.0 ) {
	gammaSpectrumWeighted5[ sampleLayer ]->fill( std::log10( kinEnergy/MeV ), 
						     1.0/(p*std::abs(kinEnergy)) );
      }
    }
  } else if ( particleDef == G4PionPlus::PionPlusDefinition() ||
	      particleDef == G4PionMinus::PionMinusDefinition() ||
	      particleDef == G4PionZero::PionZeroDefinition() ) {
    pionSpectrum1[ sampleLayer ]->fill( kinEnergy/GeV );
    pionSpectrum2[ sampleLayer ]->fill( kinEnergy/GeV );
    pionSpectrum3[ sampleLayer ]->fill( kinEnergy/GeV );
    pionSpectrum4[ sampleLayer ]->fill( kinEnergy/GeV );
    if ( kinEnergy > 0.0 ) {
      pionSpectrum5[ sampleLayer ]->fill( std::log10( kinEnergy/MeV ) );
    }
    if ( p > 0.0 ) {
      pionSpectrumWeighted1[ sampleLayer ]->fill( kinEnergy/GeV, 1.0/p );
      pionSpectrumWeighted2[ sampleLayer ]->fill( kinEnergy/GeV, 1.0/p );
      pionSpectrumWeighted3[ sampleLayer ]->fill( kinEnergy/GeV, 1.0/p );
      pionSpectrumWeighted4[ sampleLayer ]->fill( kinEnergy/GeV, 1.0/p );
      if ( kinEnergy > 0.0 ) {
	pionSpectrumWeighted5[ sampleLayer ]->fill( std::log10( kinEnergy/MeV ),
						    1.0/(p*std::abs(kinEnergy)) );
      }
    }
    // Now separate between pi+ and pi-, because they have, 
    // in general, different spectra.
    if ( particleDef == G4PionPlus::PionPlusDefinition() ) { 
      pionPlusSpectrum1[ sampleLayer ]->fill( kinEnergy/GeV );
      pionPlusSpectrum2[ sampleLayer ]->fill( kinEnergy/GeV );
      pionPlusSpectrum3[ sampleLayer ]->fill( kinEnergy/GeV );
      pionPlusSpectrum4[ sampleLayer ]->fill( kinEnergy/GeV );
      if ( kinEnergy > 0.0 ) {
	pionPlusSpectrum5[ sampleLayer ]->fill( std::log10( kinEnergy/MeV ) );
      }
      if ( p > 0.0 ) {
	pionPlusSpectrumWeighted1[ sampleLayer ]->fill( kinEnergy/GeV, 1.0/p );
	pionPlusSpectrumWeighted2[ sampleLayer ]->fill( kinEnergy/GeV, 1.0/p );
	pionPlusSpectrumWeighted3[ sampleLayer ]->fill( kinEnergy/GeV, 1.0/p );
	pionPlusSpectrumWeighted4[ sampleLayer ]->fill( kinEnergy/GeV, 1.0/p );
	if ( kinEnergy > 0.0 ) {
	  pionPlusSpectrumWeighted5[ sampleLayer ]->fill( std::log10( kinEnergy/MeV ),
							  1.0/(p*std::abs(kinEnergy)) );
	}
      }
    }
    if ( particleDef == G4PionMinus::PionMinusDefinition() ) { 
      pionMinusSpectrum1[ sampleLayer ]->fill( kinEnergy/GeV );
      pionMinusSpectrum2[ sampleLayer ]->fill( kinEnergy/GeV );
      pionMinusSpectrum3[ sampleLayer ]->fill( kinEnergy/GeV );
      pionMinusSpectrum4[ sampleLayer ]->fill( kinEnergy/GeV );
      if ( kinEnergy > 0.0 ) {
	//pionMinusSpectrum5[ sampleLayer ]->fill( std::log10( kinEnergy/MeV ) );
      }
      if ( p > 0.0 ) {
	pionMinusSpectrumWeighted1[ sampleLayer ]->fill( kinEnergy/GeV, 1.0/p );
	pionMinusSpectrumWeighted2[ sampleLayer ]->fill( kinEnergy/GeV, 1.0/p );
	pionMinusSpectrumWeighted3[ sampleLayer ]->fill( kinEnergy/GeV, 1.0/p );
	pionMinusSpectrumWeighted4[ sampleLayer ]->fill( kinEnergy/GeV, 1.0/p );
	if ( kinEnergy > 0.0 ) {
	  pionMinusSpectrumWeighted5[ sampleLayer ]->fill( std::log10( kinEnergy/MeV ),
							   1.0/(p*std::abs(kinEnergy)) );
	}
      }
    }
  } else if ( particleDef == G4Proton::ProtonDefinition() ) {
    protonSpectrum1[ sampleLayer ]->fill( kinEnergy/GeV );
    protonSpectrum2[ sampleLayer ]->fill( kinEnergy/GeV );
    protonSpectrum3[ sampleLayer ]->fill( kinEnergy/GeV );
    protonSpectrum4[ sampleLayer ]->fill( kinEnergy/GeV );
    if ( kinEnergy > 0.0 ) {
      protonSpectrum5[ sampleLayer ]->fill( std::log10( kinEnergy/MeV ) );  
    }
    if ( p > 0.0 ) {
      protonSpectrumWeighted1[ sampleLayer ]->fill( kinEnergy/GeV, 1.0/p );
      protonSpectrumWeighted2[ sampleLayer ]->fill( kinEnergy/GeV, 1.0/p );
      protonSpectrumWeighted3[ sampleLayer ]->fill( kinEnergy/GeV, 1.0/p );
      protonSpectrumWeighted4[ sampleLayer ]->fill( kinEnergy/GeV, 1.0/p );
      if ( kinEnergy > 0.0 ) {
	protonSpectrumWeighted5[ sampleLayer ]->fill( std::log10( kinEnergy/MeV ),
						      1.0/(p*std::abs(kinEnergy)) );
      }
    }
  } else if ( particleDef == G4Neutron::NeutronDefinition() ) {
    neutronSpectrum1[ sampleLayer ]->fill( kinEnergy/GeV );
    neutronSpectrum2[ sampleLayer ]->fill( kinEnergy/GeV );
    neutronSpectrum3[ sampleLayer ]->fill( kinEnergy/GeV );
    neutronSpectrum4[ sampleLayer ]->fill( kinEnergy/GeV );
    if ( kinEnergy > 0.0 ) {
      neutronSpectrum5[ sampleLayer ]->fill( std::log10( kinEnergy/MeV ) );
    }
    if ( p > 0.0 ) {
      neutronSpectrumWeighted1[ sampleLayer ]->fill( kinEnergy/GeV, 1.0/p );
      neutronSpectrumWeighted2[ sampleLayer ]->fill( kinEnergy/GeV, 1.0/p );
      neutronSpectrumWeighted3[ sampleLayer ]->fill( kinEnergy/GeV, 1.0/p );
      neutronSpectrumWeighted4[ sampleLayer ]->fill( kinEnergy/GeV, 1.0/p );
      if ( kinEnergy > 0.0 ) {
	neutronSpectrumWeighted5[ sampleLayer ]->fill( std::log10( kinEnergy/MeV ),
						       1.0/(p*std::abs(kinEnergy)) );
      }
    }
  }
}


void StatAccepTestAnalysis::fillNtuple( float incidentParticleId, 
					float incidentParticleEnergy, 
					float totalEnergyDepositedInActiveLayers,
					float totalEnergyDepositedInCalorimeter ) {
  primaryParticleId = static_cast< G4int >( incidentParticleId );
  beamEnergy = incidentParticleEnergy;
  if ( totalEnergyDepositedInCalorimeter - beamEnergy > 0.001*MeV ) {
    G4cout << "\t ***ENERGY-NON-CONSERVATION*** " 
	   << totalEnergyDepositedInCalorimeter << " MeV" << G4endl;
    countEnergyNonConservation++;
  }
  if ( totalEnergyDepositedInCalorimeter > maxEdepTot ) {
    maxEdepTot = totalEnergyDepositedInCalorimeter;
  }
  if (tuple) {
    //G4cout << " StatAccepTestAnalysis::fillNtuple : DEBUG Info " << G4endl
    //       << "\t incidentParticleId = " << incidentParticleId << G4endl
    //       << "\t incidentParticleEnergy = " << incidentParticleEnergy << G4endl
    //       << "\t totalEnergyDepositedInActiveLayers = " 
    //       << totalEnergyDepositedInActiveLayers << G4endl
    //       << "\t totalEnergyDepositedInCalorimeter = " 
    //       << totalEnergyDepositedInCalorimeter << G4endl;       // ***DEBUG***

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
    //G4cout << " StatAccepTestAnalysis::fillNtuple : DEBUG Info " << G4endl
    //       << "\t Longitudinal profile: layer = " << layer
    //       << "   energy = " << longitudinalProfile[ layer ] / MeV 
    //       << " MeV " << G4endl;                                 //***DEBUG***
    longitudinalProfile[ layer ] = 0.0;
  }
  for ( int ir = 0; ir < numberOfRadiusBins; ir++ ) {
    //G4cout << " StatAccepTestAnalysis::fillNtuple : DEBUG Info " << G4endl
    //       << "\t Transverse profile: iBinRadius = " << ir / mm
    //       << " mm   energy = " << transverseProfile[ ir ] / MeV 
    //       << " MeV " << G4endl;                                 //***DEBUG***
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
      //G4cout << " tree commit ,  at event=" << numberOfEvents-1 
      //      << G4endl; //***DEBUG***
    }
  }

  // Store information of the visible energy, for later computing
  // of the energy resolution.
  vecEvis.push_back( totalEnergyDepositedInActiveLayers );

}


void StatAccepTestAnalysis::
fillShowerProfile( G4int replica, const G4double radius, 
		   const G4double edep, const G4int particlePDG ) {

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

  // Consider now the separate contribution due to the following particles:
  //    -  electron (e-  and  e+    together)
  //    -  muons    (mu- and  mu+   together)
  //    -  pions    (pi- and  pi+   together)
  //    -  kaons    (k-  and  k+    together)
  //    -  protons  (p   and  pbar  together)
  //    -  pdg0     (all particles with PDG code = 0)
  if ( particlePDG == G4Electron::ElectronDefinition()->GetPDGEncoding()  ||
       particlePDG == G4Positron::PositronDefinition()->GetPDGEncoding() ) {
    sumEdepAct_electron += edep;
    sumEdepTot_electron += edep;
    longitudinalProfile_electron[ readoutLayer ] += edep;
    transverseProfile_electron[ iBinRadius ] += edep;
  } else if ( particlePDG == G4MuonMinus::MuonMinusDefinition()->GetPDGEncoding()  ||
	      particlePDG == G4MuonPlus::MuonPlusDefinition()->GetPDGEncoding() ) {
    sumEdepAct_muon += edep;
    sumEdepTot_muon += edep;
    longitudinalProfile_muon[ readoutLayer ] += edep;
    transverseProfile_muon[ iBinRadius ] += edep;
  } else if ( particlePDG == G4PionPlus::PionPlusDefinition()->GetPDGEncoding()  ||
	      particlePDG == G4PionMinus::PionMinusDefinition()->GetPDGEncoding() ) {
    sumEdepAct_pion += edep;
    sumEdepTot_pion += edep;
    longitudinalProfile_pion[ readoutLayer ] += edep;
    transverseProfile_pion[ iBinRadius ] += edep;
  } else if ( particlePDG == G4KaonMinus::KaonMinusDefinition()->GetPDGEncoding()  ||
	      particlePDG == G4KaonPlus::KaonPlusDefinition()->GetPDGEncoding() ) {
    sumEdepAct_kaon += edep;
    sumEdepTot_kaon += edep;
    longitudinalProfile_kaon[ readoutLayer ] += edep;
    transverseProfile_kaon[ iBinRadius ] += edep;
  } else if ( particlePDG == G4Proton::ProtonDefinition()->GetPDGEncoding()  ||
	      particlePDG == G4AntiProton::AntiProtonDefinition()->GetPDGEncoding() ) {
    sumEdepAct_proton += edep;
    sumEdepTot_proton += edep;
    longitudinalProfile_proton[ readoutLayer ] += edep;
    transverseProfile_proton[ iBinRadius ] += edep;
  } else if ( particlePDG == 0 ) {
    sumEdepAct_pdg0 += edep;
    sumEdepTot_pdg0 += edep;
    longitudinalProfile_pdg0[ readoutLayer ] += edep;
    transverseProfile_pdg0[ iBinRadius ] += edep;
  }

}


void StatAccepTestAnalysis::infoStep( const G4Step* aStep ) {
  classifyParticle( false , aStep->GetTrack()->GetDefinition() );

  // 2D plots on Step Energy vs. Step Length.
  if ( isHistogramOn && is2DHistogramStepLvsEOn ) {
    if ( aStep->GetTrack()->GetVolume()->GetName() == "physiActive" ) {
      h2stepEvsL_active->fill( aStep->GetStepLength() / mm, 
			       aStep->GetTotalEnergyDeposit() / MeV );
    } else if ( aStep->GetTrack()->GetVolume()->GetName() == "physiAbsorber" ) {
      h2stepEvsL_absorber->fill( aStep->GetStepLength() / mm, 
				 aStep->GetTotalEnergyDeposit() / MeV );
    }
    if ( aStep->GetTrack()->GetDefinition() == G4Electron::ElectronDefinition() || 
	 aStep->GetTrack()->GetDefinition() == G4Positron::PositronDefinition() ) {
      if ( aStep->GetTrack()->GetVolume()->GetName() == "physiActive" ) {
	h2stepEvsL_electron_active->fill( aStep->GetStepLength() / mm, 
					  aStep->GetTotalEnergyDeposit() / MeV );
      } else if ( aStep->GetTrack()->GetVolume()->GetName() == "physiAbsorber" ) {
	h2stepEvsL_electron_absorber->fill( aStep->GetStepLength() / mm,
					    aStep->GetTotalEnergyDeposit() / MeV );
      }
    } else if ( aStep->GetTrack()->GetDefinition() == G4Gamma::GammaDefinition() ) {
      if ( aStep->GetTrack()->GetVolume()->GetName() == "physiActive" ) {
	h2stepEvsL_gamma_active->fill( aStep->GetStepLength() / mm, 
				       aStep->GetTotalEnergyDeposit() / MeV );
      } else if ( aStep->GetTrack()->GetVolume()->GetName() == "physiAbsorber" ) {
	h2stepEvsL_gamma_absorber->fill( aStep->GetStepLength() / mm, 
					 aStep->GetTotalEnergyDeposit() / MeV );
      }
    } else if ( aStep->GetTrack()->GetDefinition() == 
		G4MuonMinus::MuonMinusDefinition()  ||  
		aStep->GetTrack()->GetDefinition() == 
		G4MuonPlus::MuonPlusDefinition() ) {
      if ( aStep->GetTrack()->GetVolume()->GetName() == "physiActive" ) {
	h2stepEvsL_muon_active->fill( aStep->GetStepLength() / mm, 
				      aStep->GetTotalEnergyDeposit() / MeV );
      } else if ( aStep->GetTrack()->GetVolume()->GetName() == "physiAbsorber" ) {
	h2stepEvsL_muon_absorber->fill( aStep->GetStepLength() / mm, 
					aStep->GetTotalEnergyDeposit() / MeV );
      }
    } else if ( aStep->GetTrack()->GetDefinition() == 
		G4PionPlus::PionPlusDefinition()  ||
		aStep->GetTrack()->GetDefinition() == 
		G4PionMinus::PionMinusDefinition() ) {
      if ( aStep->GetTrack()->GetVolume()->GetName() == "physiActive" ) {
	h2stepEvsL_pionCharged_active->fill( aStep->GetStepLength() / mm, 
					     aStep->GetTotalEnergyDeposit() / MeV );
      } else if ( aStep->GetTrack()->GetVolume()->GetName() == "physiAbsorber" ) {
	h2stepEvsL_pionCharged_absorber->fill( aStep->GetStepLength() / mm, 
					       aStep->GetTotalEnergyDeposit() / MeV );
      }
    } else if ( aStep->GetTrack()->GetDefinition() == G4Proton::ProtonDefinition() ) {
      if ( aStep->GetTrack()->GetVolume()->GetName() == "physiActive" ) {
	h2stepEvsL_proton_active->fill( aStep->GetStepLength() / mm, 
					aStep->GetTotalEnergyDeposit() / MeV );
      } else if ( aStep->GetTrack()->GetVolume()->GetName() == "physiAbsorber" ) {
	h2stepEvsL_proton_absorber->fill( aStep->GetStepLength() / mm,
					  aStep->GetTotalEnergyDeposit() / MeV );
      }
    } else if ( aStep->GetTrack()->GetDefinition() == G4Neutron::NeutronDefinition() ) {
      if ( aStep->GetTrack()->GetVolume()->GetName() == "physiActive" ) {
	h2stepEvsL_neutron_active->fill( aStep->GetStepLength() / mm,
					 aStep->GetTotalEnergyDeposit() / MeV );
      } else if ( aStep->GetTrack()->GetVolume()->GetName() == "physiAbsorber" ) {
	h2stepEvsL_neutron_absorber->fill( aStep->GetStepLength() / mm,
					   aStep->GetTotalEnergyDeposit() / MeV );
      }    
    }
  }

  // Update the information on the energy deposition in the absorber
  // for the following particles (in the case of the energy deposition
  // in active layer, it has already been considered in the method
  // StatAccepTestAnalysis::fillShowerProfile):
  //    -  electron (e-  and  e+    together)
  //    -  muons    (mu- and  mu+   together)
  //    -  pions    (pi- and  pi+   together)
  //    -  kaons    (k-  and  k+    together)
  //    -  protons  (p   and  pbar  together)
  //    -  pdg0     (all particles with PDG code = 0)
  if ( aStep->GetTrack()->GetVolume()->GetName() == "physiAbsorber" ) {
    G4double edep = aStep->GetTotalEnergyDeposit() * aStep->GetTrack()->GetWeight();
    if ( aStep->GetTrack()->GetDefinition() == 
	 G4Electron::ElectronDefinition()  ||
	 aStep->GetTrack()->GetDefinition() == 
	 G4Positron::PositronDefinition() ) {
      sumEdepTot_electron += edep;
    } else if ( aStep->GetTrack()->GetDefinition() == 
		G4MuonMinus::MuonMinusDefinition()  ||
		aStep->GetTrack()->GetDefinition() == 
		G4MuonPlus::MuonPlusDefinition() ) {
      sumEdepTot_muon += edep;
    } else if ( aStep->GetTrack()->GetDefinition() == 
		G4PionPlus::PionPlusDefinition()  ||
		aStep->GetTrack()->GetDefinition() == 
		G4PionMinus::PionMinusDefinition() ) {
      sumEdepTot_pion += edep;
    } else if ( aStep->GetTrack()->GetDefinition() == 
		G4KaonMinus::KaonMinusDefinition()  ||
		aStep->GetTrack()->GetDefinition() == 
		G4KaonPlus::KaonPlusDefinition() ) {
      sumEdepTot_kaon += edep;
    } else if ( aStep->GetTrack()->GetDefinition() == 
		G4Proton::ProtonDefinition()  ||
		aStep->GetTrack()->GetDefinition() == 
		G4AntiProton::AntiProtonDefinition() ) {
      sumEdepTot_proton += edep;
    } else if ( aStep->GetTrack()->GetDefinition()->GetPDGEncoding() == 0 ) {
      sumEdepTot_pdg0 += edep;
    }
  }

}


void StatAccepTestAnalysis::infoTrack( const G4Track* aTrack ) {

  if ( aTrack->GetTrackStatus() == fStopAndKill ) {
    //G4cout << "\t --- Info Track when fStopAndKill --- " << G4endl
    //	     << "\t TrackID = " << aTrack->GetTrackID() 
    //       << "\t Name = " << aTrack->GetDefinition()->GetParticleName() << G4endl
    //       << "\t Volume = " << aTrack->GetVolume()->GetName() 
    //       << "\t Material = " << aTrack->GetMaterial()->GetName() << G4endl
    //       << "\t Vertex (origin) Volume = " 
    //       << aTrack->GetLogicalVolumeAtVertex()->GetName() << G4endl
    //       << "\t Ekin = " << aTrack->GetKineticEnergy() << " MeV "
    //       << "\t Length track = " << aTrack->GetTrackLength() 
    //       << G4endl;  //***DEBUG*** 
    //if ( ! ( aTrack->GetVolume()->GetName() == "expHall" ||
    //	     aTrack->GetVolume()->GetName() == "physiAbsorber" ||
    //	     aTrack->GetVolume()->GetName() == "physiActive" ) ) {
    //  G4cout << " ***STRANGE VOLUME *** : " << aTrack->GetVolume()->GetName() << G4endl;
    //}
  
    // To avoid bias in the track length due to the big world volume
    // (which can affect significantly the track length of neutrons)
    // we consider only those tracks that are fully contained inside
    // the calorimeter, i.e. created and terminated inside it.
    if ( aTrack->GetLogicalVolumeAtVertex()->GetName() != "expHall"  &&
	 aTrack->GetVolume()->GetName() != "expHall" ) {
      G4double trackLength = aTrack->GetTrackLength();
      if ( aTrack->GetDefinition() == G4Electron::ElectronDefinition() || 
	   aTrack->GetDefinition() == G4Positron::PositronDefinition() ) {
	electronTrackLength += trackLength;
	electronTrackLength2 += trackLength * trackLength;
      } else if ( aTrack->GetDefinition() == G4Gamma::GammaDefinition() ) {
	gammaTrackLength += trackLength;
	gammaTrackLength2 += trackLength * trackLength;
      } else if ( aTrack->GetDefinition() == G4MuonMinus::MuonMinusDefinition() ||  
		  aTrack->GetDefinition() == G4MuonPlus::MuonPlusDefinition() ) {
	muonTrackLength += trackLength;
	muonTrackLength2 += trackLength * trackLength;
      } else if ( aTrack->GetDefinition() == G4PionPlus::PionPlusDefinition() ||
		  aTrack->GetDefinition() == G4PionMinus::PionMinusDefinition() ) {
	pionChargedTrackLength += trackLength;
	pionChargedTrackLength2 += trackLength * trackLength;
      } else if ( aTrack->GetDefinition() == G4PionZero::PionZeroDefinition() ) {
	pion0TrackLength += trackLength;
	pion0TrackLength2 += trackLength * trackLength;
      } else if ( aTrack->GetDefinition() == G4Proton::ProtonDefinition() ) {
	protonTrackLength += trackLength;
	protonTrackLength2 += trackLength * trackLength;
      } else if ( aTrack->GetDefinition() == G4Neutron::NeutronDefinition() ) {
	neutronTrackLength += trackLength;
	neutronTrackLength2 += trackLength * trackLength;
      }
    } else if ( aTrack->GetVolume()->GetName() == "expHall" ) {
      kinEnergyExiting += aTrack->GetKineticEnergy();
      kinEnergyExiting2 += aTrack->GetKineticEnergy() * aTrack->GetKineticEnergy();
      numExiting++;
      //G4cout << " Exiting particle: " 
      //       << aTrack->GetDefinition()->GetParticleName() 
      //       << "\t" << kinEnergyExiting << " MeV" << G4endl
      //       << "\t Production: " 
      //       << aTrack->GetCreatorProcess()->GetProcessName() << G4endl
      //       << "\t" << aTrack->GetVertexPosition() << " mm"
      //       << "\t" << aTrack->GetLogicalVolumeAtVertex()->GetName() << G4endl
      //       << "\t" << aTrack->GetVertexKineticEnergy() << " MeV"
      //       << "\t vertex Ekin = " << aTrack->GetVertexKineticEnergy() 
      //       << " MeV" << G4endl;  //***DEBUG***
      if ( aTrack->GetDefinition() == G4Gamma::GammaDefinition() ) {
	kinEnergyExitingGammas += aTrack->GetKineticEnergy();
	numExitingGammas++;
      } else if ( aTrack->GetDefinition() == G4Neutron::NeutronDefinition() ) {
	kinEnergyExitingNeutrons += aTrack->GetKineticEnergy();
	numExitingNeutrons++;
      } else if ( aTrack->GetDefinition() == G4NeutrinoE::NeutrinoEDefinition()  ||
		  aTrack->GetDefinition() == G4NeutrinoMu::NeutrinoMuDefinition()  ||
		  aTrack->GetDefinition() == G4NeutrinoTau::NeutrinoTauDefinition()  ||
		  aTrack->GetDefinition() == G4AntiNeutrinoE::AntiNeutrinoEDefinition()  ||
		  aTrack->GetDefinition() == G4AntiNeutrinoMu::AntiNeutrinoMuDefinition()  ||
		  aTrack->GetDefinition() == G4AntiNeutrinoTau::AntiNeutrinoTauDefinition() ) {
	kinEnergyExitingNeutrinos += aTrack->GetKineticEnergy();
	numExitingNeutrinos++;
      } else if ( aTrack->GetDefinition() == G4MuonMinus::MuonMinusDefinition() ||  
		  aTrack->GetDefinition() == G4MuonPlus::MuonPlusDefinition() ) {
	kinEnergyExitingMuons += aTrack->GetKineticEnergy();
	numExitingMuons++;        
      } else if ( aTrack->GetDefinition() == G4Electron::ElectronDefinition() ||  
		  aTrack->GetDefinition() == G4Positron::PositronDefinition() ) {
	kinEnergyExitingElectrons += aTrack->GetKineticEnergy();
	numExitingElectrons++;        
      } else {
	kinEnergyExitingOthers += aTrack->GetKineticEnergy();
	numExitingOthers++;
      }
    }
  } else {
    classifyParticle( true , aTrack->GetDefinition() );
  }
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



void StatAccepTestAnalysis::endOfEvent() {
  // This method is useful to update the "squared" event variables
  // which are used at the end of the Run to compute the statistical
  // uncertainties of various quantities.
  // Notice that only the quantities that are meaningful on an event
  // by event basis, like the number of steps, the number of tracks,
  // the exiting kinetic energy, the number of exiting particles;
  // this does not apply to the track length.

  static G4double numStep_previous = 0.0;
  static G4double numStepPositive_previous = 0.0;
  static G4double numStepNeutral_previous = 0.0;
  static G4double numStepNegative_previous = 0.0;
  static G4double numStepPDGCodeZero_previous = 0.0;
  static G4double numStepPDGCodeUnrecognized_previous = 0.0;
  static G4double numStepEM_previous = 0.0;
  static G4double numStepEWK_previous = 0.0;
  static G4double numStepHAD_previous = 0.0;
  static G4double numStepMeson_previous = 0.0;
  static G4double numStepBaryon_previous = 0.0;
  static G4double numStepMesonLight_previous = 0.0;
  static G4double numStepBaryonLight_previous = 0.0;
  static G4double numStepMesonStrange_previous = 0.0;
  static G4double numStepBaryonStrange_previous = 0.0;
  static G4double numStepMesonHeavy_previous = 0.0;
  static G4double numStepBaryonHeavy_previous = 0.0;
  static G4double numStepElectron_previous = 0.0;
  static G4double numStepGamma_previous = 0.0;
  static G4double numStepPositron_previous = 0.0;
  static G4double numStepMuMinus_previous = 0.0;
  static G4double numStepMuPlus_previous = 0.0;
  static G4double numStepTauMinus_previous = 0.0;
  static G4double numStepTauPlus_previous = 0.0;
  static G4double numStepNeutrino_previous = 0.0;
  static G4double numStepPiPlus_previous = 0.0;
  static G4double numStepPi0_previous = 0.0;
  static G4double numStepPiMinus_previous = 0.0;
  static G4double numStepKPlus_previous = 0.0;
  static G4double numStepKNeutral_previous = 0.0;
  static G4double numStepKMinus_previous = 0.0;
  static G4double numStepProton_previous = 0.0;
  static G4double numStepAntiProton_previous = 0.0;
  static G4double numStepNeutron_previous = 0.0;
  static G4double numStepAntiNeutron_previous = 0.0;

  numStep2 +=
    ( numStep - numStep_previous ) *
    ( numStep - numStep_previous );
  numStepPositive2 +=
    ( numStepPositive - numStepPositive_previous ) *
    ( numStepPositive - numStepPositive_previous );
  numStepNeutral2 +=
    ( numStepNeutral - numStepNeutral_previous ) *
    ( numStepNeutral - numStepNeutral_previous );
  numStepNegative2 +=
    ( numStepNegative - numStepNegative_previous ) *
    ( numStepNegative - numStepNegative_previous );
  numStepPDGCodeZero2 +=
    ( numStepPDGCodeZero - numStepPDGCodeZero_previous ) *
    ( numStepPDGCodeZero - numStepPDGCodeZero_previous );
  numStepPDGCodeUnrecognized2 +=
    ( numStepPDGCodeUnrecognized - numStepPDGCodeUnrecognized_previous ) *
    ( numStepPDGCodeUnrecognized - numStepPDGCodeUnrecognized_previous );
  numStepEM2 +=
    ( numStepEM - numStepEM_previous ) *
    ( numStepEM - numStepEM_previous );
  numStepEWK2 +=
    ( numStepEWK - numStepEWK_previous ) *
    ( numStepEWK - numStepEWK_previous );
  numStepHAD2 +=
    ( numStepHAD - numStepHAD_previous ) *
    ( numStepHAD - numStepHAD_previous );
  numStepMeson2 +=
    ( numStepMeson - numStepMeson_previous ) *
    ( numStepMeson - numStepMeson_previous );
  numStepBaryon2 +=
    ( numStepBaryon - numStepBaryon_previous ) *
    ( numStepBaryon - numStepBaryon_previous );
  numStepMesonLight2 +=
    ( numStepMesonLight - numStepMesonLight_previous ) *
    ( numStepMesonLight - numStepMesonLight_previous );
  numStepBaryonLight2 +=
    ( numStepBaryonLight - numStepBaryonLight_previous ) *
    ( numStepBaryonLight - numStepBaryonLight_previous );
  numStepMesonStrange2 +=
    ( numStepMesonStrange - numStepMesonStrange_previous ) *
    ( numStepMesonStrange - numStepMesonStrange_previous );
  numStepBaryonStrange2 +=
    ( numStepBaryonStrange - numStepBaryonStrange_previous ) *
    ( numStepBaryonStrange - numStepBaryonStrange_previous );
  numStepMesonHeavy2 +=
    ( numStepMesonHeavy - numStepMesonHeavy_previous ) *
    ( numStepMesonHeavy - numStepMesonHeavy_previous );
  numStepBaryonHeavy2 +=
    ( numStepBaryonHeavy - numStepBaryonHeavy_previous ) *
    ( numStepBaryonHeavy - numStepBaryonHeavy_previous );
  numStepElectron2 +=
    ( numStepElectron - numStepElectron_previous ) *
    ( numStepElectron - numStepElectron_previous );
  numStepGamma2 +=
    ( numStepGamma - numStepGamma_previous ) *
    ( numStepGamma - numStepGamma_previous );
  numStepPositron2 +=
    ( numStepPositron - numStepPositron_previous ) *
    ( numStepPositron - numStepPositron_previous );
  numStepMuMinus2 +=
    ( numStepMuMinus - numStepMuMinus_previous ) *
    ( numStepMuMinus - numStepMuMinus_previous );
  numStepMuPlus2 +=
    ( numStepMuPlus - numStepMuPlus_previous ) *
    ( numStepMuPlus - numStepMuPlus_previous );
  numStepTauMinus2 +=
    ( numStepTauMinus - numStepTauMinus_previous ) *
    ( numStepTauMinus - numStepTauMinus_previous );
  numStepTauPlus2 +=
    ( numStepTauPlus - numStepTauPlus_previous ) *
    ( numStepTauPlus - numStepTauPlus_previous );
  numStepNeutrino2 +=
    ( numStepNeutrino - numStepNeutrino_previous ) *
    ( numStepNeutrino - numStepNeutrino_previous );
  numStepPiPlus2 +=
    ( numStepPiPlus - numStepPiPlus_previous ) *
    ( numStepPiPlus - numStepPiPlus_previous );
  numStepPi02 +=
    ( numStepPi0 - numStepPi0_previous ) *
    ( numStepPi0 - numStepPi0_previous );
  numStepPiMinus2 +=
    ( numStepPiMinus - numStepPiMinus_previous ) *
    ( numStepPiMinus - numStepPiMinus_previous );
  numStepKPlus2 +=
    ( numStepKPlus - numStepKPlus_previous ) *
    ( numStepKPlus - numStepKPlus_previous );
  numStepKNeutral2 +=
    ( numStepKNeutral - numStepKNeutral_previous ) *
    ( numStepKNeutral - numStepKNeutral_previous );
  numStepKMinus2 +=
    ( numStepKMinus - numStepKMinus_previous ) *
    ( numStepKMinus - numStepKMinus_previous );
  numStepProton2 +=
    ( numStepProton - numStepProton_previous ) *
    ( numStepProton - numStepProton_previous );
  numStepAntiProton2 +=
    ( numStepAntiProton - numStepAntiProton_previous ) *
    ( numStepAntiProton - numStepAntiProton_previous );
  numStepNeutron2 +=
    ( numStepNeutron - numStepNeutron_previous ) *
    ( numStepNeutron - numStepNeutron_previous );
  numStepAntiNeutron2 +=
    ( numStepAntiNeutron - numStepAntiNeutron_previous ) *
    ( numStepAntiNeutron - numStepAntiNeutron_previous );

  numStep_previous = numStep;
  numStepPositive_previous = numStepPositive;
  numStepNeutral_previous = numStepNeutral;
  numStepNegative_previous = numStepNegative;
  numStepPDGCodeZero_previous = numStepPDGCodeZero;
  numStepPDGCodeUnrecognized_previous = numStepPDGCodeUnrecognized;
  numStepEM_previous = numStepEM;
  numStepEWK_previous = numStepEWK;
  numStepHAD_previous = numStepHAD;
  numStepMeson_previous = numStepMeson;
  numStepBaryon_previous = numStepBaryon;
  numStepMesonLight_previous = numStepMesonLight;
  numStepBaryonLight_previous = numStepBaryonLight;
  numStepMesonStrange_previous = numStepMesonStrange;
  numStepBaryonStrange_previous = numStepBaryonStrange;
  numStepMesonHeavy_previous = numStepMesonHeavy;
  numStepBaryonHeavy_previous = numStepBaryonHeavy;
  numStepElectron_previous = numStepElectron;
  numStepGamma_previous = numStepGamma;
  numStepPositron_previous = numStepPositron;
  numStepMuMinus_previous = numStepMuMinus;
  numStepMuPlus_previous = numStepMuPlus;
  numStepTauMinus_previous = numStepTauMinus;
  numStepTauPlus_previous = numStepTauPlus;
  numStepNeutrino_previous = numStepNeutrino;
  numStepPiPlus_previous = numStepPiPlus;
  numStepPi0_previous = numStepPi0;
  numStepPiMinus_previous = numStepPiMinus;
  numStepKPlus_previous = numStepKPlus;
  numStepKNeutral_previous = numStepKNeutral;
  numStepKMinus_previous = numStepKMinus;
  numStepProton_previous = numStepProton;
  numStepAntiProton_previous = numStepAntiProton;
  numStepNeutron_previous = numStepNeutron;
  numStepAntiNeutron_previous = numStepAntiNeutron;

  static G4double numTrack_previous = 0.0;
  static G4double numTrackPositive_previous = 0.0;
  static G4double numTrackNeutral_previous = 0.0;
  static G4double numTrackNegative_previous = 0.0;
  static G4double numTrackPDGCodeZero_previous = 0.0;
  static G4double numTrackPDGCodeUnrecognized_previous = 0.0;
  static G4double numTrackEM_previous = 0.0;
  static G4double numTrackEWK_previous = 0.0;
  static G4double numTrackHAD_previous = 0.0;
  static G4double numTrackMeson_previous = 0.0;
  static G4double numTrackBaryon_previous = 0.0;
  static G4double numTrackMesonLight_previous = 0.0;
  static G4double numTrackBaryonLight_previous = 0.0;
  static G4double numTrackMesonStrange_previous = 0.0;
  static G4double numTrackBaryonStrange_previous = 0.0;
  static G4double numTrackMesonHeavy_previous = 0.0;
  static G4double numTrackBaryonHeavy_previous = 0.0;
  static G4double numTrackElectron_previous = 0.0;
  static G4double numTrackGamma_previous = 0.0;
  static G4double numTrackPositron_previous = 0.0;
  static G4double numTrackMuMinus_previous = 0.0;
  static G4double numTrackMuPlus_previous = 0.0;
  static G4double numTrackTauMinus_previous = 0.0;
  static G4double numTrackTauPlus_previous = 0.0;
  static G4double numTrackNeutrino_previous = 0.0;
  static G4double numTrackPiPlus_previous = 0.0;
  static G4double numTrackPi0_previous = 0.0;
  static G4double numTrackPiMinus_previous = 0.0;
  static G4double numTrackKPlus_previous = 0.0;
  static G4double numTrackKNeutral_previous = 0.0;
  static G4double numTrackKMinus_previous = 0.0;
  static G4double numTrackProton_previous = 0.0;
  static G4double numTrackAntiProton_previous = 0.0;
  static G4double numTrackNeutron_previous = 0.0;
  static G4double numTrackAntiNeutron_previous = 0.0;

  numTrack2 +=
    ( numTrack - numTrack_previous ) *
    ( numTrack - numTrack_previous );
  numTrackPositive2 +=
    ( numTrackPositive - numTrackPositive_previous ) *
    ( numTrackPositive - numTrackPositive_previous );
  numTrackNeutral2 +=
    ( numTrackNeutral - numTrackNeutral_previous ) *
    ( numTrackNeutral - numTrackNeutral_previous );
  numTrackNegative2 +=
    ( numTrackNegative - numTrackNegative_previous ) *
    ( numTrackNegative - numTrackNegative_previous );
  numTrackPDGCodeZero2 +=
    ( numTrackPDGCodeZero - numTrackPDGCodeZero_previous ) *
    ( numTrackPDGCodeZero - numTrackPDGCodeZero_previous );
  numTrackPDGCodeUnrecognized2 +=
    ( numTrackPDGCodeUnrecognized - numTrackPDGCodeUnrecognized_previous ) *
    ( numTrackPDGCodeUnrecognized - numTrackPDGCodeUnrecognized_previous );
  numTrackEM2 +=
    ( numTrackEM - numTrackEM_previous ) *
    ( numTrackEM - numTrackEM_previous );
  numTrackEWK2 +=
    ( numTrackEWK - numTrackEWK_previous ) *
    ( numTrackEWK - numTrackEWK_previous );
  numTrackHAD2 +=
    ( numTrackHAD - numTrackHAD_previous ) *
    ( numTrackHAD - numTrackHAD_previous );
  numTrackMeson2 +=
    ( numTrackMeson - numTrackMeson_previous ) *
    ( numTrackMeson - numTrackMeson_previous );
  numTrackBaryon2 +=
    ( numTrackBaryon - numTrackBaryon_previous ) *
    ( numTrackBaryon - numTrackBaryon_previous );
  numTrackMesonLight2 +=
    ( numTrackMesonLight - numTrackMesonLight_previous ) *
    ( numTrackMesonLight - numTrackMesonLight_previous );
  numTrackBaryonLight2 +=
    ( numTrackBaryonLight - numTrackBaryonLight_previous ) *
    ( numTrackBaryonLight - numTrackBaryonLight_previous );
  numTrackMesonStrange2 +=
    ( numTrackMesonStrange - numTrackMesonStrange_previous ) *
    ( numTrackMesonStrange - numTrackMesonStrange_previous );
  numTrackBaryonStrange2 +=
    ( numTrackBaryonStrange - numTrackBaryonStrange_previous ) *
    ( numTrackBaryonStrange - numTrackBaryonStrange_previous );
  numTrackMesonHeavy2 +=
    ( numTrackMesonHeavy - numTrackMesonHeavy_previous ) *
    ( numTrackMesonHeavy - numTrackMesonHeavy_previous );
  numTrackBaryonHeavy2 +=
    ( numTrackBaryonHeavy - numTrackBaryonHeavy_previous ) *
    ( numTrackBaryonHeavy - numTrackBaryonHeavy_previous );
  numTrackElectron2 +=
    ( numTrackElectron - numTrackElectron_previous ) *
    ( numTrackElectron - numTrackElectron_previous );
  numTrackGamma2 +=
    ( numTrackGamma - numTrackGamma_previous ) *
    ( numTrackGamma - numTrackGamma_previous );
  numTrackPositron2 +=
    ( numTrackPositron - numTrackPositron_previous ) *
    ( numTrackPositron - numTrackPositron_previous );
  numTrackMuMinus2 +=
    ( numTrackMuMinus - numTrackMuMinus_previous ) *
    ( numTrackMuMinus - numTrackMuMinus_previous );
  numTrackMuPlus2 +=
    ( numTrackMuPlus - numTrackMuPlus_previous ) *
    ( numTrackMuPlus - numTrackMuPlus_previous );
  numTrackTauMinus2 +=
    ( numTrackTauMinus - numTrackTauMinus_previous ) *
    ( numTrackTauMinus - numTrackTauMinus_previous );
  numTrackTauPlus2 +=
    ( numTrackTauPlus - numTrackTauPlus_previous ) *
    ( numTrackTauPlus - numTrackTauPlus_previous );
  numTrackNeutrino2 +=
    ( numTrackNeutrino - numTrackNeutrino_previous ) *
    ( numTrackNeutrino - numTrackNeutrino_previous );
  numTrackPiPlus2 +=
    ( numTrackPiPlus - numTrackPiPlus_previous ) *
    ( numTrackPiPlus - numTrackPiPlus_previous );
  numTrackPi02 +=
    ( numTrackPi0 - numTrackPi0_previous ) *
    ( numTrackPi0 - numTrackPi0_previous );
  numTrackPiMinus2 +=
    ( numTrackPiMinus - numTrackPiMinus_previous ) *
    ( numTrackPiMinus - numTrackPiMinus_previous );
  numTrackKPlus2 +=
    ( numTrackKPlus - numTrackKPlus_previous ) *
    ( numTrackKPlus - numTrackKPlus_previous );
  numTrackKNeutral2 +=
    ( numTrackKNeutral - numTrackKNeutral_previous ) *
    ( numTrackKNeutral - numTrackKNeutral_previous );
  numTrackKMinus2 +=
    ( numTrackKMinus - numTrackKMinus_previous ) *
    ( numTrackKMinus - numTrackKMinus_previous );
  numTrackProton2 +=
    ( numTrackProton - numTrackProton_previous ) *
    ( numTrackProton - numTrackProton_previous );
  numTrackAntiProton2 +=
    ( numTrackAntiProton - numTrackAntiProton_previous ) *
    ( numTrackAntiProton - numTrackAntiProton_previous );
  numTrackNeutron2 +=
    ( numTrackNeutron - numTrackNeutron_previous ) *
    ( numTrackNeutron - numTrackNeutron_previous );
  numTrackAntiNeutron2 +=
    ( numTrackAntiNeutron - numTrackAntiNeutron_previous ) *
    ( numTrackAntiNeutron - numTrackAntiNeutron_previous );

  numTrack_previous = numTrack;
  numTrackPositive_previous = numTrackPositive;
  numTrackNeutral_previous = numTrackNeutral;
  numTrackNegative_previous = numTrackNegative;
  numTrackPDGCodeZero_previous = numTrackPDGCodeZero;
  numTrackPDGCodeUnrecognized_previous = numTrackPDGCodeUnrecognized;
  numTrackEM_previous = numTrackEM;
  numTrackEWK_previous = numTrackEWK;
  numTrackHAD_previous = numTrackHAD;
  numTrackMeson_previous = numTrackMeson;
  numTrackBaryon_previous = numTrackBaryon;
  numTrackMesonLight_previous = numTrackMesonLight;
  numTrackBaryonLight_previous = numTrackBaryonLight;
  numTrackMesonStrange_previous = numTrackMesonStrange;
  numTrackBaryonStrange_previous = numTrackBaryonStrange;
  numTrackMesonHeavy_previous = numTrackMesonHeavy;
  numTrackBaryonHeavy_previous = numTrackBaryonHeavy;
  numTrackElectron_previous = numTrackElectron;
  numTrackGamma_previous = numTrackGamma;
  numTrackPositron_previous = numTrackPositron;
  numTrackMuMinus_previous = numTrackMuMinus;
  numTrackMuPlus_previous = numTrackMuPlus;
  numTrackTauMinus_previous = numTrackTauMinus;
  numTrackTauPlus_previous = numTrackTauPlus;
  numTrackNeutrino_previous = numTrackNeutrino;
  numTrackPiPlus_previous = numTrackPiPlus;
  numTrackPi0_previous = numTrackPi0;
  numTrackPiMinus_previous = numTrackPiMinus;
  numTrackKPlus_previous = numTrackKPlus;
  numTrackKNeutral_previous = numTrackKNeutral;
  numTrackKMinus_previous = numTrackKMinus;
  numTrackProton_previous = numTrackProton;
  numTrackAntiProton_previous = numTrackAntiProton;
  numTrackNeutron_previous = numTrackNeutron;
  numTrackAntiNeutron_previous = numTrackAntiNeutron;

  static G4double kinEnergyExiting_previous = 0.0;
  static G4double kinEnergyExitingGammas_previous = 0.0; 
  static G4double kinEnergyExitingNeutrons_previous = 0.0;
  static G4double kinEnergyExitingNeutrinos_previous = 0.0;
  static G4double kinEnergyExitingMuons_previous = 0.0;
  static G4double kinEnergyExitingElectrons_previous = 0.0;
  static G4double kinEnergyExitingOthers_previous = 0.0;
  static G4double numExiting_previous = 0.0;
  static G4double numExitingGammas_previous = 0.0;
  static G4double numExitingNeutrons_previous = 0.0;
  static G4double numExitingNeutrinos_previous = 0.0;
  static G4double numExitingMuons_previous = 0.0;
  static G4double numExitingElectrons_previous = 0.0;
  static G4double numExitingOthers_previous = 0.0;
  
  kinEnergyExiting2 += 
    ( kinEnergyExiting - kinEnergyExiting_previous ) *
    ( kinEnergyExiting - kinEnergyExiting_previous );
  kinEnergyExitingGammas2 += 
    ( kinEnergyExitingGammas - kinEnergyExitingGammas_previous ) *
    ( kinEnergyExitingGammas - kinEnergyExitingGammas_previous );
  kinEnergyExitingNeutrons2 += 
    ( kinEnergyExitingNeutrons - kinEnergyExitingNeutrons_previous ) *
    ( kinEnergyExitingNeutrons - kinEnergyExitingNeutrons_previous );
  kinEnergyExitingNeutrinos2 += 
    ( kinEnergyExitingNeutrinos - kinEnergyExitingNeutrinos_previous ) *
    ( kinEnergyExitingNeutrinos - kinEnergyExitingNeutrinos_previous );
  kinEnergyExitingMuons2 += 
    ( kinEnergyExitingMuons - kinEnergyExitingMuons_previous ) *
    ( kinEnergyExitingMuons - kinEnergyExitingMuons_previous );
  kinEnergyExitingElectrons2 += 
    ( kinEnergyExitingElectrons - kinEnergyExitingElectrons_previous ) *
    ( kinEnergyExitingElectrons - kinEnergyExitingElectrons_previous );
  kinEnergyExitingOthers2 += 
    ( kinEnergyExitingOthers - kinEnergyExitingOthers_previous ) *
    ( kinEnergyExitingOthers - kinEnergyExitingOthers_previous );
  numExiting2 += 
    ( numExiting - numExiting_previous ) *
    ( numExiting - numExiting_previous );
  numExitingGammas2 += 
    ( numExitingGammas - numExitingGammas_previous ) *
    ( numExitingGammas - numExitingGammas_previous );
  numExitingNeutrons2 += 
    ( numExitingNeutrons - numExitingNeutrons_previous ) *
    ( numExitingNeutrons - numExitingNeutrons_previous );
  numExitingNeutrinos2 += 
    ( numExitingNeutrinos - numExitingNeutrinos_previous ) *
    ( numExitingNeutrinos - numExitingNeutrinos_previous );
  numExitingMuons2 += 
    ( numExitingMuons - numExitingMuons_previous ) *
    ( numExitingMuons - numExitingMuons_previous );
  numExitingElectrons2 += 
    ( numExitingElectrons - numExitingElectrons_previous ) *
    ( numExitingElectrons - numExitingElectrons_previous );
  numExitingOthers2 += 
    ( numExitingOthers - numExitingOthers_previous ) *
    ( numExitingOthers - numExitingOthers_previous );

  kinEnergyExiting_previous = kinEnergyExiting;
  kinEnergyExitingGammas_previous = kinEnergyExitingGammas;
  kinEnergyExitingNeutrons_previous =  kinEnergyExitingNeutrons;
  kinEnergyExitingNeutrinos_previous = kinEnergyExitingNeutrinos;
  kinEnergyExitingMuons_previous = kinEnergyExitingMuons;
  kinEnergyExitingElectrons_previous = kinEnergyExitingElectrons;
  kinEnergyExitingOthers_previous = kinEnergyExitingOthers;
  numExiting_previous = numExiting;
  numExitingGammas_previous = numExitingGammas;
  numExitingNeutrons_previous = numExitingNeutrons;
  numExitingNeutrinos_previous = numExitingNeutrinos;
  numExitingMuons_previous = numExitingMuons;
  numExitingElectrons_previous = numExitingElectrons;
  numExitingOthers_previous = numExitingOthers;

  // Apply the same trick for the energy deposits and shower profiles
  // of the following groups of particles:
  // e-  and  e+  together
  static G4double sumEdepAct_electron_previous = 0.0;
  sumEdepAct_electron2 += ( sumEdepAct_electron - sumEdepAct_electron_previous ) *
                          ( sumEdepAct_electron - sumEdepAct_electron_previous );
  sumEdepAct_electron_previous = sumEdepAct_electron;
  static G4double sumEdepTot_electron_previous = 0.0;
  sumEdepTot_electron2 += ( sumEdepTot_electron - sumEdepTot_electron_previous ) *
                          ( sumEdepTot_electron - sumEdepTot_electron_previous );
  sumEdepTot_electron_previous = sumEdepTot_electron;
  for ( int iLayer = 0; iLayer < numberOfReadoutLayers; iLayer++ ) {
    sumL_electron[ iLayer ]  += longitudinalProfile_electron[ iLayer ];
    sumL_electron2[ iLayer ] += longitudinalProfile_electron[ iLayer ] * 
                                longitudinalProfile_electron[ iLayer ];  
    longitudinalProfile_electron[ iLayer ] = 0.0;  // Reset it for the next event.
  }
  for ( int iBinR = 0; iBinR < numberOfRadiusBins; iBinR++ ) {
    sumR_electron[ iBinR ]  += transverseProfile_electron[ iBinR ];
    sumR_electron2[ iBinR ] += transverseProfile_electron[ iBinR ] * 
                               transverseProfile_electron[ iBinR ];  
    transverseProfile_electron[ iBinR ] = 0.0;  // Reset it for the next event.
  }

  // mu-  and  mu+  together
  static G4double sumEdepAct_muon_previous = 0.0;
  sumEdepAct_muon2 += ( sumEdepAct_muon - sumEdepAct_muon_previous ) *
                      ( sumEdepAct_muon - sumEdepAct_muon_previous );
  sumEdepAct_muon_previous = sumEdepAct_muon;
  static G4double sumEdepTot_muon_previous = 0.0;
  sumEdepTot_muon2 += ( sumEdepTot_muon - sumEdepTot_muon_previous ) *
                      ( sumEdepTot_muon - sumEdepTot_muon_previous );
  sumEdepTot_muon_previous = sumEdepTot_muon;
  for ( int iLayer = 0; iLayer < numberOfReadoutLayers; iLayer++ ) {
    sumL_muon[ iLayer ]  += longitudinalProfile_muon[ iLayer ];
    sumL_muon2[ iLayer ] += longitudinalProfile_muon[ iLayer ] * 
                            longitudinalProfile_muon[ iLayer ];  
    longitudinalProfile_muon[ iLayer ] = 0.0;  // Reset it for the next event.
  }
  for ( int iBinR = 0; iBinR < numberOfRadiusBins; iBinR++ ) {
    sumR_muon[ iBinR ]  += transverseProfile_muon[ iBinR ];
    sumR_muon2[ iBinR ] += transverseProfile_muon[ iBinR ] * 
                           transverseProfile_muon[ iBinR ];  
    transverseProfile_muon[ iBinR ] = 0.0;  // Reset it for the next event.
  }

  // pi+  and  pi-  together
  static G4double sumEdepAct_pion_previous = 0.0;
  sumEdepAct_pion2 += ( sumEdepAct_pion - sumEdepAct_pion_previous ) *
                      ( sumEdepAct_pion - sumEdepAct_pion_previous );
  sumEdepAct_pion_previous = sumEdepAct_pion;
  static G4double sumEdepTot_pion_previous = 0.0;
  sumEdepTot_pion2 += ( sumEdepTot_pion - sumEdepTot_pion_previous ) *
                      ( sumEdepTot_pion - sumEdepTot_pion_previous );
  sumEdepTot_pion_previous = sumEdepTot_pion;
  for ( int iLayer = 0; iLayer < numberOfReadoutLayers; iLayer++ ) {
    sumL_pion[ iLayer ]  += longitudinalProfile_pion[ iLayer ];
    sumL_pion2[ iLayer ] += longitudinalProfile_pion[ iLayer ] * 
                            longitudinalProfile_pion[ iLayer ];  
    longitudinalProfile_pion[ iLayer ] = 0.0;  // Reset it for the next event.
  }
  for ( int iBinR = 0; iBinR < numberOfRadiusBins; iBinR++ ) {
    sumR_pion[ iBinR ]  += transverseProfile_pion[ iBinR ];
    sumR_pion2[ iBinR ] += transverseProfile_pion[ iBinR ] * 
                           transverseProfile_pion[ iBinR ];  
    transverseProfile_pion[ iBinR ] = 0.0;  // Reset it for the next event.
  }

  // k+  and  k-  together
  static G4double sumEdepAct_kaon_previous = 0.0;
  sumEdepAct_kaon2 += ( sumEdepAct_kaon - sumEdepAct_kaon_previous ) *
                      ( sumEdepAct_kaon - sumEdepAct_kaon_previous );
  sumEdepAct_kaon_previous = sumEdepAct_kaon;
  static G4double sumEdepTot_kaon_previous = 0.0;
  sumEdepTot_kaon2 += ( sumEdepTot_kaon - sumEdepTot_kaon_previous ) *
                      ( sumEdepTot_kaon - sumEdepTot_kaon_previous );
  sumEdepTot_kaon_previous = sumEdepTot_kaon;
  for ( int iLayer = 0; iLayer < numberOfReadoutLayers; iLayer++ ) {
    sumL_kaon[ iLayer ]  += longitudinalProfile_kaon[ iLayer ];
    sumL_kaon2[ iLayer ] += longitudinalProfile_kaon[ iLayer ] * 
                            longitudinalProfile_kaon[ iLayer ];  
    longitudinalProfile_kaon[ iLayer ] = 0.0;  // Reset it for the next event.
  }
  for ( int iBinR = 0; iBinR < numberOfRadiusBins; iBinR++ ) {
    sumR_kaon[ iBinR ]  += transverseProfile_kaon[ iBinR ];
    sumR_kaon2[ iBinR ] += transverseProfile_kaon[ iBinR ] * 
                           transverseProfile_kaon[ iBinR ];  
    transverseProfile_kaon[ iBinR ] = 0.0;  // Reset it for the next event.
  }

  // p  and  pbar  together
  static G4double sumEdepAct_proton_previous = 0.0;
  sumEdepAct_proton2 += ( sumEdepAct_proton - sumEdepAct_proton_previous ) *
                        ( sumEdepAct_proton - sumEdepAct_proton_previous );
  sumEdepAct_proton_previous = sumEdepAct_proton;
  static G4double sumEdepTot_proton_previous = 0.0;
  sumEdepTot_proton2 += ( sumEdepTot_proton - sumEdepTot_proton_previous ) *
                        ( sumEdepTot_proton - sumEdepTot_proton_previous );
  sumEdepTot_proton_previous = sumEdepTot_proton;
  for ( int iLayer = 0; iLayer < numberOfReadoutLayers; iLayer++ ) {
    sumL_proton[ iLayer ]  += longitudinalProfile_proton[ iLayer ];
    sumL_proton2[ iLayer ] += longitudinalProfile_proton[ iLayer ] * 
                              longitudinalProfile_proton[ iLayer ];  
    longitudinalProfile_proton[ iLayer ] = 0.0;  // Reset it for the next event.
  }
  for ( int iBinR = 0; iBinR < numberOfRadiusBins; iBinR++ ) {
    sumR_proton[ iBinR ]  += transverseProfile_proton[ iBinR ];
    sumR_proton2[ iBinR ] += transverseProfile_proton[ iBinR ] * 
                             transverseProfile_proton[ iBinR ];  
    transverseProfile_proton[ iBinR ] = 0.0;  // Reset it for the next event.
  }

  // all particles with PDG code = 0
  static G4double sumEdepAct_pdg0_previous = 0.0;
  sumEdepAct_pdg02 += ( sumEdepAct_pdg0 - sumEdepAct_pdg0_previous ) *
                      ( sumEdepAct_pdg0 - sumEdepAct_pdg0_previous );
  sumEdepAct_pdg0_previous = sumEdepAct_pdg0;
  static G4double sumEdepTot_pdg0_previous = 0.0;
  sumEdepTot_pdg02 += ( sumEdepTot_pdg0 - sumEdepTot_pdg0_previous ) *
                      ( sumEdepTot_pdg0 - sumEdepTot_pdg0_previous );
  sumEdepTot_pdg0_previous = sumEdepTot_pdg0;
  for ( int iLayer = 0; iLayer < numberOfReadoutLayers; iLayer++ ) {
    sumL_pdg0[ iLayer ]  += longitudinalProfile_pdg0[ iLayer ];
    sumL_pdg02[ iLayer ] += longitudinalProfile_pdg0[ iLayer ] * 
                            longitudinalProfile_pdg0[ iLayer ];  
    longitudinalProfile_pdg0[ iLayer ] = 0.0;  // Reset it for the next event.
  }
  for ( int iBinR = 0; iBinR < numberOfRadiusBins; iBinR++ ) {
    sumR_pdg0[ iBinR ]  += transverseProfile_pdg0[ iBinR ];
    sumR_pdg02[ iBinR ] += transverseProfile_pdg0[ iBinR ] * 
                           transverseProfile_pdg0[ iBinR ];  
    transverseProfile_pdg0[ iBinR ] = 0.0;  // Reset it for the next event.
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
  G4cout << " Primary particle PDG Id = " << primaryParticleId << G4endl;
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
  if ( maxEdepTot - beamEnergy > 0.001*MeV ) {
    G4cout << " maxEdepTot  [MeV] = " << maxEdepTot 
           << "\t ***ENERGY-NON-CONSERVATION*** " << G4endl;
    G4cout << "\t number of event with energy non conservation = " 
	   << countEnergyNonConservation << G4endl;
  } else {
    G4cout << " maxEdepTot  [MeV] = " << maxEdepTot << G4endl;
  }

  G4double n = static_cast< G4double >( numberOfEvents );
  if ( n <= 1.0 ) {
    n = 2.0;         // To avoid division by zero in  sigma
  }
  //G4cout << " n=" << n << G4endl;                             //***DEBUG***
  G4double sum, sum2, mu, sigma, mu_sigma;
  sum  = sumEdepAct;
  sum2 = sumEdepAct2;
  mu       = sum / n;
  sigma    = std::sqrt( std::abs( ( sum2 - sum*sum/n ) ) / (n - 1.0) );
  mu_sigma = sigma / std::sqrt( n );
  G4cout << " Average <E> [MeV] deposited in all active layers = " 
         << mu << " +/- " << mu_sigma << G4endl;
  G4double mu_Evis = mu;             // For later usage.
  G4double mu_Evis_sigma = mu_sigma; //  "    "     "

  sum  = sumEdepTot;
  sum2 = sumEdepTot2;
  mu       = sum / n;
  sigma    = std::sqrt( std::abs( ( sum2 - sum*sum/n ) ) / (n - 1.0) );
  mu_sigma = sigma / std::sqrt( n );
  G4cout << " Average <E> [MeV] deposited in the whole calorimeter = " 
         << mu << " +/- " << mu_sigma << G4endl;
  G4double mu_Etot = mu;             // For later usage.
  G4double mu_Etot_sigma = mu_sigma; //  "    "     "

  G4double fractionLongitudinal1stQuarter = 0.0;
  G4double fractionLongitudinal2ndQuarter = 0.0;
  G4double fractionLongitudinal3rdQuarter = 0.0;
  G4double fractionLongitudinal4thQuarter = 0.0;
  G4double fractionLongitudinal1stQuarter_sigma = 0.0;
  G4double fractionLongitudinal2ndQuarter_sigma = 0.0;
  G4double fractionLongitudinal3rdQuarter_sigma = 0.0;
  G4double fractionLongitudinal4thQuarter_sigma = 0.0;
  G4cout << " Average <E> [MeV] in each Layer " << G4endl; 
  for ( int iLayer = 0; iLayer < numberOfReadoutLayers; iLayer++ ) {
    sum  = sumL[ iLayer ];
    sum2 = sumL2[ iLayer ];
    mu       = sum / n;
    sigma    = std::sqrt( std::abs( ( sum2 - sum*sum/n ) ) / (n - 1.0) );
    mu_sigma = sigma / std::sqrt( n );
    //if ( mu > 1.0E-06 ) {
    G4cout << "\t layer = " << iLayer << "\t <E> = " 
	   << mu << " +/- " << mu_sigma << G4endl;
    //}
    if ( longitudinalProfileHisto ) {
      //***LOOKHERE*** :  mu_sigma  should be set as error.
      longitudinalProfileHisto->fill( 1.0*iLayer, mu );
    }
    if ( iLayer < numberOfReadoutLayers/4 ) {
      fractionLongitudinal1stQuarter += mu;
      fractionLongitudinal1stQuarter_sigma += mu_sigma * mu_sigma;
    } else if ( iLayer < 2*numberOfReadoutLayers/4 ) {
      fractionLongitudinal2ndQuarter += mu;
      fractionLongitudinal2ndQuarter_sigma += mu_sigma * mu_sigma;
    } else if ( iLayer < 3*numberOfReadoutLayers/4 ) {
      fractionLongitudinal3rdQuarter += mu;
      fractionLongitudinal3rdQuarter_sigma += mu_sigma * mu_sigma;
    } else {
      fractionLongitudinal4thQuarter += mu;
      fractionLongitudinal4thQuarter_sigma += mu_sigma * mu_sigma;
    }
  }
  if ( mu_Evis > 1.0E-06 ) {
    G4cout << "  sumL_1 = " << fractionLongitudinal1stQuarter << " +/- "
           << std::sqrt( fractionLongitudinal1stQuarter_sigma ) << G4endl
           << "  sumL_2 = " << fractionLongitudinal2ndQuarter << " +/- "
           << std::sqrt( fractionLongitudinal2ndQuarter_sigma ) << G4endl
           << "  sumL_3 = " << fractionLongitudinal3rdQuarter << " +/- "
           << std::sqrt( fractionLongitudinal3rdQuarter_sigma ) << G4endl
           << "  sumL_4 = " << fractionLongitudinal4thQuarter << " +/- "
           << std::sqrt( fractionLongitudinal4thQuarter_sigma ) << G4endl;
    if ( fractionLongitudinal1stQuarter > 1.0E-06 ) { 
      fractionLongitudinal1stQuarter_sigma /= 
	( fractionLongitudinal1stQuarter * fractionLongitudinal1stQuarter );
    }
    fractionLongitudinal1stQuarter /= mu_Evis;
    fractionLongitudinal1stQuarter_sigma = fractionLongitudinal1stQuarter *
      std::sqrt( fractionLongitudinal1stQuarter_sigma +
		 ( mu_Evis_sigma / mu_Evis ) * ( mu_Evis_sigma / mu_Evis ) );

    if ( fractionLongitudinal2ndQuarter > 1.0E-06 ) { 
      fractionLongitudinal2ndQuarter_sigma /= 
	( fractionLongitudinal2ndQuarter * fractionLongitudinal2ndQuarter );
    }
    fractionLongitudinal2ndQuarter /= mu_Evis;
    fractionLongitudinal2ndQuarter_sigma = fractionLongitudinal2ndQuarter *
      std::sqrt( fractionLongitudinal2ndQuarter_sigma +
		 ( mu_Evis_sigma / mu_Evis ) * ( mu_Evis_sigma / mu_Evis ) );

    if ( fractionLongitudinal3rdQuarter > 1.0E-06 ) { 
      fractionLongitudinal3rdQuarter_sigma /= 
	( fractionLongitudinal3rdQuarter * fractionLongitudinal3rdQuarter );
    }
    fractionLongitudinal3rdQuarter /= mu_Evis;
    fractionLongitudinal3rdQuarter_sigma = fractionLongitudinal3rdQuarter *
      std::sqrt( fractionLongitudinal3rdQuarter_sigma +
		 ( mu_Evis_sigma / mu_Evis ) * ( mu_Evis_sigma / mu_Evis ) );

    if ( fractionLongitudinal4thQuarter > 1.0E-06 ) { 
      fractionLongitudinal4thQuarter_sigma /= 
	( fractionLongitudinal4thQuarter * fractionLongitudinal4thQuarter );
    }
    fractionLongitudinal4thQuarter /= mu_Evis;
    fractionLongitudinal4thQuarter_sigma = fractionLongitudinal4thQuarter *
      std::sqrt( fractionLongitudinal4thQuarter_sigma +
		 ( mu_Evis_sigma / mu_Evis ) * ( mu_Evis_sigma / mu_Evis ) );
  }
  G4cout << " longitudinal fraction in the 1st quarter = "
         << fractionLongitudinal1stQuarter*100.0 << " +/- "
         << fractionLongitudinal1stQuarter_sigma*100.0 << " %" << std::endl
	 << "                              2nd         = "
         << fractionLongitudinal2ndQuarter*100.0 << " +/- "
         << fractionLongitudinal2ndQuarter_sigma*100.0 << " %" << std::endl
	 << "                              3rd         = "
         << fractionLongitudinal3rdQuarter*100.0 << " +/- "
         << fractionLongitudinal3rdQuarter_sigma*100.0 << " %" << std::endl
	 << "                              4th         = "
         << fractionLongitudinal4thQuarter*100.0 << " +/- "
         << fractionLongitudinal4thQuarter_sigma*100.0 << " %" << std::endl;

  G4double fractionTransverse1stThird = 0.0;
  G4double fractionTransverse2ndThird = 0.0;
  G4double fractionTransverse3rdThird = 0.0;
  G4double fractionTransverse1stThird_sigma = 0.0;
  G4double fractionTransverse2ndThird_sigma = 0.0;
  G4double fractionTransverse3rdThird_sigma = 0.0;
  // // // std::vector< G4double > rmsTransverseProfile;   //***TEMPORARY WORK-AROUND***
  G4cout << " Average <E> [MeV] in each Radius bin " << G4endl; 
  for ( int iBinR = 0; iBinR < numberOfRadiusBins; iBinR++ ) {
    sum  = sumR[ iBinR ];
    sum2 = sumR2[ iBinR ];
    mu       = sum / n;
    sigma    = std::sqrt( std::abs( ( sum2 - sum*sum/n ) ) / (n - 1.0) );
    mu_sigma = sigma / std::sqrt( n );
    //if ( mu > 1.0E-06 ) {
    G4cout << "\t iBinR = " << iBinR << "\t <E> = " 
	   << mu << " +/- " << mu_sigma << G4endl;
    //}
    if ( transverseProfileHisto ) {
      //***LOOKHERE*** :  mu_sigma  should be set as error.
      transverseProfileHisto->fill( 1.0*iBinR, mu );
    }
    // // // rmsTransverseProfile.push_back( mu_sigma );   //***TEMPORARY WORK-AROUND***
    if ( iBinR < numberOfRadiusBins/3 ) {
      fractionTransverse1stThird += mu;
      fractionTransverse1stThird_sigma += mu_sigma * mu_sigma;
    } else if ( iBinR < 2*numberOfRadiusBins/3 ) {
      fractionTransverse2ndThird += mu;
      fractionTransverse2ndThird_sigma += mu_sigma * mu_sigma;
    } else {
      fractionTransverse3rdThird += mu;
      fractionTransverse3rdThird_sigma += mu_sigma * mu_sigma;
    }
  }  
  if ( mu_Evis > 1.0E-06 ) {
    G4cout << "  sumR_1 = " << fractionTransverse1stThird << " +/- "
           << std::sqrt( fractionTransverse1stThird_sigma ) << G4endl
           << "  sumR_2 = " << fractionTransverse2ndThird << " +/- "
           << std::sqrt( fractionTransverse2ndThird_sigma ) << G4endl
           << "  sumR_3 = " << fractionTransverse3rdThird << " +/- "
           << std::sqrt( fractionTransverse3rdThird_sigma ) << G4endl;
    if ( fractionTransverse1stThird > 1.0E-06 ) {
      fractionTransverse1stThird_sigma /=
	( fractionTransverse1stThird * fractionTransverse1stThird );
    } 
    fractionTransverse1stThird /= mu_Evis;
    fractionTransverse1stThird_sigma = fractionTransverse1stThird *
      std::sqrt( fractionTransverse1stThird_sigma +  
		 ( mu_Evis_sigma / mu_Evis ) * ( mu_Evis_sigma / mu_Evis ) );

    if ( fractionTransverse2ndThird > 1.0E-06 ) {
      fractionTransverse2ndThird_sigma /=
	( fractionTransverse2ndThird * fractionTransverse2ndThird );
    } 
    fractionTransverse2ndThird /= mu_Evis;
    fractionTransverse2ndThird_sigma = fractionTransverse2ndThird *
      std::sqrt( fractionTransverse2ndThird_sigma +  
		 ( mu_Evis_sigma / mu_Evis ) * ( mu_Evis_sigma / mu_Evis ) );

    if ( fractionTransverse3rdThird > 1.0E-06 ) {
      fractionTransverse3rdThird_sigma /=
	( fractionTransverse3rdThird * fractionTransverse3rdThird );
    } 
    fractionTransverse3rdThird /= mu_Evis;
    fractionTransverse3rdThird_sigma = fractionTransverse3rdThird *
      std::sqrt( fractionTransverse3rdThird_sigma +  
		 ( mu_Evis_sigma / mu_Evis ) * ( mu_Evis_sigma / mu_Evis ) );
  }
  G4cout << " transverse fraction in the 1st third = "
         << fractionTransverse1stThird*100.0 << " +/- "
         << fractionTransverse1stThird_sigma*100.0 << " %" << std::endl
	 << "                            2nd       = "
         << fractionTransverse2ndThird*100.0 << " +/- " 
         << fractionTransverse2ndThird_sigma*100.0 << " %" << std::endl
	 << "                            3rd       = "
         << fractionTransverse3rdThird*100.0 << " +/- " 
         << fractionTransverse3rdThird_sigma*100.0 << " %" << std::endl;

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
  width_Evis *= std::sqrt( 3.141592654/2.0 ) / n ;
  G4double width_Evis_sigma = width_Evis / std::sqrt( 2.0*(n - 1) );
  G4double energyResolution = 0.0;
  G4double energyResolution_sigma = 0.0;
  G4double samplingFraction = 0.0;
  G4double samplingFraction_sigma = 0.0;
  if ( mu_Evis > 1.0E-06) {
    energyResolution = width_Evis / mu_Evis;
    if ( width_Evis > 1.0E-06 ) {
      energyResolution_sigma = energyResolution *
	std::sqrt( ( width_Evis_sigma * width_Evis_sigma ) / ( width_Evis * width_Evis ) 
		   +
		   ( mu_Evis_sigma * mu_Evis_sigma ) / ( mu_Evis + mu_Evis ) );
    }
    samplingFraction = mu_Evis / beamEnergy;
    samplingFraction_sigma = samplingFraction * mu_Evis_sigma / mu_Evis;
  }

  G4cout << " Visible energy information [MeV] " << G4endl
	 << "\t mu_Evis    = " << mu_Evis << " +/- " << mu_Evis_sigma << G4endl
         << "\t sigma_Evis = " << width_Evis << " +/- " << width_Evis_sigma << G4endl
         << "\t energy resolution = " << energyResolution << " +/- "  
	 << energyResolution_sigma << G4endl
         << "\t sampling fraction = " << samplingFraction << " +/- " 
	 << samplingFraction_sigma << G4endl;

  if ( tree ) tree->commit();

  // Print information on the different particle contributions to
  // the visible energy, the total energy, and the shower shapes.

  G4cout << G4endl << " Contributions of the main particle types [MeV] " << G4endl;
  for ( int iCase = 0; iCase < 6; iCase++ ) {

    std::string caseName = "";
    G4double sumVis = 0.0;
    G4double sumVis2 = 0.0;
    G4double sumTot = 0.0;
    G4double sumTot2 = 0.0;
    std::vector< G4double > vecSumL;
    std::vector< G4double > vecSumL2;
    std::vector< G4double > vecSumR;
    std::vector< G4double > vecSumR2;

    switch ( iCase ) {
    case 0 : {
      caseName = "electron";
      sumVis = sumEdepAct_electron;
      sumVis2 = sumEdepAct_electron2;
      sumTot = sumEdepTot_electron;
      sumTot2 = sumEdepTot_electron2;
      for ( int iLayer = 0; iLayer < numberOfReadoutLayers; iLayer++ ) {
	vecSumL.push_back( sumL_electron[ iLayer ] );
	vecSumL2.push_back( sumL_electron2[ iLayer ] );
      }
      for ( int iBinR = 0; iBinR < numberOfRadiusBins; iBinR++ ) {
	vecSumR.push_back( sumR_electron[ iBinR ] );
	vecSumR2.push_back( sumR_electron2[ iBinR ] );
      }
      break;
    }
    case 1 : {
      caseName = "muon";
      sumVis = sumEdepAct_muon;
      sumVis2 = sumEdepAct_muon2;
      sumTot = sumEdepTot_muon;
      sumTot2 = sumEdepTot_muon2;
      for ( int iLayer = 0; iLayer < numberOfReadoutLayers; iLayer++ ) {
	vecSumL.push_back( sumL_muon[ iLayer ] );
	vecSumL2.push_back( sumL_muon2[ iLayer ] );
      }
      for ( int iBinR = 0; iBinR < numberOfRadiusBins; iBinR++ ) {
	vecSumR.push_back( sumR_muon[ iBinR ] );
	vecSumR2.push_back( sumR_muon2[ iBinR ] );
      }
      break;
    }
    case 2 : {
      caseName = "pion";
      sumVis = sumEdepAct_pion;
      sumVis2 = sumEdepAct_pion2;
      sumTot = sumEdepTot_pion;
      sumTot2 = sumEdepTot_pion2;
      for ( int iLayer = 0; iLayer < numberOfReadoutLayers; iLayer++ ) {
	vecSumL.push_back( sumL_pion[ iLayer ] );
	vecSumL2.push_back( sumL_pion2[ iLayer ] );
      }
      for ( int iBinR = 0; iBinR < numberOfRadiusBins; iBinR++ ) {
	vecSumR.push_back( sumR_pion[ iBinR ] );
	vecSumR2.push_back( sumR_pion2[ iBinR ] );
      }
      break;
    }
    case 3 : {
      caseName = "kaon";
      sumVis = sumEdepAct_kaon;
      sumVis2 = sumEdepAct_kaon2;
      sumTot = sumEdepTot_kaon;
      sumTot2 = sumEdepTot_kaon2;
      for ( int iLayer = 0; iLayer < numberOfReadoutLayers; iLayer++ ) {
	vecSumL.push_back( sumL_kaon[ iLayer ] );
	vecSumL2.push_back( sumL_kaon2[ iLayer ] );
      }
      for ( int iBinR = 0; iBinR < numberOfRadiusBins; iBinR++ ) {
	vecSumR.push_back( sumR_kaon[ iBinR ] );
	vecSumR2.push_back( sumR_kaon2[ iBinR ] );
      }
      break;
    }
    case 4 : {
      caseName = "proton";
      sumVis = sumEdepAct_proton;
      sumVis2 = sumEdepAct_proton2;
      sumTot = sumEdepTot_proton;
      sumTot2 = sumEdepTot_proton2;
      for ( int iLayer = 0; iLayer < numberOfReadoutLayers; iLayer++ ) {
	vecSumL.push_back( sumL_proton[ iLayer ] );
	vecSumL2.push_back( sumL_proton2[ iLayer ] );
      }
      for ( int iBinR = 0; iBinR < numberOfRadiusBins; iBinR++ ) {
	vecSumR.push_back( sumR_proton[ iBinR ] );
	vecSumR2.push_back( sumR_proton2[ iBinR ] );
      }
      break;
    }
    case 5 : {
      caseName = "pdg0";
      sumVis = sumEdepAct_pdg0;
      sumVis2 = sumEdepAct_pdg02;
      sumTot = sumEdepTot_pdg0;
      sumTot2 = sumEdepTot_pdg02;
      for ( int iLayer = 0; iLayer < numberOfReadoutLayers; iLayer++ ) {
	vecSumL.push_back( sumL_pdg0[ iLayer ] );
	vecSumL2.push_back( sumL_pdg02[ iLayer ] );
      }
      for ( int iBinR = 0; iBinR < numberOfRadiusBins; iBinR++ ) {
	vecSumR.push_back( sumR_pdg0[ iBinR ] );
	vecSumR2.push_back( sumR_pdg02[ iBinR ] );
      }
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
    if ( mu_Evis > 1.0E-06 ) {
      f = mu / mu_Evis;
      if ( mu > 1.0E-06 ) {
	f_sigma = f * std::sqrt( ( (mu_sigma*mu_sigma) / (mu*mu) ) +
				 ( (mu_Evis_sigma*mu_Evis_sigma) / (mu_Evis*mu_Evis) ) );
      }
    }
    G4cout << "\t \t <E_vis> = " << mu << " +/- " << mu_sigma 
	   << "  ( " << 100.0*f << " +/- " << 100.0*f_sigma << " % )" << G4endl;
    G4double particle_mu_Evis = mu;             // For later usage.
    G4double particle_mu_Evis_sigma = mu_sigma; //  "    "     "

    sum  = sumTot;
    sum2 = sumTot2;
    mu       = sum / n;
    sigma    = std::sqrt( std::abs( sum2 - sum*sum/n ) / (n - 1.0) );
    mu_sigma = sigma / std::sqrt( n );
    f = 0.0;
    f_sigma = 0.0;
    if ( mu_Etot > 1.0E-06 ) {
      f = mu / mu_Etot;
      if ( mu > 1.0E-06 ) {
	f_sigma = f * std::sqrt( ( (mu_sigma*mu_sigma) / (mu*mu) ) +
				 ( (mu_Etot_sigma*mu_Etot_sigma) / (mu_Etot*mu_Etot) ) );
      }
    }
    G4cout << "\t \t <E_tot> = " << mu << " +/- " << mu_sigma 
	   << "  ( " << 100.0*f << " +/- " << 100.0*f_sigma << " % )" << G4endl;

    G4double fLongitudinal1stQuarter = 0.0;
    G4double fLongitudinal2ndQuarter = 0.0;
    G4double fLongitudinal3rdQuarter = 0.0;
    G4double fLongitudinal4thQuarter = 0.0;
    G4double fLongitudinal1stQuarter_sigma = 0.0;
    G4double fLongitudinal2ndQuarter_sigma = 0.0;
    G4double fLongitudinal3rdQuarter_sigma = 0.0;
    G4double fLongitudinal4thQuarter_sigma = 0.0;
    G4cout << "\t \t Average <E> [MeV] in each Layer " << G4endl; 
    for ( int iLayer = 0; iLayer < numberOfReadoutLayers; iLayer++ ) {
      sum  = vecSumL[ iLayer ];
      sum2 = vecSumL2[ iLayer ];
      mu       = sum / n;
      sigma    = std::sqrt( std::abs( ( sum2 - sum*sum/n ) ) / (n - 1.0) );
      mu_sigma = sigma / std::sqrt( n );
      G4cout << "\t \t \t layer = " << iLayer << "\t <E> = " 
	     << mu << " +/- " << mu_sigma << G4endl;
      if ( iLayer < numberOfReadoutLayers/4 ) {
	fLongitudinal1stQuarter += mu;
	fLongitudinal1stQuarter_sigma += mu_sigma * mu_sigma;
      } else if ( iLayer < 2*numberOfReadoutLayers/4 ) {
	fLongitudinal2ndQuarter += mu;
	fLongitudinal2ndQuarter_sigma += mu_sigma * mu_sigma;
      } else if ( iLayer < 3*numberOfReadoutLayers/4 ) {
	fLongitudinal3rdQuarter += mu;
	fLongitudinal3rdQuarter_sigma += mu_sigma * mu_sigma;
      } else {
	fLongitudinal4thQuarter += mu;
	fLongitudinal4thQuarter_sigma += mu_sigma * mu_sigma;
      }
    }
    if ( particle_mu_Evis > 1.0E-06 ) {
      G4cout << "\t \t  sumL_1 = " << fLongitudinal1stQuarter << " +/- "
	     << std::sqrt( fLongitudinal1stQuarter_sigma ) << G4endl
	     << "\t \t  sumL_2 = " << fLongitudinal2ndQuarter << " +/- "
	     << std::sqrt( fLongitudinal2ndQuarter_sigma ) << G4endl
	     << "\t \t  sumL_3 = " << fLongitudinal3rdQuarter << " +/- "
	     << std::sqrt( fLongitudinal3rdQuarter_sigma ) << G4endl
	     << "\t \t  sumL_4 = " << fLongitudinal4thQuarter << " +/- "
	     << std::sqrt( fLongitudinal4thQuarter_sigma ) << G4endl;
      if ( fLongitudinal1stQuarter > 1.0E-06 ) { 
	fLongitudinal1stQuarter_sigma /= 
	  ( fLongitudinal1stQuarter * fLongitudinal1stQuarter );
      }
      fLongitudinal1stQuarter /= particle_mu_Evis;
      fLongitudinal1stQuarter_sigma = fLongitudinal1stQuarter *
	std::sqrt( fLongitudinal1stQuarter_sigma +
		   ( particle_mu_Evis_sigma / particle_mu_Evis ) * 
		   ( particle_mu_Evis_sigma / particle_mu_Evis ) );
      
      if ( fLongitudinal2ndQuarter > 1.0E-06 ) { 
	fLongitudinal2ndQuarter_sigma /= 
	  ( fLongitudinal2ndQuarter * fLongitudinal2ndQuarter );
      }
      fLongitudinal2ndQuarter /= particle_mu_Evis;
      fLongitudinal2ndQuarter_sigma = fLongitudinal2ndQuarter *
	std::sqrt( fLongitudinal2ndQuarter_sigma +
		   ( particle_mu_Evis_sigma / particle_mu_Evis ) * 
		   ( particle_mu_Evis_sigma / particle_mu_Evis ) );

      if ( fLongitudinal3rdQuarter > 1.0E-06 ) { 
	fLongitudinal3rdQuarter_sigma /= 
	  ( fLongitudinal3rdQuarter * fLongitudinal3rdQuarter );
      }
      fLongitudinal3rdQuarter /= particle_mu_Evis;
      fLongitudinal3rdQuarter_sigma = fLongitudinal3rdQuarter *
	std::sqrt( fLongitudinal3rdQuarter_sigma +
		   ( particle_mu_Evis_sigma / particle_mu_Evis ) * 
		   ( particle_mu_Evis_sigma / particle_mu_Evis ) );
      
      if ( fLongitudinal4thQuarter > 1.0E-06 ) { 
	fLongitudinal4thQuarter_sigma /= 
	  ( fLongitudinal4thQuarter * fLongitudinal4thQuarter );
      }
      fLongitudinal4thQuarter /= particle_mu_Evis;
      fLongitudinal4thQuarter_sigma = fLongitudinal4thQuarter *
	std::sqrt( fLongitudinal4thQuarter_sigma +
		   ( particle_mu_Evis_sigma / particle_mu_Evis ) * 
		   ( particle_mu_Evis_sigma / particle_mu_Evis ) );
    }
    G4cout << "\t \t longitudinal fraction in the 1st quarter = "
	   << fLongitudinal1stQuarter*100.0 << " +/- "
	   << fLongitudinal1stQuarter_sigma*100.0 << " %" << std::endl
	   << "\t \t                              2nd         = "
	   << fLongitudinal2ndQuarter*100.0 << " +/- "
	   << fLongitudinal2ndQuarter_sigma*100.0 << " %" << std::endl
	   << "\t \t                              3rd         = "
	   << fLongitudinal3rdQuarter*100.0 << " +/- "
	   << fLongitudinal3rdQuarter_sigma*100.0 << " %" << std::endl
	   << "\t \t                              4th         = "
	   << fLongitudinal4thQuarter*100.0 << " +/- "
	   << fLongitudinal4thQuarter_sigma*100.0 << " %" << std::endl;
    
    G4double fTransverse1stThird = 0.0;
    G4double fTransverse2ndThird = 0.0;
    G4double fTransverse3rdThird = 0.0;
    G4double fTransverse1stThird_sigma = 0.0;
    G4double fTransverse2ndThird_sigma = 0.0;
    G4double fTransverse3rdThird_sigma = 0.0;
    G4cout << "\t \t Average <E> [MeV] in each Radius bin " << G4endl; 
    for ( int iBinR = 0; iBinR < numberOfRadiusBins; iBinR++ ) {
      sum  = vecSumR[ iBinR ];
      sum2 = vecSumR2[ iBinR ];
      mu       = sum / n;
      sigma    = std::sqrt( std::abs( ( sum2 - sum*sum/n ) ) / (n - 1.0) );
      mu_sigma = sigma / std::sqrt( n );
      G4cout << "\t \t \t iBinR = " << iBinR << "\t <E> = " 
	     << mu << " +/- " << mu_sigma << G4endl;
      if ( iBinR < numberOfRadiusBins/3 ) {
	fTransverse1stThird += mu;
	fTransverse1stThird_sigma += mu_sigma * mu_sigma;
      } else if ( iBinR < 2*numberOfRadiusBins/3 ) {
	fTransverse2ndThird += mu;
	fTransverse2ndThird_sigma += mu_sigma * mu_sigma;
      } else {
	fTransverse3rdThird += mu;
	fTransverse3rdThird_sigma += mu_sigma * mu_sigma;
      }
    }  
    if ( particle_mu_Evis > 1.0E-06 ) {
      G4cout << "\t \t  sumR_1 = " << fTransverse1stThird << " +/- "
	     << std::sqrt( fTransverse1stThird_sigma ) << G4endl
	     << "\t \t  sumR_2 = " << fTransverse2ndThird << " +/- "
	     << std::sqrt( fTransverse2ndThird_sigma ) << G4endl
	     << "\t \t  sumR_3 = " << fTransverse3rdThird << " +/- "
	     << std::sqrt( fTransverse3rdThird_sigma ) << G4endl;
      if ( fTransverse1stThird > 1.0E-06 ) {
	fTransverse1stThird_sigma /=
	  ( fTransverse1stThird * fTransverse1stThird );
      } 
      fTransverse1stThird /= particle_mu_Evis;
      fTransverse1stThird_sigma = fTransverse1stThird *
	std::sqrt( fTransverse1stThird_sigma +  
		   ( particle_mu_Evis_sigma / particle_mu_Evis ) * 
		   ( particle_mu_Evis_sigma / particle_mu_Evis ) );

      if ( fTransverse2ndThird > 1.0E-06 ) {
	fTransverse2ndThird_sigma /=
	  ( fTransverse2ndThird * fTransverse2ndThird );
      } 
      fTransverse2ndThird /= particle_mu_Evis;
      fTransverse2ndThird_sigma = fTransverse2ndThird *
	std::sqrt( fTransverse2ndThird_sigma +  
		   ( particle_mu_Evis_sigma / particle_mu_Evis ) * 
		   ( particle_mu_Evis_sigma / particle_mu_Evis ) );
      
      if ( fTransverse3rdThird > 1.0E-06 ) {
	fTransverse3rdThird_sigma /=
	  ( fTransverse3rdThird * fTransverse3rdThird );
      } 
      fTransverse3rdThird /= particle_mu_Evis;
      fTransverse3rdThird_sigma = fTransverse3rdThird *
	std::sqrt( fTransverse3rdThird_sigma +  
		   ( particle_mu_Evis_sigma / particle_mu_Evis ) * 
		   ( particle_mu_Evis_sigma / particle_mu_Evis ) );
    }
    G4cout << "\t \t transverse fraction in the 1st third = "
	   << fTransverse1stThird*100.0 << " +/- "
	   << fTransverse1stThird_sigma*100.0 << " %" << std::endl
	   << "\t \t                            2nd       = "
	   << fTransverse2ndThird*100.0 << " +/- " 
	   << fTransverse2ndThird_sigma*100.0 << " %" << std::endl
	   << "\t \t                            3rd       = "
	   << fTransverse3rdThird*100.0 << " +/- " 
	   << fTransverse3rdThird_sigma*100.0 << " %" << std::endl;
    
  } // End of the loop over the particle types.

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
	  sum = numStep;
	  sum2 = numStep2;
	} else {
	  sum = numTrack;
	  sum2 = numTrack2;
	}
	break;
      }
      case 1 : {
	caseName += "# positives";
	if ( i == 0 ) {
	  sum = numStepPositive;
	  sum2 = numStepPositive2;
	} else {
	  sum = numTrackPositive;
	  sum2 = numTrackPositive2;
	}
	break;
      }
      case 2 : {
	caseName += "# neutrals";
	if ( i == 0 ) {
	  sum = numStepNeutral;
	  sum2 = numStepNeutral2;
	} else {
	  sum = numTrackNeutral;
	  sum2 = numTrackNeutral2;
	}
	break;
      }
      case 3 : {
	caseName += "# negatives";
	if ( i == 0 ) {
	  sum = numStepNegative;
	  sum2 = numStepNegative2;
	} else {
	  sum = numTrackNegative;
	  sum2 = numTrackNegative2;
	}
	break;
      }
      case 4 : {
	caseName += "# particles with 0 PDG code";
	if ( i == 0 ) {
	  sum = numStepPDGCodeZero;
	  sum2 = numStepPDGCodeZero2;
	} else {
	  sum = numTrackPDGCodeZero;
	  sum2 = numTrackPDGCodeZero2;
	}
	break;
      }
      case 5 : {
	caseName += "# particles with Unrecognized PDG code";
	if ( i == 0 ) {
	  sum = numStepPDGCodeUnrecognized;
	  sum2 = numStepPDGCodeUnrecognized2;
	} else {
	  sum = numTrackPDGCodeUnrecognized;
	  sum2 = numTrackPDGCodeUnrecognized2;
	}
	break;
      }
      case 6 : {
	caseName += "# electromagnetic (e+ , e- , gammas)";
	if ( i == 0 ) {
	  sum = numStepEM;
	  sum2 = numStepEM2;
	} else {
	  sum = numTrackEM;
	  sum2 = numTrackEM2;
	}
	break;
      }
      case 7 : {
	caseName += "# electroweak (mu+, mu-, tau+, tau-, neutrinos)";
	if ( i == 0 ) {
	  sum = numStepEWK;
	  sum2 = numStepEWK2;
	} else {
	  sum = numTrackEWK;
	  sum2 = numTrackEWK2;
	}
	break;
      }
      case 8 : {
	caseName += "# hadrons";
	if ( i == 0 ) {
	  sum = numStepHAD;
	  sum2 = numStepHAD2;
	} else {
	  sum = numTrackHAD;
	  sum2 = numTrackHAD2;
	}
	break;
      }
      case 9 : {
	caseName += "# mesons";
	if ( i == 0 ) {
	  sum = numStepMeson;
	  sum2 = numStepMeson2;
	} else {
	  sum = numTrackMeson;
	  sum2 = numTrackMeson2;
	}
	break;
      }
      case 10 : {
	caseName += "# baryons";
	if ( i == 0 ) {
	  sum = numStepBaryon;
	  sum2 = numStepBaryon2;
	} else {
	  sum = numTrackBaryon;
	  sum2 = numTrackBaryon2;
	}
	break;
      }
      case 11 : {
	caseName += "# light mesons (u/ubar/d/dbar)";
	if ( i == 0 ) {
	  sum = numStepMesonLight;
	  sum2 = numStepMesonLight2;
	} else {
	  sum = numTrackMesonLight;
	  sum2 = numTrackMesonLight2;
	}
	break;
      }
      case 12 : {
	caseName += "# light baryons (u/ubar/d/dbar)";
	if ( i == 0 ) {
	  sum = numStepBaryonLight;
	  sum2 = numStepBaryonLight2;
	} else {
	  sum = numTrackBaryonLight;
	  sum2 = numTrackBaryonLight2;
	}
	break;
      }
      case 13 : {
	caseName += "# strange (s/sbar) mesons";
	if ( i == 0 ) {
	  sum = numStepMesonStrange;
	  sum2 = numStepMesonStrange2;
	} else {
	  sum = numTrackMesonStrange;
	  sum2 = numTrackMesonStrange2;
	}
	break;
      }
      case 14 : {
	caseName += "# strange (s/sbar) baryons";
	if ( i == 0 ) {
	  sum = numStepBaryonStrange;
	  sum2 = numStepBaryonStrange2;
	} else {
	  sum = numTrackBaryonStrange;
	  sum2 = numTrackBaryonStrange2;
	}
	break;
      }
      case 15 : {
	caseName += "# heavy (c/cbar or b/bbar) mesons";
	if ( i == 0 ) {
	  sum = numStepMesonHeavy;
	  sum2 = numStepMesonHeavy2;
	} else {
	  sum = numTrackMesonHeavy;
	  sum2 = numTrackMesonHeavy2;
	}
	break;
      }
      case 16 : {
	caseName += "# heavy (c/cbar or b/bbar) baryons";
	if ( i == 0 ) {
	  sum = numStepBaryonHeavy;
	  sum2 = numStepBaryonHeavy2;
	} else {
	  sum = numTrackBaryonHeavy;
	  sum2 = numTrackBaryonHeavy2;
	}
	break;
      }
      case 17 : {
	caseName += "# electrons";
	if ( i == 0 ) {
	  sum = numStepElectron;
	  sum2 = numStepElectron2;
	} else {
	  sum = numTrackElectron;
	  sum2 = numTrackElectron2;
	}
	break;
      }
      case 18 : {
	caseName += "# gammas";
	if ( i == 0 ) {
	  sum = numStepGamma;
	  sum2 = numStepGamma2;
	} else {
	  sum = numTrackGamma;
	  sum2 = numTrackGamma2;
	}
	break;
      }
      case 19 : {
	caseName += "# positrons";
	if ( i == 0 ) {
	  sum = numStepPositron;
	  sum2 = numStepPositron2;
	} else {
	  sum = numTrackPositron;
	  sum2 = numTrackPositron2;
	}
	break;
      }
      case 20 : {
	caseName += "# mu-";
	if ( i == 0 ) {
	  sum = numStepMuMinus;
	  sum2 = numStepMuMinus2;
	} else {
	  sum = numTrackMuMinus;
	  sum2 = numTrackMuMinus2;
	}
	break;
      }
      case 21 : {
	caseName += "# mu+";
	if ( i == 0 ) {
	  sum = numStepMuPlus;
	  sum2 = numStepMuPlus2;
	} else {
	  sum = numTrackMuPlus;
	  sum2 = numTrackMuPlus2;
	}
	break;
      }
      case 22 : {
	caseName += "# tau-";
	if ( i == 0 ) {
	  sum = numStepTauMinus;
	  sum2 = numStepTauMinus2;
	} else {
	  sum = numTrackTauMinus;
	  sum2 = numTrackTauMinus2;
	}
	break;
      }
      case 23 : {
	caseName += "# tau+";
	if ( i == 0 ) {
	  sum = numStepTauPlus;
	  sum2 = numStepTauPlus2;
	} else {
	  sum = numTrackTauPlus;
	  sum2 = numTrackTauPlus2;
	}
	break;
      }
      case 24 : {
	caseName += "# neutrinos";
	if ( i == 0 ) {
	  sum = numStepNeutrino;
	  sum2 = numStepNeutrino2;
	} else {
	  sum = numTrackNeutrino;
	  sum2 = numTrackNeutrino2;
	}
	break;
      }
      case 25 : {
	caseName += "# pi+";
	if ( i == 0 ) {
	  sum = numStepPiPlus;
	  sum2 = numStepPiPlus2;
	} else {
	  sum = numTrackPiPlus;
	  sum2 = numTrackPiPlus2;
	}
	break;
      }
      case 26 : {
	caseName += "# pi0";
	if ( i == 0 ) {
	  sum = numStepPi0;
	  sum2 = numStepPi02;
	} else {
	  sum = numTrackPi0;
	  sum2 = numTrackPi02;
	}
	break;
      }
      case 27 : {
	caseName += "# pi-";
	if ( i == 0 ) {
	  sum = numStepPiMinus;
	  sum2 = numStepPiMinus2;
	} else {
	  sum = numTrackPiMinus;
	  sum2 = numTrackPiMinus2;
	}
	break;
      }
      case 28 : {
	caseName += "# K+";
	if ( i == 0 ) {
	  sum = numStepKPlus;
	  sum2 = numStepKPlus2;
	} else {
	  sum = numTrackKPlus;
	  sum2 = numTrackKPlus2;
	}
	break;
      }
      case 29 : {
	caseName += "# K-neutral (K0/K0bar or K0_S/K0_L)";
	if ( i == 0 ) {
	  sum = numStepKNeutral;
	  sum2 = numStepKNeutral2;
	} else {
	  sum = numTrackKNeutral;
	  sum2 = numTrackKNeutral2;
	}
	break;
      }
      case 30 : {
	caseName += "# K-";
	if ( i == 0 ) {
	  sum = numStepKMinus;
	  sum2 = numStepKMinus2;
	} else {
	  sum = numTrackKMinus;
	  sum2 = numTrackKMinus2;
	}
	break;
      }
      case 31 : {
	caseName += "# protons";
	if ( i == 0 ) {
	  sum = numStepProton;
	  sum2 = numStepProton2;
	} else {
	  sum = numTrackProton;
	  sum2 = numTrackProton2;
	}
	break;
      }
      case 32 : {
	caseName += "# anti-protons";
	if ( i == 0 ) {
	  sum = numStepAntiProton;
	  sum2 = numStepAntiProton2;
	} else {
	  sum = numTrackAntiProton;
	  sum2 = numTrackAntiProton2;
	}
	break;
      }
      case 33 : {
	caseName += "# neutrons";
	if ( i == 0 ) {
	  sum = numStepNeutron;
	  sum2 = numStepNeutron2;
	} else {
	  sum = numTrackNeutron;
	  sum2 = numTrackNeutron2;
	}
	break;
      }
      case 34 : {
	caseName += "# anti-neutrons";
	if ( i == 0 ) {
	  sum = numStepAntiNeutron;
	  sum2 = numStepAntiNeutron2;
	} else {
	  sum = numTrackAntiNeutron;
	  sum2 = numTrackAntiNeutron2;
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
      mu = sum / n;
      sigma = std::sqrt( std::abs( ( sum2 - sum*sum/n ) ) / (n - 1.0) );
      mu_sigma = sigma / std::sqrt( n );
      G4cout << "\t" << caseName << " = " << mu  
	     << " +/- " << mu_sigma << G4endl;
    }
  }

  // Print information about track length:
  G4cout << G4endl << " Average track LENGTH [mm] " << G4endl;
  G4double nn = 0.0;
  for ( int iCase = 0; iCase < 7; iCase++ ) {
    if ( iCase == 1  ||  iCase == 5 ) continue; // Exclude muons and pi0s
    std::string caseName;
    switch ( iCase ) {
    case 0 : {
      nn = numTrackElectron + numTrackPositron;
      sum  = electronTrackLength;
      sum2 = electronTrackLength2; 
      caseName = "electron/positron";
      break;
    }
    case 1 : {
      nn = numTrackMuMinus + numTrackMuPlus;
      sum  = muonTrackLength;
      sum2 = muonTrackLength2; 
      caseName = "muon-/muon+";
      break;
    }
    case 2 : {
      nn = numTrackPiPlus + numTrackPiMinus;
      sum  = pionChargedTrackLength;
      sum2 = pionChargedTrackLength2; 
      caseName = "pion-/pion+";
      break;
    }
    case 3 : {
      nn = numTrackProton;
      sum  = protonTrackLength;
      sum2 = protonTrackLength2; 
      caseName = "proton";
      break;
    }
    case 4 : {
      nn = numTrackGamma;
      sum  = gammaTrackLength;
      sum2 = gammaTrackLength2; 
      caseName = "gamma";
      break;
    }
    case 5 : {
      nn = numTrackPi0;
      sum  = pion0TrackLength;
      sum2 = pion0TrackLength2; 
      caseName = "pion0";
      break;
    }
    case 6 : {
      nn = numTrackNeutron;
      sum  = neutronTrackLength;
      sum2 = neutronTrackLength2; 
      caseName = "neutron";
      break;
    }
    }
    mu = mu_sigma = 0.0;
    if ( nn > 0.0 ) {
      mu = sum / nn;
      if ( nn <= 1.0 ) nn = 2.0;
      sigma = std::sqrt( std::abs( ( sum2 - sum*sum/nn ) ) / (nn - 1.0) );
      mu_sigma = sigma / std::sqrt( nn );
    }
    G4cout << "\t" << caseName << " : " << mu << " +/- " << mu_sigma << G4endl;
  }

  // Print information about step length and number of steps
  G4cout << G4endl << " Average number of steps and step LENGTH [mm]" << G4endl;
  for ( int iCase = 0; iCase < 7; iCase++ ) {
    if ( iCase == 1  ||  iCase == 5 ) continue; // Exclude muons and pi0s
    std::string caseName;
    G4double numerator_a = 0.0;
    G4double numerator_b = 0.0;
    G4double numerator_a2 = 0.0;
    G4double numerator_b2 = 0.0;
    G4double denominator_a = 0.0;
    G4double denominator_b = 0.0;
    G4double denominator_a2 = 0.0;
    G4double denominator_b2 = 0.0;
    G4double trackLength = 0.0;
    G4double trackLength2 = 0.0;
    switch ( iCase ) {
    case 0 : {
      numerator_a = numStepElectron;
      numerator_b = numStepPositron;
      denominator_a = numTrackElectron;
      denominator_b = numTrackPositron;
      numerator_a2 = numStepElectron2;
      numerator_b2 = numStepPositron2;
      denominator_a2 = numTrackElectron2;
      denominator_b2 = numTrackPositron2;
      trackLength = electronTrackLength;
      trackLength2 = electronTrackLength2;
      caseName = "electron/positron";
      break;
    }
    case 1 : {
      numerator_a = numStepMuMinus;
      numerator_b = numStepMuPlus;
      denominator_a = numTrackMuMinus;
      denominator_b = numTrackMuPlus;
      numerator_a2 = numStepMuMinus2;
      numerator_b2 = numStepMuPlus2;
      denominator_a2 = numTrackMuMinus2;
      denominator_b2 = numTrackMuPlus2;
      trackLength = muonTrackLength;
      trackLength2 = muonTrackLength2;
      caseName = "muon-/muon+";
     break;
    }
    case 2 : {
      numerator_a = numStepPiMinus;
      numerator_b = numStepPiPlus;
      denominator_a = numTrackPiMinus;
      denominator_b = numTrackPiPlus;
      numerator_a2 = numStepPiMinus2;
      numerator_b2 = numStepPiPlus2;
      denominator_a2 = numTrackPiMinus2;
      denominator_b2 = numTrackPiPlus2;
      trackLength = pionChargedTrackLength;
      trackLength2 = pionChargedTrackLength2;
      caseName = "pion-/pion+";
      break;
    }
    case 3 : {
      numerator_a = numStepProton;
      numerator_b = 0.0;
      denominator_a = numTrackProton;
      denominator_b = 0.0;
      numerator_a2 = numStepProton2;
      numerator_b2 = 0.0;
      denominator_a2 = numTrackProton2;
      denominator_b2 = 0.0;
      trackLength = protonTrackLength;
      trackLength2 = protonTrackLength2;
      caseName = "proton";
      break;
    }
    case 4 : {
      numerator_a = numStepGamma;
      numerator_b = 0.0;
      denominator_a = numTrackGamma;
      denominator_b = 0.0;
      numerator_a2 = numStepGamma2;
      numerator_b2 = 0.0;
      denominator_a2 = numTrackGamma2;
      denominator_b2 = 0.0;
      trackLength = gammaTrackLength;
      trackLength2 = gammaTrackLength2;
      caseName = "gamma";
      break;
    }
    case 5 : {
      numerator_a = numStepPi0;
      numerator_b = 0.0;
      denominator_a = numTrackPi0;
      denominator_b = 0.0;
      numerator_a2 = numStepPi02;
      numerator_b2 = 0.0;
      denominator_a2 = numTrackPi02;
      denominator_b2 = 0.0;
      trackLength = pion0TrackLength;
      trackLength2 = pion0TrackLength2;
      caseName = "pion0";
      break;
    }
    case 6 : {
      numerator_a = numStepNeutron;
      numerator_b = 0.0;
      denominator_a = numTrackNeutron;
      denominator_b = 0.0;
      numerator_a2 = numStepNeutron2;
      numerator_b2 = 0.0;
      denominator_a2 = numTrackNeutron2;
      denominator_b2 = 0.0;
      trackLength = neutronTrackLength;
      trackLength2 = neutronTrackLength2;
      caseName = "neutron";
      break;
    }
    }
    G4double sigma2_numerator_a = 
      ( std::abs( ( numerator_a2 - numerator_a*numerator_a/n ) ) / (n - 1.0) ) / n ;
    G4double sigma2_numerator_b = 
      ( std::abs( ( numerator_b2 - numerator_b*numerator_b/n ) ) / (n - 1.0) ) / n ;
    G4double sigma2_denominator_a = 
      ( std::abs( ( denominator_a2 - denominator_a*denominator_a/n ) ) / (n - 1.0) ) / n;
    G4double sigma2_denominator_b = 
      ( std::abs( ( denominator_b2 - denominator_b*denominator_b/n ) ) / (n - 1.0) ) / n;
    numerator_a /= n;
    numerator_b /= n;
    nn = denominator_a + denominator_b;
    if ( nn <= 1.0 ) nn = 2;
    denominator_a /= n;
    denominator_b /= n;
    G4double mu_numSteps = 0.0;
    G4double mu_numSteps_sigma = 0.0;
    if ( ( denominator_a + denominator_b ) > 0.0 ) {
      mu_numSteps = ( numerator_a + numerator_b ) / ( denominator_a + denominator_b );
      if ( ( numerator_a + numerator_b ) > 0.0 ) {
	mu_numSteps_sigma = mu_numSteps * 
	  std::sqrt( ( sigma2_numerator_a + sigma2_numerator_b ) / 
		     ( ( numerator_a + numerator_b ) * 
		       ( numerator_a + numerator_b ) )
		     +
		     ( sigma2_denominator_a + sigma2_denominator_b ) / 
		     ( ( denominator_a + denominator_b ) * 
		       ( denominator_a + denominator_b ) ) );
      }
    }
    G4cout << "\t" << caseName << " : numSteps = " 
	   << mu_numSteps << " +/- " << mu_numSteps_sigma << G4endl;

    G4double sigma2_trackLength =  
      ( std::abs( ( trackLength2 - trackLength*trackLength/nn ) ) / (nn - 1.0) ) / nn;
    trackLength /= nn;
    G4double mu_stepLength = 0.0;
    G4double mu_stepLength_sigma = 0.0;
    if ( mu_numSteps > 0.0 ) {
      mu_stepLength = trackLength / mu_numSteps;
      mu_stepLength_sigma = mu_stepLength * 
	std::sqrt( sigma2_trackLength / 
		   ( trackLength * trackLength ) 
		   +
		   ( mu_numSteps_sigma * mu_numSteps_sigma ) / 
		   ( mu_numSteps * mu_numSteps ) );
    }
    G4cout << "\t \t stepLength = " << mu_stepLength << " +/- " 
	   << mu_stepLength_sigma << G4endl;
  }

  // Print information about exiting kinetic energy.
  sum  = kinEnergyExiting; 
  sum2 = kinEnergyExiting2;
  mu = sum / n;
  sigma = std::sqrt( std::abs( ( sum2 - sum*sum/n ) ) / (n - 1.0) );
  mu_sigma = sigma / std::sqrt( n );
  G4cout << G4endl << " Average exiting Kinetic Energy = " 
         << mu << " +/- " << mu_sigma << " MeV " << G4endl;
  if ( mu > 1.0E-06 ) {
    G4double mu_tot = mu;
    G4double mu_tot_sigma = mu_sigma;
    for ( int iCase = 0; iCase < 6; iCase++ ) {
      std::string caseName;
      switch ( iCase ) {
      case 0 : {
	sum  = kinEnergyExitingGammas; 
	sum2 = kinEnergyExitingGammas2; 
        caseName = "Gammas";
	break;
      }
      case 1 : {
	sum  = kinEnergyExitingNeutrons; 
	sum2 = kinEnergyExitingNeutrons2; 
        caseName = "Neutrons";
	break;
      }
      case 2 : {
	sum  = kinEnergyExitingNeutrinos; 
	sum2 = kinEnergyExitingNeutrinos2; 
        caseName = "Neutrinos";
	break;
      }
      case 3 : {
	sum  = kinEnergyExitingMuons; 
	sum2 = kinEnergyExitingMuons2; 
        caseName = "Muons";
	break;
      }
      case 4 : {
	sum  = kinEnergyExitingElectrons; 
	sum2 = kinEnergyExitingElectrons2; 
        caseName = "Electrons";
	break;
      }
      case 5 : {
	sum  = kinEnergyExitingOthers; 
	sum2 = kinEnergyExitingOthers2; 
        caseName = "Others";
	break;
      }
      }
      mu = sum / n;
      sigma = std::sqrt( std::abs( ( sum2 - sum*sum/n ) ) / (n - 1.0) );
      mu_sigma = sigma / std::sqrt( n );
      G4double ratio = mu / mu_tot ;
      G4double ratio_sigma = 0.0;
      if ( mu > 1.0E-06 ) {
	ratio_sigma = ratio * 
	  std::sqrt( ( mu_sigma / mu ) * ( mu_sigma / mu ) +
		     ( mu_tot_sigma / mu_tot ) * ( mu_tot_sigma / mu_tot ) );
      }
      G4cout << "\t fraction due to " << caseName << " = " 
	     << 100.0 * ratio << " +/- " << 100.0 * ratio_sigma << " %" << G4endl;
    }

    sum  = numExiting;
    sum2 = numExiting2;
    mu = sum / n;
    sigma = std::sqrt( std::abs( ( sum2 - sum*sum/n ) ) / (n - 1.0) );
    mu_sigma = sigma / std::sqrt( n );
    G4cout << "\t number of exiting particles = " 
	   << mu << " +/- " << mu_sigma << G4endl;
    if ( mu > 1.0E-06 ) {
      mu_tot = mu;
      mu_tot_sigma = mu_sigma;
      for ( int iCase = 0; iCase < 6; iCase++ ) {
	std::string caseName;
	switch ( iCase ) {
	case 0 : {
	  sum  = numExitingGammas; 
	  sum2 = numExitingGammas2; 
	  caseName = "Gammas";
	  break;
	}
	case 1 : {
	  sum  = numExitingNeutrons; 
	  sum2 = numExitingNeutrons2; 
	  caseName = "Neutrons";
	  break;
	}
	case 2 : {
	  sum  = numExitingNeutrinos; 
	  sum2 = numExitingNeutrinos2; 
	  caseName = "Neutrinos";
	  break;
	}
	case 3 : {
	  sum  = numExitingMuons; 
	  sum2 = numExitingMuons2; 
	  caseName = "Muons";
	  break;
	}
	case 4 : {
	  sum  = numExitingElectrons; 
	  sum2 = numExitingElectrons2; 
	  caseName = "Electrons";
	  break;
	}
	case 5 : {
	  sum  = numExitingOthers; 
	  sum2 = numExitingOthers2; 
	  caseName = "Others";
	  break;
	}
	}
	mu = sum / n;
	sigma = std::sqrt( std::abs( ( sum2 - sum*sum/n ) ) / (n - 1.0) );
	mu_sigma = sigma / std::sqrt( n );
	G4double ratio = mu / mu_tot ;
	G4double ratio_sigma = 0.0;
	if ( mu > 1.0E-06 ) {
	  ratio_sigma = ratio * 
	    std::sqrt( ( mu_sigma / mu ) * ( mu_sigma / mu ) +
		       ( mu_tot_sigma / mu_tot ) * ( mu_tot_sigma / mu_tot ) );
	}
	G4cout << "\t number of exiting " << caseName << " = " 
	       << mu << " +/- " << mu_sigma << "  ("
	       << 100.0 * ratio 
	//     << " +/- " << 100.0 * ratio_sigma 
	       << " %)" << G4endl;
      }
    }
  }

}


void StatAccepTestAnalysis::setIsHistogramOn( const bool choice ) {
  isHistogramOn = choice;
  G4cout << "  --->  StatAccepTestAnalysis::isHistogramOn = " << isHistogramOn 
	 << "  <---" << G4endl << G4endl; 
}


