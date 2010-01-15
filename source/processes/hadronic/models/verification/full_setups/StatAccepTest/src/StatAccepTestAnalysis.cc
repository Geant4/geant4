#include "StatAccepTestAnalysis.hh"
#include <string>
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4VProcess.hh"
#include "G4ProcessType.hh"
#include "G4HadronicProcessType.hh"

#ifdef G4ANALYSIS_USEROOT
#include <TTree.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2D.h>
#endif

//***LOOKHERE***
bool StatAccepTestAnalysis::isHistogramOn = true;
bool StatAccepTestAnalysis::is2DHistogramStepLvsEOn = true;
bool StatAccepTestAnalysis::isHistogramSpectrumUnweightedOn = true;
bool StatAccepTestAnalysis::isHistogramSpectrumWeightedOn = true;
bool StatAccepTestAnalysis::isCountingProcessesOn = false;  
bool StatAccepTestAnalysis::isMapParticleNamesOn = false;  
bool StatAccepTestAnalysis::isMapInfoAboutTrackOn = false;  

G4double StatAccepTestAnalysis::infParticleEkin_electron = -999.9;
G4double StatAccepTestAnalysis::supParticleEkin_electron = 1.0E+30; 
G4double StatAccepTestAnalysis::infParticleEkin_muon     = -999.9;     
G4double StatAccepTestAnalysis::supParticleEkin_muon     = 1.0E+30;     
G4double StatAccepTestAnalysis::infParticleEkin_pion     = -999.9;     
G4double StatAccepTestAnalysis::supParticleEkin_pion     = 1.0E+30; 
G4double StatAccepTestAnalysis::infParticleEkin_kaon     = -999.9;     
G4double StatAccepTestAnalysis::supParticleEkin_kaon     = 1.0E+30; 
G4double StatAccepTestAnalysis::infParticleEkin_proton   = -999.9;   
G4double StatAccepTestAnalysis::supParticleEkin_proton   = 1.0E+30;
G4double StatAccepTestAnalysis::infParticleEkin_nuclei   = -999.9;   
G4double StatAccepTestAnalysis::supParticleEkin_nuclei   = 1.0E+30;

G4double StatAccepTestAnalysis::infTimeWindow = -999.9;  // in [nanosec]
G4double StatAccepTestAnalysis::supTimeWindow = 1.0E+30; // in [nanosec]

StatAccepTestAnalysis* StatAccepTestAnalysis::instance = 0;


StatAccepTestAnalysis::StatAccepTestAnalysis() : 
  numberOfEvents( 0 ), numberOfReplicas( 0 ), 
  numberOfReadoutLayers( 0 ), numberOfActiveLayersPerReadoutLayer( 1 ), 
  numberOfRadiusBins( 0 ), radiusBin( 0.0 )
  //
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
#ifdef G4ANALYSIS_USEROOT
  , outFile( "ntuple.root", "recreate" ), tree( 0 )
  , plongitudinalProfile( 0 ), ptransverseProfile( 0 )
  , h_longitudinalProfile( 0 ), h_transverseProfile( 0 )
#endif
{
#ifdef G4ANALYSIS_USEROOT
    if ( outFile.IsZombie() ) {
      G4cerr << "*ERROR* Cannot open ntuple.root file." 
	     << "        Ntuple and histograms will not be saved." << G4endl;
    }
    tree = new TTree( "SimplifiedCalorimeter", "SimplifiedCalorimeter Tree" );
    tree->Branch( "ID",       &primaryParticleId,     "id/I" );
    tree->Branch( "E",        &beamEnergy,            "E/D" );
    tree->Branch( "nLayers",  &numberOfReadoutLayers, "nLayers/I" );
    tree->Branch( "nBinR",    &numberOfRadiusBins,    "nBinR/I" );
    tree->Branch( "EDEP_ACT", &totEDepActiveLayer,    "EDEP_ACT/F" );
    tree->Branch( "EDEP_CAL", &totEDepCalorimeter,    "EDEP_CAL/F");
    // This needs an explanation: longitudinal and transverse are 
    // of type  std::vector< G4double >  therefore we need to create
    // the branch using the syntax for objects: 
    //  TTree::Branch( const char** barnchName, 
    //                 const char** className, 
    //                 void ** object, ... ) 
    // so (see the third argument) we need a pointer to the 
    // pointer of the class.
    plongitudinalProfile = &longitudinalProfile;
    tree->Branch( "L", "std::vector<G4double>", &plongitudinalProfile );
    ptransverseProfile = &transverseProfile;
    tree->Branch( "R", "std::vector<G4double>", &ptransverseProfile );
#endif // G4ANALYSIS_USEROOT
}


StatAccepTestAnalysis::~StatAccepTestAnalysis() {}


#ifdef G4ANALYSIS_USEROOT
std::pair< StatAccepTestAnalysis::spectra_t, std::string >
StatAccepTestAnalysis::getKey( const G4int& pId, 
			       const SpectraERange& etype, 
			       const unsigned int& layer, 
			       bool combinedPions ) const {
  std::string p, t;
  if ( pId == electronId  ||  pId == positronId ) { 
    p = "e"; 
    t = "e^{-} e^{+}";
  } else if ( combinedPions  &&  ( pId == pionPlusId   || 
				   pId == pionMinusId  || 
				   pId == pionZeroId ) ) { 
    p = "pi"; 
    t = "#pi^{+-0}";
  } else if ( pId == protonId  ||  pId==antiProtonId ) { 
    p = "p"; 
    t = "proton";
  } else if ( pId == neutronId  ||  pId == antiNeutronId ) { 
    p = "n"; 
    t = "neutron";
  } else if (  pId == gammaId ) { 
    p = "gamma"; 
    t = "#gamma";
  } else if ( pId == pionPlusId ) { 
    p = "piPlus"; 
    t = "#pi^{+}"; 
  } else if ( pId == pionMinusId ) { 
    p = "piMinus"; 
    t = "#pi^{-}"; 
  } else {
    return std::pair< spectra_t, std::string >( spectra_t( "" ), "" );
  }
  std::map< SpectraERange, std::string > m;
  m[ E100 ] = "E100";
  m[ E10 ]  = "E10";
  m[ E1 ]   = "E1";
  m[ E01 ]  = "E01";
  m[ ELOG ] = "ELOG";
  char buf[30];
  sprintf( buf, "h_%s_%s_%d", p.c_str(), m[ etype ].c_str(), layer );
  char title[256];
  sprintf( title, "Spectrum %s , Layer %d ", t.c_str(), layer );
  return std::make_pair< spectra_t, std::string >( spectra_t( buf ), std::string( title ) );
}
#endif // G4ANALYSIS_USEROOT


void StatAccepTestAnalysis::close() {
#ifdef G4ANALYSIS_USEROOT
#ifdef G4VERBOSE
  if ( tree ) {
    tree->Print();
  }
#endif
  int nbytes = 0;
  if ( ! outFile.IsZombie() ) {
    // We used the AutoSave feature to save the ntuple tree periodically.
    // We thus use option kWriteDelete to ask ROOT to delete from the TFile
    // all previous (if any) copies of the tree from the file.
    // We also use (the slow) "delete after write" option to increase security
    // otherwise, if, for any reason, Write() would fail we would lose the tree.
    nbytes = outFile.Write( 0, TObject::kWriteDelete );
  }
  if ( tree ) { // The tree will be deleted by calling TFile::Close()
    tree = NULL;
  }
  // Actually histograms will be saved in the file, thus they will be deleted
  // automatically when the file is closed, but this is harmless anyway
  // (the TFile will recognize the histos have been already deleted).
  // Note that if the TH1::AddDirectory(false) method has been called
  // histograms are not added to the file and thus not deleted automatically.
  // This will clean also them (provided that have been booked with the
  // book<>(...) method
  histoList.Delete();
  if ( ! outFile.IsZombie() ) { // Safely close the file
    G4cout << "Wrote file " << outFile.GetName() << " (" 
	   << nbytes/(1.e+6) << "MB)" << G4endl;
    outFile.Close();
  }
  G4cout << "Spectra do not include the following particles:" << G4endl
         << "=====================" << G4endl;
  std::for_each( unknownCounter.begin(), unknownCounter.end(), printUnknownMap() );
  G4cout << "=====================" << G4endl;
#endif // G4ANALYSIS_USEROOT
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

#ifdef G4ANALYSIS_USEROOT
  tree->Reset();
#endif

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

  vecEvis_electron.clear();
  vecEvis_muon.clear();
  vecEvis_pion.clear();
  vecEvis_kaon.clear();
  vecEvis_proton.clear();
  vecEvis_nuclei.clear();

  vecEvis_no_electron.clear();
  vecEvis_no_muon.clear();
  vecEvis_no_pion.clear();
  vecEvis_no_kaon.clear();
  vecEvis_no_proton.clear();
  vecEvis_no_nuclei.clear();

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

  sumEdepAct_nuclei  = 0.0;
  sumEdepAct_nuclei2 = 0.0;
  sumEdepTot_nuclei  = 0.0;
  sumEdepTot_nuclei2 = 0.0;
  sumL_nuclei.clear();
  sumL_nuclei2.clear();
  longitudinalProfile_nuclei.clear();
  for ( int layer = 0; layer < numberOfReadoutLayers; layer++ ) {
    longitudinalProfile_nuclei.push_back( 0.0 );
    sumL_nuclei.push_back( 0.0 );
    sumL_nuclei2.push_back( 0.0 );
  }
  sumR_nuclei.clear();
  sumR_nuclei2.clear();
  transverseProfile_nuclei.clear();
  for ( int ir = 0; ir < numberOfRadiusBins; ir++ ) {
    transverseProfile_nuclei.push_back( 0.0 );
    sumR_nuclei.push_back( 0.0 );
    sumR_nuclei2.push_back( 0.0 );
  } 

  sumEdepAct_no_electron  = 0.0;
  sumEdepAct_no_electron2 = 0.0;
  sumEdepAct_no_muon  = 0.0;
  sumEdepAct_no_muon2 = 0.0;
  sumEdepAct_no_pion  = 0.0;
  sumEdepAct_no_pion2 = 0.0;
  sumEdepAct_no_kaon  = 0.0;
  sumEdepAct_no_kaon2 = 0.0;
  sumEdepAct_no_proton  = 0.0;
  sumEdepAct_no_proton2 = 0.0;
  sumEdepAct_no_nuclei  = 0.0;
  sumEdepAct_no_nuclei2 = 0.0;

  numStep = 0.0;
  numStepPositive = numStepNeutral = numStepNegative = 0.0;
  numStepNucleus = numStepPDGCodeUnrecognized = 0.0;
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
  numStepNucleus2 = numStepPDGCodeUnrecognized2 = 0.0;
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
  numTrackNucleus = numTrackPDGCodeUnrecognized = 0.0;
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
  numTrackNucleus2 = numTrackPDGCodeUnrecognized2 = 0.0;
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
  neutronTrackLength = 0.0;

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

  numInelasticProcesses = 0;
  numProtonInelasticProcesses = 0;
  numAntiProtonInelasticProcesses = 0;
  numNeutronInelasticProcesses = 0;
  numAntiNeutronInelasticProcesses = 0;
  numPionPlusInelasticProcesses = 0;
  numPionMinusInelasticProcesses = 0;
  numKaonPlusInelasticProcesses = 0;
  numKaonMinusInelasticProcesses = 0;
  numKaonNeutralInelasticProcesses = 0;
  numLambdaInelasticProcesses = 0;
  numAntiLambdaInelasticProcesses = 0;
  numSigmaMinusInelasticProcesses = 0;
  numMuonMinusInelasticProcesses = 0;
  numOtherInelasticProcesses = 0;
  numStoppingAtRestProcesses = 0;
  numAntiProtonStoppingAtRestProcesses = 0;
  numNeutronStoppingAtRestProcesses = 0;
  numAntiNeutronStoppingAtRestProcesses = 0;
  numPionMinusStoppingAtRestProcesses = 0;
  numKaonMinusStoppingAtRestProcesses = 0;
  numKaonNeutralStoppingAtRestProcesses = 0;
  numLambdaStoppingAtRestProcesses = 0;
  numAntiLambdaStoppingAtRestProcesses = 0;
  numSigmaMinusStoppingAtRestProcesses = 0;
  numMuonMinusStoppingAtRestProcesses = 0;
  numOtherStoppingAtRestProcesses = 0;

  mapParticleNames.clear();

  eventTimeSet.clear();

  mapInfoAboutTrack.clear();
  mapInfoAboutVertex.clear();

  suggestedTrackPosition = mapInfoAboutTrack.begin();
  suggestedVertexPosition = mapInfoAboutVertex.begin();

  countNumberOfTracks = countNumberOfVertices = 
    countNumberOfElectromagneticVertices = 
    countNumberOfPhotoleptonHadronVertices =
    countNumberOfDecayVertices =
    countNumberOfHadronicVertices = 
    countNumberOfHadronElasticVertices_notFromNeutrons =
    countNumberOfHadronElasticVertices_fromNeutrons =
    countNumberOfHadronInelasticVertices_notFromNeutrons =
    countNumberOfHadronInelasticVertices_fromNeutrons =
    countNumberOfCaptureVertices_notFromNeutrons =
    countNumberOfCaptureVertices_fromNeutrons =
    countNumberOfFissionVertices_notFromNeutrons =
    countNumberOfFissionVertices_fromNeutrons =
    countNumberOfAtRestVertices_notFromNeutrons =
    countNumberOfAtRestVertices_fromNeutrons =
    countNumberOfChargeExchangeVertices_notFromNeutrons =
    countNumberOfChargeExchangeVertices_fromNeutrons = 0;

  sumEvis_em = sumEtot_em = sumEvis_p = sumEtot_p = sumEvis_pi = sumEtot_pi =
    sumEvis_ion = sumEtot_ion = 0.0;
  sumEvis_from1stInterac_pi0 = sumEtot_from1stInterac_pi0 = 
    sumEvis_from1stInterac_pip = sumEtot_from1stInterac_pip =
    sumEvis_from1stInterac_pim = sumEtot_from1stInterac_pim =
    sumEvis_from1stInterac_p = sumEtot_from1stInterac_p =
    sumEvis_from1stInterac_n = sumEtot_from1stInterac_n =
    sumEvis_from1stInterac_lightion = sumEtot_from1stInterac_lightion = 0.0;
  sumEvis_closest_pi0 = sumEtot_closest_pi0 =
    sumEvis_closest_pip = sumEtot_closest_pip =
    sumEvis_closest_pim = sumEtot_closest_pim =
    sumEvis_closest_p = sumEtot_closest_p =
    sumEvis_closest_n = sumEtot_closest_n = 
    sumEvis_closest_lightion = sumEtot_closest_lightion = 0.0;
  sumEvis_em_from_pi0 = sumEtot_em_from_pi0 =
    sumEvis_em_from_pip = sumEtot_em_from_pip =
    sumEvis_em_from_pim = sumEtot_em_from_pim =
    sumEvis_em_from_p = sumEtot_em_from_p =
    sumEvis_em_from_n = sumEtot_em_from_n =
    sumEvis_em_from_lightion = sumEtot_em_from_lightion = 0.0;

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

#ifdef G4ANALYSIS_USEROOT
  if ( outFile.IsZombie() == kFALSE ) outFile.cd(); // Save histograms in file.
  histoList.Delete(); // This recursively cleans the list calling delete 
                      // for each element of the list.
  h_longitudinalProfile = book1D< TH1D >( "hLongProfile", 
					  "Longitudinal Shower Profile",
					  numberOfReadoutLayers,
					  -0.5,
					  1.0*numberOfReadoutLayers - 0.5,
					  "Layer",
					  "<E> (MeV)" );
  h_transverseProfile =   book1D< TH1D >( "hTransvProfile",
					  "Transverse Shower Profile",
					  numberOfRadiusBins,
					  -0.5, 
					  1.0*numberOfRadiusBins - 0.5,
					  "Radial Layer",
					  "<E> (MeV)" );
  const Int_t nbinx = 1000;
  const Float_t minx = -4.0;
  const Float_t maxx = 6.0;
  if ( ! outFile.IsZombie() ) {
    TDirectory* dir = outFile.mkdir( "Spectra" );
    if ( dir ) {
      dir->cd();
    } else {
      outFile.cd( "Spectra" );
    }
  }
  h_Spectrum[ "gamma" ]     = book1D< TH1F >( "hGammaSpectrum",
					      "Kinetic Energy for #gamma",
					      nbinx, minx, maxx,
					      "log10(E_{kin}) (MeV)", "Entries" );
  h_Spectrum[ "neutron" ]   = book1D< TH1F >( "hNeutronSpectrum",
					      "Kinetic Energy for neutron",
					      nbinx, minx, maxx,
					      "log10(E_{kin}) (MeV)", "Entries" );
  h_Spectrum[ "proton" ]    = book1D< TH1F >( "hProtonSpectrum",
					      "Kinetic Energy for proton",
					      nbinx, minx, maxx,
					      "log10(E_{kin}) (MeV)", "Entries" );
  h_Spectrum[ "pionZero" ]  = book1D< TH1F >( "hPionZeroSpectrum",
					      "Kinetic Energy for #pi^{0}",
					      nbinx, minx, maxx,
					      "log10(E_{kin}) (MeV)", "Entries" );
  h_Spectrum[ "pionMinus" ] = book1D< TH1F >( "hPionMinusSpectrum",
					      "Kinetic Energy for #pi^{-}",
					      nbinx, minx, maxx,
					      "log10(E_{kin}) (MeV)", "Entries" );
  h_Spectrum[ "pionPlus" ]  = book1D< TH1F >( "hPionPlusSpectrum",
					      "Kinetic Energy for #pi^{+}",
					      nbinx, minx, maxx,
					      "log10(E_{kin}) (MeV)", "Entries" );
  // Other histograms are optional.
  if ( isHistogramOn ) {
    if ( is2DHistogramStepLvsEOn ) {
      if ( ! outFile.IsZombie() ) {
	TDirectory* dir = outFile.mkdir( "StepEdepVsLength" );
	if ( dir ) {
	  dir->cd();
	} else {
	  outFile.cd( "StepEdepVsLength" );
	}
      }
      // Active
      h_stepEvsL[ "active" ] = 
	book2D< TH2D >( "hStepEVsL_active",
			"Step: Energy vs. Length,active",
			100, 0.0, 10.0, 
			100, 0.0, 10.0,
			"L (mm)", "E_{dep} (MeV)", "Entries" );
      h_stepEvsL[ "electron_active" ] = 
	book2D< TH2D >( "hStepEVsL_electron_active",
			"Step: Energy vs. Length,active, ELECTRON",
			100, 0.0, 10.0,
			100, 0.0, 10.0,
			"L (mm)", "E_{dep} (MeV)", "Entries" );
      h_stepEvsL[ "muon_active" ] = 
	book2D< TH2D >( "hStepEVsL_muon_active",
			"Step: Energy vs. Length,active, MUON",
			100, 0.0, 10.0,
			100, 0.0, 10.0,
			"L (mm)", "E_{dep} (MeV)", "Entries" );
      h_stepEvsL[ "pioncharged_active" ] =
	book2D< TH2D >( "hStepEVsL_pioncharged_active",
			"Step: Energy vs. Length,active, PION charged",
			100, 0.0, 10.0,
			100, 0.0, 10.0,
			"L (mm)", "E_{dep} (MeV)", "Entries" );
      h_stepEvsL[ "proton_active" ] =
	book2D< TH2D >( "hStepEVsL_proton_active",
			"Step: Energy vs. Length,active, PROTON",
			100, 0.0, 10.0, 
			100, 0.0, 10.0,
			"L (mm)", "E_{dep} (MeV)","Entries" );
      h_stepEvsL[ "gamma_active" ] = 
	book2D< TH2D >( "hStepEVsL_gamma_active",
			"Step: Energy vs. Length,active, GAMMA",
			100, 0.0, 10.0, 
			100, 0.0, 10.0,
			"L (mm)", "E_{dep} (MeV)", "Entries" );
      h_stepEvsL[ "neutron_active" ] =
	book2D< TH2D >( "hStepEVsL_neutron_active",
			"Step: Energy vs. Length,active, NEUTRON",
			100, 0.0, 10.0,
			100, 0.0, 10.0,
			"L (mm)", "E_{dep} (MeV)", "Entries" );
      h_stepEvsL[ "other_active" ] =
	book2D< TH2D >( "hStepEVsL_other_active",
			"Step: Energy vs. Length,active, OTHER",
			100, 0.0, 10.0,
			100, 0.0, 10.0,
			"L (mm)", "E_{dep} (MeV)", "Entries" );
      // Absorber
      h_stepEvsL[ "absorber" ] =
	book2D< TH2D >( "hStepEVsL_absorber",
			"Step: Energy vs. Length,absorber",
			100, 0.0, 10.0,
			100, 0.0, 10.0,
			"L (mm)", "E_{dep} (MeV)", "Entries" );
      h_stepEvsL[ "electron_absorber" ] = 
	book2D< TH2D >( "hStepEVsL_electron_absorber",
			"Step: Energy vs. Length,absorber, ELECTRON",
			100, 0.0, 10.0,
			100, 0.0, 10.0,
			"L (mm)", "E_{dep} (MeV)", "Entries" );
      h_stepEvsL[ "muon_absorber" ] = 
	book2D< TH2D >( "hStepEVsL_muon_absorber",
			"Step: Energy vs. Length,absorber, MUON",
			100, 0.0, 10.0,
			100, 0.0, 10.0,
			"L (mm)", "E_{dep} (MeV)", "Entries" );
      h_stepEvsL[ "pioncharged_absorber" ] =
	book2D< TH2D >( "hStepEVsL_pioncharged_absorber",
			"Step: Energy vs. Length,absorber, PION charged",
			100, 0.0, 10.0,
			100, 0.0, 10.0,
			"L (mm)", "E_{dep} (MeV)", "Entries" );
      h_stepEvsL[ "proton_absorber" ] =
	book2D< TH2D >( "hStepEVsL_proton_absorber",
			"Step: Energy vs. Length,absorber, PROTON",
			100, 0.0, 10.0,
			100, 0.0, 10.0,
			"L (mm)", "E_{dep} (MeV)", "Entries" );
      h_stepEvsL[ "gamma_absorber" ] = 
	book2D< TH2D >( "hStepEVsL_gamma_absorber", 
			"Step: Energy vs. Length,absorber, GAMMA",
			100, 0.0, 10.0,
			100, 0.0, 10.0,
			"L (mm)", "E_{dep} (MeV)", "Entries" );
      h_stepEvsL[ "neutron_absorber" ] = 
	book2D< TH2D >( "hStepEVsL_neutron_absorber",
			"Step: Energy vs. Length,absorber, NEUTRON",
			100, 0.0, 10.0,
			100, 0.0, 10.0,
			"L (mm)", "E_{dep} (MeV)", "Entries" );
      h_stepEvsL[ "other_absorber" ] =
	book2D< TH2D >( "hStepEVsL_other_absorber",
			"Step: Energy vs. Length,absorber, OTHER",
			100, 0.0, 10.0, 
			100, 0.0, 10.0,
			"L (mm)", "E_{dep} (MeV)", "Entries" );
    }
    // Kinetic spectra of some particles in some planes (normal
    // to the beam direction).
    // We evaluate the kinetic spectra in 10 active layers,
    // uniformely spread along the calorimeter.
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
    // The name of the histogram and the key in the container map codes 
    // the particle type, the energy binning and the layer (0-9)
    if ( isHistogramSpectrumUnweightedOn  ||  isHistogramSpectrumWeightedOn ) {
      if ( ! outFile.IsZombie() ) {
	TDirectory* dir = outFile.mkdir( "SpectraDetailed" );
	if ( dir ) {
	  dir->cd();
	} else {
	  outFile.cd( "SpectraDetailed" );
	}
      }
      G4int parts[ 6 ] = { electronId , pionMinusId , pionPlusId , 
			   protonId , neutronId , gammaId };
      SpectraERange allE[ 5 ] = { E100 , E10 , E1 , E01 , ELOG };
      const char* xtitle[ 5 ] = { "E_{kin} (GeV)" , "E_{kin} (GeV)" , "E_{kin} (GeV)" ,
				  "E_{kin} (GeV)" , "log10(E_{kin} (GeV))" };
      const char* ytitle[ 5 ] = { "#frac{dN}{dE_{kin}} (GeV^{-1})" , 
				  "#frac{dN}{dE_{kin}} (GeV^{-1})", 
				  "#frac{dN}{dE_{kin}} (GeV^{-1})" , 
				  "#frac{dN}{dE_{kin}} (GeV^{-1})" , 
				  "#frac{dN}{d(log10(E_{kin}))} (log10(GeV)^{-1})" };
      const char* ytitleW[ 5 ] = { "#frac{1}{P}*#frac{dN}{dE_{kin}} (GeV^{-2})" , 
				   "#frac{1}{P}*#frac{dN}{dE_{kin}} (GeV^{-2})", 
				   "#frac{1}{P}*#frac{dN}{dE_{kin}} (GeV^{-2})" , 
				   "#frac{1}{P}*#frac{dN}{dE_{kin}} (GeV^{-2})", 
				   "#frac{1}{E_{kin}P}*#frac{dN}{d(log10(E_{kin}))} ((GeV)^{-2}*log10(GeV)^{-1})" };
      float eminx[ 5 ] = { -100, -10, -1, -0.1, -4.0 };
      float emaxx[ 5 ] = { 100, 10, 1, 0.1, 6 };
      int nbinsx = 1000;
      for ( unsigned int layerIndex = 0 ; layerIndex < 10 ; ++layerIndex ) {
	for ( unsigned int et = 0 ; et < 5 ; ++et ) {
	  for ( unsigned ipart = 0 ; ipart < 6 ; ++ipart ) {
	    std::pair< spectra_t, std::string > key = 
	      getKey( parts[ ipart ], allE[ et ], layerIndex, false );
	    // First consider pions separately
	    if ( key.first == "" ) {
	      G4cerr << "Non valid key for: particle ID:" 
		     << parts[ ipart ] << " E=" << allE[ et ] 
		     << "and layer:" << layerIndex << G4endl;
	    } else {
	      if ( isHistogramSpectrumUnweightedOn ) {
		h_SpectraFlux[ key.first ] = 
		  book1D< TH1F >( key.first.c_str(), key.second.c_str(),
				  nbinsx, eminx[ et ], emaxx[ et ],
				  xtitle[ et ], ytitle[ et ] );
	      }
	      if ( isHistogramSpectrumWeightedOn ) {
		h_SpectraFluxWeight[ key.first ] = 
		  book1D< TH1F >( std::string( key.first + "W" ).c_str(), 
				  std::string( key.second + " (Weighted)" ).c_str(),
				  nbinsx, eminx[ et ], emaxx[ et ],
				  xtitle[ et ], ytitleW[ et ] );
	      }
	    }
	  } // End of loop on different particles types, 
            // we still need to handle the special case of "all pions together"
	  std::pair< spectra_t, std::string > key = 
	    getKey( pionMinusId, allE[ et ], layerIndex, true ); 
	  // Key for combined pions, just take one pion id
	  // G4cout << key.first << " " << key.second << G4endl;
	  if ( key.first == "" ) {
	    G4cerr << "Non valid key for: particle ID:" << pionMinusId
		   << " E=" << allE[ et ] << "and layer:" << layerIndex << G4endl;
	  } else {
	    if ( isHistogramSpectrumUnweightedOn )
	      h_SpectraFlux[ key.first ] = 
		book1D< TH1F >( key.first.c_str(), key.second.c_str(),
				nbinsx, eminx[ et ], emaxx[ et ],
				xtitle[ et ], ytitle[ et ] );
	    if ( isHistogramSpectrumWeightedOn )
	      h_SpectraFluxWeight[ key.first ] = 
		book1D< TH1F >( std::string( key.first + "W" ).c_str(), 
				std::string( key.second + " (Weighted)" ).c_str(),
				nbinsx, eminx[ et ], emaxx[ et ],
				xtitle[ et ], ytitleW[ et ] );
	  }
	}
      }
      
    }
  } // isHistogrammOn
#endif // G4ANALYSIS_USEROOT
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
#ifdef G4ANALYSIS_USEROOT
  if ( ! isHistogramOn  || 
       ( ! isHistogramSpectrumUnweightedOn  &&  
	 ! isHistogramSpectrumWeightedOn ) ) return;
  if ( numberOfReplicas < 10  ||
       layerNumber % (numberOfReplicas/10) != 0  ||
       layerNumber / (numberOfReplicas/10) >= 10 ) return;

  // To improve CPU performance, get particle properties only once 
  // and save them on variables.
  G4double particleMass = particleDef->GetPDGMass();
  G4int particleId = particleDef->GetPDGEncoding();

  G4int sampleLayer = layerNumber / (numberOfReplicas/10); 
  G4double p = std::abs( kinEnergy ); // Momentum for the weighting...
  if ( particleMass/MeV > 1.0e-06 ) {
    p = std::sqrt( std::abs( kinEnergy ) * 
		   ( std::abs( kinEnergy ) + 2.0 * particleMass ) );
    //G4cout << "  m=" << particleMass / MeV
    //       << "  Ekin=" << std::abs( kinEnergy ) / MeV 
    //	     <<  "  p=" << p / MeV << " MeV" << G4endl;       //***DEBUG***
  }
  G4double p_inGeV = p / GeV;
  //float weight = ( eidx < 4 )  ?  1.0/p_inGeV  :  1.0 / (p*std::abs(kinEnergy));
  static const float cw = 4*pi*std::log(10.0);
  float w_lin = 1.0 / (4*pi*p_inGeV);
  float w_log = 1.0 / (cw*p_inGeV*std::abs(kinEnergy));
  SpectraERange etype[ 5 ] = { E100 , E10 , E1 , E01 , ELOG };
  // We need to find the key of the histogram for all Energy types for particleId.
  // The complication is that we need to fill two histograms for pions: 
  // the "all together" and the pi+/pi- separately.
  for ( int eidx = 0 ; eidx <5 ; ++eidx ) {
    if ( eidx == 4  &&  kinEnergy <= 0 ) continue; 
    // Cannot fill Log histogram if energy is <=0 .
    spectra_t key = getKey( particleId, etype[ eidx ], sampleLayer, true ).first;
    // First all pions together
    histoSpectraMapIt_t uWe = h_SpectraFlux.find( key );
    histoSpectraMapIt_t Wei = h_SpectraFluxWeight.find( key );
    if ( uWe == h_SpectraFlux.end() ||  Wei == h_SpectraFluxWeight.end() ) {
      // Log this information once, it is a particle for which no histogram 
      // is defined (neutrino, Kaon for example).
      if ( eidx == 0  &&  sampleLayer == 0 ) {
	char buf[100];
	sprintf( buf, "particleId=%d", particleId );
	std::string theId( buf );
	unknownCounterIt_t it = unknownCounter.find( theId );
	if ( it == unknownCounter.end() ) {  // Element not present
	  unknownCounter[ theId ] = 
	    std::make_pair< unsigned int, G4double >( 1, std::abs( kinEnergy/MeV) );
	} else {
	  std::pair< unsigned int, G4double >& element = (*it).second;
	  element.first  += 1; // Increase counter
	  element.second += std::abs( kinEnergy/MeV );
	}
      }
      continue;
    }
    float weight = ( eidx < 4 )  ?  w_lin :  w_log ; 
    // Correct weight for lorentz invariant: EdN/d3p
    TH1* huWe = (*uWe).second;
    TH1* hWei = (*Wei).second;
    if ( isHistogramSpectrumUnweightedOn ) huWe->Fill( kinEnergy/GeV ); 
    if ( isHistogramSpectrumWeightedOn )   hWei->Fill( kinEnergy/GeV, weight );
    // In case the particle is a pi+ or a pi- fill two additional sets 
    // of spectra for pi+ and pi- alone.
    if ( particleId == pionMinusId  ||  particleId == pionPlusId ) {
      key = getKey( particleId, etype[ eidx ], sampleLayer, false ).first;
      // All pions together.
      uWe = h_SpectraFlux.find( key );
      Wei = h_SpectraFluxWeight.find( key );
      assert ( uWe != h_SpectraFlux.end() );
      assert ( Wei != h_SpectraFluxWeight.end() );
      huWe = (*uWe).second;
      hWei = (*Wei).second;
      if ( isHistogramSpectrumUnweightedOn ) huWe->Fill( kinEnergy/GeV );
      if ( isHistogramSpectrumWeightedOn )   hWei->Fill( kinEnergy/GeV, weight );
    }
  }
#endif // G4ANALYSIS_USEROOT
}


void StatAccepTestAnalysis::fillNtuple( float incidentParticleId, 
					float incidentParticleEnergy, 
					float totalEnergyDepositedInActiveLayers,
					float totalEnergyDepositedInCalorimeter ) {
  primaryParticleId = static_cast< G4int >( incidentParticleId );
  beamEnergy = incidentParticleEnergy;

  // For meson primary beam particles, consider the total energy.
  if ( primaryParticleId == -211  ||  primaryParticleId == +211 ) {  // pi- or pi+
    beamEnergy += G4PionPlus::PionPlusDefinition()->GetPDGMass();
  } else if ( primaryParticleId == -321  ||  primaryParticleId == +321 ) {  // K- or K+
    beamEnergy += G4KaonPlus::KaonPlusDefinition()->GetPDGMass();
  } else if ( primaryParticleId == 130 ) { // K0L
    beamEnergy += G4KaonZeroLong::KaonZeroLongDefinition()->GetPDGMass(); 
  } else if ( primaryParticleId == 310 ) { // K0S
    beamEnergy += G4KaonZeroShort::KaonZeroShortDefinition()->GetPDGMass(); 
  }

  if ( totalEnergyDepositedInCalorimeter - beamEnergy > 0.001*MeV ) {
    G4cout << "\t ***ENERGY-NON-CONSERVATION*** " 
	   << totalEnergyDepositedInCalorimeter / MeV << " MeV" << G4endl;
    countEnergyNonConservation++;
  }
  if ( totalEnergyDepositedInCalorimeter > maxEdepTot ) {
    maxEdepTot = totalEnergyDepositedInCalorimeter;
  }

  //G4cout << " StatAccepTestAnalysis::fillNtuple : DEBUG Info " << G4endl
  //       << "\t incidentParticleId = " << incidentParticleId << G4endl
  //       << "\t incidentParticleEnergy = " << incidentParticleEnergy << G4endl
  //       << "\t totalEnergyDepositedInActiveLayers = " 
  //       << totalEnergyDepositedInActiveLayers << G4endl
  //       << "\t totalEnergyDepositedInCalorimeter = " 
  //       << totalEnergyDepositedInCalorimeter << G4endl;       //***DEBUG***

#ifdef G4ANALYSIS_USEROOT
  if ( tree ) {
    // Fill copy variables with data.
    totEDepActiveLayer = (Float_t) totalEnergyDepositedInActiveLayers;
    totEDepCalorimeter = (Float_t) totalEnergyDepositedInCalorimeter;
    // Verified with a cout that the content of the produced ntuple is identical
    // to the value of the longitudinalProfile variable.
    //for ( int i= 0; i < numberOfRadiusBins; ++i ) 
    //  G4cout << transverseProfile[ i ] << G4endl;
    tree->Fill();
    
    // This method is called at each event, so it is useful to commit
    // the tree from time to time, for instance every 1000 events, in
    // such a way that it is possible to see the ntuple while the job 
    // is running, for instance for a quick check that it makes sense. 
    // Or, if the job crashes after many events, we have anyhow some data
    // already stored in the ntuple to be checked.
    if ( numberOfEvents > 0  and  numberOfEvents % 1000 == 0 ) {
      tree->AutoSave( "SaveSelf" ); 
      // SaveSelf allows another process to concurrently open the root file.
    }
  }
#endif // G4ANALYSIS_USEROOT

  sumEdepAct  += totalEnergyDepositedInActiveLayers;
  sumEdepAct2 += 
    totalEnergyDepositedInActiveLayers * totalEnergyDepositedInActiveLayers;
  sumEdepTot  += totalEnergyDepositedInCalorimeter;
  sumEdepTot2 += 
    totalEnergyDepositedInCalorimeter * totalEnergyDepositedInCalorimeter;
  for ( int iLayer = 0; iLayer < numberOfReadoutLayers; iLayer++ ) {
    //G4cout << " StatAccepTestAnalysis::fillNtuple : DEBUG Info " << G4endl
    //       << "\t Longitudinal profile: layer = " << iLayer
    //       << "   energy = " << longitudinalProfile[ iLayer ] / MeV 
    //       << " MeV " << G4endl;                                 //***DEBUG***
    sumL[ iLayer ]  += longitudinalProfile[ iLayer ];
    sumL2[ iLayer ] += longitudinalProfile[ iLayer ] * longitudinalProfile[ iLayer ];  
    longitudinalProfile[ iLayer ] = 0.0;
  }
  for ( int iBinR = 0; iBinR < numberOfRadiusBins; iBinR++ ) {
    //G4cout << " StatAccepTestAnalysis::fillNtuple : DEBUG Info " << G4endl
    //       << "\t Transverse profile: iBinRadius = " << iBinR 
    //       << "   energy = " << transverseProfile[ iBinR ] / MeV 
    //       << " MeV " << G4endl;                                 //***DEBUG***
    sumR[ iBinR ]  += transverseProfile[ iBinR ];
    sumR2[ iBinR ] += transverseProfile[ iBinR ] * transverseProfile[ iBinR ];  
    transverseProfile[ iBinR ] = 0.0;
  }

  numberOfEvents++;

  // Store information of the visible energy, for later computing
  // of the energy resolution.
  vecEvis.push_back( totalEnergyDepositedInActiveLayers );

}


void StatAccepTestAnalysis::
fillShowerProfile( G4int replica, const G4double radius, const G4double edep,
                   const G4int particlePDG, const G4double particleEkin ) {

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
  //       << " \t edep = " << edep / MeV << " MeV "  << G4endl;  //***DEBUG***

  // Consider now the separate contribution due to the following particles:
  //    -  electron (e-  and  e+    together)
  //    -  muons    (mu- and  mu+   together)
  //    -  pions    (pi- and  pi+   together)
  //    -  kaons    (k-  and  k+    together)
  //    -  protons  (p   and  pbar  together)
  //    -  nuclei   (all particles with PDG code = 0 , or with
  //                 the new PDG code with 10-digits; also
  //                 neutrons are included, as an effective way
  //                 to take into account the recoil of nuclei
  //                 below a certain threshold)

  // 02-Apr-2007 : I have tried to transform the following long  if  
  //               statement into a  switch  statement in which 
  // "electronId", "positronId", etc. appear as "case label", but 
  // I got a strange compilation error saying that 
  //      "case label does not reduce to an integer constant"
  // in spite of being defined as integer constants. 
  // The errors do not disappear if static constants are used instead
  // of simple constant integers. So I gave up!

  //***HERE***

  if ( particlePDG == electronId  ||  particlePDG == positronId ) {
    if ( particleEkin > infParticleEkin_electron  &&  
	 particleEkin < supParticleEkin_electron ) {
      sumEdepAct_electron += edep;
      sumEdepTot_electron += edep;
      longitudinalProfile_electron[ readoutLayer ] += edep;
      transverseProfile_electron[ iBinRadius ] += edep;
    }
  } else if ( particlePDG == muonMinusId  ||  particlePDG == muonPlusId ) {
    if ( particleEkin > infParticleEkin_muon  &&  
	 particleEkin < supParticleEkin_muon ) {
      sumEdepAct_muon += edep;
      sumEdepTot_muon += edep;
      longitudinalProfile_muon[ readoutLayer ] += edep;
      transverseProfile_muon[ iBinRadius ] += edep;
    }
  } else if ( particlePDG == pionPlusId  ||  particlePDG == pionMinusId ) {
    if ( particleEkin > infParticleEkin_pion  &&  
	 particleEkin < supParticleEkin_pion ) {
      sumEdepAct_pion += edep;
      sumEdepTot_pion += edep;
      longitudinalProfile_pion[ readoutLayer ] += edep;
      transverseProfile_pion[ iBinRadius ] += edep;
    }
  } else if ( particlePDG == kaonMinusId  ||  particlePDG == kaonPlusId ) {
    if ( particleEkin > infParticleEkin_kaon  &&  
	 particleEkin < supParticleEkin_kaon ) {
      sumEdepAct_kaon += edep;
      sumEdepTot_kaon += edep;
      longitudinalProfile_kaon[ readoutLayer ] += edep;
      transverseProfile_kaon[ iBinRadius ] += edep;
    }
  } else if ( particlePDG == protonId  ||  particlePDG == antiProtonId ) {
    if ( particleEkin > infParticleEkin_proton  &&  
	 particleEkin < supParticleEkin_proton ) {
      sumEdepAct_proton += edep;
      sumEdepTot_proton += edep;
      longitudinalProfile_proton[ readoutLayer ] += edep;
      transverseProfile_proton[ iBinRadius ] += edep;
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
    if ( particleEkin > infParticleEkin_nuclei  &&  
	 particleEkin < supParticleEkin_nuclei ) {
      sumEdepAct_nuclei += edep;
      sumEdepTot_nuclei += edep;
      longitudinalProfile_nuclei[ readoutLayer ] += edep;
      transverseProfile_nuclei[ iBinRadius ] += edep;
    }
  }

}


void StatAccepTestAnalysis::infoStep( const G4Step* aStep ) {

  classifyParticle( false , aStep->GetTrack()->GetDefinition() );

  // To improve CPU performance, get particle properties only once 
  // and save them on variables.
  G4int particleId = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();
  G4String particleName = aStep->GetTrack()->GetDefinition()->GetParticleName();
  G4String volumeName = aStep->GetTrack()->GetVolume()->GetName();
  G4double stepLength = aStep->GetStepLength();
  G4double stepEnergyDeposit = aStep->GetTotalEnergyDeposit();
  G4double stepWeight = aStep->GetTrack()->GetWeight();

  //G4cout << "L=" << stepLength/mm << " Edep=" << stepEnergyDeposit/MeV << G4endl;

  // 02-Apr-2007 : I have tried to transform the following long  if  
  //               statement into a  switch  statement in which 
  // "electronId", "positronId", etc. appear as "case label", but 
  // I got a strange compilation error saying that 
  //      "case label does not reduce to an integer constant"
  // in spite of being defined as integer constants. 
  // The errors do not disappear if static constants are used instead
  // of simple constant integers. So I gave up!

  // Count the number of inelastic and stoppingAtRest processes.
  if ( StatAccepTestAnalysis::isCountingProcessesOn ) { 
    G4String processStr = "UserLimit";
    if ( aStep->GetPostStepPoint()->GetProcessDefinedStep() ) {
      processStr = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
    }  
    //G4cout << " Process = " << processStr << G4endl; //***DEBUG***
    if ( processStr.find( "Inelastic" ) != G4String::npos  ||
	 processStr.find( "AbsorptionAtRest" ) != G4String::npos  ||
	 processStr.find( "CaptureAtRest" ) != G4String::npos ) {
      //G4cout << " Process = " << processStr 
      //       << "\t particle = " << particleId << "\t" << particleName
      //       << G4endl;  //***DEBUG***
      if ( processStr.find( "Inelastic" ) != G4String::npos ) {
	numInelasticProcesses++;
	if ( particleId == protonId ) {
          numProtonInelasticProcesses++;
	} else if ( particleId == antiProtonId ) {
          numAntiProtonInelasticProcesses++;
	} else if ( particleId == neutronId ) {
          numNeutronInelasticProcesses++;
	} else if ( particleId == antiNeutronId ) {
          numAntiNeutronInelasticProcesses++;
	} else if ( particleId == pionPlusId ) {
          numPionPlusInelasticProcesses++;
	} else if ( particleId == pionMinusId ) {
          numPionMinusInelasticProcesses++;
	} else if ( particleId == kaonPlusId ) {
          numKaonPlusInelasticProcesses++;
	} else if ( particleId == kaonMinusId ) {
          numKaonMinusInelasticProcesses++;
	} else if ( particleId == kaonZeroId  ||  particleId == antiKaonZeroId  ||
                    particleId == kaonLongId  ||  particleId == kaonShortId ) {
          numKaonNeutralInelasticProcesses++;
	} else if ( particleId == 3122 ) {
	  numLambdaInelasticProcesses++;
	} else if ( particleId == -3122 ) {
	  numAntiLambdaInelasticProcesses++;
 	} else if ( particleId == 3112 ) {
	  numSigmaMinusInelasticProcesses++;
	} else if ( particleId == muonMinusId ) {
          numMuonMinusInelasticProcesses++;
        } else {
	  numOtherInelasticProcesses++;
	  //G4cout << " OtherInelasticProcess = " << processStr 
	  //	   << "\t particle = " <<  particleId 
	  //	   << "\t" << particleName
	  //	   << G4endl;  //***DEBUG***
        }
      } else {
	numStoppingAtRestProcesses++;
	if ( particleId == antiProtonId ) {
          numAntiProtonStoppingAtRestProcesses++;
	} else if ( particleId == neutronId ) {
          numNeutronStoppingAtRestProcesses++;
	} else if ( particleId == antiNeutronId ) {
          numAntiNeutronStoppingAtRestProcesses++;
	} else if ( particleId == pionMinusId ) {
          numPionMinusStoppingAtRestProcesses++;
	} else if ( particleId == kaonMinusId ) {
          numKaonMinusStoppingAtRestProcesses++;
	} else if ( particleId == kaonZeroId  ||  particleId == antiKaonZeroId  ||
                    particleId == kaonLongId  ||  particleId == kaonShortId ) {
          numKaonNeutralStoppingAtRestProcesses++;
	} else if ( particleId == 3122 ) {
	  numLambdaStoppingAtRestProcesses++;
	} else if ( particleId == -3122 ) {
	  numAntiLambdaStoppingAtRestProcesses++;
 	} else if ( particleId == 3112 ) {
	  numSigmaMinusStoppingAtRestProcesses++;
	} else if ( particleId == muonMinusId ) {
          numMuonMinusStoppingAtRestProcesses++;
        } else {
	  numOtherStoppingAtRestProcesses++;
	  G4cout << " OtherAtRestProcess = " << processStr 
		 << "\t particle = " << particleId << "\t" << particleName
		 << G4endl;  //***DEBUG***
        }
      }
    }  
  }

#ifdef G4ANALYSIS_USEROOT
  // 2D plots on Step Energy vs. Step Length.
  if ( isHistogramOn  &&  is2DHistogramStepLvsEOn ) {
    if ( volumeName == "physiActive" ) {
    	h_stepEvsL[ "active" ]->Fill( stepLength / mm, stepEnergyDeposit / MeV );
    } else if ( volumeName == "physiAbsorber" ) {
    	h_stepEvsL[ "absorber" ]->Fill( stepLength / mm, stepEnergyDeposit / MeV );
    }
    // 02-Apr-2007 : I have tried to transform the following long  if
    //               statement into a  switch  statement in which
    // "electronId", "positronId", etc. appear as "case label", but
    // I got a strange compilation error saying that
    //      "case label does not reduce to an integer constant"
    // in spite of being defined as integer constants.
    // The errors do not disappear if static constants are used instead
    // of simple constant integers. So I gave up!
    std::string postfix = "";
    if ( volumeName == "physiActive" ) {
      postfix="_active";
    } else if ( volumeName == "physiAbsorber" ) {
      postfix="_absorber";
    }
    if ( postfix != "" ) { // We still have to deal with the case the step is in expHall.
      std::string partname = "";
      if ( particleId == electronId  ||  particleId == positronId ) { 
	partname = "electron"; 
      } else if ( particleId == gammaId ) { 
	partname = "gamma";
      } else if ( particleId == muonMinusId  ||  particleId == muonPlusId ) { 
	partname = "muon"; 
      } else if ( particleId == pionPlusId  ||  particleId == pionMinusId ) { 
	partname = "pioncharged"; 
      } else if ( particleId == protonId  ||  particleId==antiProtonId ) { 
	partname = "proton"; 
      } else if ( particleId == neutronId  ||  particleId==antiNeutronId ) { 
	partname = "neutron"; 
      } else { 
	partname = "other"; 
      }
      std::string id = partname + postfix;
      histoMapConstIt_t it = h_stepEvsL.find( id );
      if ( it == h_stepEvsL.end() ) {
	G4cerr << "Error, cannot find histogram \"Step Edep Vs Length\" for particle id:"
	       << id << "(" << partname << ")" << " - Edep=" 
	       << stepEnergyDeposit/MeV << "(MeV) L=" <<stepLength/mm 
	       << "(mm)" << G4endl;
      } else {
	it->second->Fill( stepLength/mm, stepEnergyDeposit/MeV );
      }
    }
  }
#endif // G4ANALYSIS_USEROOT

  // Update the information on the energy deposition in the absorber
  // for the following particles (in the case of the energy deposition
  // in active layer, it has already been considered in the method
  // StatAccepTestAnalysis::fillShowerProfile):
  //    -  electron (e-  and  e+    together)
  //    -  muons    (mu- and  mu+   together)
  //    -  pions    (pi- and  pi+   together)
  //    -  kaons    (k-  and  k+    together)
  //    -  protons  (p   and  pbar  together)
  //    -  nuclei   (all particles with PDG code = 0 , or with
  //                 the new PDG code with 10-digits; also
  //                 neutrons are included, as an effective way
  //                 to take into account the recoil of nuclei
  //                 below a certain threshold)
  if ( volumeName == "physiAbsorber" ) {
    G4double edep = stepEnergyDeposit * stepWeight;
    if ( particleId == electronId  ||  particleId == positronId ) {
      sumEdepTot_electron += edep;
    } else if ( particleId == muonMinusId  ||  particleId == muonPlusId ) {
      sumEdepTot_muon += edep;
    } else if ( particleId == pionPlusId  ||  particleId == pionMinusId ) {
      sumEdepTot_pion += edep;
    } else if ( particleId == kaonMinusId  ||  particleId == kaonPlusId ) {
      sumEdepTot_kaon += edep;
    } else if ( particleId == protonId  ||  particleId == antiProtonId ) {
      sumEdepTot_proton += edep;
    } else if ( particleId == 0   ||
	        particleId / 1000000000 >= 1   ||
		particleId == neutronId ) {
      // Before 8.1.ref04, the PDG code for nuclei was 0; 
      // from 8.1.ref04, the new PDG code convention for nuclei has been
      // supported: in this schema, the PDG code of nuclei is characterized
      // by 10 digits (all other particles have lower codes).
      // Starting with 8.1, neutrons can also deposit energy, as an 
      // effective way to take into account the recoil of nuclei 
      // below a certain threshold.
      sumEdepTot_nuclei += edep;
    }
  }

  if ( StatAccepTestAnalysis::isMapInfoAboutTrackOn ) { 
    G4int trackId = aStep->GetTrack()->GetTrackID();
    if ( mapInfoAboutTrack.find( trackId ) != mapInfoAboutTrack.end() ) {
      G4double edep = stepEnergyDeposit * stepWeight;
      if ( volumeName == "physiActive" ) {
	mapInfoAboutTrack.find( trackId )->second.theTrackVisibleEdep += edep; 
	mapInfoAboutTrack.find( trackId )->second.theTrackDepositedEtot += edep; 	
      } 
      if ( volumeName == "physiAbsorber" ) {
	mapInfoAboutTrack.find( trackId )->second.theTrackDepositedEtot += edep;
      }
    } else {
      G4cout << "*** StatAccepTestAnalysis::infoStep : WARNING ***" << G4endl
	     << "\t  trackId=" << trackId << " NOT found in mapInfoAboutTrack !"
	     << G4endl;
    }
  }

}


void StatAccepTestAnalysis::infoTrack( const G4Track* aTrack ) {

  // To improve CPU performance, get particle properties only once 
  // and save them on variables.  
  G4int particleId = aTrack->GetDefinition()->GetPDGEncoding();
  G4String particleName = aTrack->GetDefinition()->GetParticleName();
  G4String volumeName = aTrack->GetVolume()->GetName();
  G4String volumeAtVertexName = aTrack->GetLogicalVolumeAtVertex()->GetName();
  G4double trackEkin = aTrack->GetKineticEnergy();
  G4double trackLength = aTrack->GetTrackLength();

  if ( aTrack->GetTrackStatus() == fStopAndKill ) {   
    //G4cout << "\t --- Info Track when fStopAndKill --- " << G4endl
    //	     << "\t TrackID = " << aTrack->GetTrackID() 
    //       << "\t Name = " << particleName << G4endl
    //       << "\t Volume = " << volumeName 
    //       << "\t Material = " << aTrack->GetMaterial()->GetName() << G4endl
    //       << "\t Vertex (origin) Volume = " << volumeAtVertexName << G4endl
    //       << "\t Ekin = " << trackEkin / MeV << " MeV "
    //       << "\t Length track = " << trackLength / mm << " mm " 
    //       << G4endl;  //***DEBUG*** 
    //if ( ! ( volumeName == "expHall"  ||
    //         volumeName == "physiAbsorber" ||
    //	       volumeName == "physiActive" ) ) {
    //  G4cout << " ***STRANGE VOLUME *** : " << volumeName << G4endl;
    //}
  
    // 02-Apr-2007 : I have tried to transform the following long  if  
    //               statement into a  switch  statement in which 
    // "electronId", "positronId", etc. appear as "case label", but 
    // I got a strange compilation error saying that 
    //      "case label does not reduce to an integer constant"
    // in spite of being defined as integer constants. 
    // The errors do not disappear if static constants are used instead
    // of simple constant integers. So I gave up!

    // To avoid bias in the track length due to the big world volume
    // (which can affect significantly the track length of neutrons)
    // we consider only those tracks that are fully contained inside
    // the calorimeter, i.e. created and terminated inside it.
    if ( volumeAtVertexName != "expHall"  &&  volumeName != "expHall" ) {
      if ( particleId == electronId || particleId == positronId ) {
	electronTrackLength += trackLength;
	electronTrackLength2 += trackLength * trackLength;
      } else if ( particleId == gammaId ) {
	gammaTrackLength += trackLength;
	gammaTrackLength2 += trackLength * trackLength;
      } else if ( particleId == muonMinusId  ||  particleId == muonPlusId ) {
	muonTrackLength += trackLength;
	muonTrackLength2 += trackLength * trackLength;
      } else if ( particleId == pionPlusId  ||  particleId == pionMinusId ) {
	pionChargedTrackLength += trackLength;
	pionChargedTrackLength2 += trackLength * trackLength;
      } else if ( particleId == pionZeroId ) {
	pion0TrackLength += trackLength;
	pion0TrackLength2 += trackLength * trackLength;
      } else if ( particleId == protonId ) {
	protonTrackLength += trackLength;
	protonTrackLength2 += trackLength * trackLength;
      } else if ( particleId == neutronId ) {
	neutronTrackLength += trackLength;
	neutronTrackLength2 += trackLength * trackLength;
      }
    } else if ( volumeName == "expHall" ) {
      kinEnergyExiting += trackEkin;
      kinEnergyExiting2 += trackEkin * trackEkin;
      numExiting++;
      //G4cout << " Exiting particle: " << particleName 
      //       << "\t" << kinEnergyExiting / MeV << " MeV" << G4endl
      //       << "\t Production: " 
      //       << aTrack->GetCreatorProcess()->GetProcessName() << G4endl
      //       << "\t" << aTrack->GetVertexPosition() / mm << " mm"
      //       << "\t" << volumeAtVertexName << G4endl
      //       << "\t vertex Ekin = " << aTrack->GetVertexKineticEnergy() / MeV 
      //       << " MeV" << G4endl;  //***DEBUG***
      if ( particleId == gammaId ) {
	kinEnergyExitingGammas += trackEkin;
	numExitingGammas++;
      } else if ( particleId == neutronId ) {
	kinEnergyExitingNeutrons += trackEkin;
	numExitingNeutrons++;
      } else if ( particleId == eNeutrinoId  ||  particleId == antiENeutrinoId  ||
		  particleId == muNeutrinoId  ||  particleId == antiMuNeutrinoId  ||
		  particleId == tauNeutrinoId  ||  particleId == antiTauNeutrinoId ) {
	kinEnergyExitingNeutrinos += trackEkin;
	numExitingNeutrinos++;
      } else if ( particleId == muonMinusId  ||  particleId == muonPlusId ) {
	kinEnergyExitingMuons += trackEkin;
	numExitingMuons++;        
      } else if ( particleId == electronId  ||  particleId == positronId ) {
	kinEnergyExitingElectrons += trackEkin;
	numExitingElectrons++;        
      } else {
	kinEnergyExitingOthers += trackEkin;
	numExitingOthers++;
      }
    }
  } else if ( aTrack->GetCurrentStepNumber() == 0 ) {  // When a track is created.

    classifyParticle( true , aTrack->GetDefinition() );

    if ( StatAccepTestAnalysis::isMapParticleNamesOn  &&  
	 particleId != gammaId  &&  
         particleId != electronId  &&  particleId != positronId ) {
      if ( mapParticleNames.find( particleName ) == mapParticleNames.end() ) {
	mapParticleNames.insert( std::pair< std::string, int >( particleName, 1 ) );
      } else {
	mapParticleNames.find( particleName )->second += 1;
      }
    }
#ifdef G4ANALYSIS_USEROOT
    if ( trackEkin < 1.0*eV ) {
      // This is to avoid problems with the logarithm.
      trackEkin = 1.0*eV;
    }
    if ( particleId == gammaId ) {
    	h_Spectrum[ "gamma" ]->Fill( std::log10( trackEkin/MeV ) );
    } else if ( particleId == neutronId ) {
    	h_Spectrum[ "neutron" ]->Fill( std::log10( trackEkin/MeV ) );
    } else if ( particleId == protonId ) {
    	h_Spectrum[ "proton" ]->Fill( std::log10( trackEkin/MeV ) );
    } else if ( particleId == pionZeroId ) {
    	h_Spectrum[ "pionZero" ]->Fill( std::log10( trackEkin/MeV ) );
    } else if ( particleId == pionPlusId ) {
    	h_Spectrum[ "pionPlus" ]->Fill( std::log10( trackEkin/MeV ) );
    } else if ( particleId == pionMinusId ) {
    	h_Spectrum[ "pionMinus" ]->Fill( std::log10( trackEkin/MeV ) );
    }
#endif // G4ANALYSIS_USEROOT


    if ( StatAccepTestAnalysis::isMapInfoAboutTrackOn ) { 

      G4int trackId = aTrack->GetTrackID();
      G4int parentId = aTrack->GetParentID();
      const G4VProcess* creatorProcess = aTrack->GetCreatorProcess();

      if ( parentId > 0  &&  creatorProcess ) {

        // Find the closest vertex with the same parent and same process type.
	G4double distance2 = 9.9E+99; // A very large value.
	std::map< G4int, structInfoAboutVertex >::const_iterator found_citvtx = 
	  mapInfoAboutVertex.end();
	for ( std::map< G4int, structInfoAboutVertex >::const_iterator 
		citvtx = mapInfoAboutVertex.begin();
	      citvtx != mapInfoAboutVertex.end() ; ++citvtx ) {
	  if ( citvtx->second.theCreatorTrackID == parentId  
	       &&
	       citvtx->second.theCreatorProcessType == 
               creatorProcess->GetProcessType()              
	       &&
	       citvtx->second.theCreatorProcessSubType == 
	       creatorProcess->GetProcessSubType()          
	     ) {
	    G4double d2 = 
	      ( aTrack->GetPosition().x() - citvtx->second.theVertexPosition_x ) * 
	      ( aTrack->GetPosition().x() - citvtx->second.theVertexPosition_x )
	      +
	      ( aTrack->GetPosition().y() - citvtx->second.theVertexPosition_y ) * 
	      ( aTrack->GetPosition().y() - citvtx->second.theVertexPosition_y )
	      +
	      ( aTrack->GetPosition().z() - citvtx->second.theVertexPosition_z ) * 
	      ( aTrack->GetPosition().z() - citvtx->second.theVertexPosition_z );
	    if ( d2 < distance2 ) {
	      distance2 = d2;
	      found_citvtx = citvtx;
            }
          }
        }
        // Create a new vertex either if a vertex with the same
        // parent and process type does not exist, or if it exists
        // but it is too displaced, by more than 1 nm .
        if ( found_citvtx == mapInfoAboutVertex.end()  ||
	     distance2 > 1.0E-12 ) {
	  
	  G4ProcessType creatorProcessType = creatorProcess->GetProcessType();
	  G4int creatorProcessSubType = creatorProcess->GetProcessSubType();
	  G4String creatorProcessName = creatorProcess->GetProcessName();
	  G4String volumeName = aTrack->GetVolume()->GetName();
	  G4double const x = aTrack->GetPosition().x(); 
	  G4double const y = aTrack->GetPosition().y(); 
	  G4double const z = aTrack->GetPosition().z(); 
	  G4double const t = aTrack->GetGlobalTime(); // Time since the start of 
                                                      // the event [nsec].

	  structInfoAboutVertex aStructInfoAboutVertex;

	  aStructInfoAboutVertex.theCreatorTrackID = parentId;
	  aStructInfoAboutVertex.theCreatorProcessType = 
	    static_cast< G4int >( creatorProcessType );
	  aStructInfoAboutVertex.theCreatorProcessSubType = creatorProcessSubType;
	  aStructInfoAboutVertex.theCreatorProcessName = creatorProcessName;
	  aStructInfoAboutVertex.theVertexVolumeName = volumeName;
	  aStructInfoAboutVertex.theVertexPosition_x = x;
	  aStructInfoAboutVertex.theVertexPosition_y = y;
	  aStructInfoAboutVertex.theVertexPosition_z = z;
	  aStructInfoAboutVertex.theVertexTime = t;

	  G4int numberOfVertices = mapInfoAboutVertex.size();

	  suggestedVertexPosition = 
	    mapInfoAboutVertex.insert( suggestedVertexPosition,
	  			       std::pair< G4int, structInfoAboutVertex >
	  			       ( numberOfVertices + 1, aStructInfoAboutVertex ) );
	}
      }

      structInfoAboutTrack aStructInfoAboutTrack;

      aStructInfoAboutTrack.theParticleName = particleName;
      aStructInfoAboutTrack.theParticlePDGCode = particleId;
      aStructInfoAboutTrack.theTrackEkinAtCreation = trackEkin;
      aStructInfoAboutTrack.theTrackVisibleEdep = 0.0;
      aStructInfoAboutTrack.theTrackDepositedEtot = 0.0;
      aStructInfoAboutTrack.theParentTrackID = parentId;

      // The starting vertex for the track is the closest vertex 
      // originated by parentId, and with the same process type.
      G4double distance2 = 9.9E+99; // A very large value.
      std::map< G4int, structInfoAboutVertex >::const_iterator found_citvtx = 
	mapInfoAboutVertex.end();
      for ( std::map< G4int, structInfoAboutVertex >::const_iterator 
	      citvtx = mapInfoAboutVertex.begin();
	    citvtx != mapInfoAboutVertex.end() ; ++citvtx ) {
	if ( citvtx->second.theCreatorTrackID == parentId  
	     &&
	     citvtx->second.theCreatorProcessType == 
	     creatorProcess->GetProcessType()              
	     &&
	     citvtx->second.theCreatorProcessSubType == 
	     creatorProcess->GetProcessSubType()          
	   ) {
	  G4double d2 = 
	    ( aTrack->GetPosition().x() - citvtx->second.theVertexPosition_x ) * 
	    ( aTrack->GetPosition().x() - citvtx->second.theVertexPosition_x )
	    +
	    ( aTrack->GetPosition().y() - citvtx->second.theVertexPosition_y ) * 
	    ( aTrack->GetPosition().y() - citvtx->second.theVertexPosition_y )
	    +
	    ( aTrack->GetPosition().z() - citvtx->second.theVertexPosition_z ) * 
	    ( aTrack->GetPosition().z() - citvtx->second.theVertexPosition_z );
	  if ( d2 < distance2 ) {
	    distance2 = d2;
	    found_citvtx = citvtx;
	  }
	}
      }
      if ( found_citvtx != mapInfoAboutVertex.end() ) {
	aStructInfoAboutTrack.theStartingVertexID = found_citvtx->first;	
      } else {  // It should happens only for the primary track.
	aStructInfoAboutTrack.theStartingVertexID = 0;
      }

      // To find the closest hadronic relative, we need to look
      // at the parent track, see if it produces a hadronic
      // interaction; if this is the case, then the closest
      // hadronic relative is the parent and we are done;
      // else, we have to repeat the same procedure for the 
      // grand-parent, then eventually for the grand-grand-parent,
      // and so on, until either we find that one of these produces
      // a hadronic interaction, or we scan the whole map, reaching
      // the primary.
      bool isSearchCompleted = false;
      bool isFound = false;
      G4int currentParent = parentId;
      G4int currentParentVertex = aStructInfoAboutTrack.theStartingVertexID;
      while ( ! isSearchCompleted  &&  ! isFound ) {
	std::map< G4int, structInfoAboutVertex >::const_iterator citvtx_bis = 
	  mapInfoAboutVertex.find( currentParentVertex );
	if ( citvtx_bis != mapInfoAboutVertex.end() ) {
	  G4int subtype = citvtx_bis->second.theCreatorProcessSubType;
	  if ( 
	      //subtype == 111  ||     // G4HadronicProcessType::fHadronElastic
	      subtype == 121  ||     // G4HadronicProcessType::fHadronInelastic
	      subtype == 131  ||     // G4HadronicProcessType::fCapture
	      subtype == 141  ||     // G4HadronicProcessType::fFission
	      subtype == 151  ||     // G4HadronicProcessType::fHadronAtRest
	      subtype == 161         // G4HadronicProcessType::fChargeExchange 
	      ) {
	    isFound = true;
          } else {
	    std::map< G4int, structInfoAboutTrack >::const_iterator cittrk =
	      mapInfoAboutTrack.find( currentParent );
	    if ( cittrk != mapInfoAboutTrack.end() ) {
	      currentParent = cittrk->second.theParentTrackID;
	      currentParentVertex = cittrk->second.theStartingVertexID;
            } else {
	      isSearchCompleted = true;
            }
          }
	} else {
	  isSearchCompleted = true;
	}
      } // end of while loop
      if ( isFound ) {
	aStructInfoAboutTrack.theClosestHadronicRelativeID = currentParent;
        aStructInfoAboutTrack.theClosestHadronicVertexID = currentParentVertex;
      } else {  // It should happens only for the primary track and
                // eventual ionization electrons or bremsstrahlung
	        // gammas produced by the primary track before the
	        // first inelastic hadronic interaction.
	aStructInfoAboutTrack.theClosestHadronicRelativeID = 0;
        aStructInfoAboutTrack.theClosestHadronicVertexID = 0;
      }

      suggestedTrackPosition = 
	mapInfoAboutTrack.insert( suggestedTrackPosition,
				  std::pair< G4int, structInfoAboutTrack >
				  ( trackId , aStructInfoAboutTrack ) );

    } // end of if ( StatAccepTestAnalysis::isMapInfoAboutTrackOn )

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

  // 02-Apr-2007 : I have tried to transform the following long  if  
  //               statement into a  switch  statement in which 
  // "electronId", "positronId", etc. appear as "case label", but 
  // I got a strange compilation error saying that 
  //      "case label does not reduce to an integer constant"
  // in spite of being defined as integer constants. 
  // The errors do not disappear if static constants are used instead
  // of simple constant integers. So I gave up!

  if ( id == 0  ||  id / 1000000000 >= 1 ) {
    // Before 8.1.ref04, the PDG code for nuclei was 0; 
    // from 8.1.ref04, the new PDG code convention for nuclei has been
    // supported: in this schema, the PDG code of nuclei is characterized
    // by 10 digits (all other particles have lower codes).
    if ( isTrack ) {
      numTrackNucleus++;
    } else {
      numStepNucleus++;
    }
  } else if ( id == gammaId ) {
    if ( isTrack ) {
      numTrackNeutral++;
      numTrackEM++;
      numTrackGamma++;
    } else {
      numStepNeutral++;
      numStepEM++;
      numStepGamma++;
    }
  } else if ( id == electronId ) {
    if ( isTrack ) {
      numTrackNegative++;
      numTrackEM++;
      numTrackElectron++;
    } else {
      numStepNegative++;
      numStepEM++;
      numStepElectron++;
    }
  } else if ( id == positronId ) {
    if ( isTrack ) {
      numTrackPositive++;
      numTrackEM++;
      numTrackPositron++;
    } else {
      numStepPositive++;
      numStepEM++;
      numStepPositron++;
    }
  } else if ( id == muonMinusId ) {
    if ( isTrack ) {
      numTrackNegative++;
      numTrackEWK++;
      numTrackMuMinus++;
    } else {
      numStepNegative++;
      numStepEWK++;
      numStepMuMinus++;
    } 
  } else if ( id == muonPlusId ) {
    if ( isTrack ) {
      numTrackPositive++;
      numTrackEWK++;
      numTrackMuPlus++;
    } else {
      numStepPositive++;
      numStepEWK++;
      numStepMuPlus++;
    }
  } else if ( id == tauMinusId ) {
    if ( isTrack ) {
      numTrackNegative++;
      numTrackEWK++;
      numTrackTauMinus++;
    } else {
      numStepNegative++;
      numStepEWK++; 
      numStepTauMinus++;
    } 
  } else if ( id == tauPlusId ) {
    if ( isTrack ) {
      numTrackPositive++;
      numTrackEWK++;
      numTrackTauPlus++;
    } else {
      numStepPositive++;
      numStepEWK++; 
      numStepTauPlus++;
    } 
  } else if ( id == eNeutrinoId  ||  id == antiENeutrinoId  ||
	      id == muNeutrinoId  ||  id == antiMuNeutrinoId  ||
	      id == tauNeutrinoId  ||  id == antiTauNeutrinoId ) {
    if ( isTrack ) {
      numTrackNeutral++;
      numTrackEWK++;
      numTrackNeutrino++;
    } else {
      numStepNeutral++;
      numStepEWK++;
      numStepNeutrino++;
    }
  } else if ( id == pionMinusId ) {
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
  } else if ( id == pionZeroId ) {
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
  } else if ( id == pionPlusId ) {
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
  } else if ( id == kaonMinusId ) {
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
  } else if ( id == kaonZeroId  ||  id == antiKaonZeroId  || 
	      id == kaonLongId  ||  id == kaonShortId ) {
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
  } else if ( id == kaonPlusId ) {
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
  } else if ( id == neutronId ) {
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
  } else if ( id == antiNeutronId ) {
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
  } else if ( id == protonId ) {
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
  } else if ( id == antiProtonId ) {
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


void StatAccepTestAnalysis::endOfEvent( const G4double timeEventInSec ) {
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
  static G4double numStepNucleus_previous = 0.0;
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
  numStepNucleus2 +=
    ( numStepNucleus - numStepNucleus_previous ) *
    ( numStepNucleus - numStepNucleus_previous );
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
  numStepNucleus_previous = numStepNucleus;
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
  static G4double numTrackNucleus_previous = 0.0;
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
  numTrackNucleus2 +=
    ( numTrackNucleus - numTrackNucleus_previous ) *
    ( numTrackNucleus - numTrackNucleus_previous );
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
  numTrackNucleus_previous = numTrackNucleus;
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
  G4double totalEnergyDepositedInActiveLayers_electron = 
    sumEdepAct_electron - sumEdepAct_electron_previous;
  vecEvis_electron.push_back( totalEnergyDepositedInActiveLayers_electron );
  sumEdepAct_electron2 += totalEnergyDepositedInActiveLayers_electron *
                          totalEnergyDepositedInActiveLayers_electron;
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
  G4double totalEnergyDepositedInActiveLayers_muon = 
    sumEdepAct_muon - sumEdepAct_muon_previous;
  vecEvis_muon.push_back( totalEnergyDepositedInActiveLayers_muon );
  sumEdepAct_muon2 += totalEnergyDepositedInActiveLayers_muon *
                          totalEnergyDepositedInActiveLayers_muon;
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
  G4double totalEnergyDepositedInActiveLayers_pion = 
    sumEdepAct_pion - sumEdepAct_pion_previous;
  vecEvis_pion.push_back( totalEnergyDepositedInActiveLayers_pion );
  sumEdepAct_pion2 += totalEnergyDepositedInActiveLayers_pion *
                          totalEnergyDepositedInActiveLayers_pion;
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
  G4double totalEnergyDepositedInActiveLayers_kaon = 
    sumEdepAct_kaon - sumEdepAct_kaon_previous;
  vecEvis_kaon.push_back( totalEnergyDepositedInActiveLayers_kaon );
  sumEdepAct_kaon2 += totalEnergyDepositedInActiveLayers_kaon *
                          totalEnergyDepositedInActiveLayers_kaon;
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
  G4double totalEnergyDepositedInActiveLayers_proton = 
    sumEdepAct_proton - sumEdepAct_proton_previous;
  vecEvis_proton.push_back( totalEnergyDepositedInActiveLayers_proton );
  sumEdepAct_proton2 += totalEnergyDepositedInActiveLayers_proton *
                          totalEnergyDepositedInActiveLayers_proton;
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

  // Nuclei  and  neutrons  together
  static G4double sumEdepAct_nuclei_previous = 0.0;
  G4double totalEnergyDepositedInActiveLayers_nuclei = 
    sumEdepAct_nuclei - sumEdepAct_nuclei_previous;
  vecEvis_nuclei.push_back( totalEnergyDepositedInActiveLayers_nuclei );
  sumEdepAct_nuclei2 += totalEnergyDepositedInActiveLayers_nuclei *
                          totalEnergyDepositedInActiveLayers_nuclei;
  sumEdepAct_nuclei_previous = sumEdepAct_nuclei;
  static G4double sumEdepTot_nuclei_previous = 0.0;
  sumEdepTot_nuclei2 += ( sumEdepTot_nuclei - sumEdepTot_nuclei_previous ) *
                        ( sumEdepTot_nuclei - sumEdepTot_nuclei_previous );
  sumEdepTot_nuclei_previous = sumEdepTot_nuclei;
  for ( int iLayer = 0; iLayer < numberOfReadoutLayers; iLayer++ ) {
    sumL_nuclei[ iLayer ]  += longitudinalProfile_nuclei[ iLayer ];
    sumL_nuclei2[ iLayer ] += longitudinalProfile_nuclei[ iLayer ] * 
                              longitudinalProfile_nuclei[ iLayer ];  
    longitudinalProfile_nuclei[ iLayer ] = 0.0;  // Reset it for the next event.
  }
  for ( int iBinR = 0; iBinR < numberOfRadiusBins; iBinR++ ) {
    sumR_nuclei[ iBinR ]  += transverseProfile_nuclei[ iBinR ];
    sumR_nuclei2[ iBinR ] += transverseProfile_nuclei[ iBinR ] * 
                             transverseProfile_nuclei[ iBinR ];  
    transverseProfile_nuclei[ iBinR ] = 0.0;  // Reset it for the next event.
  }

  G4double totalEnergyDepositedInActiveLayers =     
    totalEnergyDepositedInActiveLayers_electron +
    totalEnergyDepositedInActiveLayers_muon +
    totalEnergyDepositedInActiveLayers_pion +
    totalEnergyDepositedInActiveLayers_kaon +
    totalEnergyDepositedInActiveLayers_proton +
    totalEnergyDepositedInActiveLayers_nuclei;

  // All but  e-  and  e+
  G4double totalEnergyDepositedInActiveLayers_no_electron =
    totalEnergyDepositedInActiveLayers - totalEnergyDepositedInActiveLayers_electron;
  vecEvis_no_electron.push_back( totalEnergyDepositedInActiveLayers_no_electron );
  sumEdepAct_no_electron += totalEnergyDepositedInActiveLayers_no_electron;
  sumEdepAct_no_electron2 += ( totalEnergyDepositedInActiveLayers_no_electron * 
			       totalEnergyDepositedInActiveLayers_no_electron ); 

  // All but  mu-  and  mu+
  G4double totalEnergyDepositedInActiveLayers_no_muon =
    totalEnergyDepositedInActiveLayers - totalEnergyDepositedInActiveLayers_muon;
  vecEvis_no_muon.push_back( totalEnergyDepositedInActiveLayers_no_muon );
  sumEdepAct_no_muon  += totalEnergyDepositedInActiveLayers_no_muon;
  sumEdepAct_no_muon2 += ( totalEnergyDepositedInActiveLayers_no_muon * 
			   totalEnergyDepositedInActiveLayers_no_muon ); 

  // All but  pi+  and  pi-
  G4double totalEnergyDepositedInActiveLayers_no_pion =
    totalEnergyDepositedInActiveLayers - totalEnergyDepositedInActiveLayers_pion;
  vecEvis_no_pion.push_back( totalEnergyDepositedInActiveLayers_no_pion );
  sumEdepAct_no_pion += totalEnergyDepositedInActiveLayers_no_pion;
  sumEdepAct_no_pion2 += ( totalEnergyDepositedInActiveLayers_no_pion * 
			   totalEnergyDepositedInActiveLayers_no_pion ); 

  // All but  k+  and  k-
  G4double totalEnergyDepositedInActiveLayers_no_kaon =
    totalEnergyDepositedInActiveLayers - totalEnergyDepositedInActiveLayers_kaon;
  vecEvis_no_kaon.push_back( totalEnergyDepositedInActiveLayers_no_kaon );
  sumEdepAct_no_kaon += totalEnergyDepositedInActiveLayers_no_kaon;
  sumEdepAct_no_kaon2 += ( totalEnergyDepositedInActiveLayers_no_kaon * 
			   totalEnergyDepositedInActiveLayers_no_kaon ); 

  // All but  p  and  pbar
  G4double totalEnergyDepositedInActiveLayers_no_proton =
    totalEnergyDepositedInActiveLayers - totalEnergyDepositedInActiveLayers_proton;
  vecEvis_no_proton.push_back( totalEnergyDepositedInActiveLayers_no_proton );
  sumEdepAct_no_proton += totalEnergyDepositedInActiveLayers_no_proton;
  sumEdepAct_no_proton2 += ( totalEnergyDepositedInActiveLayers_no_proton * 
			     totalEnergyDepositedInActiveLayers_no_proton ); 

  // All but  Nuclei  and  neutrons
  G4double totalEnergyDepositedInActiveLayers_no_nuclei =
    totalEnergyDepositedInActiveLayers - totalEnergyDepositedInActiveLayers_nuclei;
  vecEvis_no_nuclei.push_back( totalEnergyDepositedInActiveLayers_no_nuclei );
  sumEdepAct_no_nuclei += totalEnergyDepositedInActiveLayers_no_nuclei;
  sumEdepAct_no_nuclei2 += ( totalEnergyDepositedInActiveLayers_no_nuclei * 
			     totalEnergyDepositedInActiveLayers_no_nuclei ); 

  // Keep the time per event in a multiset.
  eventTimeSet.insert( timeEventInSec );

  if ( StatAccepTestAnalysis::isMapInfoAboutTrackOn ) { 
    //analysisTrackAndVertices_0();  //***DEBUG*** : write details event information
    analysisTrackAndVertices_1();
    analysisTrackAndVertices_2();
    analysisTrackAndVertices_3() ;
    analysisTrackAndVertices_4() ;

    mapInfoAboutTrack.clear();
    mapInfoAboutVertex.clear();

    suggestedTrackPosition = mapInfoAboutTrack.begin();
    suggestedVertexPosition = mapInfoAboutVertex.begin();
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
  G4cout << " Beam energy [MeV] = " << beamEnergy / MeV << G4endl;

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
    G4cout << " maxEdepTot  [MeV] = " << maxEdepTot / MeV 
           << "\t ***ENERGY-NON-CONSERVATION*** " << G4endl;
    G4cout << "\t number of event with energy non conservation = " 
	   << countEnergyNonConservation << G4endl;
  } else {
    G4cout << " maxEdepTot  [MeV] = " << maxEdepTot / MeV << G4endl;
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

#ifdef G4ANALYSIS_USEROOT
    h_longitudinalProfile->SetBinContent( h_longitudinalProfile->FindBin( iLayer ), mu );
    h_longitudinalProfile->SetBinError( h_longitudinalProfile->FindBin( iLayer ), mu_sigma );
#endif

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
#ifdef G4ANALYSIS_USEROOT
    h_transverseProfile->SetBinContent( h_transverseProfile->FindBin( iBinR ), mu );
    h_transverseProfile->SetBinError( h_transverseProfile->FindBin( iBinR ), mu_sigma );
#endif
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

  // Print information on the different particle contributions to
  // the visible energy, the total energy, and the shower shapes.

  G4cout << G4endl << " Contributions of the main particle types [MeV] " << G4endl;
  for ( int iCase = 0; iCase < 6; iCase++ ) {

    std::string caseName = "";
    G4double sumVis = 0.0;
    G4double sumVis2 = 0.0;
    G4double sumTot = 0.0;
    G4double sumTot2 = 0.0;
    G4double no_particle_sumVis = 0.0;
    G4double no_particle_sumVis2 = 0.0;
    std::vector< G4double > vecSumL;
    std::vector< G4double > vecSumL2;
    std::vector< G4double > vecSumR;
    std::vector< G4double > vecSumR2;
    std::vector< G4double > particle_vecEvis;
    std::vector< G4double > no_particle_vecEvis;

    switch ( iCase ) {
    case 0 : {
      caseName = "electron";
      sumVis = sumEdepAct_electron;
      sumVis2 = sumEdepAct_electron2;
      sumTot = sumEdepTot_electron;
      sumTot2 = sumEdepTot_electron2;
      no_particle_sumVis = sumEdepAct_no_electron;
      no_particle_sumVis2 = sumEdepAct_no_electron2;
      for ( int iLayer = 0; iLayer < numberOfReadoutLayers; iLayer++ ) {
	vecSumL.push_back( sumL_electron[ iLayer ] );
	vecSumL2.push_back( sumL_electron2[ iLayer ] );
      }
      for ( int iBinR = 0; iBinR < numberOfRadiusBins; iBinR++ ) {
	vecSumR.push_back( sumR_electron[ iBinR ] );
	vecSumR2.push_back( sumR_electron2[ iBinR ] );
      }
      for ( std::vector< G4double >::const_iterator cit = vecEvis_electron.begin();
	    cit != vecEvis_electron.end() ; ++cit ) {
	particle_vecEvis.push_back( *cit );
      }
      for ( std::vector< G4double >::const_iterator cit = vecEvis_no_electron.begin();
	    cit != vecEvis_no_electron.end() ; ++cit ) {
	no_particle_vecEvis.push_back( *cit );
      }
      break;
    }
    case 1 : {
      caseName = "muon";
      sumVis = sumEdepAct_muon;
      sumVis2 = sumEdepAct_muon2;
      sumTot = sumEdepTot_muon;
      sumTot2 = sumEdepTot_muon2;
      no_particle_sumVis = sumEdepAct_no_muon;
      no_particle_sumVis2 = sumEdepAct_no_muon2;
      for ( int iLayer = 0; iLayer < numberOfReadoutLayers; iLayer++ ) {
	vecSumL.push_back( sumL_muon[ iLayer ] );
	vecSumL2.push_back( sumL_muon2[ iLayer ] );
      }
      for ( int iBinR = 0; iBinR < numberOfRadiusBins; iBinR++ ) {
	vecSumR.push_back( sumR_muon[ iBinR ] );
	vecSumR2.push_back( sumR_muon2[ iBinR ] );
      }
      for ( std::vector< G4double >::const_iterator cit = vecEvis_muon.begin();
	    cit != vecEvis_muon.end() ; ++cit ) {
	particle_vecEvis.push_back( *cit );
      }
      for ( std::vector< G4double >::const_iterator cit = vecEvis_no_muon.begin();
	    cit != vecEvis_no_muon.end() ; ++cit ) {
	no_particle_vecEvis.push_back( *cit );
      }
      break;
    }
    case 2 : {
      caseName = "pion";
      sumVis = sumEdepAct_pion;
      sumVis2 = sumEdepAct_pion2;
      sumTot = sumEdepTot_pion;
      sumTot2 = sumEdepTot_pion2;
      no_particle_sumVis = sumEdepAct_no_pion;
      no_particle_sumVis2 = sumEdepAct_no_pion2;
      for ( int iLayer = 0; iLayer < numberOfReadoutLayers; iLayer++ ) {
	vecSumL.push_back( sumL_pion[ iLayer ] );
	vecSumL2.push_back( sumL_pion2[ iLayer ] );
      }
      for ( int iBinR = 0; iBinR < numberOfRadiusBins; iBinR++ ) {
	vecSumR.push_back( sumR_pion[ iBinR ] );
	vecSumR2.push_back( sumR_pion2[ iBinR ] );
      }
      for ( std::vector< G4double >::const_iterator cit = vecEvis_pion.begin();
	    cit != vecEvis_pion.end() ; ++cit ) {
	particle_vecEvis.push_back( *cit );
      }
      for ( std::vector< G4double >::const_iterator cit = vecEvis_no_pion.begin();
	    cit != vecEvis_no_pion.end() ; ++cit ) {
	no_particle_vecEvis.push_back( *cit );
      }
      break;
    }
    case 3 : {
      caseName = "kaon";
      sumVis = sumEdepAct_kaon;
      sumVis2 = sumEdepAct_kaon2;
      sumTot = sumEdepTot_kaon;
      sumTot2 = sumEdepTot_kaon2;
      no_particle_sumVis = sumEdepAct_no_kaon;
      no_particle_sumVis2 = sumEdepAct_no_kaon2;
      for ( int iLayer = 0; iLayer < numberOfReadoutLayers; iLayer++ ) {
	vecSumL.push_back( sumL_kaon[ iLayer ] );
	vecSumL2.push_back( sumL_kaon2[ iLayer ] );
      }
      for ( int iBinR = 0; iBinR < numberOfRadiusBins; iBinR++ ) {
	vecSumR.push_back( sumR_kaon[ iBinR ] );
	vecSumR2.push_back( sumR_kaon2[ iBinR ] );
      }
      for ( std::vector< G4double >::const_iterator cit = vecEvis_kaon.begin();
	    cit != vecEvis_kaon.end() ; ++cit ) {
	particle_vecEvis.push_back( *cit );
      }
      for ( std::vector< G4double >::const_iterator cit = vecEvis_no_kaon.begin();
	    cit != vecEvis_no_kaon.end() ; ++cit ) {
	no_particle_vecEvis.push_back( *cit );
      }
      break;
    }
    case 4 : {
      caseName = "proton";
      sumVis = sumEdepAct_proton;
      sumVis2 = sumEdepAct_proton2;
      sumTot = sumEdepTot_proton;
      sumTot2 = sumEdepTot_proton2;
      no_particle_sumVis = sumEdepAct_no_proton;
      no_particle_sumVis2 = sumEdepAct_no_proton2;
      for ( int iLayer = 0; iLayer < numberOfReadoutLayers; iLayer++ ) {
	vecSumL.push_back( sumL_proton[ iLayer ] );
	vecSumL2.push_back( sumL_proton2[ iLayer ] );
      }
      for ( int iBinR = 0; iBinR < numberOfRadiusBins; iBinR++ ) {
	vecSumR.push_back( sumR_proton[ iBinR ] );
	vecSumR2.push_back( sumR_proton2[ iBinR ] );
      }
      for ( std::vector< G4double >::const_iterator cit = vecEvis_proton.begin();
	    cit != vecEvis_proton.end() ; ++cit ) {
	particle_vecEvis.push_back( *cit );
      }
      for ( std::vector< G4double >::const_iterator cit = vecEvis_no_proton.begin();
	    cit != vecEvis_no_proton.end() ; ++cit ) {
	no_particle_vecEvis.push_back( *cit );
      }
      break;
    }
    case 5 : {
      caseName = "nuclei";
      sumVis = sumEdepAct_nuclei;
      sumVis2 = sumEdepAct_nuclei2;
      sumTot = sumEdepTot_nuclei;
      sumTot2 = sumEdepTot_nuclei2;
      no_particle_sumVis = sumEdepAct_no_nuclei;
      no_particle_sumVis2 = sumEdepAct_no_nuclei2;
      for ( int iLayer = 0; iLayer < numberOfReadoutLayers; iLayer++ ) {
	vecSumL.push_back( sumL_nuclei[ iLayer ] );
	vecSumL2.push_back( sumL_nuclei2[ iLayer ] );
      }
      for ( int iBinR = 0; iBinR < numberOfRadiusBins; iBinR++ ) {
	vecSumR.push_back( sumR_nuclei[ iBinR ] );
	vecSumR2.push_back( sumR_nuclei2[ iBinR ] );
      }
      for ( std::vector< G4double >::const_iterator cit = vecEvis_nuclei.begin();
	    cit != vecEvis_nuclei.end() ; ++cit ) {
	particle_vecEvis.push_back( *cit );
      }
      for ( std::vector< G4double >::const_iterator cit = vecEvis_no_nuclei.begin();
	    cit != vecEvis_no_nuclei.end() ; ++cit ) {
	no_particle_vecEvis.push_back( *cit );
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

    sum  = no_particle_sumVis;
    sum2 = no_particle_sumVis2;
    mu       = sum / n;
    sigma    = std::sqrt( std::abs( sum2 - sum*sum/n ) / (n - 1.0) );
    mu_sigma = sigma / std::sqrt( n );
    G4double no_particle_mu_Evis = mu;             // For later usage.
    G4double no_particle_mu_Evis_sigma = mu_sigma; //  "    "     "

    G4double width_Evis = 0.0;
    for ( std::vector< G4double >::const_iterator cit = particle_vecEvis.begin();
	  cit != particle_vecEvis.end() ; ++cit ) {
      width_Evis += std::abs( *cit - particle_mu_Evis );
    }
    width_Evis *= std::sqrt( 3.141592654/2.0 ) / n;
    G4double width_Evis_sigma = width_Evis / std::sqrt( 2.0*(n - 1) );
    G4double energyResolution = 0.0;
    G4double energyResolution_sigma = 0.0;
    if ( particle_mu_Evis > 1.0E-06 ) {
      energyResolution = width_Evis / particle_mu_Evis;
      if ( width_Evis > 1.0E-06 ) {
	energyResolution_sigma = energyResolution *
	  std::sqrt( ( width_Evis_sigma * width_Evis_sigma ) / 
		     ( width_Evis * width_Evis ) 
		     +
		     ( particle_mu_Evis_sigma * particle_mu_Evis_sigma ) / 
		     ( particle_mu_Evis + particle_mu_Evis ) );
      }
    }
    G4cout << "\t \t sigma_Evis = " << width_Evis << " +/- " 
	   << width_Evis_sigma << G4endl
	   << "\t \t energy resolution = " << energyResolution << " +/- "  
	   << energyResolution_sigma << G4endl;

    width_Evis = 0.0;
    for ( std::vector< G4double >::const_iterator cit = no_particle_vecEvis.begin();
	  cit != no_particle_vecEvis.end() ; ++cit ) {
      width_Evis += std::abs( *cit - no_particle_mu_Evis );
    }
    width_Evis *= std::sqrt( 3.141592654/2.0 ) / n;
    width_Evis_sigma = width_Evis / std::sqrt( 2.0*(n - 1) );
    energyResolution = 0.0;
    energyResolution_sigma = 0.0;
    if ( no_particle_mu_Evis > 1.0E-06 ) {
      energyResolution = width_Evis / no_particle_mu_Evis;
      if ( width_Evis > 1.0E-06 ) {
	energyResolution_sigma = energyResolution *
	  std::sqrt( ( width_Evis_sigma * width_Evis_sigma ) / 
		     ( width_Evis * width_Evis ) 
		     +
		     ( no_particle_mu_Evis_sigma * no_particle_mu_Evis_sigma ) / 
		     ( no_particle_mu_Evis + no_particle_mu_Evis ) );
      }
    }
    G4cout << "\t \t (All contributions without this particle type: " << G4endl
           << "\t \t   -> sigma_Evis = " << width_Evis << " +/- " 
	   << width_Evis_sigma << G4endl
	   << "\t \t   -> energy resolution = " << energyResolution << " +/- "  
	   << energyResolution_sigma << G4endl
           << "\t \t )" << G4endl;

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
	caseName += "# nuclei";
	if ( i == 0 ) {
	  sum = numStepNucleus;
	  sum2 = numStepNucleus2;
	} else {
	  sum = numTrackNucleus;
	  sum2 = numTrackNucleus2;
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
    if ( mu_numSteps > 0.0  &&  trackLength > 1.0E-06*mm ) {
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

  // Print information about number of inelastic and stoppingAtRest processes.
  if ( StatAccepTestAnalysis::isCountingProcessesOn ) { 
    G4cout << G4endl
           << " Number of inelastic and stoppingAtRest processes (total, not per event)" 
	   << G4endl
           << "\t numInelasticProcesses = " 
	   << numInelasticProcesses << G4endl
           << "\t numProtonInelasticProcesses = " 
	   << numProtonInelasticProcesses << G4endl
	   << "\t numAntiProtonInelasticProcesses = " 
	   << numAntiProtonInelasticProcesses << G4endl
	   << "\t numNeutronInelasticProcesses = " 
	   << numNeutronInelasticProcesses << G4endl
	   << "\t numAntiNeutronInelasticProcesses = " 
	   << numAntiNeutronInelasticProcesses << G4endl
	   << "\t numPionPlusInelasticProcesses = " 
	   << numPionPlusInelasticProcesses << G4endl
	   << "\t numPionMinusInelasticProcesses = " 
	   << numPionMinusInelasticProcesses << G4endl
	   << "\t numKaonPlusInelasticProcesses = " 
	   << numKaonPlusInelasticProcesses << G4endl
	   << "\t numKaonMinusInelasticProcesses = " 
	   << numKaonMinusInelasticProcesses << G4endl
	   << "\t numKaonNeutralInelasticProcesses = " 
	   << numKaonNeutralInelasticProcesses << G4endl
	   << "\t numLambdaInelasticProcesses = " 
	   << numLambdaInelasticProcesses << G4endl
	   << "\t numAntiLambdaInelasticProcesses = " 
	   << numAntiLambdaInelasticProcesses << G4endl
	   << "\t numSigmaMinusInelasticProcesses = " 
	   << numSigmaMinusInelasticProcesses << G4endl
	   << "\t numMuonMinusInelasticProcesses = " 
	   << numMuonMinusInelasticProcesses << G4endl
	   << "\t numOtherInelasticProcesses = " 
	   << numOtherInelasticProcesses << G4endl
	   << "\t numStoppingAtRestProcesses = " 
	   << numStoppingAtRestProcesses << G4endl
	   << "\t numAntiProtonStoppingAtRestProcesses = " 
	   << numAntiProtonStoppingAtRestProcesses << G4endl
	   << "\t numNeutronStoppingAtRestProcesses = " 
	   << numNeutronStoppingAtRestProcesses << G4endl
	   << "\t numAntiNeutronStoppingAtRestProcesses = " 
	   << numAntiNeutronStoppingAtRestProcesses << G4endl
	   << "\t numPionMinusStoppingAtRestProcesses = " 
	   << numPionMinusStoppingAtRestProcesses << G4endl
	   << "\t numKaonMinusStoppingAtRestProcesses = " 
	   << numKaonMinusStoppingAtRestProcesses << G4endl
	   << "\t numKaonNeutralStoppingAtRestProcesses = " 
	   << numKaonNeutralStoppingAtRestProcesses << G4endl
	   << "\t numLambdaStoppingAtRestProcesses = " 
	   << numLambdaStoppingAtRestProcesses << G4endl
	   << "\t numAntiLambdaStoppingAtRestProcesses = " 
	   << numAntiLambdaStoppingAtRestProcesses << G4endl
	   << "\t numSigmaMinusStoppingAtRestProcesses = " 
	   << numSigmaMinusStoppingAtRestProcesses << G4endl
	   << "\t numMuonMinusStoppingAtRestProcesses = " 
	   << numMuonMinusStoppingAtRestProcesses << G4endl
	   << "\t numOtherStoppingAtRestProcesses = " 
	   << numOtherStoppingAtRestProcesses 
	   << G4endl;
    if ( numInelasticProcesses != 
	 ( numProtonInelasticProcesses + numAntiProtonInelasticProcesses +
	   numNeutronInelasticProcesses + numAntiNeutronInelasticProcesses +
	   numPionPlusInelasticProcesses + numPionMinusInelasticProcesses +
	   numKaonPlusInelasticProcesses + numKaonMinusInelasticProcesses +
	   numKaonNeutralInelasticProcesses + 
           numLambdaInelasticProcesses + numAntiLambdaInelasticProcesses +
	   numSigmaMinusInelasticProcesses +
	   numMuonMinusInelasticProcesses +
           numOtherInelasticProcesses ) ) {
      G4cout << " NOT CONSERVED #inelastic processes " << G4endl;
    }
    if ( numStoppingAtRestProcesses != 
	 ( numAntiProtonStoppingAtRestProcesses +
	   numNeutronStoppingAtRestProcesses + numAntiNeutronStoppingAtRestProcesses +
	   numPionMinusStoppingAtRestProcesses +
	   numKaonMinusStoppingAtRestProcesses +
	   numKaonNeutralStoppingAtRestProcesses + 
           numLambdaStoppingAtRestProcesses + numAntiLambdaStoppingAtRestProcesses +
	   numSigmaMinusStoppingAtRestProcesses +
	   numMuonMinusStoppingAtRestProcesses +
	   numOtherStoppingAtRestProcesses ) ) {
      G4cout << " NOT CONSERVED #stoppingAtRest processes " << G4endl;
    }
  }

  // Print information about the names of particles and the total 
  // number of times they have been created (excluding e-, e+, gamma).
  if ( StatAccepTestAnalysis::isMapParticleNamesOn ) { 
    G4cout << G4endl
           << " Particle names and total (not per event) tracks created \n"
           << " (excluding e-, e+, gamma)"  
	   << G4endl;
    for ( std::map< std::string, int >::const_iterator cit = mapParticleNames.begin(); 
	  cit != mapParticleNames.end(); ++cit ) {
      G4cout << "\t particle = " << cit->first 
	     << "\t total number of tracks = " << cit->second
	     << G4endl;
    }   
  }


  if ( StatAccepTestAnalysis::isMapInfoAboutTrackOn ) { 
    G4cout << G4endl << " Analysis of interactions (Tracks & Vertices)" << G4endl;

    G4cout << "\t average number of tracks = " << countNumberOfTracks / n 
	   << G4endl
	   << "\t average number of vertices = " << countNumberOfVertices / n 
	   << G4endl
	   << "\t    electromagnetic = " << countNumberOfElectromagneticVertices / n 
	   << G4endl
	   << "\t    photolepton-hadron = " << countNumberOfPhotoleptonHadronVertices / n
	   << G4endl
	   << "\t    decay = " << countNumberOfDecayVertices / n 
	   << G4endl
	   << "\t    hadronic elastic : from neutrons = " 
	   << countNumberOfHadronElasticVertices_fromNeutrons / n 
	   << "   from others = " 
	   << countNumberOfHadronElasticVertices_notFromNeutrons / n
	   << G4endl
	   << "\t \t hadron inelastic : from neutrons = " 
	   << countNumberOfHadronInelasticVertices_fromNeutrons / n
	   << "   from others = " 
	   << countNumberOfHadronInelasticVertices_notFromNeutrons / n
	   << G4endl
	   << "\t \t capture : from neutrons = " 
	   << countNumberOfCaptureVertices_fromNeutrons / n
	   << "   from others = " 
	   << countNumberOfCaptureVertices_notFromNeutrons / n
	   << G4endl
	   << "\t \t fission : from neutrons = " 
	   << countNumberOfFissionVertices_fromNeutrons / n
	   << "   from others = " 
	   << countNumberOfFissionVertices_notFromNeutrons / n
	   << G4endl
	   << "\t \t at rest : from neutrons = " 
	   << countNumberOfAtRestVertices_fromNeutrons / n
	   << "   from others = " 
	   << countNumberOfAtRestVertices_notFromNeutrons / n
	   << G4endl
	   << "\t \t charge exchange : from neutrons = " 
	   << countNumberOfChargeExchangeVertices_fromNeutrons / n
	   << "   from others = " 
	   << countNumberOfChargeExchangeVertices_notFromNeutrons / n
	   << G4endl
	   << G4endl;

    G4cout << "\t <Evis_em>  = " << sumEvis_em  / n 
	   << "   <Etot_em>  = " << sumEtot_em  / n << "  MeV" << G4endl
           << "\t <Evis_p>   = " << sumEvis_p   / n 
	   << "   <Etot_p>   = " << sumEtot_p   / n << "  MeV" << G4endl
           << "\t <Evis_pi>  = " << sumEvis_pi  / n 
	   << "   <Etot_pi>  = " << sumEtot_pi  / n << "  MeV" << G4endl
           << "\t <Evis_ion> = " << sumEvis_ion / n 
	   << "   <Etot_ion> = " << sumEtot_ion / n << "  MeV" << G4endl
           << G4endl;

    G4cout << "\t <Evis_from1stInterac_pi0> = " << sumEvis_from1stInterac_pi0 / n
           << "  MeV \t ( f_Evis = " << 100.0 * sumEvis_from1stInterac_pi0 / sumEdepAct
           << "% )" << G4endl
	   << "\t <Etot_from1stInterac_pi0> = " << sumEtot_from1stInterac_pi0 / n
           << "  MeV \t ( f_Ebeam = " 
           << 100.0 * sumEtot_from1stInterac_pi0 / ( n * beamEnergy ) 
	   << "% )" << G4endl
           << "\t <Evis_from1stInterac_pi+> = " << sumEvis_from1stInterac_pip / n
           << "  MeV \t ( f_Evis = " << 100.0 * sumEvis_from1stInterac_pip / sumEdepAct
           << "% )" << G4endl
	   << "\t <Etot_from1stInterac_pi+> = " << sumEtot_from1stInterac_pip / n
           << "  MeV \t ( f_Ebeam = " 
           << 100.0 * sumEtot_from1stInterac_pip / ( n * beamEnergy ) 
	   << "% )" << G4endl
           << "\t <Evis_from1stInterac_pi-> = " << sumEvis_from1stInterac_pim / n
           << "  MeV \t ( f_Evis = " << 100.0 * sumEvis_from1stInterac_pim / sumEdepAct
           << "% )" << G4endl
	   << "\t <Etot_from1stInterac_pi-> = " << sumEtot_from1stInterac_pim / n
           << "  MeV \t ( f_Ebeam = " 
           << 100.0 * sumEtot_from1stInterac_pim / ( n * beamEnergy ) 
	   << "% )" << G4endl
           << "\t <Evis_from1stInterac_p> = " << sumEvis_from1stInterac_p / n
           << "  MeV \t ( f_Evis = " << 100.0 * sumEvis_from1stInterac_p / sumEdepAct
           << "% )" << G4endl
	   << "\t <Etot_from1stInterac_p> = " << sumEtot_from1stInterac_p / n
           << "  MeV \t ( f_Ebeam = " 
           << 100.0 * sumEtot_from1stInterac_p / ( n * beamEnergy ) 
	   << "% )" << G4endl
           << "\t <Evis_from1stInterac_n> = " << sumEvis_from1stInterac_n / n
           << "  MeV \t ( f_Evis = " << 100.0 * sumEvis_from1stInterac_n / sumEdepAct
           << "% )" << G4endl
	   << "\t <Etot_from1stInterac_n> = " << sumEtot_from1stInterac_n / n
           << "  MeV \t ( f_Ebeam = " 
           << 100.0 * sumEtot_from1stInterac_n / ( n * beamEnergy ) 
	   << "% )" << G4endl
           << "\t <Evis_from1stInterac_lightion> = " 
	   << sumEvis_from1stInterac_lightion / n << "  MeV \t ( f_Evis = " 
	   << 100.0 * sumEvis_from1stInterac_lightion / sumEdepAct
           << "% )" << G4endl
	   << "\t <Etot_from1stInterac_lightion> = " 
	   << sumEtot_from1stInterac_lightion / n << "  MeV \t ( f_Ebeam = " 
           << 100.0 * sumEtot_from1stInterac_lightion / ( n * beamEnergy ) 
	   << "% )" << G4endl
	   << G4endl;

    G4cout << "\t <Evis_closest_pi0> = " << sumEvis_closest_pi0 / n 
	   << " MeV \t ( f_Evis = " 
	   << 100.0 * sumEvis_closest_pi0 / sumEdepAct
	   << "% )" << G4endl
           << "\t <Etot_closest_pi0> = " << sumEtot_closest_pi0 / n 
	   << " MeV \t ( f_Ebeam = " 
	   << 100.0 * sumEtot_closest_pi0 / ( n * beamEnergy )
	   << "% )" << G4endl
           << "\t <Evis_closest_pi+> = " << sumEvis_closest_pip / n 
	   << " MeV \t ( f_Evis = " 
	   << 100.0 * sumEvis_closest_pip / sumEdepAct
	   << "% )" << G4endl
           << "\t <Etot_closest_pi+> = " << sumEtot_closest_pip / n 
	   << " MeV \t ( f_Ebeam = " 
	   << 100.0 * sumEtot_closest_pip / ( n * beamEnergy )
	   << "% )" << G4endl
           << "\t <Evis_closest_pi-> = " << sumEvis_closest_pim / n 
	   << " MeV \t ( f_Evis = " 
	   << 100.0 * sumEvis_closest_pim / sumEdepAct 
	   << "% )" << G4endl
           << "\t <Etot_closest_pi-> = " << sumEtot_closest_pim / n 
	   << " MeV \t ( f_Ebeam = " 
	   << 100.0 * sumEtot_closest_pim / ( n * beamEnergy )
	   << "% )" << G4endl
           << "\t <Evis_closest_p> = " << sumEvis_closest_p / n 
	   << " MeV \t ( f_Evis = " 
	   << 100.0 * sumEvis_closest_p / sumEdepAct 
	   << "% )" << G4endl
           << "\t <Etot_closest_p> = " << sumEtot_closest_p / n 
	   << " MeV \t ( f_Ebeam = " 
	   << 100.0 * sumEtot_closest_p / ( n * beamEnergy )
	   << "% )" << G4endl
           << "\t <Evis_closest_n> = " << sumEvis_closest_n / n 
	   << " MeV \t ( f_Evis = " 
	   << 100.0 * sumEvis_closest_n / sumEdepAct 
	   << "% )" << G4endl
           << "\t <Etot_closest_n> = " << sumEtot_closest_n / n 
	   << " MeV \t ( f_Ebeam = " 
	   << 100.0 * sumEtot_closest_n / ( n * beamEnergy )
	   << "% )" << G4endl
           << "\t <Evis_closest_lightion> = " << sumEvis_closest_lightion / n 
	   << " MeV \t ( f_Evis = " 
	   << 100.0 * sumEvis_closest_lightion / sumEdepAct 
	   << "% )" << G4endl
           << "\t <Etot_closest_lightion> = " << sumEtot_closest_lightion / n 
	   << " MeV \t ( f_Ebeam = " 
	   << 100.0 * sumEtot_closest_lightion / ( n * beamEnergy )
	   << "% )" << G4endl
	   << G4endl;

    G4cout << "\t <Evis_em_from_pi0> = " << sumEvis_em_from_pi0 / n 
	   << " MeV \t ( f_Evis_em = " << 100.0 * sumEvis_em_from_pi0 / sumEvis_em 
	   << "% )" << G4endl
           << "\t <Etot_em_from_pi0> = " << sumEtot_em_from_pi0 / n 
	   << " MeV \t ( f_Etot_em = " << 100.0 * sumEtot_em_from_pi0 / sumEtot_em
	   << "% )" << G4endl
           << "\t <Evis_em_from_pi+> = " << sumEvis_em_from_pip / n 
	   << " MeV \t ( f_Evis_em = " << 100.0 * sumEvis_em_from_pip / sumEvis_em 
	   << "% )" << G4endl
           << "\t <Etot_em_from_pi+> = " << sumEtot_em_from_pip / n 
	   << " MeV \t ( f_Etot_em = " << 100.0 * sumEtot_em_from_pip / sumEtot_em
	   << "% )" << G4endl
           << "\t <Evis_em_from_pi-> = " << sumEvis_em_from_pim / n 
	   << " MeV \t ( f_Evis_em = " << 100.0 * sumEvis_em_from_pim / sumEvis_em 
	   << "% )" << G4endl
           << "\t <Etot_em_from_pi-> = " << sumEtot_em_from_pim / n 
	   << " MeV \t ( f_Etot_em = " << 100.0 * sumEtot_em_from_pim / sumEtot_em
	   << "% )" << G4endl
           << "\t <Evis_em_from_p> = " << sumEvis_em_from_p / n 
	   << " MeV \t ( f_Evis_em = " << 100.0 * sumEvis_em_from_p / sumEvis_em 
	   << "% )" << G4endl
           << "\t <Etot_em_from_p> = " << sumEtot_em_from_p / n 
	   << " MeV \t ( f_Etot_em = " << 100.0 * sumEtot_em_from_p / sumEtot_em
	   << "% )" << G4endl
           << "\t <Evis_em_from_n> = " << sumEvis_em_from_n / n 
	   << " MeV \t ( f_Evis_em = " << 100.0 * sumEvis_em_from_n / sumEvis_em 
	   << "% )" << G4endl
           << "\t <Etot_em_from_n> = " << sumEtot_em_from_n / n 
	   << " MeV \t ( f_Etot_em = " << 100.0 * sumEtot_em_from_n / sumEtot_em
	   << "% )" << G4endl
           << "\t <Evis_em_from_lightion> = " << sumEvis_em_from_lightion / n 
	   << " MeV \t ( f_Evis_em = " << 100.0 * sumEvis_em_from_lightion / sumEvis_em 
	   << "% )" << G4endl
           << "\t <Etot_em_from_lightion> = " << sumEtot_em_from_lightion / n 
	   << " MeV \t ( f_Etot_em = " << 100.0 * sumEtot_em_from_lightion / sumEtot_em
	   << "% )" << G4endl
	   << G4endl;
  }


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


void StatAccepTestAnalysis::analysisTrackAndVertices_0() {
  // Event-by-event debugging information.
  // Be careful that it writes a lot of information, so you should 
  // use it only for few events.

  G4cout << G4endl
	 << "\t =============================================== " << G4endl
	 << "\t ===  END OF EVENT INFO ABOUT TRACK & VERTEX === " << G4endl
	 << "\t =============================================== " << G4endl
	 << G4endl;

  G4cout << G4endl
	 << "\t          *** TRACKS *** "  << G4endl << G4endl
	 << "\t mapInfoAboutTrack.size()=" << mapInfoAboutTrack.size() << G4endl
	 << G4endl;
  for ( std::map< G4int, structInfoAboutTrack >::const_iterator 
	  cit = mapInfoAboutTrack.begin(); 
	cit != mapInfoAboutTrack.end(); ++cit ) {
    G4cout << "\t -------------------------------------------------------" << G4endl
	   << "\t trackId=" << cit->first << G4endl
	   << "\t particlePDGCode=" << cit->second.theParticlePDGCode 
	   << "\t particleName=" << cit->second.theParticleName << G4endl
	   << "\t trackEkinAtCreation=" << cit->second.theTrackEkinAtCreation 
	   << " MeV" << G4endl
	   << "\t trackVisibleEdep=" << cit->second.theTrackVisibleEdep 
	   << " MeV" << G4endl
	   << "\t trackDepositedEtot=" << cit->second.theTrackDepositedEtot 
	   << " MeV" << G4endl
	   << "\t parentTrackId=" << cit->second.theParentTrackID << G4endl 
	   << "\t closestHadronicRelativeId=" 
	   << cit->second.theClosestHadronicRelativeID << G4endl
	   << "\t closestHadronicVertexId="
	   << cit->second.theClosestHadronicVertexID << G4endl
	   << "\t startingVertexId=" << cit->second.theStartingVertexID
	   << G4endl;
  }   

  G4cout << G4endl << G4endl 
	 << "\t          *** VERTICES *** "  << G4endl << G4endl
	 << "\t mapInfoAboutVertex.size()=" << mapInfoAboutVertex.size() << G4endl
	 << G4endl;
  for ( std::map< G4int, structInfoAboutVertex >::const_iterator 
	  cit = mapInfoAboutVertex.begin(); 
	cit != mapInfoAboutVertex.end(); ++cit ) {
    G4cout << "\t -------------------------------------------------------" << G4endl
	   << "\t vertexId=" << cit->first << G4endl
	   << "\t parentTrackId=" << cit->second.theCreatorTrackID << G4endl
	   << "\t creatorProcessType=" << cit->second.theCreatorProcessType << G4endl
	   << "\t creatorProcessSubType=" << cit->second.theCreatorProcessSubType 
	   << G4endl
	   << "\t creatorProcessName=" << cit->second.theCreatorProcessName << G4endl
	   << "\t vertexVolumeName=" << cit->second.theVertexVolumeName << G4endl
	   << "\t theVertexPosition=(" 
	   << cit->second.theVertexPosition_x << ","
	   << cit->second.theVertexPosition_y << "," 
	   << cit->second.theVertexPosition_z << ") mm" << G4endl
	   << "\t global time=" << cit->second.theVertexTime << " nsec" << G4endl
	   << G4endl;
  }

  G4cout << G4endl << G4endl 
	 << "\t          *** MOTHER/CHILDREN *** "  << G4endl << G4endl
	 << G4endl;
  for ( std::map< G4int, structInfoAboutVertex >::const_iterator 
	  citvtx = mapInfoAboutVertex.begin(); 
	citvtx != mapInfoAboutVertex.end(); ++citvtx ) {
    G4int mother = citvtx->second.theCreatorTrackID;
    G4cout << "\t -------------------------------------------------------" << G4endl
	   << "\t vertexId=" << citvtx->first 
	   << "\t mother=" << mother;
    std::map< G4int, structInfoAboutTrack >::const_iterator cittrk = 
      mapInfoAboutTrack.find( mother ); 
    if ( cittrk != mapInfoAboutTrack.end() ) {
      G4cout << "  " << cittrk->second.theParticleName
	     << "  vtx: " << cittrk->second.theStartingVertexID
	     << "  (closestHadronicRelative: " 
	     << cittrk->second.theClosestHadronicRelativeID;
      std::map< G4int, structInfoAboutTrack >::const_iterator cittrk_bis = 
	mapInfoAboutTrack.find( cittrk->second.theClosestHadronicRelativeID );
      if ( cittrk_bis != mapInfoAboutTrack.end() ) {
	G4cout << "  " << cittrk_bis->second.theParticleName;
      }
      G4cout << "  vtx: " << cittrk->second.theClosestHadronicVertexID
	     << ")" << G4endl;
    }
    G4cout << "\t -> children:" << G4endl;
    for ( ++cittrk ; cittrk != mapInfoAboutTrack.end() ; ++cittrk ) {
      if ( cittrk->second.theStartingVertexID == citvtx->first ) {
	G4cout << "\t \t" << cittrk->first 
	       << "  " << cittrk->second.theParticleName 
	       << "  Ekin=" << cittrk->second.theTrackEkinAtCreation
	       << "  parentId=" << cittrk->second.theParentTrackID
	       << G4endl;
      }
    }
  }
  G4cout << G4endl
	 << "\t =============================================== " << G4endl
	 << G4endl;

}


void StatAccepTestAnalysis::analysisTrackAndVertices_1() {
  // Some basic quantities are computed. 
  // Some of these can be compared with independent estimations made 
  // in other parts of this file, as a kind of sanity check.

  countNumberOfTracks += mapInfoAboutTrack.size();
  countNumberOfVertices += mapInfoAboutVertex.size();
  for ( std::map< G4int, structInfoAboutVertex >::const_iterator 
	  cit = mapInfoAboutVertex.begin(); 
	cit != mapInfoAboutVertex.end(); ++cit ) {
    G4int type = cit->second.theCreatorProcessType;
    G4int subtype = cit->second.theCreatorProcessSubType;
    G4int creator_pdgcode = 0;
    if ( mapInfoAboutTrack.find( cit->second.theCreatorTrackID ) !=  
	 mapInfoAboutTrack.end() ) {
      creator_pdgcode = mapInfoAboutTrack.find( cit->second.theCreatorTrackID )->
	second.theParticlePDGCode;
    }
    if ( type == 2 ) {            // G4ProcessType::fElectromagnetic
      countNumberOfElectromagneticVertices += 1;
    } else if ( type == 5 ) {     // G4ProcessType::fPhotolepton_hadron
      countNumberOfPhotoleptonHadronVertices += 1;
    } else if ( type == 6 ) {     // G4ProcessType::fDecay
      countNumberOfDecayVertices += 1;
    }
    if ( subtype == 111 ) {  // G4HadronicProcessType::fHadronElastic
      if ( creator_pdgcode == neutronId ) {
	countNumberOfHadronElasticVertices_fromNeutrons += 1;
      } else {
	countNumberOfHadronElasticVertices_notFromNeutrons += 1;
      }
    } else {
      if (
	  subtype == 121  ||     // G4HadronicProcessType::fHadronInelastic
	  subtype == 131  ||     // G4HadronicProcessType::fCapture
	  subtype == 141  ||     // G4HadronicProcessType::fFission
	  subtype == 151  ||     // G4HadronicProcessType::fHadronAtRest
	  subtype == 161         // G4HadronicProcessType::fChargeExchange 
	  ) {
	countNumberOfHadronicVertices++;
	if ( subtype == 121 ) {
	  if ( creator_pdgcode == neutronId ) {
	    countNumberOfHadronInelasticVertices_fromNeutrons += 1;
	  } else {
	    countNumberOfHadronInelasticVertices_notFromNeutrons += 1;
	  }
	} else if ( subtype == 131 ) {
	  if ( creator_pdgcode == neutronId ) {
	    countNumberOfCaptureVertices_fromNeutrons += 1;
	  } else {
	    countNumberOfCaptureVertices_notFromNeutrons += 1;
	  }
	} else if ( subtype == 141 ) {
	  if ( creator_pdgcode == neutronId ) {
	    countNumberOfFissionVertices_fromNeutrons += 1;
	  } else {
	    countNumberOfFissionVertices_notFromNeutrons += 1;
	  }
	} else if ( subtype == 151 ) {
	  if ( creator_pdgcode == neutronId ) {
	    countNumberOfAtRestVertices_fromNeutrons += 1;
	  } else {
	    countNumberOfAtRestVertices_notFromNeutrons += 1;
	  }
	} else if ( subtype == 161 ) {
	  if ( creator_pdgcode == neutronId ) {
	    countNumberOfChargeExchangeVertices_fromNeutrons += 1;
	  } else {
	    countNumberOfChargeExchangeVertices_notFromNeutrons += 1;
	  }
	}
      } 
    }
  }
  
  for ( std::map< G4int, structInfoAboutTrack >::const_iterator 
	  cit = mapInfoAboutTrack.begin(); 
	cit != mapInfoAboutTrack.end(); ++cit ) {
    G4int pdgcode = cit->second.theParticlePDGCode;
    G4double evis = cit->second.theTrackVisibleEdep;
    G4double etot = cit->second.theTrackDepositedEtot;
    if ( pdgcode == electronId  ||  pdgcode == positronId ) {
      sumEvis_em += evis;
      sumEtot_em += etot;
    } else if ( pdgcode == pionPlusId  ||  pdgcode == pionMinusId ) {
      sumEvis_pi += evis;
      sumEtot_pi += etot;
    } else if ( pdgcode == protonId  ||  pdgcode == antiProtonId ) {
      sumEvis_p += evis;
      sumEtot_p += etot;
    } else if ( pdgcode == 0  ||  pdgcode/1000000000 >= 1  ||  pdgcode == neutronId ) {
      sumEvis_ion += evis;
      sumEtot_ion += etot;
    }
  }
}


void StatAccepTestAnalysis::analysisTrackAndVertices_2() {
  // We try here to calculate some quantities that should help linking
  // the calorimeter observables with the quantities we have calculated
  // at model-level.
  // We consider only the first hadronic inelastic interaction,
  // and we look at the visible and total energies in the calorimeter
  // that are coming from the particles produced directly, or
  // indirectly (i.e. later generations) from that first hadronic 
  // inelastic interaction.

  for ( std::map< G4int, structInfoAboutTrack >::const_iterator 
	  cit = mapInfoAboutTrack.begin(); 
	cit != mapInfoAboutTrack.end(); ++cit ) {

    G4double evis = cit->second.theTrackVisibleEdep;
    G4double etot = cit->second.theTrackDepositedEtot;

    // For each track, we look for its ancestors until we find a 
    // daughter track of the primary particle (that always has 
    // track id = 1).
    std::map< G4int, structInfoAboutTrack >::const_iterator cit_bis = cit;
    bool isFound = false;
    G4int parentId = cit->second.theParentTrackID;
    while ( parentId > 0  &&  ! isFound ) {
      if ( parentId == 1 ) {  // primary beam particle
	isFound = true;
      } else {
	cit_bis = mapInfoAboutTrack.find( parentId );
	if ( cit_bis != mapInfoAboutTrack.end() ) {
	  parentId = cit_bis->second.theParentTrackID;
        } else {
	  G4cout << "*** StatAccepTestAnalysis::analysisTrackAndVertices_2 : WARNING ***"
		 << "\t parentId = " << parentId << " NOT found in the map!" << G4endl;
	  parentId = 0;
	}
      }
    }

    //G4cout << "\t Examining track_id = " << cit->first              //***DEBUG*** 
    //	     << "\t" << cit->second.theParticleName << G4endl; 
    if ( isFound ) {
      //G4cout << "\t \t FOUND 1st interaction relative: " << cit_bis->first 
      //       << "\t" << cit_bis->second.theParticleName << G4endl;  //***DEBUG***
      G4int pdgcode = cit_bis->second.theParticlePDGCode;
      G4String name = cit_bis->second.theParticleName;
      if ( pdgcode == pionZeroId ) {
        sumEvis_from1stInterac_pi0 += evis;
	sumEtot_from1stInterac_pi0 += etot;
      } else if ( pdgcode == pionPlusId ) {  
        sumEvis_from1stInterac_pip += evis;
	sumEtot_from1stInterac_pip += etot;
      } else if ( pdgcode == pionMinusId ) {
        sumEvis_from1stInterac_pim += evis;
	sumEtot_from1stInterac_pim += etot;
      } else if ( pdgcode == protonId ) {
        sumEvis_from1stInterac_p += evis;
	sumEtot_from1stInterac_p += etot;
      } else if ( pdgcode == neutronId ) {
        sumEvis_from1stInterac_n += evis;
	sumEtot_from1stInterac_n += etot;
      } else if ( name == "deuteron"  || 
		  name == "triton"    ||
		  name == "alpha"     ||
		  name == "He3" ) {
        sumEvis_from1stInterac_lightion += evis;
	sumEtot_from1stInterac_lightion += etot;
      }
    } // end if ( isFound )
  } // End of loop over mapInfoAboutTrack
}


void StatAccepTestAnalysis::analysisTrackAndVertices_3() {
  // We calculate the visible and total energy in the calorimeter
  // due to pi0, pi+, pi-, p, n, light ion (d, t, 3He, alpha),
  // when these particles are either directly depositing energy
  // (via ionization) or the closest hadron of the actual
  // track that deposited energy.
  // Also for non-light nuclear fragments (i.e. excluding d, t, 3He, alpha)
  // we consider the closest hadron that actually produced them
  // (often neutrons).
  // Notice that in the previous analysis (analysisTrackAndVertices_2)
  // we have classified the visible and total energy in the calorimeter
  // according to only the hadrons produced directly by the primary
  // hadronic inelastic interaction, whereas here we are considering
  // all hadrons, no matter when they have been generated.

  for ( std::map< G4int, structInfoAboutTrack >::const_iterator 
	  cit = mapInfoAboutTrack.begin(); 
	cit != mapInfoAboutTrack.end(); ++cit ) {

    G4double evis = cit->second.theTrackVisibleEdep;
    G4double etot = cit->second.theTrackDepositedEtot; 

    // For each track, we look for its ancestors until we find a 
    // track of the following type: pi0, or pi+, or pi-, or p, or n,
    // or light ion (d, t, 3He, alpha).
    bool isFound = false;
    std::map< G4int, structInfoAboutTrack >::const_iterator cit_bis = cit;
    G4int parentId = cit->second.theParentTrackID;
    G4int pdgcode = 0;
    G4String name = "";
    while ( parentId > 0  &&  ! isFound ) {
      pdgcode = cit_bis->second.theParticlePDGCode;
      name = cit_bis->second.theParticleName;
      if ( pdgcode == pionZeroId  ||
	   pdgcode == pionPlusId  ||  
	   pdgcode == pionMinusId ||
	   pdgcode == protonId    ||
	   pdgcode == neutronId   ||
	   name    == "deuteron"  || 
	   name    == "triton"    ||
	   name    == "alpha"     ||
	   name    == "He3"      ) {
	isFound = true;
      } else {
	cit_bis = mapInfoAboutTrack.find( parentId );
	if ( cit_bis != mapInfoAboutTrack.end() ) {
	  parentId = cit_bis->second.theParentTrackID;
        } else {
	  G4cout << "*** StatAccepTestAnalysis::analysisTrackAndVertices_3 : WARNING ***"
		 << "\t parentId = " << parentId << " NOT found in the map!" << G4endl;
	  parentId = 0;
	}
      }
    }

    //G4cout << "\t Examining track_id = " << cit->first              //***DEBUG*** 
    //       << "\t" << cit->second.theParticleName << G4endl; 
    if ( isFound ) {
      //G4cout << "\t \t FOUND  closest hadron: " << cit_bis->first 
      //       << "\t" << cit_bis->second.theParticleName << G4endl;  //***DEBUG***
      if ( pdgcode == pionZeroId ) {
         sumEvis_closest_pi0 += evis;
	 sumEtot_closest_pi0 += etot;
      } else if ( pdgcode == pionPlusId ) {  
         sumEvis_closest_pip += evis;
	 sumEtot_closest_pip += etot;
      } else if ( pdgcode == pionMinusId ) {
         sumEvis_closest_pim += evis;
	 sumEtot_closest_pim += etot;
      } else if ( pdgcode == protonId ) {
         sumEvis_closest_p += evis;
	 sumEtot_closest_p += etot;
      } else if ( pdgcode == neutronId ) {
         sumEvis_closest_n += evis;
	 sumEtot_closest_n += etot;
      } else if ( name == "deuteron"  || 
		  name == "triton"    ||
		  name == "alpha"     ||
		  name == "He3" ) {
         sumEvis_closest_lightion += evis;
	 sumEtot_closest_lightion += etot;
      }
    } // end if ( isFound )
  } // End of loop over mapInfoAboutTrack
}


void StatAccepTestAnalysis::analysisTrackAndVertices_4() {
  // We know that the dominat component of the visible and
  // and total energy deposited in a calorimeter is due to
  // electromagnetic, i.e. ionization of e- and e+.
  // We want to investigate here which are the closest hadrons
  // that produce those e- and e+ : we expect that the largest
  // contribution would come from pi0, but we want to see how
  // much and who produces the rest.

  for ( std::map< G4int, structInfoAboutTrack >::const_iterator 
	  cit = mapInfoAboutTrack.begin(); 
	cit != mapInfoAboutTrack.end(); ++cit ) {

    // Skip any other track but electron and positron.
    if ( cit->second.theParticlePDGCode == electronId  ||
         cit->second.theParticlePDGCode == positronId ) {
	 
      G4double evis = cit->second.theTrackVisibleEdep;
      G4double etot = cit->second.theTrackDepositedEtot; 

      // We look for the ancestors of the selected track (e- or e+)
      // until we find a track of the following type: 
      // pi0, or pi+, or pi-, or p, or n, or 
      // light ion (d, t, 3He, alpha).
      bool isFound = false;
      std::map< G4int, structInfoAboutTrack >::const_iterator cit_bis = cit;
      G4int parentId = cit->second.theParentTrackID;
      G4int pdgcode = 0;
      G4String name = "";
      while ( parentId > 0  &&  ! isFound ) {
	pdgcode = cit_bis->second.theParticlePDGCode;
	name = cit_bis->second.theParticleName;
	// The first time the if statement below is evaluated
	// it is not satisfied for sure, because the track is
	// e- or e+. However, we prefer to keep the structure
	// of the loop exactly as it was for the previous
	// analysis (analysisTrackAndVertices_3).
	if ( pdgcode == pionZeroId  ||
	     pdgcode == pionPlusId  ||  
	     pdgcode == pionMinusId ||
	     pdgcode == protonId    ||
	     pdgcode == neutronId   ||
	     name    == "deuteron"  || 
	     name    == "triton"    ||
	     name    == "alpha"     ||
	     name    == "He3"      ) {
	  isFound = true;
	} else {
	  cit_bis = mapInfoAboutTrack.find( parentId );
	  if ( cit_bis != mapInfoAboutTrack.end() ) {
	    parentId = cit_bis->second.theParentTrackID;
	  } else {
	    G4cout << "***StatAccepTestAnalysis::analysisTrackAndVertices_4: WARNING ***"
		   << "\t parentId = " << parentId << " NOT found in the map!" << G4endl;
	    parentId = 0;
	  }
	}
      }

      //G4cout << "\t Examining track_id = " << cit->first              //***DEBUG*** 
      //       << "\t" << cit->second.theParticleName << G4endl; 
      if ( isFound ) {
	//G4cout << "\t \t FOUND  em closest hadron : " << cit_bis->first 
	//       << "\t" << cit_bis->second.theParticleName << G4endl;  //***DEBUG***
	if ( pdgcode == pionZeroId ) {
	  sumEvis_em_from_pi0 += evis;
	  sumEtot_em_from_pi0 += etot;
	} else if ( pdgcode == pionPlusId ) {  
	  sumEvis_em_from_pip += evis;
	  sumEtot_em_from_pip += etot;
	} else if ( pdgcode == pionMinusId ) {
	  sumEvis_em_from_pim += evis;
	  sumEtot_em_from_pim += etot;
	} else if ( pdgcode == protonId ) {
	  sumEvis_em_from_p += evis;
	  sumEtot_em_from_p += etot;
	} else if ( pdgcode == neutronId ) {
	  sumEvis_em_from_n += evis;
	  sumEtot_em_from_n += etot;
	} else if ( name == "deuteron"  || 
		    name == "triton"    ||
		    name == "alpha"     ||
		    name == "He3" ) {
	  sumEvis_em_from_lightion += evis;
	  sumEtot_em_from_lightion += etot;
	}
      } // end if ( isFound )
    } // end if electron or positron
  } // End of loop over mapInfoAboutTrack
}
