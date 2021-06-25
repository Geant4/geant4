//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file Run.cc
/// \brief Implementation of the Run class
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Run.hh"
#include "G4SystemOfUnits.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run() : G4Run(),
             fNumEvents( 0 ),
             fPrimaryParticleId( 0 ),
             fPrimaryParticleInitialKineticEnergy( 0.0 ),
             fPrimaryParticleInitialTotalEnergy( 0.0 ),
             fPrimaryParticleInitialMomentum( 0.0 ),
             fPrimaryParticleInitialBeta( 0.0 ),
             fPrimaryParticleInitialGamma( 0.0 ),
             fPrimaryParticleInitial3Momentum( G4ThreeVector() ),
             fPrimaryParticleInitialPosition( G4ThreeVector() ),
             fToleranceEPviolations( 0.0 ),
             fToleranceDeltaDecayRadius( 0.0 ),
             fIsPreassignedDecayEnabled( true ),
             fIsBoostToLabEnabled( true ),
             fNumDecays( 0 ),
             fNumBadDecays( 0 ),
             fNumUnexpectedDecays( 0 ),
             fNumEviolations( 0 ),
             fNumPviolations( 0 ),
             fNum_mc_truth_rPos_deltaMax_above( 0 ),
             fNum_underestimated_mc_truth_rPos_delta_above( 0 ),
             fNum_overestimated_mc_truth_rPos_delta_above( 0 ),
             fNumLargeUnderestimates( 0 ),
             fNumLargeOverestimates( 0 ),
             fDecayT( 0.0 ),
             fSumDecayT( 0.0 ),
             fMinDecayT(  999999.9 ),
             fMaxDecayT( -999999.9 ),
             fDecayR_mc_truth( 0.0 ),
             fDecayR( 0.0 ),
             fSumDecayR( 0.0 ),
             fMinDecayR(  999999.9 ),
             fMaxDecayR( -999999.9 ),
             fDecayX( 0.0 ),
             fSumDecayX( 0.0 ),
             fMinDecayX(  999999.9 ),
             fMaxDecayX( -999999.9 ),
             fDecayY( 0.0 ),
             fSumDecayY( 0.0 ),
             fMinDecayY(  999999.9 ),
             fMaxDecayY( -999999.9 ),
             fDecayZ( 0.0 ),
             fSumDecayZ( 0.0 ),
             fMinDecayZ(  999999.9 ),
             fMaxDecayZ( -999999.9 ),
             fDeltaDecayR( 0.0 ),
             fSumDeltaDecayR( 0.0 ),
             fMinDeltaDecayR(  999999.9 ),
             fMaxDeltaDecayR( -999999.9 ),
             fDeflectionAngle( 0.0 ),
             fSumDeflectionAngle( 0.0 ),
             fMinDeflectionAngle(  999999.9 ),
             fMaxDeflectionAngle( -999999.9 ),
             fDeltaEkin( 0.0 ),
             fSumDeltaEkin( 0.0 ),
             fMinDeltaEkin(  999999.9 ),
             fMaxDeltaEkin( -999999.9 ),
             fDecayEkin( 0.0 ),
             fSumDecayEkin( 0.0 ),
             fMinDecayEkin(  999999.9 ),
             fMaxDecayEkin( -999999.9 ),
             fDecayPx( 0.0 ),
             fSumDecayPx( 0.0 ),
             fMinDecayPx(  999999.9 ),
             fMaxDecayPx( -999999.9 ),
             fDecayPy( 0.0 ),
             fSumDecayPy( 0.0 ),
             fMinDecayPy(  999999.9 ),
             fMaxDecayPy( -999999.9 ),
             fDecayPz( 0.0 ),
             fSumDecayPz( 0.0 ),
             fMinDecayPz(  999999.9 ),
             fMaxDecayPz( -999999.9 ),
             fDecayEtotViolation( 0.0 ),
             fSumDecayEtotViolation( 0.0 ),
             fMinDecayEtotViolation(  999999.9 ),
             fMaxDecayEtotViolation( -999999.9 ),
             fDecayPxViolation( 0.0 ),
             fSumDecayPxViolation( 0.0 ),
             fMinDecayPxViolation(  999999.9 ),
             fMaxDecayPxViolation( -999999.9 ),
             fDecayPyViolation( 0.0 ),
             fSumDecayPyViolation( 0.0 ),
             fMinDecayPyViolation(  999999.9 ),
             fMaxDecayPyViolation( -999999.9 ),
             fDecayPzViolation( 0.0 ),
             fSumDecayPzViolation( 0.0 ),
             fMinDecayPzViolation(  999999.9 ),
             fMaxDecayPzViolation( -999999.9 ),
             fMaxEkin_deltaMax( 0.0 ),
             fMaxEtot_deltaMax( 0.0 ),
             fMaxP_deltaMax( 0.0 ),
             fMaxPdir_deltaMax( 0.0 ),
             fMaxMass_deltaMax1( 0.0 ),
             fMaxMass_deltaMax2( 0.0 ),
             fSumMass_deltaMax3( 0.0 ),
             fMaxMass_deltaMax3( 0.0 ),
             fMaxBeta_deltaMax1( 0.0 ),
             fMaxBeta_deltaMax2( 0.0 ),
             fMaxGamma_deltaMax1( 0.0 ),
             fMaxGamma_deltaMax2( 0.0 ),
             fMaxGamma_deltaMax3( 0.0 ),
             fMaxT_proper_deltaMax( 0.0 ),
             fMaxT_lab_deltaMax( 0.0 ),
             fSumMc_truth_rPos_deltaMax( 0.0 ),  
             fMaxMc_truth_rPos_deltaMax( 0.0 ),
             fSumUnderestimated_mc_truth_rPos_delta( 0.0 ),
             fMinUnderestimated_mc_truth_rPos_delta( 999999.9 ),
             fSumOverestimated_mc_truth_rPos_delta( 0.0 ),
             fMaxOverestimated_mc_truth_rPos_delta( -999999.9 ),
             fSumUnderestimated_rDeltaPos( 0.0 ),
             fMinUnderestimated_rDeltaPos( 999999.9 ),
             fSumOverestimated_rDeltaPos( 0.0 ),
             fMaxOverestimated_rDeltaPos( -999999.9 ),
             fMaxFloat_rDeltaPos_deltaMax(-999999.9 ) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::~Run() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::RecordEvent( const G4Event* anEvent ) {
  // This method is called automatically by the Geant4 kernel (not by the user!)
  // at the end of each event : in MT-mode, it is called only for the Working thread
  // that handled that event.
  G4int nEvt = anEvent->GetEventID();
  if ( nEvt % 1000 == 0 ) {
    G4cout << std::setprecision(6) << " Event#=" << nEvt
           << " ; t[ns]=" << fDecayT << " ; r[mm]=" << fDecayR
           << " ; deltaR[mum]=" << ( fDecayR_mc_truth - fDecayR ) / CLHEP::micrometer
           << " ; deltaEkin[MeV]=" << fDeltaEkin
           << " ; deltaAngle(deg)=" << fDeflectionAngle << G4endl;
  }
  G4Run::RecordEvent( anEvent );  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::Merge( const G4Run* aRun ) {
  // This method is called automatically by the Geant4 kernel (not by the user!)
  // only in the case of multithreaded mode and only for Working threads.
  const Run* localRun = static_cast< const Run* >( aRun );
  fNumEvents += localRun->GetNumberOfEvent();

  fPrimaryParticleId = localRun->GetPrimaryParticleId();
  fPrimaryParticleInitialKineticEnergy = localRun->GetPrimaryParticleInitialKineticEnergy();
  fPrimaryParticleInitialTotalEnergy = localRun->GetPrimaryParticleInitialTotalEnergy();
  fPrimaryParticleInitialMomentum = localRun->GetPrimaryParticleInitialMomentum();
  fPrimaryParticleInitialBeta = localRun->GetPrimaryParticleInitialBeta();
  fPrimaryParticleInitialGamma = localRun->GetPrimaryParticleInitialGamma();
  fPrimaryParticleInitial3Momentum = localRun->GetPrimaryParticleInitial3Momentum();
  fPrimaryParticleInitialPosition = localRun->GetPrimaryParticleInitialPosition();
  fToleranceEPviolations = localRun->GetToleranceEPviolations();
  fToleranceDeltaDecayRadius = localRun->GetToleranceDeltaDecayRadius();
  fIsPreassignedDecayEnabled = localRun->GetIsPreassignedDecayEnabled();
  fIsBoostToLabEnabled = localRun->GetIsBoostToLabEnabled();
  
  fNumDecays += localRun->GetNumberDecays();
  fNumBadDecays += localRun->GetNumberBadDecays();
  fNumUnexpectedDecays += localRun->GetNumberUnexpectedDecays();
  fNumEviolations += localRun->GetNumberEviolations();
  fNumPviolations += localRun->GetNumberPviolations();
  fNum_mc_truth_rPos_deltaMax_above += localRun->GetNumber_mc_truth_rPos_deltaMax_above();
  fNum_underestimated_mc_truth_rPos_delta_above +=
    localRun->GetNumberUnderestimated_mc_truth_rPos_delta_above();
  fNum_overestimated_mc_truth_rPos_delta_above  +=
    localRun->GetNumberOverestimated_mc_truth_rPos_delta_above();
  fNumLargeUnderestimates += localRun->GetNumberLargeUnderestimates();
  fNumLargeOverestimates  += localRun->GetNumberLargeOverestimates();
  
  fSumDecayT += localRun->GetSumDecayT();
  fMinDecayT = std::min( fMinDecayT, localRun->GetMinDecayT() );
  fMaxDecayT = std::max( fMaxDecayT, localRun->GetMaxDecayT() );
  fSumDecayR += localRun->GetSumDecayR();
  fMinDecayR = std::min( fMinDecayR, localRun->GetMinDecayR() );
  fMaxDecayR = std::max( fMaxDecayR, localRun->GetMaxDecayR() );
  fSumDecayX += localRun->GetSumDecayX();
  fMinDecayX = std::min( fMinDecayX, localRun->GetMinDecayX() );
  fMaxDecayX = std::max( fMaxDecayX, localRun->GetMaxDecayX() );
  fSumDecayY += localRun->GetSumDecayY();
  fMinDecayY = std::min( fMinDecayY, localRun->GetMinDecayY() );
  fMaxDecayY = std::max( fMaxDecayY, localRun->GetMaxDecayY() );
  fSumDecayZ += localRun->GetSumDecayZ();
  fMinDecayZ = std::min( fMinDecayZ, localRun->GetMinDecayZ() );
  fMaxDecayZ = std::max( fMaxDecayZ, localRun->GetMaxDecayZ() );
  fSumDeltaDecayR += localRun->GetSumDeltaDecayR();
  fMinDeltaDecayR = std::min( fMinDeltaDecayR, localRun->GetMinDeltaDecayR() );
  fMaxDeltaDecayR = std::max( fMaxDeltaDecayR, localRun->GetMaxDeltaDecayR() );
  fSumDeflectionAngle += localRun->GetSumDeflectionAngle();
  fMinDeflectionAngle = std::min( fMinDeflectionAngle, localRun->GetMinDeflectionAngle() );
  fMaxDeflectionAngle = std::max( fMaxDeflectionAngle, localRun->GetMaxDeflectionAngle() );
  fSumDeltaEkin += localRun->GetSumDeltaEkin();
  fMinDeltaEkin = std::min( fMinDeltaEkin, localRun->GetMinDeltaEkin() );
  fMaxDeltaEkin = std::max( fMaxDeltaEkin, localRun->GetMaxDeltaEkin() );
  fSumDecayEkin += localRun->GetSumDecayEkin();
  fMinDecayEkin = std::min( fMinDecayEkin, localRun->GetMinDecayEkin() );
  fMaxDecayEkin = std::max( fMaxDecayEkin, localRun->GetMaxDecayEkin() );
  fSumDecayPx += localRun->GetSumDecayPx();
  fMinDecayPx = std::min( fMinDecayPx, localRun->GetMinDecayPx() );
  fMaxDecayPx = std::max( fMaxDecayPx, localRun->GetMaxDecayPx() );
  fSumDecayPy += localRun->GetSumDecayPy();
  fMinDecayPy = std::min( fMinDecayPy, localRun->GetMinDecayPy() );
  fMaxDecayPy = std::max( fMaxDecayPy, localRun->GetMaxDecayPy() );
  fSumDecayPz += localRun->GetSumDecayPz();
  fMinDecayPz = std::min( fMinDecayPz, localRun->GetMinDecayPz() );
  fMaxDecayPz = std::max( fMaxDecayPz, localRun->GetMaxDecayPz() );
  fSumDecayEtotViolation += localRun->GetSumDecayEtotViolation();
  fMinDecayEtotViolation = std::min( fMinDecayEtotViolation,
                                     localRun->GetMinDecayEtotViolation() );
  fMaxDecayEtotViolation = std::max( fMaxDecayEtotViolation,
                                     localRun->GetMaxDecayEtotViolation() );
  fSumDecayPxViolation += localRun->GetSumDecayPxViolation();
  fMinDecayPxViolation = std::min( fMinDecayPxViolation, localRun->GetMinDecayPxViolation() );
  fMaxDecayPxViolation = std::max( fMaxDecayPxViolation, localRun->GetMaxDecayPxViolation() );
  fSumDecayPyViolation += localRun->GetSumDecayPyViolation();
  fMinDecayPyViolation = std::min( fMinDecayPyViolation, localRun->GetMinDecayPyViolation() );
  fMaxDecayPyViolation = std::max( fMaxDecayPyViolation, localRun->GetMaxDecayPyViolation() );
  fSumDecayPzViolation += localRun->GetSumDecayPzViolation();
  fMinDecayPzViolation = std::min( fMinDecayPzViolation, localRun->GetMinDecayPzViolation() );
  fMaxDecayPzViolation = std::max( fMaxDecayPzViolation, localRun->GetMaxDecayPzViolation() );

  fMaxEkin_deltaMax  = std::max( fMaxEkin_deltaMax, localRun->GetMaxEkin_deltaMax() );
  fMaxEtot_deltaMax  = std::max( fMaxEtot_deltaMax, localRun->GetMaxEtot_deltaMax() );
  fMaxP_deltaMax     = std::max( fMaxP_deltaMax, localRun->GetMaxP_deltaMax() );
  fMaxPdir_deltaMax  = std::max( fMaxPdir_deltaMax, localRun->GetMaxPdir_deltaMax() );
  fMaxMass_deltaMax1 = std::max( fMaxMass_deltaMax1, localRun->GetMaxMass_deltaMax1() );
  fMaxMass_deltaMax2 = std::max( fMaxMass_deltaMax2, localRun->GetMaxMass_deltaMax2() );
  fMaxMass_deltaMax3 = std::max( fMaxMass_deltaMax3, localRun->GetMaxMass_deltaMax3() );
  fSumMass_deltaMax3 += localRun->GetSumMass_deltaMax3();
  fMaxBeta_deltaMax1 = std::max( fMaxBeta_deltaMax1, localRun->GetMaxBeta_deltaMax1() );
  fMaxBeta_deltaMax2 = std::max( fMaxBeta_deltaMax2, localRun->GetMaxBeta_deltaMax2() );
  fMaxGamma_deltaMax1 = std::max( fMaxGamma_deltaMax1, localRun->GetMaxGamma_deltaMax1() );
  fMaxGamma_deltaMax2 = std::max( fMaxGamma_deltaMax2, localRun->GetMaxGamma_deltaMax2() );
  fMaxGamma_deltaMax3 = std::max( fMaxGamma_deltaMax3, localRun->GetMaxGamma_deltaMax3() );
  fMaxT_proper_deltaMax = std::max( fMaxT_proper_deltaMax, localRun->GetMaxT_proper_deltaMax() );
  fMaxT_lab_deltaMax    = std::max( fMaxT_lab_deltaMax, localRun->GetMaxT_lab_deltaMax() );
  fMaxMc_truth_rPos_deltaMax = std::max( fMaxMc_truth_rPos_deltaMax,
                                         localRun->GetMaxMc_truth_rPos_deltaMax() );
  fSumMc_truth_rPos_deltaMax += localRun->GetSumMc_truth_rPos_deltaMax();

  fMinUnderestimated_mc_truth_rPos_delta =
    std::min( fMinUnderestimated_mc_truth_rPos_delta,
              localRun->GetMinUnderestimated_mc_truth_rPos_delta() );
  fSumUnderestimated_mc_truth_rPos_delta += localRun->GetSumUnderestimated_mc_truth_rPos_delta();
  fMaxOverestimated_mc_truth_rPos_delta =
    std::max( fMaxOverestimated_mc_truth_rPos_delta,
              localRun->GetMaxOverestimated_mc_truth_rPos_delta() );
  fSumOverestimated_mc_truth_rPos_delta += localRun->GetSumOverestimated_mc_truth_rPos_delta();
  fMinUnderestimated_rDeltaPos = std::min( fMinUnderestimated_rDeltaPos,
                                           localRun->GetMinUnderestimated_rDeltaPos() );
  fSumUnderestimated_rDeltaPos += localRun->GetSumUnderestimated_rDeltaPos();
  fMaxOverestimated_rDeltaPos = std::max( fMaxOverestimated_rDeltaPos,
                                          localRun->GetMaxOverestimated_rDeltaPos() );
  fSumOverestimated_rDeltaPos += localRun->GetSumOverestimated_rDeltaPos();

  fMaxFloat_rDeltaPos_deltaMax = std::max( fMaxFloat_rDeltaPos_deltaMax,
                                           localRun->GetMaxFloat_rDeltaPos_deltaMax() );
 
  G4Run::Merge( aRun );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::PrintInfo() const {
  // This method is called by RunAction::EndOfRunAction.
  // In MT-mode, only the master thread calls it.
  const G4double femtosecond = 0.001 * CLHEP::picosecond;
  const G4double NN = ( fNumDecays > 0 ? G4double( fNumDecays ) : 1.0 );
  G4cout << std::setprecision(6)
         << G4endl << G4endl
         << " ============ Run::printInfo() =============== \t RunID = " << GetRunID() << G4endl
         << " primary PDG code                     = " << fPrimaryParticleId << G4endl
         << " primary initial kinetic energy [GeV] = "
         << fPrimaryParticleInitialKineticEnergy / CLHEP::GeV << G4endl
         << " primary initial total energy [GeV]   = "
         << fPrimaryParticleInitialTotalEnergy / CLHEP::GeV << G4endl
         << " primary initial momentum [GeV]       = "
         << fPrimaryParticleInitialMomentum / CLHEP::GeV << G4endl
         << " primary initial Lorentz beta         = " << fPrimaryParticleInitialBeta << G4endl
         << " primary initial Lorentz gamma        = " << fPrimaryParticleInitialGamma << G4endl
         << " primary initial 3 momentum [GeV]     = "
         << fPrimaryParticleInitial3Momentum / CLHEP::GeV << G4endl
         << " primary initial position [mm]        = "
         << fPrimaryParticleInitialPosition << G4endl
         << " toleranceEPviolations [eV]           = "
         << fToleranceEPviolations / CLHEP::eV << G4endl
         << " toleranceDeltaDecayRadius [mum]      = "
         << fToleranceDeltaDecayRadius / CLHEP::micrometer << G4endl
         << " isPreassignedDecayEnabled            = " << fIsPreassignedDecayEnabled << G4endl
         << " isBoostToLabEnabled                  = " << fIsBoostToLabEnabled << G4endl
         << " # Events               = "
         << ( fNumEvents > 0 ? fNumEvents : GetNumberOfEvent() ) << G4endl
         << " # Decays               = " << fNumDecays << G4endl
         << "   (# Bad decays        = " << fNumBadDecays << " )" << G4endl
         << "   (# Unexpected decays = " << fNumUnexpectedDecays << " )" << G4endl
         << " # E violations         = " << fNumEviolations << G4endl
         << " # P violations         = " << fNumPviolations << G4endl
         << " decay T [ns] : min=" << fMinDecayT << "\t mean=" << fSumDecayT / NN
         << "\t max=" << fMaxDecayT << G4endl
         << " decay R [mm] : min=" << fMinDecayR << "\t mean=" << fSumDecayR / NN
         << "\t max=" << fMaxDecayR << G4endl
         << " decay X [mm] : min=" << fMinDecayX << "\t mean=" << fSumDecayX / NN
         << "\t max=" << fMaxDecayX << G4endl
         << " decay Y [mm] : min=" << fMinDecayY << "\t mean=" << fSumDecayY / NN
         << "\t max=" << fMaxDecayY << G4endl
         << " decay Z [mm] : min=" << fMinDecayZ << "\t mean=" << fSumDecayZ / NN
         << "\t max=" << fMaxDecayZ << G4endl
         << " Delta decay R [mm] : min=" << fMinDeltaDecayR
         << "\t mean=" << fSumDeltaDecayR / NN
         << "\t max=" << fMaxDeltaDecayR << G4endl
         << " deflection angle [deg] : min=" << fMinDeflectionAngle
         << "\t mean=" << fSumDeflectionAngle / NN
         << "\t max=" << fMaxDeflectionAngle << G4endl
         << " Delta Ekin [MeV] : min=" << fMinDeltaEkin
         << "\t mean=" << fSumDeltaEkin / NN
         << "\t max=" << fMaxDeltaEkin << G4endl
         << " decay Ekin [GeV] : min=" << fMinDecayEkin / CLHEP::GeV
         << "\t mean=" << fSumDecayEkin / NN / CLHEP::GeV
         << "\t max=" << fMaxDecayEkin / CLHEP::GeV << G4endl
         << " decay Px [GeV]   : min=" << fMinDecayPx / CLHEP::GeV
         << "\t mean=" << fSumDecayPx / NN / CLHEP::GeV
         << "\t max=" << fMaxDecayPx / CLHEP::GeV << G4endl
         << " decay Py [GeV]   : min=" << fMinDecayPy / CLHEP::GeV
         << "\t mean=" << fSumDecayPy / NN / CLHEP::GeV
         << "\t max=" << fMaxDecayPy / CLHEP::GeV << G4endl
         << " decay Pz [GeV]   : min=" << fMinDecayPz / CLHEP::GeV
         << "\t mean=" << fSumDecayPz / NN / CLHEP::GeV
         << "\t max=" << fMaxDecayPz / CLHEP::GeV << G4endl
         << " decay Etot violation [MeV] : min=" << fMinDecayEtotViolation
         << "\t mean=" << fSumDecayEtotViolation / NN
         << "\t max=" << fMaxDecayEtotViolation << G4endl
         << " decay Px violation [MeV]   : min=" << fMinDecayPxViolation
         << "\t mean=" << fSumDecayPxViolation / NN
         << "\t max=" << fMaxDecayPxViolation << G4endl
         << " decay Py violation [MeV]   : min=" << fMinDecayPyViolation
         << "\t mean=" << fSumDecayPyViolation / NN
         << "\t max=" << fMaxDecayPyViolation << G4endl
         << " decay Pz violation [MeV]   : min=" << fMinDecayPzViolation
         << "\t mean=" << fSumDecayPzViolation / NN
         << "\t max=" << fMaxDecayPzViolation << G4endl
         << " --- Consistency checks --- " << G4endl
         << " maxEkin_deltaMax [eV]           = " << fMaxEkin_deltaMax / CLHEP::eV << G4endl
         << " maxEtot_deltaMax [eV]           = " << fMaxEtot_deltaMax / CLHEP::eV << G4endl
         << " maxP_deltaMax [eV]              = " << fMaxP_deltaMax    / CLHEP::eV << G4endl
         << " maxPdir_deltaMax                = " << fMaxPdir_deltaMax << G4endl
         << " maxMass_deltaMax{1,2,3} [eV]    = " << fMaxMass_deltaMax1 / CLHEP::eV
         << " , " << fMaxMass_deltaMax2 / CLHEP::eV << " , " << fMaxMass_deltaMax3 / CLHEP::eV
         << " (mean=" << fSumMass_deltaMax3 / NN / CLHEP::eV << ")" << G4endl
         << " maxBeta_deltaMax{1,2}           = " << fMaxBeta_deltaMax1
         << " , " << fMaxBeta_deltaMax2 << G4endl
         << " maxGamma_deltaMax{1,2,3}        = " << fMaxGamma_deltaMax1
         << " , " << fMaxGamma_deltaMax2
         << " , " << fMaxGamma_deltaMax3 << G4endl
         << " maxT_proper_deltaMax [fs]       = " << fMaxT_proper_deltaMax / femtosecond << G4endl
         << " maxT_lab_deltaMax [fs]          = " << fMaxT_lab_deltaMax   / femtosecond << G4endl
         << " maxMc_truth_rPos_deltaMax [mum] = "
         << fMaxMc_truth_rPos_deltaMax / CLHEP::micrometer
         << " (mean=" << fSumMc_truth_rPos_deltaMax / NN / CLHEP::micrometer
         << ")\t (# above threshold = " << fNum_mc_truth_rPos_deltaMax_above << ")" << G4endl
         << " --- Extra checks --- " << G4endl
         << " minUnderestimated_mc_truth_rPos_delta [mum] = "
         << fMinUnderestimated_mc_truth_rPos_delta / CLHEP::micrometer
         << " (mean=" << fSumUnderestimated_mc_truth_rPos_delta / NN / CLHEP::micrometer
         << ")\t (#above threshold = " << fNum_underestimated_mc_truth_rPos_delta_above
         << ")" << G4endl
         << " maxOverestimated_mc_truth_rPos_delta [mum]  = "
         << fMaxOverestimated_mc_truth_rPos_delta / CLHEP::micrometer
         << " (mean=" << fSumOverestimated_mc_truth_rPos_delta / NN / CLHEP::micrometer
         << ")\t (#above threshold = " << fNum_overestimated_mc_truth_rPos_delta_above
         << ")" << G4endl
         << " minUnderestimated_rDeltaPos [mum]           = "
         << fMinUnderestimated_rDeltaPos / CLHEP::micrometer
         << " (mean=" << fSumUnderestimated_rDeltaPos / NN / CLHEP::micrometer
         << ")\t (#above threshold = " << fNumLargeUnderestimates << ")" << G4endl
         << " maxOverestimated_rDeltaPos [mum]            = "
         << fMaxOverestimated_rDeltaPos / CLHEP::micrometer
         << " (mean=" << fSumOverestimated_rDeltaPos / NN / CLHEP::micrometer
         << ")\t (#above threshold = " << fNumLargeOverestimates << ")" << G4endl
         << " --- float  instead of  double --- " << G4endl
         << " fMaxFloat_rDeltaPos_deltaMax [mum] = "
         << fMaxFloat_rDeltaPos_deltaMax / CLHEP::micrometer << G4endl
         << " ============================================= " << G4endl
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetDecayT( const G4double inputValue ) {
  fDecayT = inputValue;
  fSumDecayT += inputValue;
  fMinDecayT = std::min( fMinDecayT, inputValue );
  fMaxDecayT = std::max( fMaxDecayT, inputValue );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetDecayR_mc_truth( const G4double inputValue ) {
  fDecayR_mc_truth = inputValue;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetDecayR( const G4double inputValue ) {
  fDecayR = inputValue;
  fSumDecayR += inputValue;
  fMinDecayR = std::min( fMinDecayR, inputValue );
  fMaxDecayR = std::max( fMaxDecayR, inputValue );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetDecayX( const G4double inputValue ) {
  fDecayX = inputValue;
  fSumDecayX += inputValue;
  fMinDecayX = std::min( fMinDecayX, inputValue );
  fMaxDecayX = std::max( fMaxDecayX, inputValue );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetDecayY( const G4double inputValue ) {
  fDecayY = inputValue;
  fSumDecayY += inputValue;
  fMinDecayY = std::min( fMinDecayY, inputValue );
  fMaxDecayY = std::max( fMaxDecayY, inputValue );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetDecayZ( const G4double inputValue ) {
  fDecayZ = inputValue;
  fSumDecayZ += inputValue;
  fMinDecayZ = std::min( fMinDecayZ, inputValue );
  fMaxDecayZ = std::max( fMaxDecayZ, inputValue );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetDeltaDecayR( const G4double inputValue ) {
  fDeltaDecayR = inputValue;
  fSumDeltaDecayR += inputValue;
  fMinDeltaDecayR = std::min( fMinDeltaDecayR, inputValue );
  fMaxDeltaDecayR = std::max( fMaxDeltaDecayR, inputValue );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetDeflectionAngle( const G4double inputValue ) {
  fDeflectionAngle = inputValue;
  fSumDeflectionAngle += inputValue;
  fMinDeflectionAngle = std::min( fMinDeflectionAngle, inputValue );
  fMaxDeflectionAngle = std::max( fMaxDeflectionAngle, inputValue );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetDeltaEkin( const G4double inputValue ) {
  fDeltaEkin = inputValue;
  fSumDeltaEkin += inputValue;
  fMinDeltaEkin = std::min( fMinDeltaEkin, inputValue );
  fMaxDeltaEkin = std::max( fMaxDeltaEkin, inputValue );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetDecayEkin( const G4double inputValue ) {
  fDecayEkin = inputValue;
  fSumDecayEkin += inputValue;
  fMinDecayEkin = std::min( fMinDecayEkin, inputValue );
  fMaxDecayEkin = std::max( fMaxDecayEkin, inputValue );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetDecayPx( const G4double inputValue ) {
  fDecayPx = inputValue;
  fSumDecayPx += inputValue;
  fMinDecayPx = std::min( fMinDecayPx, inputValue );
  fMaxDecayPx = std::max( fMaxDecayPx, inputValue );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetDecayPy( const G4double inputValue ) {
  fDecayPy = inputValue;
  fSumDecayPy += inputValue;
  fMinDecayPy = std::min( fMinDecayPy, inputValue );
  fMaxDecayPy = std::max( fMaxDecayPy, inputValue );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetDecayPz( const G4double inputValue ) {
  fDecayPz = inputValue;
  fSumDecayPz += inputValue;
  fMinDecayPz = std::min( fMinDecayPz, inputValue );
  fMaxDecayPz = std::max( fMaxDecayPz, inputValue );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetDecayEtotViolation( const G4double inputValue ) {
  fDecayEtotViolation = inputValue;
  fSumDecayEtotViolation += inputValue;
  fMinDecayEtotViolation = std::min( fMinDecayEtotViolation, inputValue );
  fMaxDecayEtotViolation = std::max( fMaxDecayEtotViolation, inputValue );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetDecayPxViolation( const G4double inputValue ) {
  fDecayPxViolation = inputValue;
  fSumDecayPxViolation += inputValue;
  fMinDecayPxViolation = std::min( fMinDecayPxViolation, inputValue );
  fMaxDecayPxViolation = std::max( fMaxDecayPxViolation, inputValue );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetDecayPyViolation( const G4double inputValue ) {
  fDecayPyViolation = inputValue;
  fSumDecayPyViolation += inputValue;
  fMinDecayPyViolation = std::min( fMinDecayPyViolation, inputValue );
  fMaxDecayPyViolation = std::max( fMaxDecayPyViolation, inputValue );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetDecayPzViolation( const G4double inputValue ) {
  fDecayPzViolation = inputValue;
  fSumDecayPzViolation += inputValue;
  fMinDecayPzViolation = std::min( fMinDecayPzViolation, inputValue );
  fMaxDecayPzViolation = std::max( fMaxDecayPzViolation, inputValue );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetMaxEkin_deltaMax( const G4double inputValue ) {
  fMaxEkin_deltaMax = std::max( fMaxEkin_deltaMax, inputValue );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetMaxEtot_deltaMax( const G4double inputValue ) {
  fMaxEtot_deltaMax = std::max( fMaxEtot_deltaMax, inputValue );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetMaxP_deltaMax( const G4double inputValue ) {
  fMaxP_deltaMax = std::max( fMaxP_deltaMax, inputValue );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetMaxPdir_deltaMax( const G4double inputValue ) {
  fMaxPdir_deltaMax = std::max( fMaxPdir_deltaMax, inputValue );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetMaxMass_deltaMax1( const G4double inputValue ) {
  fMaxMass_deltaMax1 = std::max( fMaxMass_deltaMax1, inputValue );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetMaxMass_deltaMax2( const G4double inputValue ) {
  fMaxMass_deltaMax2 = std::max( fMaxMass_deltaMax2, inputValue );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetMaxMass_deltaMax3( const G4double inputValue ) {
  fSumMass_deltaMax3 += inputValue;
  fMaxMass_deltaMax3 = std::max( fMaxMass_deltaMax3, inputValue );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetMaxBeta_deltaMax1( const G4double inputValue ) {
  fMaxBeta_deltaMax1 = std::max( fMaxBeta_deltaMax1, inputValue );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetMaxBeta_deltaMax2( const G4double inputValue ) {
  fMaxBeta_deltaMax2 = std::max( fMaxBeta_deltaMax2, inputValue );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetMaxGamma_deltaMax1( const G4double inputValue ) {
  fMaxGamma_deltaMax1 = std::max( fMaxGamma_deltaMax1, inputValue );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetMaxGamma_deltaMax2( const G4double inputValue ) {
  fMaxGamma_deltaMax2 = std::max( fMaxGamma_deltaMax2, inputValue );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetMaxGamma_deltaMax3( const G4double inputValue ) {
  fMaxGamma_deltaMax3 = std::max( fMaxGamma_deltaMax3, inputValue );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetMaxT_proper_deltaMax( const G4double inputValue ) {
  fMaxT_proper_deltaMax = std::max( fMaxT_proper_deltaMax, inputValue );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetMaxT_lab_deltaMax( const G4double inputValue ) {
  fMaxT_lab_deltaMax = std::max( fMaxT_lab_deltaMax, inputValue );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetMaxMc_truth_rPos_deltaMax( const G4double inputValue ) {
  fSumMc_truth_rPos_deltaMax += inputValue;
  fMaxMc_truth_rPos_deltaMax = std::max( fMaxMc_truth_rPos_deltaMax, inputValue );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetMinUnderestimated_mc_truth_rPos_delta( const G4double inputValue ) {
  fSumUnderestimated_mc_truth_rPos_delta += inputValue;
  fMinUnderestimated_mc_truth_rPos_delta = std::min( fMinUnderestimated_mc_truth_rPos_delta,
                                                     inputValue );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetMaxOverestimated_mc_truth_rPos_delta( const G4double inputValue ) {
  fSumOverestimated_mc_truth_rPos_delta += inputValue;
  fMaxOverestimated_mc_truth_rPos_delta = std::max( fMaxOverestimated_mc_truth_rPos_delta,
                                                    inputValue );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetMinUnderestimated_rDeltaPos( const G4double inputValue ) {
  fSumUnderestimated_rDeltaPos += inputValue;
  fMinUnderestimated_rDeltaPos = std::min( fMinUnderestimated_rDeltaPos, inputValue );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetMaxOverestimated_rDeltaPos( const G4double inputValue ) {
  fSumOverestimated_rDeltaPos += inputValue;
  fMaxOverestimated_rDeltaPos = std::max( fMaxOverestimated_rDeltaPos, inputValue );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetMaxFloat_rDeltaPos_deltaMax( const G4double inputValue ) {
  fMaxFloat_rDeltaPos_deltaMax = std::max( fMaxFloat_rDeltaPos_deltaMax, inputValue );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
