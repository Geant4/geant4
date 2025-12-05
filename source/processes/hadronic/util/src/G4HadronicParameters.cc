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
//
//---------------------------------------------------------------------------
//
// ClassName:      G4HadronicParameters
//
// Author:         2018 Alberto Ribon
//
// Description:    Singleton to keep global hadronic parameters.
//
// Modified:
//
//----------------------------------------------------------------------------

#include "G4HadronicParameters.hh"
#include <CLHEP/Units/PhysicalConstants.h>
#include "G4ApplicationState.hh"
#include "G4StateManager.hh"
#include "G4HadronicParametersMessenger.hh"
#include "G4Threading.hh"
#include "G4AutoLock.hh"
#include "G4UnitsTable.hh"

G4HadronicParameters* G4HadronicParameters::sInstance = nullptr;

namespace {
  G4Mutex paramMutex = G4MUTEX_INITIALIZER;
}


G4HadronicParameters* G4HadronicParameters::Instance() {
  if ( sInstance == nullptr ) {
    G4AutoLock l(&paramMutex);
    if ( sInstance == nullptr ) {
      static G4HadronicParameters theHadronicParametersObject;
      sInstance = &theHadronicParametersObject;
    }
    l.unlock();
  }
  return sInstance;
}


G4HadronicParameters::~G4HadronicParameters() {
  delete fMessenger;
}


G4HadronicParameters::G4HadronicParameters() {
  fMaxEnergy = 100.0*CLHEP::TeV;
  fMinEnergyTransitionFTF_Cascade = 3.0*CLHEP::GeV;
  fMaxEnergyTransitionFTF_Cascade = 6.0*CLHEP::GeV;
  fMinEnergyTransitionQGS_FTF = 12.0*CLHEP::GeV;
  fMaxEnergyTransitionQGS_FTF = 25.0*CLHEP::GeV;
  fMinEnergyINCLXX_Pbar = 0.0*CLHEP::GeV;
  fMaxEnergyINCLXX_Pbar = 10.0*CLHEP::GeV;
  fEnergyThresholdForHeavyHadrons = 1.1*CLHEP::GeV;
  fMessenger = new G4HadronicParametersMessenger( this );

  // read environment variables
  fReportLevel = G4GetEnv<G4int>("G4Hadronic_epReportLevel", 0);
  const char* ep1 = std::getenv("G4Hadronic_epCheckRelativeLevel");
  if(nullptr != ep1) { fRelativeDiff = std::strtod(ep1, 0); }
  const char* ep2 = std::getenv("G4Hadronic_epCheckAbsoluteLevel");
  if(nullptr != ep2) { fAbsoluteDiff = std::strtod(ep2, 0); }
  const char* v = G4FindDataDir("G4PARTICLEXSDATA");
  if(nullptr != v) {
    fDirPARTICLEXS = G4String(v);
  } else {
    if(1 < fVerboseLevel) {
      G4ExceptionDescription ed;
      ed << "Environment variable G4PARTICLEXSDATA is not defined or " 
         << " it is pointing out to not existing directory";
      G4Exception("G4LevelReader::LevelManager(..)","had014",
		  JustWarning, ed, "Check file path");
    }
  }
  const char* x = std::getenv("G4PhysListDocDir");
  if(nullptr != x) { fPhysListDocDir = G4String(x); }
  const char* y = std::getenv("G4PhysListName");
  if(nullptr != y) { fPhysListName = G4String(y); }
  const char* z = std::getenv("BINARY_CASCADE_DEBUG");
  if(nullptr != z) { fBinaryDebug = true; }
}


G4bool G4HadronicParameters::IsLocked() const {
  return ( ! G4Threading::IsMasterThread() ||
           G4StateManager::GetStateManager()->GetCurrentState() != G4State_PreInit );
}


void G4HadronicParameters::StreamInfo( std::ostream& os ) const {
  G4long prec = os.precision(5);

  // Lambda function to convert boolean to "true"/"false" string
  auto boolToString = [](G4bool value) -> const char* {
    return value ? "true" : "false";
  };

  os << "=======================================================================" << "\n";
  os << "======                 Hadronic Physics Parameters             ========" << "\n";
  os << "=======================================================================" << "\n";
  os << "Maximum energy for hadronic physics                 "
     << G4BestUnit(fMaxEnergy, "Energy") << "\n";
  os << "Energy threshold for heavy hadrons                  "
     << G4BestUnit(fEnergyThresholdForHeavyHadrons, "Energy") << "\n";
  os << "Neutron kinetic energy threshold for SVT algorithm  ";
  if (fNeutronEkinThresholdForSVT < 0.0) {
    os <<"not set" << "\n";
  } else {
    os << G4BestUnit(fNeutronEkinThresholdForSVT, "Energy") << "\n";
  }
  os << "Time threshold for radioactive decays               ";
  if (fTimeThresholdForRadioactiveDecays < 0.0) {
    os << "not set" << "\n";
  } else {
    os << G4BestUnit(fTimeThresholdForRadioactiveDecays, "Time") << "\n";
  }

  os << "=======================================================================" << "\n";
  os << "======                 Model Transition Regions                ========" << "\n";
  os << "=======================================================================" << "\n";
  os << "FTF to Cascade transition region                    "
     << G4BestUnit(fMinEnergyTransitionFTF_Cascade, "Energy") << " - "
     << G4BestUnit(fMaxEnergyTransitionFTF_Cascade, "Energy") << "\n";
  os << "QGS to FTF transition region                        "
     << G4BestUnit(fMinEnergyTransitionQGS_FTF, "Energy") << " - "
     << G4BestUnit(fMaxEnergyTransitionQGS_FTF, "Energy") << "\n";
  os << "INCLXX antiproton model energy range                "
     << G4BestUnit(fMinEnergyINCLXX_Pbar, "Energy") << " - "
     << G4BestUnit(fMaxEnergyINCLXX_Pbar, "Energy") << "\n";

  os << "=======================================================================" << "\n";
  os << "======                 Cross Section Factors                   ========" << "\n";
  os << "=======================================================================" << "\n";
  os << "Apply cross section factors                         " << boolToString(fApplyFactorXS) << "\n";
  os << "Nucleon inelastic cross section factor              " << fXSFactorNucleonInelastic << "\n";
  os << "Nucleon elastic cross section factor                " << fXSFactorNucleonElastic << "\n";
  os << "Pion inelastic cross section factor                 " << fXSFactorPionInelastic << "\n";
  os << "Pion elastic cross section factor                   " << fXSFactorPionElastic << "\n";
  os << "Hadron inelastic cross section factor               " << fXSFactorHadronInelastic << "\n";
  os << "Hadron elastic cross section factor                 " << fXSFactorHadronElastic << "\n";
  os << "EM cross section factor                             " << fXSFactorEM << "\n";
  os << "For XS neutron cross sections use detailed R-data   " << fUseRFilesForXS << "\n";

  os << "=======================================================================" << "\n";
  os << "======                 Process Control Parameters              ========" << "\n";
  os << "=======================================================================" << "\n";
  os << "Enable integral method for inelastic cross sections " << boolToString(fEnableIntegralInelasticXS) << "\n";
  os << "Enable integral method for elastic cross sections   " << boolToString(fEnableIntegralElasticXS) << "\n";
  os << "Enable diffraction dissociation for B > 10          " << boolToString(fEnableDiffDissociationForBGreater10) << "\n";
  os << "Enable neutron general process                      " << boolToString(fNeutronGeneral) << "\n";
  os << "Enable NUDEX gamma de-excitation                    " << boolToString(fEnableNUDEX) << "\n";
  os << "Enable coherent charge exchange                     " << boolToString(fChargeExchange) << "\n";

  os << "=======================================================================" << "\n";
  os << "======                 Particle Production Control             ========" << "\n";
  os << "=======================================================================" << "\n";
  os << "Enable B/C particles                                " << boolToString(fEnableBC) << "\n";
  os << "Enable hyper-nuclei                                 " << boolToString(fEnableHyperNuclei) << "\n";
  os << "Enable cosmic ray coalescence                       " << boolToString(fEnableCRCoalescence) << "\n";

  os << "=======================================================================" << "\n";
  os << "======                 Model Control Parameters                ========" << "\n";
  os << "=======================================================================" << "\n";
  os << "PT table type for URR neutrons                      ";
  // A ternary operation can't be used here as it leads to a C2445 error on Windows with C++20 and newer
  if (fTypeTablePT.empty()) {
    os << "not set" << "\n";
  } else {
    os << fTypeTablePT << "\n";
  }
  os << "Bertini angular emissions as in G4 11.2             " << boolToString(fBertiniAngularEmissionsAs11_2) << "\n";
  os << "Bertini nuclei model as in G4 11.2                  " << boolToString(fBertiniNucleiModelAs11_2) << "\n";
  os << "Bertini overall behavior as in G4 11.2              " << boolToString(IsBertiniAs11_2()) << "\n";

  os << "=======================================================================" << "\n";
  os << "======                 Debugging Options                       ========" << "\n";
  os << "=======================================================================" << "\n";
  os << "Verbose level                                       " << fVerboseLevel << "\n";
  os << "Binary cascade debug                                " << boolToString(fBinaryDebug) << "\n";
  os << "Environment reporting level                         " << fReportLevel << "\n";
  if (fRelativeDiff < DBL_MAX) {
    os << "Environment relative difference level                " << fRelativeDiff << "\n";
  }
  if (fAbsoluteDiff < DBL_MAX) {
    os << "Environment absolute difference level                " << fAbsoluteDiff << "\n";
  }

  os << "=======================================================================" << G4endl;
  os.precision(prec);
}


void G4HadronicParameters::Dump() const {
  StreamInfo(G4cout);
}


void G4HadronicParameters::SetMaxEnergy( const G4double val ) {
  if ( ! IsLocked()  &&  val > 0.0 ) { 
    fMaxEnergy = val;
  }
}


void G4HadronicParameters::SetMinEnergyTransitionFTF_Cascade( const G4double val ) {
  if ( ! IsLocked()  &&  val > 0.0 ) { 
    fMinEnergyTransitionFTF_Cascade = val;
  }
}


void G4HadronicParameters::SetMaxEnergyTransitionFTF_Cascade( const G4double val ) {
  if ( ! IsLocked()  &&  val > fMinEnergyTransitionFTF_Cascade ) { 
    fMaxEnergyTransitionFTF_Cascade = val;
  }
}


void G4HadronicParameters::SetMinEnergyTransitionQGS_FTF( const G4double val ) {
  if ( ! IsLocked()  &&  val > 0.0 ) { 
    fMinEnergyTransitionQGS_FTF = val;
  }
}


void G4HadronicParameters::SetMaxEnergyTransitionQGS_FTF( const G4double val ) {
  if ( ! IsLocked()  &&  val > fMinEnergyTransitionQGS_FTF ) { 
    fMaxEnergyTransitionQGS_FTF = val;
  }
}


void G4HadronicParameters::SetMinEnergyINCLXX_Pbar( const G4double val ) {
  if ( ! IsLocked()  &&  val >= 0.0 ) { 
    fMinEnergyINCLXX_Pbar = val;
  }
}


void G4HadronicParameters::SetMaxEnergyINCLXX_Pbar( const G4double val ) {
  if ( ! IsLocked()  &&  val > fMinEnergyINCLXX_Pbar ) { 
    fMaxEnergyINCLXX_Pbar = val;
  }
}


void G4HadronicParameters::SetEnableBCParticles( G4bool val ) {
  if ( ! IsLocked() ) fEnableBC = val;
}


void G4HadronicParameters::SetEnableHyperNuclei( G4bool val ) {
  if ( ! IsLocked() ) fEnableHyperNuclei = val;
}


void G4HadronicParameters::SetVerboseLevel( const G4int val ) {
  if ( ! IsLocked()  &&  val >= 0 ) fVerboseLevel = val;
}


void G4HadronicParameters::SetEnergyThresholdForHeavyHadrons( G4double val ) {
  if ( ! IsLocked()  &&  val >= 0 && val < 5*CLHEP::GeV ) {
    fEnergyThresholdForHeavyHadrons = val;
  }
}


void G4HadronicParameters::SetXSFactorNucleonInelastic( G4double val ) {
  if ( ! IsLocked()  &&  std::abs(val - 1.0) < fXSFactorLimit ) {
    fXSFactorNucleonInelastic = val;
  }
}


void G4HadronicParameters::SetXSFactorNucleonElastic( G4double val ) {
  if ( ! IsLocked()  &&  std::abs(val - 1.0) < fXSFactorLimit ) {
    fXSFactorNucleonElastic = val;
  }
}


void G4HadronicParameters::SetXSFactorPionInelastic( G4double val ) {
  if ( ! IsLocked()  &&  std::abs(val - 1.0) < fXSFactorLimit ) {
    fXSFactorPionInelastic = val;
  }
}


void G4HadronicParameters::SetXSFactorPionElastic( G4double val ) {
  if ( ! IsLocked()  &&  std::abs(val - 1.0) < fXSFactorLimit ) {
    fXSFactorPionElastic = val;
  }
}


void G4HadronicParameters::SetXSFactorHadronInelastic( G4double val ) {
  if ( ! IsLocked()  &&  std::abs(val - 1.0) < fXSFactorLimit ) {
    fXSFactorHadronInelastic = val;
  }
}


void G4HadronicParameters::SetXSFactorHadronElastic( G4double val ) {
  if ( ! IsLocked()  &&  std::abs(val - 1.0) < fXSFactorLimit ) {
    fXSFactorHadronElastic = val;
  }
}


void G4HadronicParameters::SetXSFactorEM( G4double val ) {
  if ( ! IsLocked()  &&  std::abs(val - 1.0) < fXSFactorLimit ) {
    fXSFactorEM = val;
  }
}


void G4HadronicParameters::SetNeutronKineticEnergyThresholdForSVT( const G4double val ) {
  // This setting works only after initialization (i.e. for G4State_Idle, 
  // whereas it does not work for G4State_PreInit).
  if ( G4Threading::IsMasterThread()  &&  val > 0.0 ) { 
    fNeutronEkinThresholdForSVT = val;
  }
}


void G4HadronicParameters::SetTimeThresholdForRadioactiveDecay( const G4double val ) {
  // This setting works only before initialization 
  // (else, if used after initialization, it will be ignored).
  if ( G4Threading::IsMasterThread()  &&  val > 0.0 ) { 
    fTimeThresholdForRadioactiveDecays = val;
  }
}


void G4HadronicParameters::SetApplyFactorXS( G4bool val ) {
  if ( ! IsLocked() ) fApplyFactorXS = val; 
}


void G4HadronicParameters::SetEnableCRCoalescence( G4bool val ) {
  if ( ! IsLocked() ) fEnableCRCoalescence = val;
}


void G4HadronicParameters::SetEnableIntegralInelasticXS( G4bool val ) {
  if ( ! IsLocked() ) fEnableIntegralInelasticXS = val;
}


void G4HadronicParameters::SetEnableIntegralElasticXS( G4bool val ) {
  if ( ! IsLocked() ) fEnableIntegralElasticXS = val;
}


void G4HadronicParameters::SetEnableDiffDissociationForBGreater10( G4bool val ) {
  if ( ! IsLocked() ) fEnableDiffDissociationForBGreater10 = val;
}


void G4HadronicParameters::SetEnableNeutronGeneralProcess( G4bool val ) {
  if ( ! IsLocked() ) fNeutronGeneral = val;
} 


void G4HadronicParameters::SetEnableNUDEX( G4bool val ) {
  if ( ! IsLocked() ) fEnableNUDEX = val;
} 


void G4HadronicParameters::SetTypeTablePT( const G4String& typeTablePT ) {
  if ( ! IsLocked() ) fTypeTablePT = typeTablePT;
}


void G4HadronicParameters::SetEnableCoherentChargeExchange( G4bool val ) {
  if ( ! IsLocked() )  fChargeExchange = val;
}


void G4HadronicParameters::SetBertiniAs11_2( G4bool val ) {
  if ( ! IsLocked() ) {
    fBertiniAngularEmissionsAs11_2 = val;
    fBertiniNucleiModelAs11_2 = val;
  }
}


void G4HadronicParameters::SetBertiniAngularEmissionsAs11_2( G4bool val ) {
  if ( ! IsLocked() ) fBertiniAngularEmissionsAs11_2 = val;
}


void G4HadronicParameters::SetBertiniNucleiModelAs11_2( G4bool val ) {
  if ( ! IsLocked() ) fBertiniNucleiModelAs11_2 = val;
}


void G4HadronicParameters::SetUseRFilesForXS( G4bool val ) {
  if ( ! IsLocked() ) fUseRFilesForXS = val;
}
