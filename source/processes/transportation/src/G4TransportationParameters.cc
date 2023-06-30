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
// Author:  J. Apostolakis Nov 2022 

#include "G4TransportationParameters.hh"

#include "G4ApplicationState.hh"
#include "G4StateManager.hh"
#include "G4Threading.hh"
#include "G4AutoLock.hh"

#include "G4SystemOfUnits.hh"

G4TransportationParameters* G4TransportationParameters::theInstance = nullptr;

#ifdef G4MULTITHREADED   
namespace
{
   G4Mutex transportParamsMutex = G4MUTEX_INITIALIZER;
}
#endif

//------------------------------------------------------------------------------

G4TransportationParameters::G4TransportationParameters()
{
  SetIntermediateLooperThresholds();
}

//------------------------------------------------------------------------------

G4TransportationParameters* G4TransportationParameters::Instance()
{
  if(theInstance == nullptr) {      
    G4MUTEXLOCK(&transportParamsMutex);
    if(theInstance == nullptr) {        
      static G4TransportationParameters manager;
      theInstance = &manager;
    }
    G4MUTEXUNLOCK(&transportParamsMutex);
  }
  return theInstance;
}

//------------------------------------------------------------------------------

G4bool G4TransportationParameters::IsLocked() const
{
  auto stateManager = G4StateManager::GetStateManager();

  auto state = stateManager->GetCurrentState();
  bool goodState = state == G4State_PreInit
                   || state == G4State_Init
                   || state == G4State_Idle ;
  return ( !G4Threading::IsMasterThread() || !goodState );
}

//------------------------------------------------------------------------------

G4bool  G4TransportationParameters::SetNumberOfTrials( G4int val )
{
   if(IsLocked()) { ReportLockError(__func__); return false; }
   fNumberOfTrials = val;
   return true;
}

//------------------------------------------------------------------------------

G4bool  G4TransportationParameters::SetWarningEnergy( double val )
{
   if(IsLocked()) { ReportLockError(__func__); return false; }
   fWarningEnergy = val;

   // Consistency check -- and trial fix
   if( fWarningEnergy > fImportantEnergy   ) {
      G4cerr << "G4TransportationParameters::GetWarningEnergy enforcing warning-E <= important-E "
             << "  resetting important energy from " << fImportantEnergy
             << " to " << val << G4endl;
      // "Enforcing Important Energy >= Warning Energy"
      fImportantEnergy = fWarningEnergy;
   }
   return true;      
}

//------------------------------------------------------------------------------

G4bool  G4TransportationParameters::SetImportantEnergy( double val )
{
   if(IsLocked()) { ReportLockError(__func__); return false; }
   fImportantEnergy = val;

   // Consistency check -- and trial fix   
   if( fImportantEnergy < fWarningEnergy ) {
      G4String mthd= G4String("G4TransportationParameters")+G4String(__func__);
      G4ExceptionDescription ed;
      ed<<"enforcing hierarchy (warning-E <= important-E): resetting important"
        <<" energy from " << fImportantEnergy << " to " << val << G4endl;
      G4Exception( mthd, "Enforcing Warning Energy <= Important Energy",
                   JustWarning, ed );

      fWarningEnergy = fImportantEnergy;
   }
   return true;
}

//------------------------------------------------------------------------------

G4bool G4TransportationParameters::SetWarningAndImportantEnergies( G4double warnE,
                                                                   G4double importE )
{
   if(IsLocked()) {
      ReportLockError(__func__);
      return false;
   }

   if( warnE <= importE )
   {
     fWarningEnergy = warnE;
     fImportantEnergy = importE;
   }
   else
   {
     fWarningEnergy = importE;
     fImportantEnergy = warnE;

     G4String mthd= G4String("G4TransportationParameters")+G4String(__func__);
     G4ExceptionDescription ed;
     ed << "To enforce hierarchy (warning-E <= important-E): "
        << " using smaller value= " <<  importE << " as Warning Energy "
        << " and larger value= " <<  warnE << " as Important Energy." << G4endl;
      G4Exception( mthd, "Enforcing Warning Energy <= Important Energy",
                   JustWarning, ed );
   }
   return true;
}

//------------------------------------------------------------------------------

void G4TransportationParameters::ReportLockError(G4String methodName, G4bool verbose) const
{
  // Report Incompatible States: GeomClosed , EventProc, (also Quit, Abort)
  G4String namesMethodClass= G4String("G4TransportationParameters") + methodName;

  auto stateManager = G4StateManager::GetStateManager();
  auto state    = stateManager->GetCurrentState();

  G4ExceptionDescription ed;
  ed << "Cannot change values of G4TransportationParameters when G4State is "
     << stateManager->GetStateString(state) << G4endl;
  ed << "Only the following Geant4 state are compatible: Pre_Init, Init and Idle." << G4endl;
  if( verbose ) {
     ed << G4endl << "Values remain as follows:" << G4endl;
     StreamInfo( ed );
  }
  G4Exception( namesMethodClass,
               "Locked, due to incompatible G4state: it not possible to change its parameters.",
               JustWarning, ed );
}

//------------------------------------------------------------------------------

void G4TransportationParameters::StreamInfo(std::ostream& os) const
{
  auto prec = os.precision(5);

  os << "Transport Parameters:  " << G4endl;
  os << "   Warning   energy = " << GetWarningEnergy()   / CLHEP::MeV << " MeV " << G4endl;
  os << "   Important energy = " << GetImportantEnergy() / CLHEP::MeV << " MeV " << G4endl;
  os << "   Number of trials = " << GetNumberOfTrials() << G4endl;
  os.precision(prec);
}

//------------------------------------------------------------------------------

void G4TransportationParameters::Dump() const
{
  G4MUTEXLOCK(&transportParamsMutex);
  StreamInfo(G4cout);
  G4MUTEXUNLOCK(&transportParamsMutex);
}

//------------------------------------------------------------------------------

G4bool G4TransportationParameters::SetHighLooperThresholds()
{
  if(IsLocked()) { return false; }
  
  // Restores the old high values -- potentially appropriate for energy-frontier
  //   HEP experiments.
  // Caution:  All tracks with E < 100 MeV that are found to loop are 
  SetWarningEnergy(    100.0 * CLHEP::MeV ); // Warn above this energy
  SetImportantEnergy(  250.0 * CLHEP::MeV ); // Extra trial above this En

  G4int maxTrials = 10;
  SetNumberOfTrials( maxTrials );
  
  return true;
}

//------------------------------------------------------------------------------

G4bool G4TransportationParameters::SetIntermediateLooperThresholds()
{
  if(IsLocked()) { return false; } // Currently must not change during loop - may relax

  // Medium values -- reasonable default for an unknown application
  fWarningEnergy=     1.0 * CLHEP::MeV; // Warn above this energy
  fImportantEnergy=  10.0 * CLHEP::MeV; // Extra trial above this En
  fNumberOfTrials= 10;
  return true;
}

//------------------------------------------------------------------------------

G4bool G4TransportationParameters::SetLowLooperThresholds() // Values for low-E applications
{
  if(IsLocked()) { return false; }
  
  // These values were the default in Geant4 10.5 - beta
  SetWarningEnergy(     1.0 * CLHEP::keV ); // Warn above this En
  SetImportantEnergy(   1.0 * CLHEP::MeV ); // Extra trials above it

  G4int maxTrials = 30;
  SetNumberOfTrials( maxTrials );

  return true;  
}

//------------------------------------------------------------------------------

G4bool G4TransportationParameters::EnableUseOfMagneticMoment(G4bool useMoment)
{
  if(IsLocked()) { return false; }
  fUseMagneticMoment= useMoment;
  return true;
}

//------------------------------------------------------------------------------


G4bool G4TransportationParameters::SetSilenceAllLooperWarnings(G4bool val)
{
  if(IsLocked()) { return false; }
  fSilenceLooperWarnings= val;
  return true;
}   
