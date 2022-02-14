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
//
//
// class G4TransportationLogger Implementation
//
// Author: J. Apostolakis, June 2018
//
// --------------------------------------------------------------------

#include <iomanip>

#include "G4SystemOfUnits.hh"
#include "G4TransportationLogger.hh"
#include "G4Track.hh"
#include "G4Step.hh"

G4TransportationLogger::G4TransportationLogger(const G4String& className,
                                                     G4int verbosity)
   : fClassName(className), fVerbose(verbosity),
     fThldWarningEnergy(0.), fThldImportantEnergy(0.), fThldTrials(0)
{
}

G4TransportationLogger::G4TransportationLogger(const char* className,
                                                           G4int verbosity)
   : fClassName(className), fVerbose(verbosity),
     fThldWarningEnergy(0.), fThldImportantEnergy(0.), fThldTrials(0)
{
}

G4TransportationLogger::~G4TransportationLogger()
{
}

// ********************************************************************
// SetThresholds
// ********************************************************************
void G4TransportationLogger::
SetThresholds( G4double newEnWarn, G4double importantEnergy,
               G4int newMaxTrials )
{
   SetThresholdWarningEnergy( newEnWarn );
   SetThresholdImportantEnergy( importantEnergy );
   SetThresholdTrials(newMaxTrials );
}

/////////////////////////////////////////////////////////////////////////////
//
void
G4TransportationLogger::ReportLooperThresholds( const char* className )
{
   G4cout << className << ":  Current values for thresholds related to "
          << " the killing of looping tracks: " << G4endl
          <<  "    Warning Energy   = " << GetThresholdWarningEnergy() / CLHEP::MeV << " MeV "
          <<  "  ( below this tracks are killed without warning ) " << G4endl
          <<  "    Important Energy = " << GetThresholdImportantEnergy() / CLHEP::MeV
          <<  "  ( above this tracks are given multiple chances ) " << G4endl
          <<  "    Extra Trials     = " << GetThresholdTrials()
          << " 'important' tracks, i.e. those above 'important' energy "
          << G4endl;
}

// ********************************************************************
// ReportLoopingTrack
// ********************************************************************
//
void G4TransportationLogger::ReportLoopingTrack( const G4Track & track,
                                                 const G4Step  & stepData,
                                                 G4int           numTrials,
                                                 G4long          noCalls,
                                                 const char* methodName
   ) const
{
  static std::atomic<unsigned int> numAdviceExcessSteps(0);
  constexpr double gramPerCm3 = gram / ( centimeter * centimeter * centimeter ) ;
  std::ostringstream msg;
  auto preStepPt= stepData.GetPreStepPoint();
  auto preStepEn= preStepPt ? preStepPt->GetKineticEnergy() / MeV : -1.0 ;
  msg << " Transportation is killing track that is looping or stuck. " << G4endl
      << "   Track is "
      << track.GetParticleDefinition()->GetParticleName()
      << " and has " << track.GetKineticEnergy() / MeV
      << " MeV energy  ( pre-Step = " << preStepEn << " ) " << G4endl;
  msg << "   momentum = " << track.GetMomentum() << " mag= " << track.GetMomentum().mag()
      << G4endl
      << "   position = " << track.GetPosition();
  auto physVolume= track.GetVolume();
  auto material= physVolume->GetLogicalVolume()->GetMaterial();
  msg << " is in volume '" << physVolume->GetName() << "', ";
  if( material )
  {
     msg << " its material is '" << material->GetName() << "'";
     msg << " with density = " << material->GetDensity() / gramPerCm3
         << " g/cm^3 ";
  }
  else
  {
     msg << " unable to obtain material information (including density.) ";
  }
  msg << G4endl;
  msg << " Total number of Steps by this track: " << track.GetCurrentStepNumber()
      << G4endl
      << " Length of this step = " << stepData.GetStepLength() / mm << " mm "
      << G4endl     
      << " Number of propagation trials = " << numTrials
      << " ( vs maximum = " << GetThresholdTrials() << " for 'important' particles ) " 
      << G4endl;

  if (noCalls)
    msg << "   ( Number of *calls* of Transport/AlongStepDoIt = " << noCalls << " )" << G4endl;

  const G4int numPrints= 5;
  if( numAdviceExcessSteps++ < numPrints )
  { 
     msg << " =============== Recommendations / advice ====================" << G4endl;             
     msg << " Recommendations to address this issue (Transport-001-ExcessSteps)" << G4endl;
     msg << " This warning is controlled by the SetThresholdWarningEnergy "
         << " method of G4Transportation.  " << G4endl
         << " Current value of 'warning' threshold= "
         << GetThresholdWarningEnergy() / MeV << " MeV " << G4endl;
     msg << " - If 'unimportant' particles (with energy low enough not to matter in your "
         << "  application, then increase its value. "   << G4endl;
     msg << " - If particles of high-enough energy to be important are being "
         << " killed, you can " << G4endl
         << "   a) Increase the trial steps using the method  SetThresholdTrials().  "
         << "  Particles above the 'important' threshold " << G4endl
         << "  will be given this many 'chances'."
         << "  The default value was 10, and the current value is " << GetThresholdTrials()
         << G4endl
         << "   b) Increase the energy which you consider 'important' (above this they are"
         << " killed only after extra trials), using the method SetThresholdImportantEnergy() " << G4endl
         << "      Note: this can incur a potentially high cost in extra simulation time " 
         << " if more tracks require very large number of integration steps . " << G4endl
         << "   c) investigate alternative integration methods " << G4endl
         << "    e.g.  Helical methods for uniform or almost uniform fields"
         << " or else higher order RK methods such as DormandPrince78 "
         << G4endl;
     msg << " This information is provided " << numPrints << " times. Current count: "
         << numAdviceExcessSteps << " / " << numPrints << G4endl;
     msg << " =============================================================" << G4endl;
  }
  const G4String fullMethodName= fClassName + "::" + methodName;
  G4Exception( fullMethodName, "Transport-001-ExcessSteps", JustWarning, msg);
}
