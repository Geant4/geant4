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
// Class Description:
//
// A utility static class, responsable for keeping common parameters
// for all transportation processes.  These parameters can be relevant only
// to some types of particles, e.g. only charged particles.
//
// It must initialized only by the master thread.
// It can be updated before the transportation processes have been instantiated or
//   after their instantiation.
// Note: It can be updated only if the state of the simulation is PreInit, Init or Idle.
// Parameters may be used in run time or at initialisation
// 
// Author:  J. Apostolakis Nov 2022 
// -------------------------------------------------------------------
//

#ifndef G4TransportationParameters_hh
#define G4TransportationParameters_hh

#include "globals.hh"

class G4TransportationParameters
{
public:
  static   G4TransportationParameters* Instance();

  ~G4TransportationParameters() = default;

  G4bool   SetNumberOfTrials(  G4int val );
  G4bool   SetWarningEnergy(   G4double val );
  G4bool   SetImportantEnergy( G4double val );

  G4int    GetNumberOfTrials() const { return fNumberOfTrials; }
  G4double GetWarningEnergy() const { return fWarningEnergy; }
  G4double GetImportantEnergy() const { return fImportantEnergy; }
  G4bool   IsMagneticMomentEnabled() const { return fUseMagneticMoment; }

  G4bool   EnableUseOfMagneticMoment(G4bool useMoment=true);
   // Whether to deflect particles with force due to magnetic moment

  G4bool   SetHighLooperThresholds(); // Shortcut method - old values (meant for HEP)
  G4bool   SetIntermediateLooperThresholds();  // Intermediate values - also used as default
  G4bool   SetLowLooperThresholds();  // Set low thresholds - for low-E applications
  G4bool   SetSilenceAllLooperWarnings(G4bool val=true);
  // return value = success or failure of setting the parameter
  
  G4bool   GetSilenceAllLooperWarnings(){ return fSilenceLooperWarnings; }

  // Probe whether the 'default' instance exists, without creating it
  static G4bool  Exists() { return theInstance != nullptr; } 

   // printing
  void StreamInfo(std::ostream& os) const;
  void Dump() const;
  friend std::ostream& operator<< (std::ostream& os, const G4TransportationParameters&);
   
private:
  G4TransportationParameters();
  //  Currently private - but potentially will open it up, to allow per-particle specialisation

  // void   Initialise();

  G4bool IsLocked() const;

  void   PrintWarning(G4ExceptionDescription& ed) const; 

private:
  static   G4TransportationParameters* theInstance;

  // STATE
  // Values for initialising 'loopers' parameters of Transport process
  G4double fWarningEnergy   =  -1.0;  //  Warn above this energy
  G4double fImportantEnergy =  -1.0;  //  Give a few trials above this E
  G4int    fNumberOfTrials  =    10;  //  Number of trials an important looper survives

  // Flags for use of gravity field(s) or fields which interact with the magnetic moment
  G4bool   fUseMagneticMoment = false;
  G4bool   fUseGravity        = false;

  // Flag to *Supress* all 'looper' warnings   
  G4bool   fSilenceLooperWarnings= false;  
};

#endif
