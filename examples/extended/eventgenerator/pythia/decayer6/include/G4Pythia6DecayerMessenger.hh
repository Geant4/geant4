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
// $Id: G4Pythia6DecayerMessenger.hh 100687 2016-10-31 11:20:33Z gcosmo $
//
/// \file eventgenerator/pythia/decayer6/include/G4Pythia6DecayerMessenger.hh
/// \brief Definition of the G4Pythia6DecayerMessenger class

#ifndef G4_PYTHIA6_DECAYER_MESSENGER_H
#define G4_PYTHIA6_DECAYER_MESSENGER_H 

#include <G4UImessenger.hh>
#include <globals.hh>

class G4Pythia6Decayer;

class G4UIdirectory;
class G4UIcmdWithAnInteger;
class G4UIcmdWithABool;

/// Messenger class that defines commands for G4Pythia6Decayer.
///
/// Implements command
/// - /pythia6Decayer/verbose [level]
/// - /pythia6Decayer/forceDecayType [decayType]

class G4Pythia6DecayerMessenger : public G4UImessenger
{
  public:

    G4Pythia6DecayerMessenger(G4Pythia6Decayer* pythia6Decayer);
    virtual ~G4Pythia6DecayerMessenger();
   
    virtual void SetNewValue(G4UIcommand* command, G4String string);
    
  private:

    /// Not implemented
    G4Pythia6DecayerMessenger();
    /// Not implemented
    G4Pythia6DecayerMessenger(const G4Pythia6DecayerMessenger& right);
    /// Not implemented
    G4Pythia6DecayerMessenger& operator=(const G4Pythia6DecayerMessenger& r);

  private:

    G4Pythia6Decayer*      fPythia6Decayer;    ///< associated class 
    G4UIdirectory*         fDirectory;         ///< command directory
    G4UIcmdWithAnInteger*  fVerboseCmd;        ///< command: verbose
    G4UIcmdWithAnInteger*  fDecayTypeCmd;      ///< command: forceDEcayeType
};

// ----------------------------------------------------------------------------

#endif
