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
// $Id$
//
//---------------------------------------------------------------------------
//
// ClassName:   G4EmQExtraPhysics
//
// Author: 16-Oct-2012 A. Ribon
//         Copied from the original G4EmExtraPhysics and renamed
//
// Modified:
//
//----------------------------------------------------------------------------
//

#ifndef G4EmQExtraPhysics_h
#define G4EmQExtraPhysics_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

#include "G4EmQMessenger.hh"

#include "G4SynchrotronRadiation.hh"
#include "G4ElectroNuclearBuilder.hh"
#include "G4MuonNuclearProcess.hh"
#include "G4MuonVDNuclearModel.hh"

class G4EmQExtraPhysics : public G4VPhysicsConstructor
{
public:

  G4EmQExtraPhysics(G4int ver = 1);

  // obsolete
  G4EmQExtraPhysics(const G4String& name);

  virtual ~G4EmQExtraPhysics();

  void ConstructParticle();
  void ConstructProcess();

  void Synch(G4String & aState);
  void GammaNuclear(G4String & aState);
  void MuonNuclear(G4String & aState);

private:

  void BuildSynch();
  void BuildGammaNuclear();
  void BuildMuonNuclear();

  G4bool wasBuilt;
  G4bool gnActivated;
  G4bool munActivated;
  G4bool synActivated;
  G4bool synchOn;
  G4bool gammNucOn;
  G4bool muNucOn;

  G4EmQMessenger*           theMessenger;
  G4SynchrotronRadiation*  theElectronSynch;
  G4SynchrotronRadiation*  thePositronSynch;
  G4ElectroNuclearBuilder* theGNPhysics;
  G4MuonNuclearProcess* muNucProcess;
  G4MuonVDNuclearModel* muNucModel;

  G4int verbose;
};

#endif





