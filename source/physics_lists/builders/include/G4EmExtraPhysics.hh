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
// $Id: G4EmExtraPhysics.hh,v 1.3 2010-06-02 17:21:29 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4EmExtraPhysics
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 10.11.2005 V.Ivanchenko edit to provide a standard
// 19.06.2006 V.Ivanchenko add mu-nuclear process
//
//----------------------------------------------------------------------------
//

#ifndef G4EmExtraPhysics_h
#define G4EmExtraPhysics_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

#include "G4EmMessenger.hh"

#include "G4SynchrotronRadiation.hh"
#include "G4ElectroNuclearBuilder.hh"
#include "G4MuNuclearInteraction.hh"

class G4EmExtraPhysics : public G4VPhysicsConstructor
{
public:

  G4EmExtraPhysics(G4int ver = 1);

  // obsolete
  G4EmExtraPhysics(const G4String& name);

  virtual ~G4EmExtraPhysics();

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

  G4EmMessenger*           theMessenger;
  G4SynchrotronRadiation*  theElectronSynch;
  G4SynchrotronRadiation*  thePositronSynch;
  G4ElectroNuclearBuilder* theGNPhysics;
  G4MuNuclearInteraction*  theMuNuc1;
  G4MuNuclearInteraction*  theMuNuc2;

  G4int verbose;
};

#endif





