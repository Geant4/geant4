//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: G4EmExtraPhysics.hh,v 1.2 2005/12/01 18:19:45 gunter Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4EmExtraPhysics
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 10.11.2005 V.Ivanchenko edit to provide a standard
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

class G4EmExtraPhysics : public G4VPhysicsConstructor
{
public:
  G4EmExtraPhysics(const G4String& name = "EM extra");
  virtual ~G4EmExtraPhysics();

  void ConstructParticle();
  void ConstructProcess();

  void Synch(G4String & aState);
  void GammaNuclear(G4String & aState);

private:

  G4bool wasActivated;
  G4bool synchOn;
  G4bool gammNucOn;

  G4EmMessenger*           theMessenger;
  G4SynchrotronRadiation*  theElectronSynch;
  G4SynchrotronRadiation*  thePositronSynch;
  G4ElectroNuclearBuilder* theGNPhysics;
};

#endif





