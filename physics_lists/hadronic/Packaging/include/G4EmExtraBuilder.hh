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
// $Id: G4EmExtraBuilder.hh,v 1.1 2005-11-11 22:56:07 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4EmExtraBuilder
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 10.11.2005 V.Ivanchenko edit to provide a standard
//
//----------------------------------------------------------------------------
//

#ifndef G4EmExtraBuilder_h
#define G4EmExtraBuilder_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

class G4EmMessenger;
class G4SynchrotronRadiation;
class G4ElectroNuclearBuilder;

class G4EmExtraBuilder : public G4VPhysicsConstructor
{
public:
  G4EmExtraBuilder(const G4String& name = "EM extra");
  virtual ~G4EmExtraBuilder();

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





