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
// ClassName:   G4QMessenger
//
// Author: 2009 M. V. Kossov
//
// Modified:
//
//----------------------------------------------------------------------------
//

#ifndef G4QMessenger_h
#define G4QMessenger_h 1

class G4EmExtraPhysics;

#include "G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

class G4QNeutrinoPhysics;
class G4QPhotoNuclearPhysics;

class G4QMessenger: public G4UImessenger
{
  G4QMessenger();

public:
  ~G4QMessenger();

  static G4QMessenger* GetPointer(); // Gives a pointer to this singletone

  void Add(G4QNeutrinoPhysics* weak);
  void Add(G4QPhotoNuclearPhysics* photo);

  void SetNewValue(G4UIcommand* aComm, G4String aS);
  G4String GetCurrentValue(G4UIcommand* aComm);

private:
  G4UIdirectory*         rootDir;

  G4UIdirectory*         weakDir;
  G4QNeutrinoPhysics*    theWeak;
  G4UIcmdWithAString* theNuElN;
  G4UIcmdWithAString* theNuMuN;
  G4UIcmdWithAString* theNuTaN;
  G4UIcmdWithADouble* biasNuNuc;

  G4UIdirectory*          photoDir;
  G4QPhotoNuclearPhysics* thePhoto;
  G4UIcmdWithAString* theSynchR;
  G4UIcmdWithADouble* minGamSR;
  G4UIcmdWithAString* theGamN;
  G4UIcmdWithAString* theEleN;
  G4UIcmdWithAString* theMuoN;
  G4UIcmdWithAString* theTauN;
  G4UIcmdWithADouble* biasPhotoN;

};

#endif
