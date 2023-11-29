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
// ClassName:   G4EmMessenger
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 09.11.2005 V.Ivanchenko edit to provide a standard
// 19.06.2006 V.Ivanchenko add mu-nuclear process
// 31.01.2018 V.Grichine, activation of neutrino-electron process
// 03.10.2018 V Grichine activation of total nneutrino-electron process 
//----------------------------------------------------------------------------
//

#ifndef G4EmMessenger_h
#define G4EmMessenger_h 1

class G4EmExtraPhysics;

#include "G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//////////////////////////////////////

class G4EmMessenger: public G4UImessenger
{
public:
  explicit G4EmMessenger(G4EmExtraPhysics* af);
  ~G4EmMessenger() override;

  void SetNewValue(G4UIcommand* aComm, G4String aS) override;

private:
  G4EmExtraPhysics*   theB;

  G4UIcmdWithABool*   theSynch;
  G4UIcmdWithABool*   theSynchAll;
  G4UIcmdWithABool*   theGN;
  G4UIcmdWithABool*   theGLENDN;
  G4UIcmdWithABool*   theEN;
  G4UIcmdWithABool*   theMUN;
  G4UIcmdWithABool*   theGMM;
  G4UIcmdWithABool*   theMMM;
  G4UIcmdWithABool*   thePMM;
  G4UIcmdWithABool*   thePH;
  G4UIcmdWithABool*   theNu;
  G4UIcmdWithABool*   theNuETX;
  G4UIcmdWithABool*   theXS;

  G4UIcmdWithADouble* theGMM1;
  G4UIcmdWithADouble* thePMM1;
  G4UIcmdWithADouble* thePH1;
  G4UIcmdWithADouble* theNuEleCcBF;
  G4UIcmdWithADouble* theNuEleNcBF;
  G4UIcmdWithADouble* theNuNucleusBF;
  G4UIcmdWithADoubleAndUnit* theGNlowe;

  G4UIcmdWithAString* theNuDN;

  G4UIdirectory*      aDir1;
  G4UIdirectory*      aDir2;
};

#endif
