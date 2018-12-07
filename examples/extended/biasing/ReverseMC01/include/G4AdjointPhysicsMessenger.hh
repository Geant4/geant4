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
/// \file biasing/ReverseMC01/include/G4AdjointPhysicsMessenger.hh
/// \brief Definition of the G4AdjointPhysicsMessenger class
//
//
//////////////////////////////////////////////////////////////
//  Class Name:        G4AdjointPhysicsMessenger
//        Author:               L. Desorgher
//        Organisation:         SpaceIT GmbH
//        Contract:        ESA contract 21435/08/NL/AT
//        Customer:             ESA/ESTEC
//////////////////////////////////////////////////////////////
// CHANGE HISTORY
//--------------
//      ChangeHistory:
//                 17-11-2009 creation by L. Desorgher
//
//-------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4AdjointPhysicsMessenger_h
#define G4AdjointPhysicsMessenger_h 1
#include "globals.hh"
#include "G4UImessenger.hh"
class G4AdjointPhysicsList;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithABool;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;
class G4UIcmdWithADouble;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4AdjointPhysicsMessenger: public G4UImessenger
{
 public:
  G4AdjointPhysicsMessenger(G4AdjointPhysicsList* );
  virtual ~G4AdjointPhysicsMessenger();
  virtual void SetNewValue(G4UIcommand*, G4String);
    
 private:
  G4AdjointPhysicsList* fPhysicsList;
  G4UIdirectory*        fPhysicsDir;
    
  //Physics Model
  G4UIcmdWithABool*  fUsepIonisationCmd;
  G4UIcmdWithABool*  fUseBremCmd;
  G4UIcmdWithABool*  fUseComptonCmd;
  G4UIcmdWithABool*  fUseMSCmd;
  G4UIcmdWithABool*  fUsePEEffectCmd;
  G4UIcmdWithABool*  fUseGammaConversionCmd;
  G4UIcmdWithABool*  fUseEgainFluctuationCmd;
  G4UIcmdWithADoubleAndUnit* fSetEminAdjModelsCmd;
  G4UIcmdWithADoubleAndUnit* fSetEmaxAdjModelsCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

