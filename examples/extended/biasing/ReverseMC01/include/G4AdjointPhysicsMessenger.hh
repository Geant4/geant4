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
// $Id: G4AdjointPhysicsMessenger.hh,v 1.1 2009-11-19 22:41:18 ldesorgh Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//////////////////////////////////////////////////////////////
//      Class Name:	G4AdjointPhysicsMessenger
//	Author:       	L. Desorgher
// 	Organisation: 	SpaceIT GmbH
//	Contract:	ESA contract 21435/08/NL/AT
// 	Customer:     	ESA/ESTEC
//////////////////////////////////////////////////////////////
// CHANGE HISTORY
//--------------
//      ChangeHistory:
//	 	17-11-2009 creation by L. Desorgher
//
//-------------------------------------------------------------
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
   ~G4AdjointPhysicsMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    G4AdjointPhysicsList* thePhysicsList;
    

    G4UIdirectory*             PhysicsDir;
    
    //Physics Model
    
    G4UIcmdWithABool*  UseIonisationCmd;
    G4UIcmdWithABool*  UsepIonisationCmd;
    G4UIcmdWithABool*  UseBremCmd;
    G4UIcmdWithABool*  UseComptonCmd;
    G4UIcmdWithABool*  UseMSCmd;
    G4UIcmdWithABool*  UsePEEffectCmd;
    G4UIcmdWithABool*  UseGammaConversionCmd;
    G4UIcmdWithABool*  UseEgainFluctuationCmd;
  
    G4UIcmdWithADoubleAndUnit* SetEminAdjModelsCmd;
    G4UIcmdWithADoubleAndUnit* SetEmaxAdjModelsCmd;

/*
#ifdef TEST_MODE    
    G4UIcmdWithADouble* SetCSBiasingFactorComptonCmd;
    G4UIcmdWithADouble* SetCSBiasingFactorBremCmd;
    G4UIcmdWithADouble* SetCSBiasingFactorIonisationCmd;
    G4UIcmdWithADouble* SetCSBiasingFactorPEeffectCmd;
 
#endif  
*/   
    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

