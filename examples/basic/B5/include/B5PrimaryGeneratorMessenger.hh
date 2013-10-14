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
/// \file B5PrimaryGeneratorMessenger.hh
/// \brief Definition of the B5PrimaryGeneratorMessenger class

#ifndef B5PrimaryGeneratorMessenger_h
#define B5PrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class B5PrimaryGeneratorAction;

class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithABool;

/// Primary generator messenger
/// 
/// It implements commands:
/// - /mydet/momentum value unit
/// - /mydet/sigmaMomentum value unit
/// - /mydet/sigmaAngle value unit
/// - /mydet/randomizePrimary true|false

class B5PrimaryGeneratorMessenger: public G4UImessenger
{
public:
    B5PrimaryGeneratorMessenger(B5PrimaryGeneratorAction* generator);
    ~B5PrimaryGeneratorMessenger();
    
    void SetNewValue(G4UIcommand * command,G4String newValues);
    G4String GetCurrentValue(G4UIcommand * command);
    
private:
    B5PrimaryGeneratorAction * fGenerator;
    
    G4UIcmdWithADoubleAndUnit*  fMomentumCmd;
    G4UIcmdWithADoubleAndUnit*  fSigmaMomCmd;
    G4UIcmdWithADoubleAndUnit*  fSigmaAngCmd;
    G4UIcmdWithABool*           fRandomCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
