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
// $Id: G4EmSaturation.hh,v 1.4 2008-03-14 14:05:57 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
#ifndef G4EmSaturation_h
#define G4EmSaturation_h 1

// -------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4EmSaturation
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 18.02.2008
//
// Modifications:
//
//
// Class Description:
//   Compution on saturation effect, which reduce visible energy 
//   deposition at the step. Default implementation takes into 
//   account Birks effect.
// 
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Step;
class G4Material;
class G4ParticleDefinition;
class G4LossTableManager;

class G4EmSaturation
{
 public: 
   G4EmSaturation();
   virtual ~G4EmSaturation();

   void Initialise();
   G4double BirksAttenuation(const G4Step*);
   
 private:
   G4double ComputeAeff(G4Material*);
   
  
 private:
   G4ParticleDefinition* proton;
   G4LossTableManager*   tableManager;
   G4bool                initialised; 
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif

