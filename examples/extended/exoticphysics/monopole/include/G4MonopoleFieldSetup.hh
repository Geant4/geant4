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
/// \file exoticphysics/monopole/include/G4MonopoleFieldSetup.hh
/// \brief Definition of the G4MonopoleFieldSetup class
//
// $Id: G4MonopoleFieldSetup.hh 104872 2017-06-23 14:19:16Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//
// ------------------------------------------------------------
//        GEANT 4  include file implementation
// ------------------------------------------------------------
//
// G4MonopoleFieldSetup is responsible for setting up a magnetic field
// and the ability to use it with two different equation of motions, 
// one for monopoles and another for the rest of the particles. 
// 
// 

// =======================================================================
// Created:  13 May 2010, B. Bozsogi
// =======================================================================

#ifndef G4MonopoleFieldSetup_H
#define G4MonopoleFieldSetup_H

#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"

class G4FieldManager;
class G4ChordFinder;
class G4Mag_UsualEqRhs;
class G4MagIntegratorStepper;
class G4MonopoleEquation;

class G4MonopoleFieldMessenger;

class G4MonopoleFieldSetup
{
public:  

  G4MonopoleFieldSetup();            
  ~G4MonopoleFieldSetup() ;     

  void InitialiseAll();    //  Set parameters and call method below
  //  void SetMagField(G4double fieldValue);
  void SetStepperAndChordFinder(G4int val);
  void SetMagField   (G4double v, bool checkIfAlreadyDefined=false);

  //  static G4MonopoleFieldSetup* GetMonopoleFieldSetup();
  double GetZmagFieldValue() const { return fZmagFieldValue; }
  
private:


  G4FieldManager*         GetGlobalFieldManager() ;   // static 

  G4FieldManager*         fFieldManager ;
  G4ChordFinder*          fChordFinder ;
  G4ChordFinder*          fUsualChordFinder ;
  G4ChordFinder*          fMonopoleChordFinder ;
  G4Mag_UsualEqRhs*       fEquation ; 
  G4MonopoleEquation*     fMonopoleEquation ;

  G4MagneticField*        fMagneticField ; 

  G4MagIntegratorStepper* fStepper ;
  G4MagIntegratorStepper* fMonopoleStepper ;

  G4double                fMinStep ;
  G4double                fZmagFieldValue;

  //  G4int                   fStepperIndex;

  //  static G4MonopoleFieldSetup*  fMonopoleFieldSetup;  
  G4MonopoleFieldMessenger*     fMonopoleFieldMessenger;
 
};

#endif
