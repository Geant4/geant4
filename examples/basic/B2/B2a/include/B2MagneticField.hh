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
/// \file B2MagneticField.hh
/// \brief Definition of the B2MagneticField class

#ifndef B2MagneticField_h
#define B2MagneticField_h 1

#include "G4UniformMagField.hh"

class G4FieldManager;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

///  A class for control of the Magnetic Field of the detector.
///  The field is assumed to be uniform.
 
class B2MagneticField : public G4UniformMagField
{
  public:
    B2MagneticField();              //  A zero field
    B2MagneticField(G4ThreeVector); //  The value of the field
    virtual ~B2MagneticField();  
      
    //Set the field (fieldValue,0,0)
    void SetMagFieldValue(G4double fieldValue);
    void SetMagFieldValue(G4ThreeVector fieldVector);
      
  protected:
    // Find the global Field Manager
    G4FieldManager* GetGlobalFieldManager(); // static 
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
