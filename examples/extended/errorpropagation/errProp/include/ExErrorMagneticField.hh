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
/// \file errProp/include/ExErrorMagneticField.hh
/// \brief Definition of the ExErrorMagneticField class
//

#ifndef ExErrorMagneticField_H
#define ExErrorMagneticField_H

#include "G4UniformMagField.hh"

class G4FieldManager;

/// Magnetic field class
///
/// A uniform 1 kilogauss field along the Z axis
///
/// History:
/// Created:   May 2007
/// \author   P. Arce 
//------------------------------------------------------------------------

class ExErrorMagneticField: public G4UniformMagField
{
  public:
  
   ExErrorMagneticField(G4ThreeVector);  //  The value of the field
   ExErrorMagneticField();               //  A zero field
      
   //Set the field (0,0,fieldValue)
   void SetFieldValue(G4double fieldValue);
   void SetFieldValue(G4ThreeVector fieldVector);

protected:

  // Find the global Field Manager
  G4FieldManager* GetGlobalFieldManager();   // static 
};

#endif
