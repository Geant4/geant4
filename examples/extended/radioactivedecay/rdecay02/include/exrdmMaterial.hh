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
/// \file radioactivedecay/rdecay02/include/exrdmMaterial.hh
/// \brief Definition of the exrdmMaterial class
//
#ifndef exrdmMaterial_HH
#define exrdmMaterial_HH
////////////////////////////////////////////////////////////////////////////////
#include "G4Material.hh"
#include "CLHEP/Units/PhysicalConstants.h"
#include <vector>

class exrdmMaterialMessenger;
////////////////////////////////////////////////////////////////////////////////
//
class exrdmMaterial
{
public:

  exrdmMaterial ();
  ~exrdmMaterial ();

public:

  void  AddMaterial (G4String, G4String, G4double, G4String, 
                     G4double tem = CLHEP::STP_Temperature, 
                     G4double pres = CLHEP::STP_Pressure);
  G4Material* GetMaterial (G4int i)  {return fMaterial[i];};
  G4Material* GetMaterial (G4String name)
    {return G4Material::GetMaterial(name);} ;
  G4int GetMaterialIndex (G4String);
  G4int GetNbOfMaterial () {return fMaterial.size();};
  void  DeleteMaterial (G4int);
  void  DeleteMaterial (G4String);

  void  ListMaterial();

private:

  exrdmMaterialMessenger         *fMaterialMessenger;

  std::vector<G4Material*>   fMaterial;
  std::vector<G4Element*>    fElement;
  std::vector<G4Isotope*>    fIsotope;

private:
  static const G4String        fELU[110];
  static const G4String        fELL[110];
  static const G4String        fEUU[110];
  static const G4double        fA[110];
       
};
////////////////////////////////////////////////////////////////////////////////
#endif
