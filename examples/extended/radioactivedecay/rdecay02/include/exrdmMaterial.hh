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
#ifndef exrdmMaterial_HH
#define exrdmMaterial_HH
////////////////////////////////////////////////////////////////////////////////
#include "G4Material.hh"
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
		     G4double tem = STP_Temperature, 
		     G4double pres = STP_Pressure);
  G4Material* GetMaterial (G4int i)  {return Material[i];};
  G4Material* GetMaterial (G4String name)
    {return G4Material::GetMaterial(name);} ;
  G4int GetMaterialIndex (G4String);
  G4int GetNbOfMaterial () {return Material.size();};
  void  DeleteMaterial (G4int);
  void  DeleteMaterial (G4String);

  void  ListMaterial();

private:

  exrdmMaterialMessenger         *materialMessenger;

  std::vector<G4Material*>   Material;
  std::vector<G4Element*>    Element;
  std::vector<G4Isotope*>    Isotope;

private:
  static const G4String        ELU[110];
  static const G4String        ELL[110];
  static const G4String        EUU[110];
  static const G4double        A[110];
       
};
////////////////////////////////////////////////////////////////////////////////
#endif
