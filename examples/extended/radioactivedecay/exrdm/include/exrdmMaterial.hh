//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
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
