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
//

//---------------------------------------------------------------------------
//
// ClassName:   G4ExtendedMaterial
//
// Description: Contains extended material properties
//
// Class description:
//
// Is used to define the additional material information. This class
// contains a map of G4VMaterialExtension associated with an integer key.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4EXTENDEDMATERIAL_HH
#define G4EXTENDEDMATERIAL_HH 1

#include "G4Material.hh"
#include <unordered_map>
#include <memory>
#include "G4VMaterialExtension.hh"//Needed for hash defintion


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


// A map for material extensions based on the hash of the name.
// Extensions are owned by the map
using G4MaterialExtensionMap=std::unordered_map<G4String, //KEY
						std::unique_ptr<G4VMaterialExtension>,//VALUE
						G4MaterialExtensionHash>;//HASHING FUNCTOR

class G4ExtendedMaterial : public G4Material
{
public:  // with description
  //
  // Constructor to create an extended material from the base-class G4Material
  //
  G4ExtendedMaterial(const G4String& name,                              //its name
                     const G4Material* baseMaterial);			//base material

  //
  // Constructor to create an extended material from single element
  //
  G4ExtendedMaterial(const G4String& name,				//its name
                   G4double  z, 				//atomic number
                   G4double  a,					//mass of mole
                   G4double  density, 				//density
                   G4State   state    = kStateUndefined,	//solid,gas
                   G4double  temp     = NTP_Temperature,	//temperature
                   G4double  pressure = CLHEP::STP_Pressure);	//pressure

  //
  // Constructor to create an extended material from a combination of elements
  // and/or materials subsequently added via AddElement and/or AddMaterial
  //
  G4ExtendedMaterial(const G4String& name,				//its name
                   G4double  density, 				//density
                   G4int     nComponents,			//nbOfComponents
                   G4State   state    = kStateUndefined,	//solid,gas
                   G4double  temp     = NTP_Temperature,	//temperature
                   G4double  pressure = CLHEP::STP_Pressure);	//pressure

  //
  // Constructor to create an extended material from the base extended material
  //
  G4ExtendedMaterial(const G4String& name,				//its name
                   G4double  density, 				//density
             const G4ExtendedMaterial* baseMaterial,			//base material
                   G4State   state    = kStateUndefined,	//solid,gas
                   G4double  temp     = NTP_Temperature,	//temperature
                   G4double  pressure = CLHEP::STP_Pressure);	//pressure

public:
 ~G4ExtendedMaterial() override = default;

private:
  G4MaterialExtensionMap fExtensionMap;
public:  // with description
  //
  // register G4VMaterialExtension
  // This class owns extensions. Register with:
  // RegisterExtension(std::unique_ptr<MyExtension>(new MyExtension("name")));
  // or:
  // RegisteerExtension(std::make_unique<MyExtension>("name"));
  void RegisterExtension(std::unique_ptr<G4VMaterialExtension> extension);
  //
  // retrieve G4VMaterialExtension, null pointer is returned if model is not available
  G4VMaterialExtension* RetrieveExtension(const G4String& name);

  inline G4int GetNumberOfExtensions() const
  { return G4int(fExtensionMap.size()); } 

  // Retrieve iterators, proxyes to c++ methods. These are const for read-only
  // access. Use Register/RetreiveExtension to modify map
  G4MaterialExtensionMap::const_iterator begin() const { return fExtensionMap.begin(); }
  G4MaterialExtensionMap::const_iterator cbegin() const { return fExtensionMap.cbegin(); }
  G4MaterialExtensionMap::const_iterator end() const { return fExtensionMap.end(); }
  G4MaterialExtensionMap::const_iterator cend() const { return fExtensionMap.cend(); }

public:
 G4bool IsExtended() const override;
 void Print(std::ostream& flux) const;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
