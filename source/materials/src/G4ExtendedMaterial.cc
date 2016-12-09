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
// $Id: G4ExtendedMaterial.cc 96792 2016-05-09 09:18:43Z vnivanch $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
// 18-10-07, move definition of material index to InitialisePointers (V.Ivanchenko) 
// 13-08-08, do not use fixed size arrays (V.Ivanchenko)
// 26-10-11, new scheme for G4Exception  (mma)
// 13-04-12, map<G4Material*,G4double> fMatComponents, filled in AddMaterial()
// 21-04-12, fMassOfMolecule, computed for AtomsCount (mma)
// 
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4ExtendedMaterial.hh"
#include "G4VMaterialExtension.hh"
#include "G4PhysicsModelCatalog.hh"

// Constructor to create an extended material from the base-class G4Material

G4ExtendedMaterial::G4ExtendedMaterial(const G4String& name,
                     const G4Material* baseMaterial)
  : G4Material(name,baseMaterial->GetDensity(),baseMaterial,
               baseMaterial->GetState(),baseMaterial->GetTemperature(),
               baseMaterial->GetPressure())
{;}

// Constructor to create an extended material from single element

G4ExtendedMaterial::G4ExtendedMaterial(const G4String& name, G4double z,
                       G4double a, G4double density, 
                       G4State state, G4double temp, G4double pressure)
  : G4Material(name,z,a,density,state,temp,pressure)
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Constructor to create an extended material from a combination of elements
// (elements and/or materials)  added with AddElement or AddMaterial

G4ExtendedMaterial::G4ExtendedMaterial(const G4String& name, G4double density,
                       G4int nComponents,
                       G4State state, G4double temp, G4double pressure)
  : G4Material(name,density,nComponents,state,temp,pressure)
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Constructor to create an extended material from the base extended material

G4ExtendedMaterial::G4ExtendedMaterial(const G4String& name, G4double density,
                       const G4ExtendedMaterial* bmat,
                       G4State state, G4double temp, G4double pressure)
  : G4Material(name,density,bmat,state,temp,pressure)
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


// register G4VMaterialExtension

void G4ExtendedMaterial::RegisterExtension(std::unique_ptr<G4VMaterialExtension> extension)
{
  auto iter = fExtensionMap.find(extension->GetName());
  if(iter!=fExtensionMap.end())
  { 
      G4ExceptionDescription msg;
      msg << "G4ExtendedMaterial <"<<GetName()<<"> already has extension for "
	  << extension->GetName()
	  << ". Extension is replaced.";
    G4Exception("G4ExtendedMaterial::RegisterExtension(...)","MatExt001",JustWarning,msg);
  }
  fExtensionMap.insert(std::make_pair(extension->GetName(),std::move(extension)));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// retrieve G4VMaterialExtension, null pointer is returned if model is not available

G4VMaterialExtension* G4ExtendedMaterial::RetrieveExtension(const G4String& name)
{
  const auto iter = fExtensionMap.find(name);
  if(iter!=fExtensionMap.end())
  { return iter->second.get(); }
  else
  {
      G4ExceptionDescription msg;
      msg << "G4ExtendedMAterial <"<<GetName()<<"> cannot find extension for "
	  << name;
      G4Exception("G4ExtendedMaterial::RetreiveExtension(...)","MatExt002",JustWarning,msg);
      return nullptr;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4ExtendedMaterial::IsExtended() const
{ return true; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ExtendedMaterial::Print(std::ostream& flux) const
{ 
  flux << "\n Registered material extensions :\n";
  auto iter = fExtensionMap.begin();
  for(;iter!=fExtensionMap.end();iter++)
  { flux << "     " << iter->first << "\n"; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
