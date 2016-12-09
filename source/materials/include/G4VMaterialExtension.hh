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
// $Id: G4VMaterialExtension.hh 96776 2016-05-06 14:11:00Z vnivanch $
//

//---------------------------------------------------------------------------
//
// ClassName:   G4VMaterialExtension
//
// Description: Contains extended material properties
//
// Class description:
//
// Is used to define the additional material information. This class
// is an abstract base class to be extended for each extension.
// Object(s) of G4VMaterialExtension should be registered to
// G4ExtendedMaterial.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4VMATERIALEXTENSION_HH
#define G4VMATERIALEXTENSION_HH 1

#include "globals.hh"
#include "G4ios.hh"
#include <functional>
#include <string>
#include "G4String.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//Hash of extension material name
using G4MaterialExtensionHash=std::hash<std::string>;

class G4VMaterialExtension
{
public:  // with description
  //
  // Base class constructor
  //
  G4VMaterialExtension(const G4String& name) :
    fName(name),fHash(G4MaterialExtensionHash{}(name)) {}
  virtual ~G4VMaterialExtension() {;}

  virtual void Print() const = 0;

  // Return the hash value of this extension
  const std::size_t& GetHash() const { return fHash;}
  // Return the extension name
  const G4String& GetName() const { return fName; }
protected:
  // Name of the extension
  const G4String& fName;
  // Hash value of the name.
  // Calculated at initialization time
  const std::size_t fHash;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
