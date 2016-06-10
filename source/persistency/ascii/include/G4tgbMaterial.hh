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
// $Id: G4tgbMaterial.hh 66363 2012-12-18 09:12:54Z gcosmo $
//
//
// class G4tgbMaterial
//
// Class description:
//
// Transient class of a material; builds a G4Material.

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#ifndef G4tgbMaterial_h
#define G4tgbMaterial_h

#include "globals.hh"

#include <vector>
#include <string>
#include <iostream>

#include "G4tgrMaterial.hh"
#include "G4Material.hh"

class G4tgbMaterial
{
  friend std::ostream& operator<<(std::ostream&, const G4tgbMaterial&);
  
  public:  // with description

    G4tgbMaterial();
    virtual ~G4tgbMaterial();

    G4tgbMaterial( G4tgrMaterial* tgr );

    virtual G4Material* BuildG4Material() = 0;

    const G4String& GetName() const
    {
      return theTgrMate->GetName();
    }

    G4double GetDensity() const
    {
      return theTgrMate->GetDensity();
    }

    G4int GetNumberOfMaterials() const
    {
      return theTgrMate->GetNumberOfComponents();
    }

    G4double GetA() const
    {
      return theTgrMate->GetA();
    }

    G4double GetZ() const
    {
      return theTgrMate->GetZ();
    }

    const G4String& GetType() const
    {
      return theTgrMate->GetType();
    }  


    G4tgrMaterial* GetTgrMate() const
    {
      return theTgrMate;
    }

  protected:

    G4tgrMaterial* theTgrMate;
    G4Material* theG4Mate;
};

#endif
