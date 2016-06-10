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
// $Id: G4tgrMaterialSimple.hh 66363 2012-12-18 09:12:54Z gcosmo $
//
//
// class G4tgrMaterialSimple
//
// Class description:
//
// Class to represent a simple material, i.e. made of a single element.

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#ifndef G4tgrMaterialSimple_h
#define G4tgrMaterialSimple_h

#include "globals.hh"

#include <vector>
#include <iostream>

#include "G4tgrMaterial.hh"

class G4tgrMaterialSimple : public G4tgrMaterial
{
  public:  // with description
 
    G4tgrMaterialSimple(const G4String& matType,
                        const std::vector<G4String>& wl);
      // Fill the data interpreting the list of words read 'wl'

    friend std::ostream& operator<<(std::ostream&, const G4tgrMaterialSimple&); 

    G4tgrMaterialSimple& operator= (const G4tgrMaterialSimple&); 

    G4double GetA() const { return theA; }
    G4double GetZ() const { return theZ; }

    const G4String& GetComponent(G4int i) const;
    G4double GetFraction(G4int i);

  public:  // without description

    G4tgrMaterialSimple();
   ~G4tgrMaterialSimple();

  private:

    G4String name;
    G4double    theA;
    G4double    theZ;
};

#endif
