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
// $Id: G4tgrSolidBoolean.hh 68052 2013-03-13 14:38:53Z gcosmo $
//
//
// class G4tgrSolidBoolean
//
// Class description:
//
// A G4tgrSolidBoolean is an G4tgrSolid object with a solid made out of
// the Boolean operation of two solids. The type of operation can be:
// Union, Substraction, Intersection.

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#ifndef G4tgrSolidBoolean_h
#define G4tgrSolidBoolean_h

#include "globals.hh"

#include <vector>

#include "G4tgrSolid.hh"

class G4tgrSolidBoolean : public G4tgrSolid
{
  public:  // with description

    G4tgrSolidBoolean(const std::vector<G4String>& wl);
   ~G4tgrSolidBoolean();

    friend std::ostream& operator<<(std::ostream&, const G4tgrSolidBoolean&);

    // Accessors

    inline const G4tgrSolid* GetSolid(G4int ii) const;
    const G4String& GetRelativeRotMatName() const;
    G4ThreeVector GetRelativePlace() const;

  private:   

    // Solid types (Box, Tube, etc) of the solids composing the Boolean solid

    std::vector< std::vector<G4double>* > theSolidParams; 
      // Vectors of parameters. 

    G4String theRelativeRotMatName;
    G4ThreeVector theRelativePlace;
      // Relative placement and rotation of solid 2 w.r.t. solid 1

    std::vector<const G4tgrSolid*> theSolids;
      // The two G4tgrSolid's that combine to make this one

};

inline const G4tgrSolid* G4tgrSolidBoolean::GetSolid( G4int ii ) const
{  
  if((ii != 0) && (ii != 1))
  {
    std::ostringstream message;
    message << "Only two G4tgrSolids (0,1) possible ! Asking for... "
            << ii;
    G4Exception("G4tgrSolidBoolean::GetSolid()", "InvalidInput",
                FatalException, message); 
  }
  return theSolids[ii];
}

#endif
