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
// G4tgrSolidScaled
//
// Class description:
//
// A G4tgrSolidScaled is a G4tgrSolid object with a solid made by 
// scaling another solid(oSolid) by a G4Scale3D vector.
//
// Author: P.Heidary, AEOI - November 2021
// --------------------------------------------------------------------

#ifndef G4tgrSolidScaled_hh
#define G4tgrSolidScaled_hh 1

#include <vector>

#include "globals.hh"
#include "G4tgrSolid.hh"
#include "G4Transform3D.hh"

class G4tgrSolidScaled : public G4tgrSolid
{
  public:

    G4tgrSolidScaled(const std::vector<G4String>& wl);
    ~G4tgrSolidScaled();

    friend std::ostream& operator<<(std::ostream&, const G4tgrSolidScaled&);

    // Accessors
    const G4Scale3D GetScale3d() const { return scale3d; }
    const G4tgrSolid* GetOrigSolid() const { return origSolid; }

  private:

    // Scale vector
    G4Scale3D scale3d;

    // Original solid
    G4tgrSolid* origSolid;
};

#endif
