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
// G4tgrSolidMultiUnion
//
// Class description:
//
// A G4tgrSolidMultiUnion is a G4tgrSolid object with a solid made by
// union of more than two(nSolid) solids with respective rotations
// (theRotMat) and positions(thePosition).
//
// Author: P.Heidary, AEOI - November 2021
// --------------------------------------------------------------------
#ifndef G4tgrSolidMultiUnion_hh
#define G4tgrSolidMultiUnion_hh 1

#include <vector>

#include "globals.hh"
#include "G4tgrSolid.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"

class G4tgrSolidMultiUnion : public G4tgrSolid
{
  public:

    G4tgrSolidMultiUnion(const std::vector<G4String>& wl);
    ~G4tgrSolidMultiUnion();

    friend std::ostream& operator<<(std::ostream&, const G4tgrSolidMultiUnion&);

    // Accessors
    inline const G4tgrSolid* GetSolid(G4int ii) const;
    inline const G4Transform3D GetTransformation(G4int ii) const;
    inline G4int GetNSolid() const { return nSolid; };

  private:

    // Placement and rotation of current solid
    G4String          theRotMatName;
    G4ThreeVector     thePosition;
    G4RotationMatrix* theRotMat;
    G4Transform3D     tr1;
    G4int             nSolid;

    // Vectors of parameters.
    std::vector<std::vector<G4double>*> theSolidParams;

    // Array to store Solids and Transformations.
    std::vector<G4Transform3D>     theTransformations;
    std::vector<const G4tgrSolid*> theSolids;
};

inline const G4tgrSolid* G4tgrSolidMultiUnion::GetSolid(G4int ii) const
{
  if(ii > nSolid)
  {
    std::ostringstream message;
    message << "Only " << nSolid + 1 << " G4tgrSolids are available! " << 
    " Asking for... " << ii + 1;
    G4Exception("G4tgrSolidMultiUnion::GetSolid()", "InvalidInput", FatalException,
                message);
  }
  return theSolids[ii];
}

inline const G4Transform3D G4tgrSolidMultiUnion::GetTransformation(G4int ii) const
{
  if(ii > nSolid)
  {
    std::ostringstream message;
    message << "Only " << nSolid + 1 << " G4tgrSolids are available! " << 
    " Asking for... " << ii + 1;
    G4Exception("G4tgrSolidMultiUnion::GetSolid()", "InvalidInput", FatalException,
                message);
  }
  return theTransformations[ii];
}

#endif
