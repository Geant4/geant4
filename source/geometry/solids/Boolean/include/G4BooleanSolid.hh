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
// G4BooleanSolid
//
// Class description:
//
// Abstract base class for solids created by boolean operations
// between other solids.

// 10.09.98 V.Grichine - created
// --------------------------------------------------------------------
#ifndef G4BOOLEANSOLID_HH
#define G4BOOLEANSOLID_HH

#include "G4DisplacedSolid.hh"

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"

class HepPolyhedronProcessor;

class G4BooleanSolid : public G4VSolid
{
  public:
 
    G4BooleanSolid( const G4String& pName,
                          G4VSolid* pSolidA ,
                          G4VSolid* pSolidB   );

    G4BooleanSolid( const G4String& pName,
                          G4VSolid* pSolidA ,
                          G4VSolid* pSolidB,
                          G4RotationMatrix* rotMatrix,
                    const G4ThreeVector& transVector    );

    G4BooleanSolid( const G4String& pName,
                          G4VSolid* pSolidA ,
                          G4VSolid* pSolidB , 
                    const G4Transform3D& transform   );

    virtual ~G4BooleanSolid();

    virtual const G4VSolid* GetConstituentSolid(G4int no) const;
    virtual       G4VSolid* GetConstituentSolid(G4int no);
      // If Solid is made up from a Boolean operation of two solids,
      // return the corresponding solid (for no=0 and 1).
      // If the solid is not a "Boolean", return 0.

    virtual G4double GetCubicVolume();
    inline G4double GetSurfaceArea();

    virtual G4GeometryType  GetEntityType() const;
    virtual G4Polyhedron* GetPolyhedron () const;

    std::ostream& StreamInfo(std::ostream& os) const;

    inline G4int GetCubVolStatistics() const;
    inline G4double GetCubVolEpsilon() const;
    inline void SetCubVolStatistics(G4int st);
    inline void SetCubVolEpsilon(G4double ep);
   
    inline G4int GetAreaStatistics() const;
    inline G4double GetAreaAccuracy() const;
    inline void SetAreaStatistics(G4int st);
    inline void SetAreaAccuracy(G4double ep);
   
    G4ThreeVector GetPointOnSurface() const;

    G4BooleanSolid(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4BooleanSolid(const G4BooleanSolid& rhs);
    G4BooleanSolid& operator=(const G4BooleanSolid& rhs);
      // Copy constructor and assignment operator.

  protected:
  
    void GetListOfPrimitives(std::vector<std::pair<G4VSolid *,G4Transform3D>>&,
                             const G4Transform3D&) const;
      // Get list of constituent primitives of the solid and their placements.

    G4Polyhedron* StackPolyhedron(HepPolyhedronProcessor&,
                                  const G4VSolid*) const;
      // Stack polyhedra for processing. Return top polyhedron.

  protected:
  
    G4VSolid* fPtrSolidA = nullptr;
    G4VSolid* fPtrSolidB = nullptr;

    G4double fCubicVolume = -1.0;
      // Stored value of fCubicVolume 

  private:

    G4int    fStatistics = 1000000;
    G4double fCubVolEpsilon = 0.001;
    G4double fAreaAccuracy = -1;
    G4double fSurfaceArea = -1.0;

    mutable G4bool fRebuildPolyhedron = false;
    mutable G4Polyhedron* fpPolyhedron = nullptr;

    mutable std::vector<std::pair<G4VSolid *,G4Transform3D>> fPrimitives;
    mutable G4double fPrimitivesSurfaceArea = 0.0;

    G4bool  createdDisplacedSolid = false;
      // If & only if this object created it, it must delete it
} ;

#include "G4BooleanSolid.icc"

#endif
