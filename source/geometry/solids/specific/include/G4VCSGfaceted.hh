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
//
// $Id: G4VCSGfaceted.hh,v 1.7 2002-10-28 11:47:51 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class header file
//
//
// G4VCSGfaceted
//
// Class description:
//
//   Virtual class defining CSG-like type shape that is built entire
//   of G4CSGface faces.

// Author:
//   David C. Williams (davidw@scipp.ucsc.edu)
// --------------------------------------------------------------------

#ifndef G4VCSGfaceted_hh
#define G4VCSGfaceted_hh

#include "G4VSolid.hh"

class G4VCSGface;
class G4VisExtent;

class G4VCSGfaceted : public G4VSolid 
{
  public:  // with description

  G4VCSGfaceted( G4String name );
  virtual ~G4VCSGfaceted();
  
  G4VCSGfaceted( const G4VCSGfaceted &source );
  const G4VCSGfaceted &operator=( const G4VCSGfaceted &source );
  
  virtual G4bool CalculateExtent( const EAxis pAxis,
                                  const G4VoxelLimits& pVoxelLimit,
                                  const G4AffineTransform& pTransform,
                                        G4double& pmin, G4double& pmax ) const;
  
  virtual EInside Inside( const G4ThreeVector& p ) const;

  virtual G4ThreeVector SurfaceNormal( const G4ThreeVector& p ) const;

  virtual G4double DistanceToIn( const G4ThreeVector& p,
                                 const G4ThreeVector& v ) const;
  virtual G4double DistanceToIn( const G4ThreeVector& p ) const;
  virtual G4double DistanceToOut( const G4ThreeVector& p,
                                  const G4ThreeVector& v,
                                  const G4bool calcNorm=false,
                                        G4bool *validNorm=0,
                                        G4ThreeVector *n=0 ) const;
  virtual G4double DistanceToOut( const G4ThreeVector& p ) const;

  virtual G4GeometryType GetEntityType() const;

  virtual G4std::ostream& StreamInfo(G4std::ostream& os) const;

  virtual G4Polyhedron* CreatePolyhedron() const = 0;

  virtual void DescribeYourselfTo( G4VGraphicsScene& scene ) const;

  virtual G4VisExtent GetExtent() const;

  protected:  // without description

  G4int    numFace;
  G4VCSGface **faces;

  virtual G4double DistanceTo( const G4ThreeVector &p,
                               const G4bool outgoing ) const;

  void CopyStuff( const G4VCSGfaceted &source );
  void DeleteStuff();
};

#endif
