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
// G4ClippablePolygon
//
// Class description:
//
// Declaration of a utility class of a polygon that can be
// clipped by a voxel.

// Author: David C. Williams (davidw@scipp.ucsc.edu)
// --------------------------------------------------------------------
#ifndef G4CLIPPABLEPOLYGON_HH
#define G4CLIPPABLEPOLYGON_HH

#include <vector>

#include "G4Types.hh"
#include "geomdefs.hh"
#include "G4ThreeVector.hh"

class G4AffineTransform;
class G4VoxelLimits;

class G4ClippablePolygon
{
  typedef std::vector<G4ThreeVector> G4ThreeVectorList;

  public:  // with description

  G4ClippablePolygon();
  virtual ~G4ClippablePolygon();
    // Constructor & virtual destructor.
  
  virtual void AddVertexInOrder( const G4ThreeVector vertex );
  virtual void ClearAllVertices();
  
  inline void SetNormal( const G4ThreeVector& newNormal );
  inline const G4ThreeVector GetNormal() const;
  
  virtual G4bool Clip( const G4VoxelLimits& voxelLimit );

  virtual G4bool PartialClip( const G4VoxelLimits& voxelLimit,
                              const EAxis IgnoreMe );
    // Clip, while ignoring the indicated axis.

  virtual void ClipAlongOneAxis( const G4VoxelLimits& voxelLimit,
                                 const EAxis axis );
    // Clip along just one axis, as specified in voxelLimit.

  virtual G4bool GetExtent( const EAxis axis, 
                                  G4double& min, G4double& max ) const;

  virtual const G4ThreeVector* GetMinPoint( const EAxis axis ) const;
    // Returns pointer to minimum point along the specified axis.
    // Take care! Do not use pointer after destroying parent polygon.

  virtual const G4ThreeVector* GetMaxPoint( const EAxis axis ) const;
    // Returns pointer to maximum point along the specified axis.
    // Take care! Do not use pointer after destroying parent polygon.

  inline std::size_t GetNumVertices() const;
  inline G4bool Empty() const;
  
  virtual G4bool InFrontOf( const G4ClippablePolygon& other, EAxis axis ) const;
    // Decide if the polygon is in "front" of another when
    // viewed along the specified axis. For our purposes here,
    // it is sufficient to use the minimum extent of the
    // polygon along the axis to determine this.

  virtual G4bool BehindOf( const G4ClippablePolygon& other, EAxis axis ) const;
    // Decide if this polygon is behind another.
    // Remarks in method "InFrontOf" are valid here too.

  virtual G4bool GetPlanerExtent( const G4ThreeVector& pointOnPlane, 
                                  const G4ThreeVector& planeNormal,
                                        G4double& min, G4double& max ) const;
    // Get min/max distance in or out of a plane.

  protected:  // with description

  void ClipToSimpleLimits( G4ThreeVectorList& pPolygon,
                           G4ThreeVectorList& outputPolygon,
                     const G4VoxelLimits& pVoxelLimit  );
    // pVoxelLimits must be only limited along one axis, and either
    // the maximum along the axis must be +kInfinity, or the minimum
    // -kInfinity

  protected:

  G4ThreeVectorList vertices;
  G4ThreeVector normal;
  G4double kCarTolerance;
};

#include "G4ClippablePolygon.icc"

#endif
