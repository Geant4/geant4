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
// $Id: G4PolyhedraSide.hh 99392 2016-09-20 13:33:55Z gcosmo $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class header file
//
//
// G4PolyhedraSide
//
// Class description:
//
//   Class implementing a face that represents one segmented side
//   of a polyhedra:
//
//   G4PolyhedraSide( const G4PolyhedraSideRZ *prevRZ,
//                    const G4PolyhedraSideRZ *tail,
//                    const G4PolyhedraSideRZ *head,
//                    const G4PolyhedraSideRZ *nextRZ,
//                          G4int    numSide,
//                          G4double phiStart, G4double phiTotal, 
//                          G4bool phiIsOpen,  G4bool isAllBehind=false )
//
//   Values for r1,z1 and r2,z2 should be specified in clockwise
//   order in (r,z).

// Author: 
//   David C. Williams (davidw@scipp.ucsc.edu)
// --------------------------------------------------------------------

#ifndef G4PolyhedraSide_hh
#define G4PolyhedraSide_hh

#include "G4VCSGface.hh"

class G4IntersectingCone;

struct G4PolyhedraSideRZ
{
  G4double r, z;  // start of vector
};

// ----------------------------------------------------------------------------
// MT-specific utility code 

#include "G4GeomSplitter.hh"

// The class G4PhSideData is introduced to encapsulate the
// fields of the class G4PolyhedraSide that may not be read-only.
//
class G4PhSideData
{
  public:
    void initialize()
    {
      fPhi.first = G4ThreeVector(0,0,0);
      fPhi.second= 0.0;
    }

    std::pair<G4ThreeVector, G4double> fPhi;  // Cached value for phi

};

// The type G4PhSideManager is introduced to encapsulate the methods used
// by both the master thread and worker threads to allocate memory space
// for the fields encapsulated by the class G4PhSideData.
//
typedef G4GeomSplitter<G4PhSideData> G4PhSideManager;

//
// ----------------------------------------------------------------------------

class G4PolyhedraSide : public G4VCSGface
{

  public:  // with description

    G4PolyhedraSide( const G4PolyhedraSideRZ *prevRZ,
                     const G4PolyhedraSideRZ *tail,
                     const G4PolyhedraSideRZ *head,
                     const G4PolyhedraSideRZ *nextRZ,
                           G4int    numSide,
                           G4double phiStart, G4double phiTotal, 
                           G4bool phiIsOpen,  G4bool isAllBehind=false );
    virtual ~G4PolyhedraSide();
  
    G4PolyhedraSide( const G4PolyhedraSide &source );
    G4PolyhedraSide& operator=( const G4PolyhedraSide &source );
  
    G4bool Intersect( const G4ThreeVector &p, const G4ThreeVector &v,  
                            G4bool outgoing, G4double surfTolerance,
                            G4double &distance, G4double &distFromSurface,
                            G4ThreeVector &normal, G4bool &allBehind );

    G4double Distance( const G4ThreeVector &p, G4bool outgoing );
  
    EInside Inside( const G4ThreeVector &p, G4double tolerance, 
                          G4double *bestDistance );
  
    G4ThreeVector Normal( const G4ThreeVector &p,  G4double *bestDistance );

    G4double Extent( const G4ThreeVector axis );
  
    void CalculateExtent( const EAxis axis, 
                          const G4VoxelLimits &voxelLimit,
                          const G4AffineTransform &tranform,
                                G4SolidExtentList &extentList );

    G4VCSGface *Clone() { return new G4PolyhedraSide( *this ); }

  public:  // without description

    // Methods used for GetPointOnSurface()

    G4double SurfaceTriangle( G4ThreeVector p1,
                              G4ThreeVector p2,
                              G4ThreeVector p3,
                              G4ThreeVector *p4 );
    G4ThreeVector GetPointOnPlane( G4ThreeVector p0, G4ThreeVector p1, 
                                   G4ThreeVector p2, G4ThreeVector p3,
                                   G4double *Area );
    G4double SurfaceArea();
    G4ThreeVector GetPointOnFace();  

  public:  // without description

    G4PolyhedraSide(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    inline G4int GetInstanceID() const  { return instanceID; }
      // Returns the instance ID.

    static const G4PhSideManager& GetSubInstanceManager();
      // Returns the private data instance manager.

    //
    // A couple internal data structures
    //
    struct sG4PolyhedraSideVec;         // Secret recipe for allowing
    friend struct sG4PolyhedraSideVec;  // protected nested structures

    typedef struct sG4PolyhedraSideEdge
    {
      G4ThreeVector  normal;       // Unit normal to this edge
      G4ThreeVector  corner[2];    // The two corners of this phi edge
      G4ThreeVector  cornNorm[2];  // The normals of these corners
    } G4PolyhedraSideEdge;
  
    typedef struct sG4PolyhedraSideVec
    {
      G4ThreeVector  normal,   // Normal (point out of the shape)
                     center,   // Point in center of side
                     surfPhi,  // Unit vector on surface pointing along phi
                     surfRZ;   // Unit vector on surface pointing along R/Z
      G4PolyhedraSideEdge *edges[2];  // The phi boundary edges to this side 
                                      //     [0]=low phi [1]=high phi
      G4ThreeVector  edgeNorm[2];     // RZ edge normals [i] at {r[i],z[i]}
    } G4PolyhedraSideVec;

  protected:

    G4bool IntersectSidePlane( const G4ThreeVector &p, const G4ThreeVector &v,
                               const G4PolyhedraSideVec& vec,
                                     G4double normSign, 
                                     G4double surfTolerance,
                                     G4double &distance,
                                     G4double &distFromSurface );

    G4int LineHitsSegments( const G4ThreeVector &p,
                            const G4ThreeVector &v,
                                  G4int *i1, G4int *i2 );

    G4int ClosestPhiSegment( G4double phi );
  
    G4int PhiSegment( G4double phi );

    G4double GetPhi( const G4ThreeVector& p );

    G4double DistanceToOneSide( const G4ThreeVector &p,
                                const G4PolyhedraSideVec &vec,
                                      G4double *normDist );

    G4double DistanceAway( const G4ThreeVector &p,
                           const G4PolyhedraSideVec &vec,
                                 G4double *normDist );
             
    void CopyStuff( const G4PolyhedraSide &source );

  protected:

    G4int   numSide;      // Number sides
    G4double r[2], z[2];  // r, z parameters, in specified order
    G4double startPhi,    // Start phi (0 to 2pi), if phiIsOpen
             deltaPhi,    // Delta phi (0 to 2pi), if phiIsOpen
             endPhi;      // End phi (>startPhi), if phiIsOpen
    G4bool   phiIsOpen;   // True if there is a phi slice
    G4bool   allBehind;   // True if the entire solid is "behind" this face
  
    G4IntersectingCone  *cone;  // Our intersecting cone
  
    G4PolyhedraSideVec  *vecs;    // Vector set for each facet of our face
    G4PolyhedraSideEdge *edges;   // The edges belong to vecs
    G4double    lenRZ,      // RZ length of each side
                lenPhi[2];  // Phi dimensions of each side
    G4double    edgeNorm;   // Normal in RZ/Phi space to each side

  private:

    G4double kCarTolerance;  // Geometrical surface thickness
    G4double fSurfaceArea;   // Surface Area 

    G4int instanceID;
      // This field is used as instance ID.
    G4GEOM_DLL static G4PhSideManager subInstanceManager;
      // This field helps to use the class G4PhSideManager introduced above.
};

#endif
