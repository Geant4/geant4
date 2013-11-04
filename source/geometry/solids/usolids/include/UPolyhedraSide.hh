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
// $Id$
//
// 
// --------------------------------------------------------------------
// GEANT 4 class header file
//
//
// UPolyhedraSide
//
// Class description:
//
//   Class implementing a face that represents one segmented side
//   of a polyhedra:
//
//   UPolyhedraSide( const UPolyhedraSideRZ *prevRZ,
//                    const UPolyhedraSideRZ *tail,
//                    const UPolyhedraSideRZ *head,
//                    const UPolyhedraSideRZ *nextRZ,
//                          int    numSide,
//                          double phiStart, double phiTotal, 
//                          bool phiIsOpen,  bool isAllBehind=false )
//
//   Values for r1,z1 and r2,z2 should be specified in clockwise
//   order in (r,z).

// Author: 
//   David C. Williams (davidw@scipp.ucsc.edu)
// --------------------------------------------------------------------

#ifndef UPolyhedraSide_hh
#define UPolyhedraSide_hh

#include "UVCSGface.hh"

class UIntersectingCone;

struct UPolyhedraSideRZ
{
  double r, z;  // start of vector
};

class UPolyhedraSide : public UVCSGface
{

  public:  // with description

    UPolyhedraSide( const UPolyhedraSideRZ *prevRZ,
                     const UPolyhedraSideRZ *tail,
                     const UPolyhedraSideRZ *head,
                     const UPolyhedraSideRZ *nextRZ,
                           int    numSide,
                           double phiStart, double phiTotal, 
                           bool phiIsOpen,  bool isAllBehind=false );
    virtual ~UPolyhedraSide();
  
    UPolyhedraSide( const UPolyhedraSide &source );
    UPolyhedraSide& operator=( const UPolyhedraSide &source );
  
    bool Distance( const UVector3 &p, const UVector3 &v,  
                            bool outgoing, double surfTolerance,
                            double &distance, double &distFromSurface,
                            UVector3 &normal, bool &allBehind );

    double Safety( const UVector3 &p, bool outgoing );
  
    VUSolid::EnumInside Inside( const UVector3 &p, double tolerance, 
                          double *bestDistance );
  
    UVector3 Normal( const UVector3 &p,  double *bestDistance );

    double Extent( const UVector3 axis );
  
    UVCSGface *Clone() { return new UPolyhedraSide( *this ); }

  public:  // without description

    // Methods used for GetPointOnSurface()

    double SurfaceTriangle( UVector3 p1,
                              UVector3 p2,
                              UVector3 p3,
                              UVector3 *p4 );
    UVector3 GetPointOnPlane( UVector3 p0, UVector3 p1, 
                                   UVector3 p2, UVector3 p3,
                                   double *Area );
    double SurfaceArea();
    UVector3 GetPointOnFace();  

  public:  // without description

    UPolyhedraSide(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

  protected:

    //
    // A couple internal data structures
    //
    struct sUPolyhedraSideVec;         // Secret recipe for allowing
    friend struct sUPolyhedraSideVec;  // protected nested structures

    typedef struct sUPolyhedraSideEdge
    {
      UVector3  normal;       // Unit normal to this edge
      UVector3  corner[2];    // The two corners of this phi edge
      UVector3  cornNorm[2];  // The normals of these corners
    } UPolyhedraSideEdge;
  
    typedef struct sUPolyhedraSideVec
    {
      UVector3  normal,   // Normal (point out of the shape)
                     center,   // Point in center of side
                     surfPhi,  // Unit vector on surface pointing along phi
                     surfRZ;   // Unit vector on surface pointing along R/Z
      UPolyhedraSideEdge *edges[2];  // The phi boundary edges to this side 
                                      //     [0]=low phi [1]=high phi
      UVector3  edgeNorm[2];     // RZ edge normals [i] at {r[i],z[i]}
    } UPolyhedraSideVec;

    bool IntersectSidePlane( const UVector3 &p, const UVector3 &v,
                               const UPolyhedraSideVec& vec,
                                     double normSign, 
                                     double surfTolerance,
                                     double &distance,
                                     double &distFromSurface );

    int LineHitsSegments( const UVector3 &p,
                            const UVector3 &v,
                                  int *i1, int *i2 );

    int ClosestPhiSegment( double phi );
  
    int PhiSegment( double phi );

    double GetPhi( const UVector3& p );

    double DistanceToOneSide( const UVector3 &p,
                                const UPolyhedraSideVec &vec,
                                      double *normDist );

    double DistanceAway( const UVector3 &p,
                           const UPolyhedraSideVec &vec,
                                 double *normDist );
             
    void CopyStuff( const UPolyhedraSide &source );

  protected:

    int   numSide;      // Number sides
    double r[2], z[2];  // r, z parameters, in specified order
    double startPhi,    // Start phi (0 to 2pi), if phiIsOpen
             deltaPhi,    // Delta phi (0 to 2pi), if phiIsOpen
             endPhi;      // End phi (>startPhi), if phiIsOpen
    bool   phiIsOpen;   // True if there is a phi slice
    bool   allBehind;   // True if the entire solid is "behind" this face
  
    UIntersectingCone  *cone;  // Our intersecting cone
  
    UPolyhedraSideVec  *vecs;    // Vector Set for each facet of our face
    UPolyhedraSideEdge *edges;   // The edges belong to vecs
    double    lenRZ,      // RZ length of each side
                lenPhi[2];  // Phi dimensions of each side
    double    edgeNorm;   // Normal in RZ/Phi space to each side

  private:

    std::pair<UVector3, double> fPhi;  // Cached value for phi
    double kCarTolerance;  // Geometrical surface thickness
    double fSurfaceArea;   // Surface Area 
};

#endif
