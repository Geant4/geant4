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
// $Id: G4Trap.cc 101121 2016-11-07 09:18:01Z gcosmo $
//
// class G4Trap
//
// Implementation for G4Trap class
//
// History:
//
// 23.09.16 E.Tcherniaev: added Extent(pmin,pmax),
//                      use G4BoundingEnvelope for CalculateExtent(),
//                      removed CreateRotatedVertices()
// 28.04.05 V.Grichine: new SurfaceNormal according to J. Apostolakis proposal 
// 26.04.05 V.Grichine: new SurfaceNormal is default 
// 19.04.05 V.Grichine: bug fixed in G4Trap("name",G4ThreeVector[8] vp)
// 12.12.04 V.Grichine: SurfaceNormal with edges/vertices 
// 15.11.04 V.Grichine: bug fixed in G4Trap("name",G4ThreeVector[8] vp)
// 13.12.99 V.Grichine: bug fixed in DistanceToIn(p,v)
// 19.11.99 V.Grichine: kUndef was added to Eside enum
// 04.06.99 S.Giani: Fixed CalculateExtent in rotated case. 
// 08.12.97 J.Allison: Added "nominal" constructor and method SetAllParameters.
// 01.11.96 V.Grichine: Costructor for Right Angular Wedge from STEP, G4Trd/Para
// 09.09.96 V.Grichine: Final modifications before to commit
// 21.03.95 P.Kent: Modified for `tolerant' geometry
//
//////////////////////////////////////////////////////////////////////////////////// 

#include "G4Trap.hh"

#if !defined(G4GEOM_USE_UTRAP)

#include "globals.hh"

#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"
#include "G4BoundingEnvelope.hh"

#include "G4VPVParameterisation.hh"

#include "Randomize.hh"

#include "G4VGraphicsScene.hh"
#include "G4Polyhedron.hh"

using namespace CLHEP;

////////////////////////////////////////////////////////////////////////
//
// Accuracy of coplanarity

const G4double kCoplanar_Tolerance = 1E-4 ;

//////////////////////////////////////////////////////////////////////////
//
// Private enum: Not for external use 
    
enum Eside {kUndef,ks0,ks1,ks2,ks3,kPZ,kMZ};

//////////////////////////////////////////////////////////////////////////
//
// Constructor - check and set half-widths as well as angles: 
// final check of coplanarity

G4Trap::G4Trap( const G4String& pName,
                      G4double pDz,
                      G4double pTheta, G4double pPhi,
                      G4double pDy1, G4double pDx1, G4double pDx2,
                      G4double pAlp1,
                      G4double pDy2, G4double pDx3, G4double pDx4,
                      G4double pAlp2)
  : G4CSGSolid(pName)
{
  if ( pDz <= 0 || pDy1 <= 0 || pDx1 <= 0 ||
       pDx2 <= 0 || pDy2 <= 0 || pDx3 <= 0 || pDx4 <= 0 )
  {
    std::ostringstream message;
    message << "Invalid length parameters for Solid: " << GetName() << G4endl
            << "        X - "
            << pDx1 << ", " << pDx2 << ", " << pDx3 << ", " << pDx4 << G4endl
            << "          Y - " << pDy1 << ", " << pDy2 << G4endl
            << "          Z - " << pDz;
    G4Exception("G4Trap::G4Trap()", "GeomSolids0002",
                FatalException, message);
  }

  fDz=pDz;
  fTthetaCphi=std::tan(pTheta)*std::cos(pPhi);
  fTthetaSphi=std::tan(pTheta)*std::sin(pPhi);
      
  fDy1=pDy1;
  fDx1=pDx1;
  fDx2=pDx2;
  fTalpha1=std::tan(pAlp1);
     
  fDy2=pDy2;
  fDx3=pDx3;
  fDx4=pDx4;
  fTalpha2=std::tan(pAlp2);

  MakePlanes();
}

////////////////////////////////////////////////////////////////////////////
//
// Constructor - Design of trapezoid based on 8 G4ThreeVector parameters, 
// which are its vertices. Checking of planarity with preparation of 
// fPlanes[] and than calculation of other members

G4Trap::G4Trap( const G4String& pName,
                const G4ThreeVector pt[8] )
  : G4CSGSolid(pName)
{
  G4bool good;

  // Start with check of centering - the center of gravity trap line
  // should cross the origin of frame

  if (!(   pt[0].z() < 0 
        && pt[0].z() == pt[1].z() && pt[0].z() == pt[2].z()
        && pt[0].z() == pt[3].z()
        && pt[4].z() > 0 
        && pt[4].z() == pt[5].z() && pt[4].z() == pt[6].z()
        && pt[4].z() == pt[7].z()
        && std::fabs( pt[0].z() + pt[4].z() ) < kCarTolerance
        && pt[0].y() == pt[1].y() && pt[2].y() == pt[3].y()
        && pt[4].y() == pt[5].y() && pt[6].y() == pt[7].y()
        && std::fabs( pt[0].y() + pt[2].y() + pt[4].y() + pt[6].y() ) < kCarTolerance 
        && std::fabs( pt[0].x() + pt[1].x() + pt[4].x() + pt[5].x() + 
           pt[2].x() + pt[3].x() + pt[6].x() + pt[7].x() ) < kCarTolerance ) )
  {
    std::ostringstream message;
    message << "Invalid vertice coordinates for Solid: " << GetName();
    G4Exception("G4Trap::G4Trap()", "GeomSolids0002",
                FatalException, message);
  }
    
  // Bottom side with normal approx. -Y
    
  good = MakePlane(pt[0],pt[4],pt[5],pt[1],fPlanes[0]);

  if (!good)
  {
    DumpInfo();
    G4Exception("G4Trap::G4Trap()", "GeomSolids0002", FatalException,
                "Face at ~-Y not planar.");
  }

  // Top side with normal approx. +Y
    
  good = MakePlane(pt[2],pt[3],pt[7],pt[6],fPlanes[1]);

  if (!good)
  {
    std::ostringstream message;
    message << "Face at ~+Y not planar for Solid: " << GetName();
    G4Exception("G4Trap::G4Trap()", "GeomSolids0002",
                FatalException, message);
  }

  // Front side with normal approx. -X
    
  good = MakePlane(pt[0],pt[2],pt[6],pt[4],fPlanes[2]);

  if (!good)
  {
    std::ostringstream message;
    message << "Face at ~-X not planar for Solid: " << GetName();
    G4Exception("G4Trap::G4Trap()", "GeomSolids0002",
                FatalException, message);
  }

  // Back side iwth normal approx. +X
    
  good = MakePlane(pt[1],pt[5],pt[7],pt[3],fPlanes[3]);
  if (!good)
  {
    std::ostringstream message;
    message << "Face at ~+X not planar for Solid: " << GetName();
    G4Exception("G4Trap::G4Trap()", "GeomSolids0002",
                FatalException, message);
  }
  fDz = (pt[7]).z() ;
      
  fDy1     = ((pt[2]).y()-(pt[1]).y())*0.5;
  fDx1     = ((pt[1]).x()-(pt[0]).x())*0.5;
  fDx2     = ((pt[3]).x()-(pt[2]).x())*0.5;
  fTalpha1 = ((pt[2]).x()+(pt[3]).x()-(pt[1]).x()-(pt[0]).x())*0.25/fDy1;

  fDy2     = ((pt[6]).y()-(pt[5]).y())*0.5;
  fDx3     = ((pt[5]).x()-(pt[4]).x())*0.5;
  fDx4     = ((pt[7]).x()-(pt[6]).x())*0.5;
  fTalpha2 = ((pt[6]).x()+(pt[7]).x()-(pt[5]).x()-(pt[4]).x())*0.25/fDy2;

  fTthetaCphi = ((pt[4]).x()+fDy2*fTalpha2+fDx3)/fDz;
  fTthetaSphi = ((pt[4]).y()+fDy2)/fDz;
}

//////////////////////////////////////////////////////////////////////////////
//
// Constructor for Right Angular Wedge from STEP

G4Trap::G4Trap( const G4String& pName,
                      G4double pZ,
                      G4double pY,
                      G4double pX, G4double pLTX )
  : G4CSGSolid(pName)
{
  G4bool good;

  if ( pZ<=0 || pY<=0 || pX<=0 || pLTX<=0 || pLTX>pX )
  {
    std::ostringstream message;
    message << "Invalid length parameters for Solid: " << GetName();
    G4Exception("G4Trap::G4Trap()", "GeomSolids0002",
                FatalException, message);
  }

  fDz = 0.5*pZ ;
  fTthetaCphi = 0 ;
  fTthetaSphi = 0 ;

  fDy1 = 0.5*pY;
  fDx1 = 0.5*pX ;
  fDx2 = 0.5*pLTX;
  fTalpha1 =  0.5*(pLTX - pX)/pY;

  fDy2 = fDy1 ;
  fDx3 = fDx1;
  fDx4 = fDx2 ;
  fTalpha2 = fTalpha1 ;

  G4ThreeVector pt[8] ;

  pt[0]=G4ThreeVector(-fDz*fTthetaCphi-fDy1*fTalpha1-fDx1,
                      -fDz*fTthetaSphi-fDy1,-fDz);
  pt[1]=G4ThreeVector(-fDz*fTthetaCphi-fDy1*fTalpha1+fDx1,
                      -fDz*fTthetaSphi-fDy1,-fDz);
  pt[2]=G4ThreeVector(-fDz*fTthetaCphi+fDy1*fTalpha1-fDx2,
                      -fDz*fTthetaSphi+fDy1,-fDz);
  pt[3]=G4ThreeVector(-fDz*fTthetaCphi+fDy1*fTalpha1+fDx2,
                      -fDz*fTthetaSphi+fDy1,-fDz);
  pt[4]=G4ThreeVector(+fDz*fTthetaCphi-fDy2*fTalpha2-fDx3,
                      +fDz*fTthetaSphi-fDy2,+fDz);
  pt[5]=G4ThreeVector(+fDz*fTthetaCphi-fDy2*fTalpha2+fDx3,
                      +fDz*fTthetaSphi-fDy2,+fDz);
  pt[6]=G4ThreeVector(+fDz*fTthetaCphi+fDy2*fTalpha2-fDx4,
                      +fDz*fTthetaSphi+fDy2,+fDz);
  pt[7]=G4ThreeVector(+fDz*fTthetaCphi+fDy2*fTalpha2+fDx4,
                      +fDz*fTthetaSphi+fDy2,+fDz);

  // Bottom side with normal approx. -Y
  //
  good=MakePlane(pt[0],pt[4],pt[5],pt[1],fPlanes[0]);
  if (!good)
  {
    std::ostringstream message;
    message << "Face at ~-Y not planar for Solid: " << GetName();
    G4Exception("G4Trap::G4Trap()", "GeomSolids0002",
                FatalException, message);
  }

  // Top side with normal approx. +Y
  //
  good=MakePlane(pt[2],pt[3],pt[7],pt[6],fPlanes[1]);
  if (!good)
  {
    std::ostringstream message;
    message << "Face at ~+Y not planar for Solid: " << GetName();
    G4Exception("G4Trap::G4Trap()", "GeomSolids0002",
                FatalException, message);
  }

  // Front side with normal approx. -X
  //
  good=MakePlane(pt[0],pt[2],pt[6],pt[4],fPlanes[2]);
  if (!good)
  {
    std::ostringstream message;
    message << "Face at ~-X not planar for Solid: " << GetName();
    G4Exception("G4Trap::G4Trap()", "GeomSolids0002",
                FatalException, message);
  }

  // Back side iwth normal approx. +X
  //
  good=MakePlane(pt[1],pt[5],pt[7],pt[3],fPlanes[3]);
  if (!good)
  {
    std::ostringstream message;
    message << "Face at ~+X not planar for Solid: " << GetName();
    G4Exception("G4Trap::G4Trap()", "GeomSolids0002",
                FatalException, message);
  }
}

///////////////////////////////////////////////////////////////////////////////
//
// Constructor for G4Trd

G4Trap::G4Trap( const G4String& pName,
                      G4double pDx1,  G4double pDx2,
                      G4double pDy1,  G4double pDy2,
                      G4double pDz )
  : G4CSGSolid(pName)
{
  G4bool good;

  if ( pDz<=0 || pDy1<=0 || pDx1<=0 || pDx2<=0 || pDy2<=0 )
  {
    std::ostringstream message;
    message << "Invalid length parameters for Solid: " << GetName();
    G4Exception("G4Trap::G4Trap()", "GeomSolids0002",
                FatalException, message);
  }

  fDz = pDz;
  fTthetaCphi = 0 ;
  fTthetaSphi = 0 ;
      
  fDy1 = pDy1 ;
  fDx1 = pDx1 ;
  fDx2 = pDx1 ;
  fTalpha1 = 0 ;
     
  fDy2 = pDy2 ;
  fDx3 = pDx2 ;
  fDx4 = pDx2 ;
  fTalpha2 = 0 ;

  G4ThreeVector pt[8] ;
     
  pt[0]=G4ThreeVector(-fDz*fTthetaCphi-fDy1*fTalpha1-fDx1,
                      -fDz*fTthetaSphi-fDy1,-fDz);
  pt[1]=G4ThreeVector(-fDz*fTthetaCphi-fDy1*fTalpha1+fDx1,
                      -fDz*fTthetaSphi-fDy1,-fDz);
  pt[2]=G4ThreeVector(-fDz*fTthetaCphi+fDy1*fTalpha1-fDx2,
                      -fDz*fTthetaSphi+fDy1,-fDz);
  pt[3]=G4ThreeVector(-fDz*fTthetaCphi+fDy1*fTalpha1+fDx2,
                      -fDz*fTthetaSphi+fDy1,-fDz);
  pt[4]=G4ThreeVector(+fDz*fTthetaCphi-fDy2*fTalpha2-fDx3,
                      +fDz*fTthetaSphi-fDy2,+fDz);
  pt[5]=G4ThreeVector(+fDz*fTthetaCphi-fDy2*fTalpha2+fDx3,
                      +fDz*fTthetaSphi-fDy2,+fDz);
  pt[6]=G4ThreeVector(+fDz*fTthetaCphi+fDy2*fTalpha2-fDx4,
                      +fDz*fTthetaSphi+fDy2,+fDz);
  pt[7]=G4ThreeVector(+fDz*fTthetaCphi+fDy2*fTalpha2+fDx4,
                      +fDz*fTthetaSphi+fDy2,+fDz);

  // Bottom side with normal approx. -Y
  //
  good=MakePlane(pt[0],pt[4],pt[5],pt[1],fPlanes[0]);
  if (!good)
  {
    std::ostringstream message;
    message << "Face at ~-Y not planar for Solid: " << GetName();
    G4Exception("G4Trap::G4Trap()", "GeomSolids0002",
                FatalException, message);
  }

  // Top side with normal approx. +Y
  //
  good=MakePlane(pt[2],pt[3],pt[7],pt[6],fPlanes[1]);
  if (!good)
  {
    std::ostringstream message;
    message << "Face at ~+Y not planar for Solid: " << GetName();
    G4Exception("G4Trap::G4Trap()", "GeomSolids0002",
                FatalException, message);
  }

  // Front side with normal approx. -X
  //
  good=MakePlane(pt[0],pt[2],pt[6],pt[4],fPlanes[2]);
  if (!good)
  {
    std::ostringstream message;
    message << "Face at ~-X not planar for Solid: " << GetName();
    G4Exception("G4Trap::G4Trap()", "GeomSolids0002",
                FatalException, message);
  }

  // Back side iwth normal approx. +X
  //
  good=MakePlane(pt[1],pt[5],pt[7],pt[3],fPlanes[3]);
  if (!good)
  {
    std::ostringstream message;
    message << "Face at ~+X not planar for Solid: " << GetName();
    G4Exception("G4Trap::G4Trap()", "GeomSolids0002",
                FatalException, message);
  }
}

////////////////////////////////////////////////////////////////////////////
//
// Constructor for G4Para

G4Trap::G4Trap( const G4String& pName,
                      G4double pDx, G4double pDy,
                      G4double pDz,
                      G4double pAlpha,
                      G4double pTheta, G4double pPhi )
  : G4CSGSolid(pName)       
{
  G4bool good;

  if ( pDz<=0 || pDy<=0 || pDx<=0 )
  {
    std::ostringstream message;
    message << "Invalid length parameters for Solid: " << GetName();
    G4Exception("G4Trap::G4Trap()", "GeomSolids0002",
                FatalException, message);
  }

  fDz = pDz ;
  fTthetaCphi = std::tan(pTheta)*std::cos(pPhi) ;
  fTthetaSphi = std::tan(pTheta)*std::sin(pPhi) ;
     
  fDy1 = pDy ;
  fDx1 = pDx ;
  fDx2 = pDx ;
  fTalpha1 = std::tan(pAlpha) ;
    
  fDy2 = pDy ;
  fDx3 = pDx ;
  fDx4 = pDx ;
  fTalpha2 = fTalpha1 ;

  G4ThreeVector pt[8] ;
     
  pt[0]=G4ThreeVector(-fDz*fTthetaCphi-fDy1*fTalpha1-fDx1,
                      -fDz*fTthetaSphi-fDy1,-fDz);
  pt[1]=G4ThreeVector(-fDz*fTthetaCphi-fDy1*fTalpha1+fDx1,
                      -fDz*fTthetaSphi-fDy1,-fDz);
  pt[2]=G4ThreeVector(-fDz*fTthetaCphi+fDy1*fTalpha1-fDx2,
                      -fDz*fTthetaSphi+fDy1,-fDz);
  pt[3]=G4ThreeVector(-fDz*fTthetaCphi+fDy1*fTalpha1+fDx2,
                      -fDz*fTthetaSphi+fDy1,-fDz);
  pt[4]=G4ThreeVector(+fDz*fTthetaCphi-fDy2*fTalpha2-fDx3,
                      +fDz*fTthetaSphi-fDy2,+fDz);
  pt[5]=G4ThreeVector(+fDz*fTthetaCphi-fDy2*fTalpha2+fDx3,
                      +fDz*fTthetaSphi-fDy2,+fDz);
  pt[6]=G4ThreeVector(+fDz*fTthetaCphi+fDy2*fTalpha2-fDx4,
                      +fDz*fTthetaSphi+fDy2,+fDz);
  pt[7]=G4ThreeVector(+fDz*fTthetaCphi+fDy2*fTalpha2+fDx4,
                      +fDz*fTthetaSphi+fDy2,+fDz);

  // Bottom side with normal approx. -Y
  //
  good=MakePlane(pt[0],pt[4],pt[5],pt[1],fPlanes[0]);
  if (!good)
  {
    std::ostringstream message;
    message << "Face at ~-Y not planar for Solid: " << GetName();
    G4Exception("G4Trap::G4Trap()", "GeomSolids0002",
                FatalException, message);
  }

  // Top side with normal approx. +Y
  //
  good=MakePlane(pt[2],pt[3],pt[7],pt[6],fPlanes[1]);
  if (!good)
  {
    std::ostringstream message;
    message << "Face at ~+Y not planar for Solid: " << GetName();
    G4Exception("G4Trap::G4Trap()", "GeomSolids0002",
                FatalException, message);
  }

  // Front side with normal approx. -X
  //
  good=MakePlane(pt[0],pt[2],pt[6],pt[4],fPlanes[2]);
  if (!good)
  {
    std::ostringstream message;
    message << "Face at ~-X not planar for Solid: " << GetName();
    G4Exception("G4Trap::G4Trap()", "GeomSolids0002",
                FatalException, message);
  }

  // Back side iwth normal approx. +X
  //
  good=MakePlane(pt[1],pt[5],pt[7],pt[3],fPlanes[3]);
  if (!good)
  {
    std::ostringstream message;
    message << "Face at ~+X not planar for Solid: " << GetName();
    G4Exception("G4Trap::G4Trap()", "GeomSolids0002",
                FatalException, message);
  }
}

///////////////////////////////////////////////////////////////////////////
//
// Nominal constructor for G4Trap whose parameters are to be set by
// a G4VParamaterisation later.  Check and set half-widths as well as
// angles: final check of coplanarity

G4Trap::G4Trap( const G4String& pName )
  : G4CSGSolid (pName), fDz(1.), fTthetaCphi(0.), fTthetaSphi(0.),
    fDy1(1.), fDx1(1.), fDx2(1.), fTalpha1(0.),
    fDy2(1.), fDx3(1.), fDx4(1.), fTalpha2(0.)
{
  MakePlanes();
}

///////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4Trap::G4Trap( __void__& a )
  : G4CSGSolid(a), fDz(1.), fTthetaCphi(0.), fTthetaSphi(0.),
    fDy1(1.), fDx1(1.), fDx2(1.), fTalpha1(0.),
    fDy2(1.), fDx3(1.), fDx4(1.), fTalpha2(0.)
{
  MakePlanes();
}

////////////////////////////////////////////////////////////////////////
//
// Destructor

G4Trap::~G4Trap()
{
}

//////////////////////////////////////////////////////////////////////////
//
// Copy constructor

G4Trap::G4Trap(const G4Trap& rhs)
  : G4CSGSolid(rhs), fDz(rhs.fDz),
    fTthetaCphi(rhs.fTthetaCphi), fTthetaSphi(rhs.fTthetaSphi),
    fDy1(rhs.fDy1), fDx1(rhs.fDx1), fDx2(rhs.fDx2), fTalpha1(rhs.fTalpha1),
    fDy2(rhs.fDy2), fDx3(rhs.fDx3), fDx4(rhs.fDx4), fTalpha2(rhs.fTalpha2)
{
  for (size_t i=0; i<4; ++i)
  {
    fPlanes[i].a = rhs.fPlanes[i].a;
    fPlanes[i].b = rhs.fPlanes[i].b;
    fPlanes[i].c = rhs.fPlanes[i].c;
    fPlanes[i].d = rhs.fPlanes[i].d;
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator

G4Trap& G4Trap::operator = (const G4Trap& rhs) 
{
  // Check assignment to self
  //
  if (this == &rhs)  { return *this; }

  // Copy base class data
  //
  G4CSGSolid::operator=(rhs);

  // Copy data
  //
  fDz = rhs.fDz;
  fTthetaCphi = rhs.fTthetaCphi; fTthetaSphi = rhs.fTthetaSphi;
  fDy1 = rhs.fDy1; fDx1 = rhs.fDx1; fDx2 = rhs.fDx2; fTalpha1 = rhs.fTalpha1;
  fDy2 = rhs.fDy2; fDx3 = rhs.fDx3; fDx4 = rhs.fDx4; fTalpha2 = rhs.fTalpha2;
  for (size_t i=0; i<4; ++i)
  {
    fPlanes[i].a = rhs.fPlanes[i].a;
    fPlanes[i].b = rhs.fPlanes[i].b;
    fPlanes[i].c = rhs.fPlanes[i].c;
    fPlanes[i].d = rhs.fPlanes[i].d;
  }

  return *this;
}

///////////////////////////////////////////////////////////////////////
//
// Set all parameters, as for constructor - check and set half-widths
// as well as angles: final check of coplanarity

void G4Trap::SetAllParameters ( G4double pDz,
                                G4double pTheta,
                                G4double pPhi,
                                G4double pDy1,
                                G4double pDx1,
                                G4double pDx2,
                                G4double pAlp1,
                                G4double pDy2,
                                G4double pDx3,
                                G4double pDx4,
                                G4double pAlp2 )
{
  if ( pDz<=0 || pDy1<=0 || pDx1<=0 || pDx2<=0 || pDy2<=0 || pDx3<=0 || pDx4<=0 )
  {
    std::ostringstream message;
    message << "Invalid Length Parameters for Solid: " << GetName() << G4endl
            << "        X - "
            << pDx1 << ", " << pDx2 << ", " << pDx3 << ", " << pDx4 << G4endl
            << "          Y - " << pDy1 << ", " << pDy2 << G4endl
            << "          Z - " << pDz;
    G4Exception("G4Trap::SetAllParameters()", "GeomSolids0002",
                FatalException, message);
  }
  fCubicVolume= 0.;
  fSurfaceArea= 0.;
  fRebuildPolyhedron = true;
  fDz=pDz;
  fTthetaCphi=std::tan(pTheta)*std::cos(pPhi);
  fTthetaSphi=std::tan(pTheta)*std::sin(pPhi);
     
  fDy1=pDy1;
  fDx1=pDx1;
  fDx2=pDx2;
  fTalpha1=std::tan(pAlp1);
    
  fDy2=pDy2;
  fDx3=pDx3;
  fDx4=pDx4;
  fTalpha2=std::tan(pAlp2);

  MakePlanes();
}

//////////////////////////////////////////////////////////////////////////
//
// Checking of coplanarity

G4bool G4Trap::MakePlanes()
{
  G4bool good = true;

  G4ThreeVector pt[8] ;
     
  pt[0]=G4ThreeVector(-fDz*fTthetaCphi-fDy1*fTalpha1-fDx1,
                      -fDz*fTthetaSphi-fDy1,-fDz);
  pt[1]=G4ThreeVector(-fDz*fTthetaCphi-fDy1*fTalpha1+fDx1,
                      -fDz*fTthetaSphi-fDy1,-fDz);
  pt[2]=G4ThreeVector(-fDz*fTthetaCphi+fDy1*fTalpha1-fDx2,
                      -fDz*fTthetaSphi+fDy1,-fDz);
  pt[3]=G4ThreeVector(-fDz*fTthetaCphi+fDy1*fTalpha1+fDx2,
                      -fDz*fTthetaSphi+fDy1,-fDz);
  pt[4]=G4ThreeVector(+fDz*fTthetaCphi-fDy2*fTalpha2-fDx3,
                      +fDz*fTthetaSphi-fDy2,+fDz);
  pt[5]=G4ThreeVector(+fDz*fTthetaCphi-fDy2*fTalpha2+fDx3,
                      +fDz*fTthetaSphi-fDy2,+fDz);
  pt[6]=G4ThreeVector(+fDz*fTthetaCphi+fDy2*fTalpha2-fDx4,
                      +fDz*fTthetaSphi+fDy2,+fDz);
  pt[7]=G4ThreeVector(+fDz*fTthetaCphi+fDy2*fTalpha2+fDx4,
                      +fDz*fTthetaSphi+fDy2,+fDz);

  // Bottom side with normal approx. -Y
  //
  good=MakePlane(pt[0],pt[4],pt[5],pt[1],fPlanes[0]) ;
  if (!good)
  {
    std::ostringstream message;
    message << "Face at ~-Y not planar for Solid: " << GetName();
    G4Exception("G4Trap::MakePlanes()", "GeomSolids0002",
                FatalException, message);
  }

  // Top side with normal approx. +Y
  //
  good=MakePlane(pt[2],pt[3],pt[7],pt[6],fPlanes[1]);
  if (!good)
  {
    std::ostringstream message;
    message << "Face at ~+Y not planar for Solid: " << GetName();
    G4Exception("G4Trap::MakePlanes()", "GeomSolids0002",
                FatalException, message);
  }

  // Front side with normal approx. -X
  //
  good=MakePlane(pt[0],pt[2],pt[6],pt[4],fPlanes[2]);
  if (!good)
  {
    std::ostringstream message;
    message << "Face at ~-X not planar for Solid: " << GetName();
    G4Exception("G4Trap::MakePlanes()", "GeomSolids0002",
                FatalException, message);
  }
   
  // Back side iwth normal approx. +X
  //
  good = MakePlane(pt[1],pt[5],pt[7],pt[3],fPlanes[3]);
  if ( !good )
  {
    std::ostringstream message;
    message << "Face at ~+X not planar for Solid: " << GetName();
    G4Exception("G4Trap::MakePlanes()", "GeomSolids0002",
                FatalException, message);
  }

  return good;
}

//////////////////////////////////////////////////////////////////////////////
//
// Calculate the coef's of the plane p1->p2->p3->p4->p1
// where the ThreeVectors 1-4 are in anti-clockwise order when viewed from
// infront of the plane (i.e. from normal direction).
//
// Return true if the ThreeVectors are coplanar + set coef;s
//        false if ThreeVectors are not coplanar

G4bool G4Trap::MakePlane( const G4ThreeVector& p1,
                          const G4ThreeVector& p2,
                          const G4ThreeVector& p3,
                          const G4ThreeVector& p4,
                                TrapSidePlane& plane )
{
  G4double a, b, c, sd;
  G4ThreeVector v12, v13, v14, Vcross;

  G4bool good;

  v12    = p2 - p1;
  v13    = p3 - p1;
  v14    = p4 - p1;
  Vcross = v12.cross(v13);

  if (std::fabs(Vcross.dot(v14)/(Vcross.mag()*v14.mag())) > kCoplanar_Tolerance)
  {
    good = false;
  }
  else
  {
    // a,b,c correspond to the x/y/z components of the
    // normal vector to the plane
     
    // a  = (p2.y()-p1.y())*(p1.z()+p2.z())+(p3.y()-p2.y())*(p2.z()+p3.z());
    // a += (p4.y()-p3.y())*(p3.z()+p4.z())+(p1.y()-p4.y())*(p4.z()+p1.z()); // ?   
    // b  = (p2.z()-p1.z())*(p1.x()+p2.x())+(p3.z()-p2.z())*(p2.x()+p3.x());
    // b += (p4.z()-p3.z())*(p3.x()+p4.x())+(p1.z()-p4.z())*(p4.x()+p1.x()); // ?      
    // c  = (p2.x()-p1.x())*(p1.y()+p2.y())+(p3.x()-p2.x())*(p2.y()+p3.y());
    // c += (p4.x()-p3.x())*(p3.y()+p4.y())+(p1.x()-p4.x())*(p4.y()+p1.y()); // ?

    // Let create diagonals 4-2 and 3-1 than (4-2)x(3-1) provides
    // vector perpendicular to the plane directed to outside !!!
    // and a,b,c, = f(1,2,3,4) external relative to trap normal

    a = +(p4.y() - p2.y())*(p3.z() - p1.z())
        - (p3.y() - p1.y())*(p4.z() - p2.z());

    b = -(p4.x() - p2.x())*(p3.z() - p1.z())
        + (p3.x() - p1.x())*(p4.z() - p2.z());
 
    c = +(p4.x() - p2.x())*(p3.y() - p1.y())
        - (p3.x() - p1.x())*(p4.y() - p2.y());

    sd = std::sqrt( a*a + b*b + c*c ); // so now vector plane.(a,b,c) is unit 

    if( sd > 0 )
    {
      plane.a = a/sd;
      plane.b = b/sd;
      plane.c = c/sd;
    }
    else
    {
      std::ostringstream message;
      message << "Invalid parameters: norm.mod() <= 0, for Solid: "
              << GetName();
      G4Exception("G4Trap::MakePlanes()", "GeomSolids0002",
                  FatalException, message) ;
    }
    // Calculate D: p1 in in plane so D=-n.p1.Vect()
    
    plane.d = -( plane.a*p1.x() + plane.b*p1.y() + plane.c*p1.z() );

    good = true;
  }
  return good;
}

//////////////////////////////////////////////////////////////////////////////
//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.

void G4Trap::ComputeDimensions(       G4VPVParameterisation* p,
                                const G4int n,
                                const G4VPhysicalVolume* pRep )
{
  p->ComputeDimensions(*this,n,pRep);
}

////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4Trap::Extent(G4ThreeVector& pMin, G4ThreeVector& pMax) const
{
  G4double dz  = GetZHalfLength();
  G4double dx1 = GetXHalfLength1();
  G4double dx2 = GetXHalfLength2();
  G4double dx3 = GetXHalfLength3();
  G4double dx4 = GetXHalfLength4();
  G4double dy1 = GetYHalfLength1();
  G4double dy2 = GetYHalfLength2();

  G4double x0 = dz*fTthetaCphi;
  G4double x1 = dy1*GetTanAlpha1();
  G4double x2 = dy2*GetTanAlpha2();
  G4double xmin =
    std::min(
    std::min(
    std::min(-x0-x1-dx1,-x0+x1-dx2),x0-x2-dx3),x0+x2-dx4);
  G4double xmax =
    std::max(
    std::max(
    std::max(-x0-x1+dx1,-x0+x1+dx2),x0-x2+dx3),x0+x2+dx4);

  G4double y0 = dz*fTthetaSphi;
  G4double ymin = std::min(-y0-dy1,y0-dy2);
  G4double ymax = std::max(-y0+dy1,y0+dy2);

  pMin.set(xmin,ymin,-dz);
  pMax.set(xmax,ymax, dz);

  // Check correctness of the bounding box
  //
  if (pMin.x() >= pMax.x() || pMin.y() >= pMax.y() || pMin.z() >= pMax.z())
  {
    std::ostringstream message;
    message << "Bad bounding box (min >= max) for solid: "
            << GetName() << " !"
            << "\npMin = " << pMin
            << "\npMax = " << pMax;
    G4Exception("G4Trap::Extent()", "GeomMgt0001", JustWarning, message);
    DumpInfo();
  }
}

////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool G4Trap::CalculateExtent( const EAxis pAxis,
                                const G4VoxelLimits& pVoxelLimit,
                                const G4AffineTransform& pTransform,
                                      G4double& pMin, G4double& pMax) const
{
  G4ThreeVector bmin, bmax;
  G4bool exist;

  // Check bounding box (bbox)
  //
  Extent(bmin,bmax);
  G4BoundingEnvelope bbox(bmin,bmax);
#ifdef G4BBOX_EXTENT
  if (true) return bbox.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
#endif
  if (bbox.BoundingBoxVsVoxelLimits(pAxis,pVoxelLimit,pTransform,pMin,pMax))
  {
    return exist = (pMin < pMax) ? true : false;
  }

  // Set bounding envelope (benv) and calculate extent
  //
  G4double dz  = GetZHalfLength();
  G4double dx1 = GetXHalfLength1();
  G4double dx2 = GetXHalfLength2();
  G4double dx3 = GetXHalfLength3();
  G4double dx4 = GetXHalfLength4();
  G4double dy1 = GetYHalfLength1();
  G4double dy2 = GetYHalfLength2();

  G4double x0 = dz*fTthetaCphi;
  G4double x1 = dy1*GetTanAlpha1();
  G4double x2 = dy2*GetTanAlpha2();
  G4double y0 = dz*fTthetaSphi;

  G4ThreeVectorList baseA(4), baseB(4);
  baseA[0].set(-x0-x1-dx1,-y0-dy1,-dz);
  baseA[1].set(-x0-x1+dx1,-y0-dy1,-dz);
  baseA[2].set(-x0+x1+dx2,-y0+dy1,-dz);
  baseA[3].set(-x0+x1-dx2,-y0+dy1,-dz);

  baseB[0].set( x0-x2-dx3, y0-dy2, dz);
  baseB[1].set( x0-x2+dx3, y0-dy2, dz);
  baseB[2].set( x0+x2+dx4, y0+dy2, dz);
  baseB[3].set( x0+x2-dx4, y0+dy2, dz);

  std::vector<const G4ThreeVectorList *> polygons(2);
  polygons[0] = &baseA;
  polygons[1] = &baseB;

  G4BoundingEnvelope benv(bmin,bmax,polygons);
  exist = benv.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
  return exist;
}

////////////////////////////////////////////////////////////////////////
//
// Return whether point inside/outside/on surface, using tolerance

EInside G4Trap::Inside( const G4ThreeVector& p ) const
{
  EInside in;
  G4double Dist;
  G4int i;
  if ( std::fabs(p.z()) <= fDz-kCarTolerance*0.5)
  {
    in = kInside;

    for ( i = 0;i < 4;i++ )
    {
      Dist = fPlanes[i].a*p.x() + fPlanes[i].b*p.y()
            +fPlanes[i].c*p.z() + fPlanes[i].d;

      if      (Dist >  kCarTolerance*0.5)  return in = kOutside;
      else if (Dist > -kCarTolerance*0.5)         in = kSurface;
       
    }
  }
  else if (std::fabs(p.z()) <= fDz+kCarTolerance*0.5)
  {
    in = kSurface;

    for ( i = 0; i < 4; i++ )
    {
      Dist =  fPlanes[i].a*p.x() + fPlanes[i].b*p.y()
             +fPlanes[i].c*p.z() + fPlanes[i].d;

      if (Dist > kCarTolerance*0.5)        return in = kOutside;      
    }
  }
  else  in = kOutside;
  
  return in;
}

/////////////////////////////////////////////////////////////////////////////
//
// Calculate side nearest to p, and return normal
// If 2+ sides equidistant, first side's normal returned (arbitrarily)

G4ThreeVector G4Trap::SurfaceNormal( const G4ThreeVector& p ) const
{
  G4int i, noSurfaces = 0;
  G4double dist, distz, distx, disty, distmx, distmy, safe = kInfinity;
  G4double delta    = 0.5*kCarTolerance;
  G4ThreeVector norm, sumnorm(0.,0.,0.);

  for (i = 0; i < 4; i++)
  {
    dist =  std::fabs(fPlanes[i].a*p.x() + fPlanes[i].b*p.y()
          + fPlanes[i].c*p.z() + fPlanes[i].d);
    if ( dist < safe )
    {
      safe = dist;
    }
  }
  distz  = std::fabs( std::fabs( p.z() ) - fDz );

  distmy = std::fabs( fPlanes[0].a*p.x() + fPlanes[0].b*p.y()
                    + fPlanes[0].c*p.z() + fPlanes[0].d      );

  disty  = std::fabs( fPlanes[1].a*p.x() + fPlanes[1].b*p.y()
                    + fPlanes[1].c*p.z() + fPlanes[1].d      );

  distmx = std::fabs( fPlanes[2].a*p.x() + fPlanes[2].b*p.y()
                    + fPlanes[2].c*p.z() + fPlanes[2].d      );

  distx  = std::fabs( fPlanes[3].a*p.x() + fPlanes[3].b*p.y()
                    + fPlanes[3].c*p.z() + fPlanes[3].d      );

  G4ThreeVector nX  = G4ThreeVector(fPlanes[3].a,fPlanes[3].b,fPlanes[3].c);
  G4ThreeVector nmX = G4ThreeVector(fPlanes[2].a,fPlanes[2].b,fPlanes[2].c);
  G4ThreeVector nY  = G4ThreeVector(fPlanes[1].a,fPlanes[1].b,fPlanes[1].c);
  G4ThreeVector nmY = G4ThreeVector(fPlanes[0].a,fPlanes[0].b,fPlanes[0].c);
  G4ThreeVector nZ  = G4ThreeVector(0.,0.,1.0);

  if (distx <= delta)      
  {
    noSurfaces ++;
    sumnorm += nX;     
  }
  if (distmx <= delta)      
  {
    noSurfaces ++;
    sumnorm += nmX;      
  }
  if (disty <= delta)
  {
    noSurfaces ++;
    sumnorm += nY;  
  }
  if (distmy <= delta)
  {
    noSurfaces ++;
    sumnorm += nmY;  
  }
  if (distz <= delta)  
  {
    noSurfaces ++;
    if ( p.z() >= 0.)  sumnorm += nZ;
    else               sumnorm -= nZ; 
  }
  if ( noSurfaces == 0 )
  {
#ifdef G4CSGDEBUG
    G4Exception("G4Trap::SurfaceNormal(p)", "GeomSolids1002",
                JustWarning, "Point p is not on surface !?" );
#endif 
     norm = ApproxSurfaceNormal(p);
  }
  else if ( noSurfaces == 1 ) norm = sumnorm;
  else                        norm = sumnorm.unit();
  return norm;
}

////////////////////////////////////////////////////////////////////////////////////
//
// Algorithm for SurfaceNormal() following the original specification
// for points not on the surface

G4ThreeVector G4Trap::ApproxSurfaceNormal( const G4ThreeVector& p ) const
{
  G4double safe=kInfinity,Dist,safez;
  G4int i,imin=0;
  for (i=0;i<4;i++)
  {
    Dist=std::fabs(fPlanes[i].a*p.x()+fPlanes[i].b*p.y()
        +fPlanes[i].c*p.z()+fPlanes[i].d);
    if (Dist<safe)
    {
      safe=Dist;
      imin=i;
    }
  }
  safez=std::fabs(std::fabs(p.z())-fDz);
  if (safe<safez)
  {
    return G4ThreeVector(fPlanes[imin].a,fPlanes[imin].b,fPlanes[imin].c);
  }
  else
  {
    if (p.z()>0)
    {
      return G4ThreeVector(0,0,1);
    }
    else
    {
      return G4ThreeVector(0,0,-1);
    }
  }
}

////////////////////////////////////////////////////////////////////////////
//
// Calculate distance to shape from outside - return kInfinity if no intersection
//
// ALGORITHM:
// For each component, calculate pair of minimum and maximum intersection
// values for which the particle is in the extent of the shape
// - The smallest (MAX minimum) allowed distance of the pairs is intersect

G4double G4Trap::DistanceToIn( const G4ThreeVector& p,
                               const G4ThreeVector& v ) const
{

  G4double snxt;    // snxt = default return value
  G4double max,smax,smin;
  G4double pdist,Comp,vdist;
  G4int i;
  //
  // Z Intersection range
  //
  if ( v.z() > 0 )
  {
    max = fDz - p.z() ;
    if (max > 0.5*kCarTolerance)
    {
      smax = max/v.z();
      smin = (-fDz-p.z())/v.z();
    }
    else
    {
      return snxt=kInfinity;
    }
  }
  else if (v.z() < 0 )
  {
    max = - fDz - p.z() ;
    if (max < -0.5*kCarTolerance )
    {
      smax=max/v.z();
      smin=(fDz-p.z())/v.z();
    }
    else
    {
      return snxt=kInfinity;
    }
  }
  else
  {
    if (std::fabs(p.z())<fDz - 0.5*kCarTolerance) // Inside was <=fDz
    {
      smin=0;
      smax=kInfinity;
    }
    else
    {
      return snxt=kInfinity;
    }
  }

  for (i=0;i<4;i++)
  {
    pdist=fPlanes[i].a*p.x()+fPlanes[i].b*p.y()
         +fPlanes[i].c*p.z()+fPlanes[i].d;
    Comp=fPlanes[i].a*v.x()+fPlanes[i].b*v.y()+fPlanes[i].c*v.z();
    if ( pdist >= -0.5*kCarTolerance )      // was >0
    {
      //
      // Outside the plane -> this is an extent entry distance
      //
      if (Comp >= 0)   // was >0
      {
        return snxt=kInfinity ;
      }
      else 
      {
        vdist=-pdist/Comp;
        if (vdist>smin)
        {
          if (vdist<smax)
          {
            smin = vdist;
          }
          else
          {
            return snxt=kInfinity;
          }
        }
      }
    }
    else
    {
      //
      // Inside the plane -> couble  be an extent exit distance (smax)
      //
      if (Comp>0)  // Will leave extent
      {
        vdist=-pdist/Comp;
        if (vdist<smax)
        {
          if (vdist>smin)
          {
            smax=vdist;
          }
          else
          {
            return snxt=kInfinity;
          }
        }  
      }
    }
  }
  //
  // Checks in non z plane intersections ensure smin<smax
  //
  if (smin >=0 )
  {
    snxt = smin ;
  }
  else
  {
    snxt = 0 ;
  }
  return snxt;
}

///////////////////////////////////////////////////////////////////////////
//
// Calculate exact shortest distance to any boundary from outside
// This is the best fast estimation of the shortest distance to trap
// - Returns 0 is ThreeVector inside

G4double G4Trap::DistanceToIn( const G4ThreeVector& p ) const
{
  G4double safe=0.0,Dist;
  G4int i;
  safe=std::fabs(p.z())-fDz;
  for (i=0;i<4;i++)
  {
    Dist=fPlanes[i].a*p.x()+fPlanes[i].b*p.y()
        +fPlanes[i].c*p.z()+fPlanes[i].d;
    if (Dist > safe) safe=Dist;
  }
  if (safe<0) safe=0;
  return safe;  
}

/////////////////////////////////////////////////////////////////////////////////
//
// Calculate distance to surface of shape from inside
// Calculate distance to x/y/z planes - smallest is exiting distance

G4double G4Trap::DistanceToOut(const G4ThreeVector& p, const G4ThreeVector& v,
                               const G4bool calcNorm,
                                     G4bool *validNorm, G4ThreeVector *n) const
{
  Eside side = kUndef;
  G4double snxt;    // snxt = return value
  G4double pdist,Comp,vdist,max;
  //
  // Z Intersections
  //
  if (v.z()>0)
  {
    max=fDz-p.z();
    if (max>kCarTolerance/2)
    {
      snxt=max/v.z();
      side=kPZ;
    }
    else
    {
      if (calcNorm)
      {
        *validNorm=true;
        *n=G4ThreeVector(0,0,1);
      }
      return snxt=0;
    }
  }
  else if (v.z()<0)
  {
    max=-fDz-p.z();
    if (max<-kCarTolerance/2)
    {
      snxt=max/v.z();
      side=kMZ;
    }
    else
    {
      if (calcNorm)
      {
        *validNorm=true;
        *n=G4ThreeVector(0,0,-1);
      }
      return snxt=0;
    }
  }
  else
  {
    snxt=kInfinity;
  }

  //
  // Intersections with planes[0] (expanded because of setting enum)
  //
  pdist=fPlanes[0].a*p.x()+fPlanes[0].b*p.y()+fPlanes[0].c*p.z()+fPlanes[0].d;
  Comp=fPlanes[0].a*v.x()+fPlanes[0].b*v.y()+fPlanes[0].c*v.z();
  if (pdist>0)
  {
    // Outside the plane
    if (Comp>0)
    {
      // Leaving immediately
      if (calcNorm)
      {
        *validNorm=true;
        *n=G4ThreeVector(fPlanes[0].a,fPlanes[0].b,fPlanes[0].c);
      }
      return snxt=0;
    }
  }
  else if (pdist<-kCarTolerance/2)
  {
    // Inside the plane
    if (Comp>0)
    {
      // Will leave extent
      vdist=-pdist/Comp;
      if (vdist<snxt)
      {
        snxt=vdist;
        side=ks0;
      }
    }
  }
  else
  {
    // On surface
    if (Comp>0)
    {
      if (calcNorm)
      {
        *validNorm=true;
        *n=G4ThreeVector(fPlanes[0].a,fPlanes[0].b,fPlanes[0].c);
      }
      return snxt=0;
    }
  }

  //
  // Intersections with planes[1] (expanded because of setting enum)
  //
  pdist=fPlanes[1].a*p.x()+fPlanes[1].b*p.y()+fPlanes[1].c*p.z()+fPlanes[1].d;
  Comp=fPlanes[1].a*v.x()+fPlanes[1].b*v.y()+fPlanes[1].c*v.z();
  if (pdist>0)
  {
    // Outside the plane
    if (Comp>0)
    {
      // Leaving immediately
      if (calcNorm)
      {
        *validNorm=true;
        *n=G4ThreeVector(fPlanes[1].a,fPlanes[1].b,fPlanes[1].c);
      }
      return snxt=0;
    }
  }
  else if (pdist<-kCarTolerance/2)
  {
    // Inside the plane
    if (Comp>0)
    {
      // Will leave extent
      vdist=-pdist/Comp;
      if (vdist<snxt)
      {
        snxt=vdist;
        side=ks1;
      }
    }
  }
  else
  {
    // On surface
    if (Comp>0)
    {
      if (calcNorm)
      {
        *validNorm=true;
        *n=G4ThreeVector(fPlanes[1].a,fPlanes[1].b,fPlanes[1].c);
      }
      return snxt=0;
    }
  }

  //
  // Intersections with planes[2] (expanded because of setting enum)
  //
  pdist=fPlanes[2].a*p.x()+fPlanes[2].b*p.y()+fPlanes[2].c*p.z()+fPlanes[2].d;
  Comp=fPlanes[2].a*v.x()+fPlanes[2].b*v.y()+fPlanes[2].c*v.z();
  if (pdist>0)
  {
    // Outside the plane
    if (Comp>0)
    {
      // Leaving immediately
      if (calcNorm)
      {
        *validNorm=true;
        *n=G4ThreeVector(fPlanes[2].a,fPlanes[2].b,fPlanes[2].c);
      }
      return snxt=0;
    }
  }
  else if (pdist<-kCarTolerance/2)
  {
    // Inside the plane
    if (Comp>0)
    {
      // Will leave extent
      vdist=-pdist/Comp;
      if (vdist<snxt)
      {
        snxt=vdist;
        side=ks2;
      }
    }
  }
  else
  {
    // On surface
    if (Comp>0)
    {
      if (calcNorm)
      {
        *validNorm=true;
        *n=G4ThreeVector(fPlanes[2].a,fPlanes[2].b,fPlanes[2].c);
      }
      return snxt=0;
    }
  }

  //
  // Intersections with planes[3] (expanded because of setting enum)
  //
  pdist=fPlanes[3].a*p.x()+fPlanes[3].b*p.y()+fPlanes[3].c*p.z()+fPlanes[3].d;
  Comp=fPlanes[3].a*v.x()+fPlanes[3].b*v.y()+fPlanes[3].c*v.z();
  if (pdist>0)
  {
    // Outside the plane
    if (Comp>0)
    {
      // Leaving immediately
      if (calcNorm)
      {
        *validNorm=true;
        *n=G4ThreeVector(fPlanes[3].a,fPlanes[3].b,fPlanes[3].c);
      }
      return snxt=0;
    }
  }
  else if (pdist<-kCarTolerance/2)
  {
    // Inside the plane
    if (Comp>0)
    {
      // Will leave extent
      vdist=-pdist/Comp;
      if (vdist<snxt)
      {
        snxt=vdist;
        side=ks3;
      }
    }
  }
  else
  {
    // On surface
    if (Comp>0)
    {
      if (calcNorm)
      {
        *validNorm=true;
        *n=G4ThreeVector(fPlanes[3].a,fPlanes[3].b,fPlanes[3].c);
      }
      return snxt=0;
    }
  }

  // set normal
  if (calcNorm)
  {
    *validNorm=true;
    switch(side)
    {
      case ks0:
        *n=G4ThreeVector(fPlanes[0].a,fPlanes[0].b,fPlanes[0].c);
        break;
      case ks1:
        *n=G4ThreeVector(fPlanes[1].a,fPlanes[1].b,fPlanes[1].c);
        break;
      case ks2:
        *n=G4ThreeVector(fPlanes[2].a,fPlanes[2].b,fPlanes[2].c);
        break;
      case ks3:
        *n=G4ThreeVector(fPlanes[3].a,fPlanes[3].b,fPlanes[3].c);
        break;
      case kMZ:
        *n=G4ThreeVector(0,0,-1);
        break;
      case kPZ:
        *n=G4ThreeVector(0,0,1);
        break;
      default:
        G4cout << G4endl;
        DumpInfo();
        std::ostringstream message;
        G4int oldprc = message.precision(16);
        message << "Undefined side for valid surface normal to solid."
                << G4endl
                << "Position:"  << G4endl << G4endl
                << "p.x() = "   << p.x()/mm << " mm" << G4endl
                << "p.y() = "   << p.y()/mm << " mm" << G4endl
                << "p.z() = "   << p.z()/mm << " mm" << G4endl << G4endl
                << "Direction:" << G4endl << G4endl
                << "v.x() = "   << v.x() << G4endl
                << "v.y() = "   << v.y() << G4endl
                << "v.z() = "   << v.z() << G4endl << G4endl
                << "Proposed distance :" << G4endl << G4endl
                << "snxt = "    << snxt/mm << " mm" << G4endl;
        message.precision(oldprc);
        G4Exception("G4Trap::DistanceToOut(p,v,..)","GeomSolids1002",
                    JustWarning, message);
        break;
    }
  }
  return snxt;
}

//////////////////////////////////////////////////////////////////////////////
//
// Calculate exact shortest distance to any boundary from inside
// - Returns 0 is ThreeVector outside

G4double G4Trap::DistanceToOut( const G4ThreeVector& p ) const
{
  G4double safe=0.0,Dist;
  G4int i;

#ifdef G4CSGDEBUG
  if( Inside(p) == kOutside )
  {
     G4int oldprc = G4cout.precision(16) ;
     G4cout << G4endl ;
     DumpInfo();
     G4cout << "Position:"  << G4endl << G4endl ;
     G4cout << "p.x() = "   << p.x()/mm << " mm" << G4endl ;
     G4cout << "p.y() = "   << p.y()/mm << " mm" << G4endl ;
     G4cout << "p.z() = "   << p.z()/mm << " mm" << G4endl << G4endl ;
     G4cout.precision(oldprc) ;
     G4Exception("G4Trap::DistanceToOut(p)",
                 "GeomSolids1002", JustWarning, "Point p is outside !?" );
  }
#endif

  safe=fDz-std::fabs(p.z());
  if (safe<0) safe=0;
  else
  {
    for (i=0;i<4;i++)
    {
      Dist=-(fPlanes[i].a*p.x()+fPlanes[i].b*p.y()
            +fPlanes[i].c*p.z()+fPlanes[i].d);
      if (Dist<safe) safe=Dist;
    }
    if (safe<0) safe=0;
  }
  return safe;  
}

//////////////////////////////////////////////////////////////////////////
//
// GetEntityType

G4GeometryType G4Trap::GetEntityType() const
{
  return G4String("G4Trap");
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object
//
G4VSolid* G4Trap::Clone() const
{
  return new G4Trap(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Stream object contents to an output stream

std::ostream& G4Trap::StreamInfo( std::ostream& os ) const
{
  G4int oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: G4Trap\n"
     << " Parameters: \n"
     << "    half length Z: " << fDz/mm << " mm \n"
     << "    half length Y of face -fDz: " << fDy1/mm << " mm \n"
     << "    half length X of side -fDy1, face -fDz: " << fDx1/mm << " mm \n"
     << "    half length X of side +fDy1, face -fDz: " << fDx2/mm << " mm \n"
     << "    half length Y of face +fDz: " << fDy2/mm << " mm \n"
     << "    half length X of side -fDy2, face +fDz: " << fDx3/mm << " mm \n"
     << "    half length X of side +fDy2, face +fDz: " << fDx4/mm << " mm \n"
     << "    std::tan(theta)*std::cos(phi): " << fTthetaCphi/degree << " degrees \n"
     << "    std::tan(theta)*std::sin(phi): " << fTthetaSphi/degree << " degrees \n"
     << "    std::tan(alpha), -fDz: " << fTalpha1/degree << " degrees \n"
     << "    std::tan(alpha), +fDz: " << fTalpha2/degree << " degrees \n"
     << "    trap side plane equations:\n"
     << "        " << fPlanes[0].a << " X + " << fPlanes[0].b << " Y + "
                   << fPlanes[0].c << " Z + " << fPlanes[0].d << " = 0\n"
     << "        " << fPlanes[1].a << " X + " << fPlanes[1].b << " Y + "
                   << fPlanes[1].c << " Z + " << fPlanes[1].d << " = 0\n"
     << "        " << fPlanes[2].a << " X + " << fPlanes[2].b << " Y + "
                   << fPlanes[2].c << " Z + " << fPlanes[2].d << " = 0\n"
     << "        " << fPlanes[3].a << " X + " << fPlanes[3].b << " Y + "
                   << fPlanes[3].c << " Z + " << fPlanes[3].d << " = 0\n"
     << "-----------------------------------------------------------\n";
  os.precision(oldprc);

  return os;
}

/////////////////////////////////////////////////////////////////////////
//
// GetPointOnPlane
//
// Auxiliary method for Get Point on Surface

G4ThreeVector G4Trap::GetPointOnPlane(G4ThreeVector p0, G4ThreeVector p1, 
                                      G4ThreeVector p2, G4ThreeVector p3,
                                      G4double& area) const
{
  G4double lambda1, lambda2, chose, aOne, aTwo;
  G4ThreeVector t, u, v, w, Area, normal;
  
  t = p1 - p0;
  u = p2 - p1;
  v = p3 - p2;
  w = p0 - p3;

  Area = G4ThreeVector(w.y()*v.z() - w.z()*v.y(),
                       w.z()*v.x() - w.x()*v.z(),
                       w.x()*v.y() - w.y()*v.x());
  
  aOne = 0.5*Area.mag();
  
  Area = G4ThreeVector(t.y()*u.z() - t.z()*u.y(),
                       t.z()*u.x() - t.x()*u.z(),
                       t.x()*u.y() - t.y()*u.x());
  
  aTwo = 0.5*Area.mag();
  
  area = aOne + aTwo;
  
  chose = G4RandFlat::shoot(0.,aOne+aTwo);

  if( (chose>=0.) && (chose < aOne) )
  {
    lambda1 = G4RandFlat::shoot(0.,1.);
    lambda2 = G4RandFlat::shoot(0.,lambda1);
    return (p2+lambda1*v+lambda2*w);    
  }
  
  // else

  lambda1 = G4RandFlat::shoot(0.,1.);
  lambda2 = G4RandFlat::shoot(0.,lambda1);

  return (p0+lambda1*t+lambda2*u);    
}

///////////////////////////////////////////////////////////////
//
// GetPointOnSurface

G4ThreeVector G4Trap::GetPointOnSurface() const
{
  G4double aOne, aTwo, aThree, aFour, aFive, aSix, chose;
  G4ThreeVector One, Two, Three, Four, Five, Six, test;
  G4ThreeVector pt[8];
     
  pt[0] = G4ThreeVector(-fDz*fTthetaCphi-fDy1*fTalpha1-fDx1,
                        -fDz*fTthetaSphi-fDy1,-fDz);
  pt[1] = G4ThreeVector(-fDz*fTthetaCphi-fDy1*fTalpha1+fDx1,
                        -fDz*fTthetaSphi-fDy1,-fDz);
  pt[2] = G4ThreeVector(-fDz*fTthetaCphi+fDy1*fTalpha1-fDx2,
                        -fDz*fTthetaSphi+fDy1,-fDz);
  pt[3] = G4ThreeVector(-fDz*fTthetaCphi+fDy1*fTalpha1+fDx2,
                        -fDz*fTthetaSphi+fDy1,-fDz);
  pt[4] = G4ThreeVector(+fDz*fTthetaCphi-fDy2*fTalpha2-fDx3,
                        +fDz*fTthetaSphi-fDy2,+fDz);
  pt[5] = G4ThreeVector(+fDz*fTthetaCphi-fDy2*fTalpha2+fDx3,
                        +fDz*fTthetaSphi-fDy2,+fDz);
  pt[6] = G4ThreeVector(+fDz*fTthetaCphi+fDy2*fTalpha2-fDx4,
                        +fDz*fTthetaSphi+fDy2,+fDz);
  pt[7] = G4ThreeVector(+fDz*fTthetaCphi+fDy2*fTalpha2+fDx4,
                        +fDz*fTthetaSphi+fDy2,+fDz);
  
  // make sure we provide the points in a clockwise fashion

  One   = GetPointOnPlane(pt[0],pt[1],pt[3],pt[2], aOne);
  Two   = GetPointOnPlane(pt[4],pt[5],pt[7],pt[6], aTwo);
  Three = GetPointOnPlane(pt[6],pt[7],pt[3],pt[2], aThree);
  Four  = GetPointOnPlane(pt[4],pt[5],pt[1],pt[0], aFour); 
  Five  = GetPointOnPlane(pt[0],pt[2],pt[6],pt[4], aFive);
  Six   = GetPointOnPlane(pt[1],pt[3],pt[7],pt[5], aSix);
 
  chose = G4RandFlat::shoot(0.,aOne+aTwo+aThree+aFour+aFive+aSix);
  if( (chose>=0.) && (chose<aOne) )                    
    { return One; }
  else if( (chose>=aOne) && (chose<aOne+aTwo) )  
    { return Two; }
  else if( (chose>=aOne+aTwo) && (chose<aOne+aTwo+aThree) )
    { return Three; }
  else if( (chose>=aOne+aTwo+aThree) && (chose<aOne+aTwo+aThree+aFour) )
    { return Four; }
  else if( (chose>=aOne+aTwo+aThree+aFour)
        && (chose<aOne+aTwo+aThree+aFour+aFive) )
    { return Five; }
  return Six;
}

//////////////////////////////////////////////////////////////////////////
//
// Methods for visualisation

void G4Trap::DescribeYourselfTo ( G4VGraphicsScene& scene ) const
{
  scene.AddSolid (*this);
}

G4Polyhedron* G4Trap::CreatePolyhedron () const
{
  G4double phi = std::atan2(fTthetaSphi, fTthetaCphi);
  G4double alpha1 = std::atan(fTalpha1);
  G4double alpha2 = std::atan(fTalpha2);
  G4double theta = std::atan(std::sqrt(fTthetaCphi*fTthetaCphi+fTthetaSphi*fTthetaSphi));

  return new G4PolyhedronTrap(fDz, theta, phi,
                              fDy1, fDx1, fDx2, alpha1,
                              fDy2, fDx3, fDx4, alpha2);
}

#endif
