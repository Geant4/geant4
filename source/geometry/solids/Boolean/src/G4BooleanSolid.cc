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
// Implementation of the base class for solids created by Boolean 
// operations between other solids
//
// 1998.09.10 V.Grichine - created
// --------------------------------------------------------------------

#include "G4BooleanSolid.hh"
#include "G4VSolid.hh"
#include "G4DisplacedSolid.hh"
#include "G4ReflectedSolid.hh"
#include "G4ScaledSolid.hh"
#include "G4Polyhedron.hh"
#include "HepPolyhedronProcessor.h"
#include "G4QuickRand.hh"

#include "G4AutoLock.hh"

namespace
{
  G4RecursiveMutex polyhedronMutex = G4MUTEX_INITIALIZER;
}

//////////////////////////////////////////////////////////////////
//
// Constructor

G4BooleanSolid::G4BooleanSolid( const G4String& pName,
                                G4VSolid* pSolidA ,
                                G4VSolid* pSolidB )
  : G4VSolid(pName), fPtrSolidA(pSolidA), fPtrSolidB(pSolidB)
{
}

//////////////////////////////////////////////////////////////////
//
// Constructor

G4BooleanSolid::G4BooleanSolid( const G4String& pName,
                                      G4VSolid* pSolidA ,
                                      G4VSolid* pSolidB ,
                                      G4RotationMatrix* rotMatrix,
                                const G4ThreeVector& transVector )
  : G4VSolid(pName), createdDisplacedSolid(true)
{
  fPtrSolidA = pSolidA ;
  fPtrSolidB = new G4DisplacedSolid("placedB",pSolidB,rotMatrix,transVector) ;
}

//////////////////////////////////////////////////////////////////
//
// Constructor

G4BooleanSolid::G4BooleanSolid( const G4String& pName,
                                      G4VSolid* pSolidA ,
                                      G4VSolid* pSolidB ,
                                const G4Transform3D& transform )
  : G4VSolid(pName), createdDisplacedSolid(true)
{
  fPtrSolidA = pSolidA ;
  fPtrSolidB = new G4DisplacedSolid("placedB",pSolidB,transform) ;
}

///////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.

G4BooleanSolid::G4BooleanSolid( __void__& a )
  : G4VSolid(a)
{
}

///////////////////////////////////////////////////////////////
//
// Destructor deletes transformation contents of the created displaced solid

G4BooleanSolid::~G4BooleanSolid() 
{
  if(createdDisplacedSolid)
  {
    ((G4DisplacedSolid*)fPtrSolidB)->CleanTransformations();
  }
  delete fpPolyhedron; fpPolyhedron = nullptr;
}

///////////////////////////////////////////////////////////////
//
// Copy constructor

G4BooleanSolid::G4BooleanSolid(const G4BooleanSolid& rhs)
  : G4VSolid (rhs), fPtrSolidA(rhs.fPtrSolidA), fPtrSolidB(rhs.fPtrSolidB),
    fCubicVolume(rhs.fCubicVolume), fStatistics(rhs.fStatistics),
    fCubVolEpsilon(rhs.fCubVolEpsilon), fAreaAccuracy(rhs.fAreaAccuracy), 
    fSurfaceArea(rhs.fSurfaceArea), fRebuildPolyhedron(false),
    fpPolyhedron(nullptr), createdDisplacedSolid(rhs.createdDisplacedSolid)
{
  fPrimitives.resize(0); fPrimitivesSurfaceArea = 0.;
}

///////////////////////////////////////////////////////////////
//
// Assignment operator

G4BooleanSolid& G4BooleanSolid::operator = (const G4BooleanSolid& rhs) 
{
  // Check assignment to self
  //
  if (this == &rhs)  { return *this; }

  // Copy base class data
  //
  G4VSolid::operator=(rhs);

  // Copy data
  //
  fPtrSolidA= rhs.fPtrSolidA; fPtrSolidB= rhs.fPtrSolidB;
  fStatistics= rhs.fStatistics; fCubVolEpsilon= rhs.fCubVolEpsilon;
  fAreaAccuracy= rhs.fAreaAccuracy; fCubicVolume= rhs.fCubicVolume;
  fSurfaceArea= rhs.fSurfaceArea;
  createdDisplacedSolid= rhs.createdDisplacedSolid;
  fRebuildPolyhedron = false;
  delete fpPolyhedron; fpPolyhedron = nullptr;
  fPrimitives.resize(0); fPrimitivesSurfaceArea = 0.;

  return *this;
}  

///////////////////////////////////////////////////////////////
//
// If solid is made up from a Boolean operation of two solids,
// return the corresponding solid (for no=0 and 1)
// If the solid is not a "Boolean", return 0

const G4VSolid* G4BooleanSolid::GetConstituentSolid(G4int no) const
{
  const G4VSolid* subSolid = nullptr;
  if( no == 0 )  
    subSolid = fPtrSolidA;
  else if( no == 1 ) 
    subSolid = fPtrSolidB;
  else
  {
    DumpInfo();
    G4Exception("G4BooleanSolid::GetConstituentSolid()",
                "GeomSolids0002", FatalException, "Invalid solid index.");
  }
  return subSolid;
}

///////////////////////////////////////////////////////////////
//
// If solid is made up from a Boolean operation of two solids,
// return the corresponding solid (for no=0 and 1)
// If the solid is not a "Boolean", return 0

G4VSolid* G4BooleanSolid::GetConstituentSolid(G4int no)
{
  G4VSolid* subSolid = nullptr;
  if( no == 0 )  
    subSolid = fPtrSolidA;
  else if( no == 1 ) 
    subSolid = fPtrSolidB;
  else
  {
    DumpInfo();
    G4Exception("G4BooleanSolid::GetConstituentSolid()",
                "GeomSolids0002", FatalException, "Invalid solid index.");
  }
  return subSolid;
}

//////////////////////////////////////////////////////////////////////////
//
// Returns entity type

G4GeometryType G4BooleanSolid::GetEntityType() const 
{
  return G4String("G4BooleanSolid");
}

//////////////////////////////////////////////////////////////////////////
//
// Stream object contents to an output stream

std::ostream& G4BooleanSolid::StreamInfo(std::ostream& os) const
{
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for Boolean solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: " << GetEntityType() << "\n"
     << " Parameters of constituent solids: \n"
     << "===========================================================\n";
  fPtrSolidA->StreamInfo(os);
  fPtrSolidB->StreamInfo(os);
  os << "===========================================================\n";

  return os;
}

//////////////////////////////////////////////////////////////////////////
//
// Creates list of constituent primitives of and their placements

void G4BooleanSolid::GetListOfPrimitives(
       std::vector<std::pair<G4VSolid*,G4Transform3D>>& primitives,
       const G4Transform3D& curPlacement) const
{
  G4Transform3D transform;
  G4VSolid* solid;
  G4String type;

  // Repeat two times, first time for fPtrSolidA and then for fPtrSolidB
  //
  for (auto i=0; i<2; ++i)
  {
    transform = curPlacement;
    solid     = (i == 0) ? fPtrSolidA : fPtrSolidB;
    type      = solid->GetEntityType();

    // While current solid is a trasformed solid just modify transform
    //
    while (type == "G4DisplacedSolid" ||
           type == "G4ReflectedSolid" ||
           type == "G4ScaledSolid")
    {
      if (type == "G4DisplacedSolid")
      {
        transform = transform * G4Transform3D(
                    ((G4DisplacedSolid*)solid)->GetObjectRotation(),
                    ((G4DisplacedSolid*)solid)->GetObjectTranslation());
        solid     = ((G4DisplacedSolid*)solid)->GetConstituentMovedSolid();
      }
      else if (type == "G4ReflectedSolid")
      {
        transform= transform*((G4ReflectedSolid*)solid)->GetDirectTransform3D();
        solid    = ((G4ReflectedSolid*)solid)->GetConstituentMovedSolid();
      }
      else if (type == "G4ScaledSolid")
      {
        transform = transform * ((G4ScaledSolid*)solid)->GetScaleTransform();
        solid     = ((G4ScaledSolid*)solid)->GetUnscaledSolid();
      }
      type  = solid->GetEntityType();
    }

    // If current solid is a Boolean solid then continue recursion,
    // otherwise add it to the list of primitives
    //
    if (type == "G4UnionSolid"        ||
        type == "G4SubtractionSolid"  ||
        type == "G4IntersectionSolid" ||
        type == "G4BooleanSolid")
    {
      ((G4BooleanSolid *)solid)->GetListOfPrimitives(primitives,transform);
    }
    else
    {
      primitives.push_back(std::pair<G4VSolid*,G4Transform3D>(solid,transform));
    }
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Returns a point (G4ThreeVector) randomly and uniformly selected
// on the surface of the solid

G4ThreeVector G4BooleanSolid::GetPointOnSurface() const
{
  std::size_t nprims = fPrimitives.size();
  std::pair<G4VSolid *, G4Transform3D> prim;

  // Get list of primitives and find the total area of their surfaces
  //
  if (nprims == 0)
  {
    GetListOfPrimitives(fPrimitives, G4Transform3D());
    nprims = fPrimitives.size();
    fPrimitivesSurfaceArea = 0.;
    for (std::size_t i=0; i<nprims; ++i)
    {
      fPrimitivesSurfaceArea += fPrimitives[i].first->GetSurfaceArea();
    }
  }

  // Select random primitive, get random point on its surface and
  // check that the point belongs to the surface of the solid
  //
  G4ThreeVector p;
  for (std::size_t k=0; k<100000; ++k) // try 100k times
  {
     G4double rand = fPrimitivesSurfaceArea * G4QuickRand();
     G4double area = 0.;
     for (std::size_t i=0; i<nprims; ++i)
     {
       prim  = fPrimitives[i];
       area += prim.first->GetSurfaceArea();
       if (rand < area) break;
     }
     p = prim.first->GetPointOnSurface();
     p = prim.second * G4Point3D(p);
     if (Inside(p) == kSurface) return p;
  }
  std::ostringstream message;
  message << "Solid - " << GetName() << "\n"
          << "All 100k attempts to generate a point on the surface have failed!\n"
          << "The solid created may be an invalid Boolean construct!";
  G4Exception("G4BooleanSolid::GetPointOnSurface()",
              "GeomSolids1001", JustWarning, message);
  return p;
}

//////////////////////////////////////////////////////////////////////////
//
// Returns polyhedron for visualization

G4Polyhedron* G4BooleanSolid::GetPolyhedron () const
{
  if (fpPolyhedron == nullptr ||
      fRebuildPolyhedron ||
      fpPolyhedron->GetNumberOfRotationStepsAtTimeOfCreation() !=
      fpPolyhedron->GetNumberOfRotationSteps())
    {
      G4RecursiveAutoLock l(&polyhedronMutex);
      delete fpPolyhedron;
      fpPolyhedron = CreatePolyhedron();
      fRebuildPolyhedron = false;
      l.unlock();
    }
  return fpPolyhedron;
}

//////////////////////////////////////////////////////////////////////////
//
// Stacks polyhedra for processing. Returns top polyhedron.

G4Polyhedron*
G4BooleanSolid::StackPolyhedron(HepPolyhedronProcessor& processor,
                                const G4VSolid* solid) const
{
  HepPolyhedronProcessor::Operation operation;
  const G4String& type = solid->GetEntityType();
  if (type == "G4UnionSolid")
    { operation = HepPolyhedronProcessor::UNION; }
  else if (type == "G4IntersectionSolid")
    { operation = HepPolyhedronProcessor::INTERSECTION; }
  else if (type == "G4SubtractionSolid")
    { operation = HepPolyhedronProcessor::SUBTRACTION; }
  else
  {
    std::ostringstream message;
    message << "Solid - " << solid->GetName()
            << " - Unrecognised composite solid" << G4endl
            << " Returning NULL !";
    G4Exception("StackPolyhedron()", "GeomSolids1001", JustWarning, message);
    return nullptr;
  }

  G4Polyhedron* top = nullptr;
  const G4VSolid* solidA = solid->GetConstituentSolid(0);
  const G4VSolid* solidB = solid->GetConstituentSolid(1);

  if (solidA->GetConstituentSolid(0) != nullptr)
  {
    top = StackPolyhedron(processor, solidA);
  }
  else
  {
    top = solidA->GetPolyhedron();
  }
  G4Polyhedron* operand = solidB->GetPolyhedron();
  if (operand != nullptr)
  {
    processor.push_back (operation, *operand);
  }
  else
  {
    std::ostringstream message;
    message << "Solid - " << solid->GetName()
            << " - No G4Polyhedron for Boolean component";
    G4Exception("G4BooleanSolid::StackPolyhedron()",
                "GeomSolids2001", JustWarning, message);
  }

  return top;
}


//////////////////////////////////////////////////////////////////////////
//
// Estimate Cubic Volume (capacity) and store it for reuse.

G4double G4BooleanSolid::GetCubicVolume()
{
  if(fCubicVolume < 0.)
  {
    fCubicVolume = EstimateCubicVolume(fStatistics,fCubVolEpsilon);
  }
  return fCubicVolume;
}
