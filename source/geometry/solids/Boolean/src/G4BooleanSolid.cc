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
// $Id: G4BooleanSolid.cc,v 1.22 2010-05-11 09:11:24 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Implementation for the abstract base class for solids created by boolean 
// operations between other solids
//
// History:
//
// 10.09.98 V.Grichine, created
//
// --------------------------------------------------------------------

#include "G4BooleanSolid.hh"
#include "G4VSolid.hh"
#include "G4Polyhedron.hh"
#include "HepPolyhedronProcessor.h"
#include "Randomize.hh"

//////////////////////////////////////////////////////////////////
//
// Constructor

G4BooleanSolid::G4BooleanSolid( const G4String& pName,
                                G4VSolid* pSolidA ,
                                G4VSolid* pSolidB   ) :
  G4VSolid(pName), fStatistics(1000000), fCubVolEpsilon(0.001),
  fAreaAccuracy(-1.), fCubicVolume(0.), fSurfaceArea(0.),
  fpPolyhedron(0), createdDisplacedSolid(false)
{
  fPtrSolidA = pSolidA ;
  fPtrSolidB = pSolidB ;
}

//////////////////////////////////////////////////////////////////
//
// Constructor

G4BooleanSolid::G4BooleanSolid( const G4String& pName,
                                      G4VSolid* pSolidA ,
                                      G4VSolid* pSolidB ,
                                      G4RotationMatrix* rotMatrix,
                                const G4ThreeVector& transVector    ) :
  G4VSolid(pName), fStatistics(1000000), fCubVolEpsilon(0.001),
  fAreaAccuracy(-1.), fCubicVolume(0.), fSurfaceArea(0.),
  fpPolyhedron(0), createdDisplacedSolid(true)
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
                                const G4Transform3D& transform    ) :
  G4VSolid(pName), fStatistics(1000000), fCubVolEpsilon(0.001),
  fAreaAccuracy(-1.), fCubicVolume(0.), fSurfaceArea(0.),
  fpPolyhedron(0), createdDisplacedSolid(true)
{
  fPtrSolidA = pSolidA ;
  fPtrSolidB = new G4DisplacedSolid("placedB",pSolidB,transform) ;
}

///////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.

G4BooleanSolid::G4BooleanSolid( __void__& a )
  : G4VSolid(a), fPtrSolidA(0), fPtrSolidB(0),
    fStatistics(1000000), fCubVolEpsilon(0.001), 
    fAreaAccuracy(-1.), fCubicVolume(0.), fSurfaceArea(0.),
    fpPolyhedron(0), createdDisplacedSolid(false)
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
  delete fpPolyhedron;
}

///////////////////////////////////////////////////////////////
//
// If Solid is made up from a Boolean operation of two solids,
//   return the corresponding solid (for no=0 and 1)
// If the solid is not a "Boolean", return 0

const G4VSolid* G4BooleanSolid::GetConstituentSolid(G4int no) const
{
  const G4VSolid*  subSolid=0;
  if( no == 0 )  
    subSolid = fPtrSolidA;
  else if( no == 1 ) 
    subSolid = fPtrSolidB;
  else
  {
    DumpInfo();
    G4Exception("G4BooleanSolid::GetConstituentSolid()",
                "WrongArgumentValue", FatalException,
                "Invalid solid index.");
  }

  return subSolid;
}

///////////////////////////////////////////////////////////////
//
// If Solid is made up from a Boolean operation of two solids,
//   return the corresponding solid (for no=0 and 1)
// If the solid is not a "Boolean", return 0

G4VSolid* G4BooleanSolid::GetConstituentSolid(G4int no)
{
  G4VSolid*  subSolid=0;
  if( no == 0 )  
    subSolid = fPtrSolidA;
  else if( no == 1 ) 
    subSolid = fPtrSolidB;
  else
  {
    DumpInfo();
    G4Exception("G4BooleanSolid::GetConstituentSolid()",
                "WrongArgumentValue", FatalException,
                "Invalid solid index.");
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
// Returns a point (G4ThreeVector) randomly and uniformly selected
// on the solid surface
//

G4ThreeVector G4BooleanSolid::GetPointOnSurface() const
{
  G4bool condition = true;
  G4double rand;
  G4ThreeVector p;

  while(condition)
  {
    rand = G4UniformRand();

    if(rand > 0.5) { p = fPtrSolidA->GetPointOnSurface(); }
    else           { p = fPtrSolidB->GetPointOnSurface(); }

    if(Inside(p) == kSurface)  { break; }
  }
  return p;
}

//////////////////////////////////////////////////////////////////////////
//
// Returns polyhedron for visualization

G4Polyhedron* G4BooleanSolid::GetPolyhedron () const
{
  if (!fpPolyhedron ||
      fpPolyhedron->GetNumberOfRotationStepsAtTimeOfCreation() !=
      fpPolyhedron->GetNumberOfRotationSteps())
    {
      delete fpPolyhedron;
      fpPolyhedron = CreatePolyhedron();
    }
  return fpPolyhedron;
}

//////////////////////////////////////////////////////////////////////////
//
// Stack polyhedra for processing.  Return top polyhedron.

G4Polyhedron* G4BooleanSolid::StackPolyhedron
(HepPolyhedronProcessor& processor, const G4VSolid* solid) const
{
  HepPolyhedronProcessor::Operation operation;
  const G4String& type = solid->GetEntityType();
  if (type == "G4UnionSolid")
    operation = HepPolyhedronProcessor::UNION;
  else if (type == "G4IntersectionSolid")
    operation = HepPolyhedronProcessor::INTERSECTION;
  else if (type == "G4SubtractionSolid")
    operation = HepPolyhedronProcessor::SUBTRACTION;
  else {
    std::ostringstream oss;
    oss << "Solid - " << solid->GetName() <<
      " - unrecognised composite solid"
      "\n  Returning NULL !";
    G4Exception("StackPolyhedron()", "InvalidSetup",
                JustWarning, oss.str().c_str());
    return 0;
  }

  G4Polyhedron* top = 0;
  const G4VSolid* solidA = solid->GetConstituentSolid(0);
  const G4VSolid* solidB = solid->GetConstituentSolid(1);

  /***
  // Debug
  G4cout << "\nComposite: "
	 << solid->GetName() << ": " << solid->GetEntityType()
	 << G4endl;
  const G4DisplacedSolid* displacedSolidB = solidB->GetDisplacedSolidPtr();
  const G4VSolid* undisplacedSolidB = solidB;
  if (displacedSolidB) undisplacedSolidB =
    displacedSolidB->GetConstituentMovedSolid();
  G4cout << "Constituents:\n"
	 << solidA->GetName() << ": " << solidA->GetEntityType()
	 << '\n' << solidB->GetName() << ": " << solidB->GetEntityType();
  if (displacedSolidB)
    G4cout << ", " << undisplacedSolidB->GetName()
	   << ": " << undisplacedSolidB->GetEntityType();
  G4cout << G4endl;
  // End debug
  ***/

  if (solidA->GetConstituentSolid(0)) {
    top = StackPolyhedron(processor, solidA);
  } else {
    top = solidA->GetPolyhedron();

    /***
    // Debug
    G4cout << "Found top solid " << solidA->GetName()
	   << ' ' << (void*)top << G4endl;
    // End debug
    ***/
  }
  G4Polyhedron* operand = solidB->GetPolyhedron();

  /***
  // Debug
  G4cout << "Pushing " << solidB->GetName()
	 << ' ' << operation << ':' << (void*)operand << G4endl;
  G4cout << "Returning " << ' ' << (void*)top << G4endl;
  // End debug
  ***/

  processor.push_back (operation, *operand);
  return top;
}
