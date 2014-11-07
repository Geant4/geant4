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
// $Id:$
//
// 
// Implementation of G4UMultiUnion wrapper class
// --------------------------------------------------------------------

#include "G4UMultiUnion.hh"
#include "G4Polyhedron.hh"
#include "G4DisplacedSolid.hh"
#include "G4RotationMatrix.hh"

////////////////////////////////////////////////////////////////////////
//
// Constructor (generic parameters)
//
G4UMultiUnion::G4UMultiUnion(const G4String& name)
  : G4USolid(name, new UMultiUnion(name))
{ 
}


////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4UMultiUnion::G4UMultiUnion(__void__& a)
  : G4USolid(a)
{
}


//////////////////////////////////////////////////////////////////////////
//
// Destructor
//
G4UMultiUnion::~G4UMultiUnion()
{
}


//////////////////////////////////////////////////////////////////////////
//
// Copy constructor
//
G4UMultiUnion::G4UMultiUnion(const G4UMultiUnion &source)
  : G4USolid(source)
{
}


//////////////////////////////////////////////////////////////////////////
//
// Assignment operator
//
G4UMultiUnion& G4UMultiUnion::operator=(const G4UMultiUnion &source)
{
  if (this == &source) return *this;
  
  G4USolid::operator=( source );
  
  return *this;
}


//////////////////////////////////////////////////////////////////////////
//
// CreatePolyhedron
//
G4Polyhedron* G4UMultiUnion::CreatePolyhedron() const
{

  HepPolyhedronProcessor processor;
  HepPolyhedronProcessor::Operation operation = HepPolyhedronProcessor::UNION;

  G4VSolid* solidA = GetSolid(0);
  const G4Transform3D* transform0=GetTransformation(0);
  G4RotationMatrix rot0=(*transform0).getRotation();
  const  G4ThreeVector transl0 = (*transform0).getTranslation();
  G4DisplacedSolid dispSolidA("placedA",solidA,&rot0,transl0);
  delete transform0;

  G4Polyhedron* top = new G4Polyhedron(*dispSolidA.GetPolyhedron());
    
  for(G4int i=1; i<GetNumberOfSolids(); ++i)
  {
    G4VSolid* solidB = GetSolid(i);
    const G4Transform3D* transform=GetTransformation(i);
    G4RotationMatrix rot=(*transform).getRotation();
    const  G4ThreeVector transl = (*transform).getTranslation();
    G4DisplacedSolid dispSolidB("placedB",solidB,&rot,transl);
    G4Polyhedron* operand = dispSolidB.GetPolyhedron();
    processor.push_back (operation, *operand);
    delete transform;
  }
   
  if (processor.execute(*top)) { return top; }
  else { return 0; } 
}
