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
// Implementation of G4UGenericTrap wrapper class
// --------------------------------------------------------------------

#include "G4GenericTrap.hh"
#include "G4UGenericTrap.hh"
#include "G4Polyhedron.hh"
#include "G4PolyhedronArbitrary.hh"

////////////////////////////////////////////////////////////////////////
//
// Constructor (generic parameters)
//
G4UGenericTrap::G4UGenericTrap(const G4String& name, G4double halfZ,
                               const std::vector<G4TwoVector>& vertices)
  : G4USolid(name, new UGenericTrap())
{
  SetZHalfLength(halfZ);
  std::vector<UVector2> v;
  for (size_t n=0; n<vertices.size(); ++n)
  {
    v.push_back(UVector2(vertices[n].x(),vertices[n].y()));
  }
  GetShape()->SetName(name);
  GetShape()->Initialise(v);
}


////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4UGenericTrap::G4UGenericTrap(__void__& a)
  : G4USolid(a)
{
}


//////////////////////////////////////////////////////////////////////////
//
// Destructor
//
G4UGenericTrap::~G4UGenericTrap()
{
}


//////////////////////////////////////////////////////////////////////////
//
// Copy constructor
//
G4UGenericTrap::G4UGenericTrap(const G4UGenericTrap &source)
  : G4USolid(source)
{
}


//////////////////////////////////////////////////////////////////////////
//
// Assignment operator
//
G4UGenericTrap&
G4UGenericTrap::operator=(const G4UGenericTrap &source)
{
  if (this == &source) return *this;
  
  G4USolid::operator=( source );
  
  return *this;
}

//////////////////////////////////////////////////////////////////////////
//
// CreatePolyhedron()
//
G4Polyhedron* G4UGenericTrap::CreatePolyhedron() const
{
  // Approximation of Twisted Side
  // Construct extra Points, if Twisted Side
  //
  G4PolyhedronArbitrary* polyhedron;
  size_t nVertices, nFacets;
  G4double fDz = GetZHalfLength();

  G4int subdivisions=0;
  G4int i;
  if(IsTwisted())
  {
    if ( GetVisSubdivisions()!= 0 )
    {
      subdivisions=GetVisSubdivisions();
    }
    else
    {
      // Estimation of Number of Subdivisions for smooth visualisation
      //
      G4double maxTwist=0.;
      for(i=0; i<4; i++)
      {
        if(GetTwistAngle(i)>maxTwist) { maxTwist=GetTwistAngle(i); }
      }

      // Computes bounding vectors for the shape
      //
      G4double Dx,Dy;
      UVector3 minBox = GetShape()->GetMinimumBBox();
      UVector3 maxBox = GetShape()->GetMaximumBBox();
      G4ThreeVector minVec(minBox.x(), minBox.y(), minBox.z());
      G4ThreeVector maxVec(maxBox.x(), maxBox.y(), maxBox.z());
      Dx = 0.5*(maxVec.x()- minVec.y());
      Dy = 0.5*(maxVec.y()- minVec.y());
      if (Dy > Dx)  { Dx=Dy; }
    
      subdivisions=8*G4int(maxTwist/(Dx*Dx*Dx)*fDz);
      if (subdivisions<4)  { subdivisions=4; }
      if (subdivisions>30) { subdivisions=30; }
    }
  }
  G4int sub4=4*subdivisions;
  nVertices = 8+subdivisions*4;
  nFacets = 6+subdivisions*4;
  G4double cf=1./(subdivisions+1);
  polyhedron = new G4PolyhedronArbitrary (nVertices, nFacets);

  // Add Vertex
  //
  for (i=0;i<4;i++)
  {
    polyhedron->AddVertex(G4ThreeVector(GetVertex(i).x(),
                                        GetVertex(i).y(),-fDz));
  }
  for( i=0;i<subdivisions;i++)
  {
    for(G4int j=0;j<4;j++)
    {
      G4TwoVector u=GetVertex(j)+cf*(i+1)*( GetVertex(j+4)-GetVertex(j));
      polyhedron->AddVertex(G4ThreeVector(u.x(),u.y(),-fDz+cf*2*fDz*(i+1)));
    }    
  }
  for (i=4;i<8;i++)
  {
    polyhedron->AddVertex(G4ThreeVector(GetVertex(i).x(),
                                        GetVertex(i).y(),fDz));
  }

  // Add Facets
  //
  polyhedron->AddFacet(1,4,3,2);  //Z-plane
  for (i=0;i<subdivisions+1;i++)
  {
    G4int is=i*4;
    polyhedron->AddFacet(5+is,8+is,4+is,1+is);
    polyhedron->AddFacet(8+is,7+is,3+is,4+is);
    polyhedron->AddFacet(7+is,6+is,2+is,3+is);
    polyhedron->AddFacet(6+is,5+is,1+is,2+is); 
  }
  polyhedron->AddFacet(5+sub4,6+sub4,7+sub4,8+sub4);  //Z-plane

  polyhedron->SetReferences();
  polyhedron->InvertFacets();

  return (G4Polyhedron*) polyhedron;
}
