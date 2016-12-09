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
// $Id: G4PrimaryVertex.cc 99159 2016-09-07 08:11:50Z gcosmo $
//

#include "G4PrimaryVertex.hh"
#include "G4SystemOfUnits.hh"
#include "G4VUserPrimaryVertexInformation.hh"
#include "G4ios.hh"

G4ThreadLocal G4Allocator<G4PrimaryVertex> *aPrimaryVertexAllocator = 0;

G4PrimaryVertex::G4PrimaryVertex()
:X0(0.),Y0(0.),Z0(0.),T0(0.),theParticle(0),theTail(0),
 nextVertex(0),tailVertex(0),numberOfParticle(0),Weight0(1.0),userInfo(0)
{
}

G4PrimaryVertex::G4PrimaryVertex(
          G4double x0,G4double y0,G4double z0,G4double t0)
:X0(x0),Y0(y0),Z0(z0),T0(t0),theParticle(0),theTail(0),
 nextVertex(0),tailVertex(0),numberOfParticle(0),Weight0(1.0),userInfo(0)
{
}

G4PrimaryVertex::G4PrimaryVertex(G4ThreeVector xyz0,G4double t0)
:T0(t0),theParticle(0),theTail(0),
 nextVertex(0),tailVertex(0),numberOfParticle(0),Weight0(1.0),userInfo(0)
{
  X0=xyz0.x();
  Y0=xyz0.y();
  Z0=xyz0.z();
}

G4PrimaryVertex::G4PrimaryVertex(const G4PrimaryVertex & right)
:theParticle(0),theTail(0),
 nextVertex(0),tailVertex(0),userInfo(0)
{
  numberOfParticle = right.numberOfParticle;
  *this = right;
}

G4PrimaryVertex::~G4PrimaryVertex()
{
  if(theParticle != 0) {
    G4PrimaryParticle* theNext = theParticle;
    while(theNext)
    {
      G4PrimaryParticle* thisPrimary = theNext;
      theNext = thisPrimary->GetNext();
      thisPrimary->ClearNext();
      delete thisPrimary;
    }
    theParticle = 0;
  }
  if(nextVertex != 0) { 
    delete nextVertex; 
    nextVertex =0;
  }
  if(userInfo != 0) { 
    delete userInfo; 
    userInfo = 0;
  }
}

G4PrimaryVertex &  G4PrimaryVertex::operator=(const G4PrimaryVertex & right)
{  
  if (this != &right) {
    X0       = right.X0;
    Y0       = right.Y0;
    Z0       = right.Z0;
    T0       = right.T0;
    Weight0  = right.Weight0;

    numberOfParticle = 0;
    if (theParticle !=0) delete theParticle;
    theParticle =0;
    theTail     =0;
    if (right.theParticle !=0 ) {
      theParticle = new G4PrimaryParticle(*(right.theParticle));
      numberOfParticle += 1;
      theTail = theParticle;
      G4PrimaryParticle * np = theParticle->GetNext();
      while (np !=0) { // Loop checking, 09.08.2015, K.Kurashige
	numberOfParticle += 1;
	theTail = np;
	np = np->GetNext();
      }
    }
    
    if (nextVertex !=0 ) delete nextVertex;
    nextVertex = 0;
    tailVertex =0;
    if (right.nextVertex !=0 ) {
      nextVertex = new G4PrimaryVertex(*(right.nextVertex));
      tailVertex = nextVertex;
      G4PrimaryVertex* nv = nextVertex->GetNext();
      while (nv !=0) { // Loop checking, 09.08.2015, K.Kurashige
	tailVertex = nv;
	nv = nv->GetNext();
      }
    }

    // userInfo can not be copied
    userInfo = 0;
  }
  return *this; 
}

G4int G4PrimaryVertex::operator==(const G4PrimaryVertex &right) const
{ return (this==&right); }

G4int G4PrimaryVertex::operator!=(const G4PrimaryVertex &right) const
{ return (this!=&right); }

G4PrimaryParticle* G4PrimaryVertex::GetPrimary(G4int i) const
{  
  if( i >= 0 && i < numberOfParticle ) {
    G4PrimaryParticle* particle = theParticle;
    for( G4int j=0; j<i; j++ ){ 
      if( particle == 0 ) return 0;
      particle = particle->GetNext();
    }
    return particle;
  } else { 
    return 0; }
}

void G4PrimaryVertex::Print() const
{ 
  G4cout << "Vertex  ( "
	 << X0/mm  << "[mm], " 
	 << Y0/mm << "[mm], " 
	 << Z0/mm << "[mm], " 
	 << T0/ns  << "[ns] )" 
	 << " Weight " << Weight0 << G4endl;
  if(userInfo!=0) userInfo->Print();
  G4cout << "  -- Primary particles :: " 
	 << "   # of primaries =" << numberOfParticle << G4endl;
  if( theParticle != 0)  theParticle->Print();
  if (nextVertex !=0 ) {
    G4cout << "Next Vertex " << G4endl;
    nextVertex->Print();
  }
}
