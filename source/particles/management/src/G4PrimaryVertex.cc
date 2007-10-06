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
// $Id: G4PrimaryVertex.cc,v 1.4 2007-10-06 06:49:29 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4PrimaryVertex.hh"
#include "G4VUserPrimaryVertexInformation.hh"
#include "G4ios.hh"

G4Allocator<G4PrimaryVertex> aPrimaryVertexAllocator;

G4PrimaryVertex::G4PrimaryVertex()
:X0(0.),Y0(0.),Z0(0.),T0(0.),theParticle(0),theTail(0),
 nextVertex(0),tailVertex(0),numberOfParticle(0),Weight0(1.0),userInfo(0)
{;}

G4PrimaryVertex::G4PrimaryVertex(
          G4double x0,G4double y0,G4double z0,G4double t0)
:X0(x0),Y0(y0),Z0(z0),T0(t0),theParticle(0),theTail(0),
 nextVertex(0),tailVertex(0),numberOfParticle(0),Weight0(1.0),userInfo(0)
{;}

G4PrimaryVertex::G4PrimaryVertex(G4ThreeVector xyz0,G4double t0)
:T0(t0),theParticle(0),theTail(0),
 nextVertex(0),tailVertex(0),numberOfParticle(0),Weight0(1.0),userInfo(0)
{
  X0=xyz0.x();
  Y0=xyz0.y();
  Z0=xyz0.z();
}

G4PrimaryVertex::~G4PrimaryVertex()
{
  if(theParticle != 0)
  { delete theParticle; }
  if(nextVertex != 0)
  { delete nextVertex; }
  if(userInfo != 0)
  { delete userInfo; }
}

const G4PrimaryVertex & 
G4PrimaryVertex::operator=(const G4PrimaryVertex &)
{ return *this; }
G4int G4PrimaryVertex::operator==(const G4PrimaryVertex &right) const
{ return (this==&right); }
G4int G4PrimaryVertex::operator!=(const G4PrimaryVertex &right) const
{ return (this!=&right); }

void G4PrimaryVertex::Print() const
{
  G4cout << "Vertex  ( "
	 << X0/mm  << "[mm], " 
	 << Y0/mm << "[mm], " 
	 << Z0/mm << "[mm], " 
	 << T0/ns  << "[ns] )" 
       << " Weight " << Weight0 << G4endl;
  if(userInfo!=0) userInfo->Print();
  G4cout << "#### Primary particles" << G4endl;
  G4PrimaryParticle* aPrim = theParticle;
  if(aPrim != 0)
  {
    aPrim->Print();
    aPrim = aPrim->GetNext();
  }
}
