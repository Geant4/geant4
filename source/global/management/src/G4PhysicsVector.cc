// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PhysicsVector.cc,v 1.6 2001-01-09 01:19:04 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 class implementation file
//
//  G4PhysicsVector.cc
//
//  History: first implementation, based on object model of
//    02 Dec. 1995, G.Cosmo : Structure created based on object model
//    03 Mar. 1996, K.Amako : Implemented the 1st version
//    01 Jul. 1996, K.Amako : Hidden bin from the user introduced
//    12 Nov. 1998, K.Amako : A bug in GetVectorLength() fixed
//    11 Nov. 2000, H.Kurashige : use g4std/vector for dataVector and binVector
//
// --------------------------------------------------------------

#include "G4PhysicsVector.hh"

G4PhysicsVector::G4PhysicsVector()
 : edgeMin(0.), edgeMax(0.), numberOfBin(0),
   lastEnergy(0.), lastValue(0.), lastBin(0),
   ptrNextTable(0) 
{}

G4PhysicsVector::~G4PhysicsVector() 
{
  dataVector.clear();
  binVector.clear();
}

G4int G4PhysicsVector::operator==(const G4PhysicsVector &right) const
{
  return (this == (G4PhysicsVector *) &right);
}

G4int G4PhysicsVector::operator!=(const G4PhysicsVector &right) const
{
  return (this != (G4PhysicsVector *) &right);
}


G4double G4PhysicsVector::GetLowEdgeEnergy(size_t binNumber) const
{
  return binVector[binNumber];
}













