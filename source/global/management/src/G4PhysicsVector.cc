// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PhysicsVector.cc,v 1.8 2001-02-02 16:23:37 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 class implementation file
//
//  G4PhysicsVector.cc
//
//  History:
//    02 Dec. 1995, G.Cosmo : Structure created based on object model
//    03 Mar. 1996, K.Amako : Implemented the 1st version
//    01 Jul. 1996, K.Amako : Hidden bin from the user introduced
//    12 Nov. 1998, K.Amako : A bug in GetVectorLength() fixed
//    11 Nov. 2000, H.Kurashige : use STL vector for dataVector and binVector
//    18 Jan. 2001, H.Kurashige : remove ptrNextTable
// --------------------------------------------------------------

#include "G4PhysicsVector.hh"

G4PhysicsVector::G4PhysicsVector()
 : edgeMin(0.), edgeMax(0.), numberOfBin(0),
   lastEnergy(0.), lastValue(0.), lastBin(0)
{}

G4PhysicsVector::~G4PhysicsVector() 
{
  dataVector.clear();
  binVector.clear();
}

G4PhysicsVector::G4PhysicsVector(const G4PhysicsVector& right)
 : edgeMin(right.edgeMin), edgeMax(right.edgeMax),
   numberOfBin(right.numberOfBin), lastEnergy(right.lastEnergy),
   lastValue(right.lastValue), lastBin(right.lastBin),
   dataVector(right.dataVector), binVector(right.binVector),
   comment(right.comment)
{}

G4PhysicsVector&
G4PhysicsVector::operator=(const G4PhysicsVector& right)
{
  if (&right==this) return *this;

  edgeMin = right.edgeMin;
  edgeMax = right.edgeMax;
  numberOfBin = right.numberOfBin;
  lastEnergy = right.lastEnergy;
  lastValue = right.lastValue;
  lastBin = right.lastBin;
  dataVector = right.dataVector;
  binVector = right.binVector;
  comment = right.comment;

  return *this;
}

G4int G4PhysicsVector::operator==(const G4PhysicsVector &right) const
{
  return (this == &right);
}

G4int G4PhysicsVector::operator!=(const G4PhysicsVector &right) const
{
  return (this != &right);
}


G4double G4PhysicsVector::GetLowEdgeEnergy(size_t binNumber) const
{
  return binVector[binNumber];
}
