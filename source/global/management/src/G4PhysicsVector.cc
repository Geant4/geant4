// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PhysicsVector.cc,v 1.2 1999-11-11 10:47:30 gunter Exp $
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
//
// --------------------------------------------------------------

#include "G4PhysicsVector.hh"


G4int G4PhysicsVector::operator==(const G4PhysicsVector &right) const
{
  return (this == (G4PhysicsVector *) &right);
}

G4int G4PhysicsVector::operator!=(const G4PhysicsVector &right) const
{
  return (this != (G4PhysicsVector *) &right);
}

void G4PhysicsVector::PutValue(size_t binNumber, G4double theValue)
{
  dataVector(binNumber) = theValue;

  // Fill the bin which is hidden to user with theValue. This is to 
  // handle correctly when Energy=theEmax in getValue.
  if(binNumber=numberOfBin-1) {
    dataVector(binNumber+1) = theValue;
  }                                 
}

G4double G4PhysicsVector::GetLowEdgeEnergy(size_t binNumber) const
{
  return binVector(binNumber);
}

size_t G4PhysicsVector::GetVectorLength() const
{
  return numberOfBin;
}

G4bool G4PhysicsVector::IsFilledVectorExist() const
{
  G4bool status=false;

  if(numberOfBin > 0) status=true;
  return status;
}

void G4PhysicsVector::LinkPhysicsTable(G4RWTPtrOrderedVector<G4PhysicsVector>& theTable)
{
  ptrNextTable = &theTable;
} 

G4bool G4PhysicsVector::IsLinkedTableExist() const
{
  G4bool status=false;

  if(ptrNextTable != 0) status=true;
  return status;
}

void G4PhysicsVector::PutComment(const G4String& theComment)
{
  comment = theComment;
}

G4String G4PhysicsVector::GetComment() const
{
  return comment;
}












