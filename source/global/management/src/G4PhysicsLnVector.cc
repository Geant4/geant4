// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//
// 
// --------------------------------------------------------------
//      GEANT 4 class implementation file
//
//  G4PhysicsLnVector.cc
//
//  History: first implementation, based on object model of
//    02 Dec. 1995, G.Cosmo : Structure created based on object model
//    15 Feb. 1996, K.Amako : Implemented the 1st version
//    01 Jul. 1996, K.Amako : Hidden bin from the user introduced
//    26 Sep. 1996, K.Amako : Constructor with only 'bin size' added.
//
// --------------------------------------------------------------

#include "G4PhysicsLnVector.hh"


G4PhysicsLnVector::G4PhysicsLnVector()
{
  ptrNextTable = 0;
  edgeMin = 0.0;
  edgeMax = 0.0;
  numberOfBin = 0;
}


G4PhysicsLnVector::G4PhysicsLnVector(size_t theNbin)
{

  // Add extra one bin (hidden to user) to handle correctly when 
  // Energy=theEmax in getValue. 
  dataVector.resize(theNbin+1);
  binVector.resize(theNbin+1); 

  ptrNextTable = 0;
  numberOfBin = theNbin;
  dBin = 0.;
  baseBin = 0.;

  edgeMin = 0.;
  edgeMax = 0.;

  lastBin = INT_MAX;
  lastEnergy = -DBL_MAX;
  lastValue = DBL_MAX;
}  


G4PhysicsLnVector::G4PhysicsLnVector(G4double theEmin, 
                                       G4double theEmax, size_t theNbin)
{

  // Add extra one bin (hidden to user) to handle correctly when 
  // Energy=theEmax in getValue. 
  dataVector.resize(theNbin+1);
  binVector.resize(theNbin+1); 

  ptrNextTable = 0;
  numberOfBin = theNbin;
  dBin = log(theEmax/theEmin) / numberOfBin;
  baseBin = log(theEmin)/dBin;

  for (G4int i=0; i<numberOfBin+1; i++) {
    binVector(i) = exp(log(theEmin)+i*dBin);
  }

  edgeMin = binVector(0);
  edgeMax = binVector(numberOfBin-1);

  lastBin = INT_MAX;
  lastEnergy = -DBL_MAX;
  lastValue = DBL_MAX;
}  


G4PhysicsLnVector::~G4PhysicsLnVector(){}






