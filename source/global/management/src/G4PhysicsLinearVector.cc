// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PhysicsLinearVector.cc,v 1.4 2001-01-09 01:19:01 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//--------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//  G4PhysicsLinearVector.cc
//
//  History: first implementation, based on object model of
//    02 Dec. 1995, G.Cosmo : Structure created based on object model
//    15 Feb. 1996, K.Amako : Implemented the 1st version
//    01 Jul. 1996, K.Amako : Cache mechanism and hidden bin from the 
//                            user introduced.
//    26 Sep. 1996, K.Amako : Constructor with only 'bin size' added.
//    11 Nov. 2000, H.Kurashige : use g4std/vector for dataVector and binVector
//
//--------------------------------------------------------------------

#include "G4PhysicsLinearVector.hh"

G4PhysicsLinearVector::G4PhysicsLinearVector()
{
  ptrNextTable = 0;
  edgeMin = 0.0;
  edgeMax = 0.0;
  numberOfBin = 0;
}

G4PhysicsLinearVector::G4PhysicsLinearVector(size_t theNbin)
{

  // Add extra one bin (hidden to user) to handle correctly when 
  // Energy=theEmax in getValue. 
  dataVector.reserve(theNbin+1);
  binVector.reserve(theNbin+1);      

  ptrNextTable = 0;
  numberOfBin = theNbin;
  dBin = 0.;
  baseBin = 0.;

  edgeMin = 0.;
  edgeMax = 0.;

  lastBin = INT_MAX;
  lastEnergy = -DBL_MAX;
  lastValue = DBL_MAX;

  for (int i=0; i<=numberOfBin; i++) {
     binVector.push_back(0.0);
     dataVector.push_back(0.0);
  }
}  

G4PhysicsLinearVector::G4PhysicsLinearVector(G4double theEmin, 
                                             G4double theEmax, size_t theNbin)
{

  // Add extra one bin (hidden to user) to handle correctly when 
  // Energy=theEmax in getValue. 
  dataVector.reserve(theNbin+1);
  binVector.reserve(theNbin+1);      

  ptrNextTable = 0;
  numberOfBin = theNbin;
  dBin = (theEmax-theEmin) / numberOfBin;
  baseBin = theEmin/dBin;

  for (G4int i=0; i<numberOfBin+1; i++) {
    binVector.push_back( theEmin + i*dBin );
    dataVector.push_back(0.0);
  }
  binVector.push_back(0.0);
  dataVector.push_back(0.0);

  edgeMin = binVector[0];
  edgeMax = binVector[numberOfBin-1];

  lastBin = INT_MAX;
  lastEnergy = -DBL_MAX;
  lastValue = DBL_MAX;
}  

G4PhysicsLinearVector::~G4PhysicsLinearVector(){}
