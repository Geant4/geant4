// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PhysicsLogVector.cc,v 1.3 2000-11-20 17:26:48 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 class implementation file
//
//  G4PhysicsLogVector.cc
//
//  History: first implementation, based on object model of
//    02 Dec. 1995, G.Cosmo : Structure created based on object model
//    15 Feb. 1996, K.Amako : Implemented the 1st version
//    01 Jul. 1996, K.Amako : Hidden bin from the user introduced
//    26 Sep. 1996, K.Amako : Constructor with only 'bin size' added.
//
// --------------------------------------------------------------

#include "G4PhysicsLogVector.hh"


G4PhysicsLogVector::G4PhysicsLogVector()
  : dBin(0.), baseBin(0.)
{}


G4PhysicsLogVector::G4PhysicsLogVector(size_t theNbin)
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


G4PhysicsLogVector::G4PhysicsLogVector(G4double theEmin, 
                                       G4double theEmax, size_t theNbin)
{

  // Add extra one bin (hidden to user) to handle correctly when 
  // Energy=theEmax in getValue. 
  dataVector.resize(theNbin+1);
  binVector.resize(theNbin+1); 

  ptrNextTable = 0;
  numberOfBin = theNbin;
  dBin = log10(theEmax/theEmin) / numberOfBin;
  baseBin = log10(theEmin)/dBin;

  for (size_t i=0; i<numberOfBin+1; i++) {
    binVector(i) = pow(10., log10(theEmin)+i*dBin);
  }

  edgeMin = binVector(0);
  edgeMax = binVector(numberOfBin-1);

  lastBin = INT_MAX;
  lastEnergy = -DBL_MAX;
  lastValue = DBL_MAX;
}  


G4PhysicsLogVector::~G4PhysicsLogVector(){}
