// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PhysicsLogVector.cc,v 1.6 2001-02-02 16:23:37 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 class implementation file
//
//  G4PhysicsLogVector.cc
//
//  History:
//    02 Dec. 1995, G.Cosmo : Structure created based on object model
//    15 Feb. 1996, K.Amako : Implemented the 1st version
//    01 Jul. 1996, K.Amako : Hidden bin from the user introduced
//    26 Sep. 1996, K.Amako : Constructor with only 'bin size' added
//    11 Nov. 2000, H.Kurashige : use STL vector for dataVector and binVector
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
  dataVector.reserve(theNbin+1);
  binVector.reserve(theNbin+1); 

  numberOfBin = theNbin;
  dBin = 0.;
  baseBin = 0.;

  edgeMin = 0.;
  edgeMax = 0.;

  lastBin = INT_MAX;
  lastEnergy = -DBL_MAX;
  lastValue = DBL_MAX;
 
  for (size_t i=0; i<=numberOfBin; i++) {
     binVector.push_back(0.0);
     dataVector.push_back(0.0);
  }
}  


G4PhysicsLogVector::G4PhysicsLogVector(G4double theEmin, 
                                       G4double theEmax, size_t theNbin)
{

  // Add extra one bin (hidden to user) to handle correctly when 
  // Energy=theEmax in getValue. 
  dataVector.reserve(theNbin+1);
  binVector.reserve(theNbin+1); 

  numberOfBin = theNbin;
  dBin = log10(theEmax/theEmin) / numberOfBin;
  baseBin = log10(theEmin)/dBin;

  for (size_t i=0; i<numberOfBin+1; i++) {
    binVector.push_back(pow(10., log10(theEmin)+i*dBin));
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


G4PhysicsLogVector::~G4PhysicsLogVector(){}
