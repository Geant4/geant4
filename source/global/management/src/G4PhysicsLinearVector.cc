// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PhysicsLinearVector.cc,v 1.7 2001-03-09 03:39:30 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//--------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//  G4PhysicsLinearVector.cc
//
//  History:
//    02 Dec. 1995, G.Cosmo : Structure created based on object model
//    15 Feb. 1996, K.Amako : Implemented the 1st version
//    01 Jul. 1996, K.Amako : Cache mechanism and hidden bin from the 
//                            user introduced.
//    26 Sep. 1996, K.Amako : Constructor with only 'bin size' added.
//    11 Nov. 2000, H.Kurashige : use STL vector for dataVector and binVector
//    9  Mar. 2001, H.Kurashige : add PhysicsVector type and Retrieve
//
//--------------------------------------------------------------------

#include "G4PhysicsLinearVector.hh"

G4PhysicsLinearVector::G4PhysicsLinearVector()
{
  edgeMin = 0.0;
  edgeMax = 0.0;
  numberOfBin = 0;
  type = T_G4PhysicsLinearVector;
}

G4PhysicsLinearVector::G4PhysicsLinearVector(size_t theNbin)
{
  type = T_G4PhysicsLinearVector;

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

G4PhysicsLinearVector::G4PhysicsLinearVector(G4double theEmin, 
                                             G4double theEmax, size_t theNbin)
{
  type = T_G4PhysicsLinearVector;

  // Add extra one bin (hidden to user) to handle correctly when 
  // Energy=theEmax in getValue. 
  dataVector.reserve(theNbin+1);
  binVector.reserve(theNbin+1);      

  numberOfBin = theNbin;
  dBin = (theEmax-theEmin) / numberOfBin;
  baseBin = theEmin/dBin;

  for (size_t i=0; i<numberOfBin+1; i++) {
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


G4bool G4PhysicsLinearVector::Retrieve(G4std::ifstream& fIn, G4bool ascii)
{
  G4bool success = G4PhysicsVector::Retrieve(fIn, ascii);
  if (success){
    G4double theEmin = binVector[0];
    dBin = binVector[1]-theEmin;
    baseBin = theEmin/dBin;
  }
  return success;
}

