// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PhysicsLnVector.cc,v 1.9 2001-03-10 04:59:15 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 class implementation file
//
//  G4PhysicsLnVector.cc
//
//  History:
//    27 Apr. 1999, M.G. Pia: Created, copying from G4PhysicsLogVector
//    11 Nov. 2000, H.Kurashige : use STL vector for dataVector and binVector
//    9  Mar. 2001, H.Kurashige : add PhysicsVector type and Retrieve
//
// --------------------------------------------------------------

#include "G4PhysicsLnVector.hh"

G4PhysicsLnVector::G4PhysicsLnVector()
  : dBin(0.), baseBin(0.)
{
  type = T_G4PhysicsLnVector;
}


G4PhysicsLnVector::G4PhysicsLnVector(size_t theNbin)
{
  type = T_G4PhysicsLnVector;

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


G4PhysicsLnVector::G4PhysicsLnVector(G4double theEmin, 
                                       G4double theEmax, size_t theNbin)
{
  type = T_G4PhysicsLnVector;

  // Add extra one bin (hidden to user) to handle correctly when 
  // Energy=theEmax in getValue. 
  dataVector.reserve(theNbin+1);
  binVector.reserve(theNbin+1); 

  numberOfBin = theNbin;
  dBin = log(theEmax/theEmin) / numberOfBin;
  baseBin = log(theEmin)/dBin;

  for (size_t i=0; i<numberOfBin+1; i++) {
    binVector.push_back(exp(log(theEmin)+i*dBin));
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


G4PhysicsLnVector::~G4PhysicsLnVector(){}

G4bool G4PhysicsLnVector::Retrieve(G4std::ifstream& fIn, G4bool ascii)
{
  G4bool success = G4PhysicsVector::Retrieve(fIn, ascii);
  if (success){
    G4double theEmin = binVector[0];
    dBin = log(binVector[1]/theEmin);
    baseBin = log(theEmin)/dBin;
  }
  return success;
}

