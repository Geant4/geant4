// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PhysicsFreeVector.cc,v 1.7 2001-03-09 03:39:30 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//--------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//  G4PhysicsFreeVector.cc
//
//  History:
//    02 Dec. 1995, G.Cosmo : Structure created based on object model
//    06 June 1996, K.Amako : The 1st version of implemented
//    01 Jul. 1996, K.Amako : Cache mechanism and hidden bin from the 
//                            user introduced
//    26 Sep. 1996, K.Amako : Constructor with only 'bin size' added
//    11 Nov. 2000, H.Kurashige : use STL vector for dataVector and binVector
//
//--------------------------------------------------------------------

#include "G4PhysicsFreeVector.hh"


G4PhysicsFreeVector::G4PhysicsFreeVector()
{
  edgeMin = 0.0;
  edgeMax = 0.0;
  numberOfBin = 0;
  type = T_G4PhysicsFreeVector;
}


G4PhysicsFreeVector::G4PhysicsFreeVector(size_t theNbin)
{
  type = T_G4PhysicsFreeVector;
  numberOfBin = theNbin;

  // Add extra one bin (hidden to user) to handle correctly when 
  // Energy=theEmax in getValue.
  dataVector.reserve(numberOfBin+1);
  binVector.reserve(numberOfBin+1);

  for (size_t i=0; i<=numberOfBin; i++) {
     binVector.push_back(0.0);
     dataVector.push_back(0.0);
  }

  edgeMin = 0.;
  edgeMax = 0.;

  lastBin = INT_MAX;
  lastEnergy = -DBL_MAX;
  lastValue = DBL_MAX;

}  


G4PhysicsFreeVector::G4PhysicsFreeVector(const G4DataVector& theBinVector, 
                                         const G4DataVector& theDataVector)
{
  type = T_G4PhysicsFreeVector;
  numberOfBin = theBinVector.size();

  // Add extra one bin (hidden to user) to handle correctly when 
  // Energy=theEmax in getValue.
  dataVector.reserve(numberOfBin+1);
  binVector.reserve(numberOfBin+1);

  for (size_t i=0; i<numberOfBin; i++) {
     binVector.push_back(theBinVector[i]);
     dataVector.push_back(theDataVector[i]);
  }

  // Put values to extra hidden bin. For 'binVector', the 'edgeMin' of the
  // extra hidden bin is assumed to have the following value. For binary
  // search, this value is completely arbitrary if it is greater than
  // the 'edgeMin' at 'numberOfBin-1'. 
  binVector.push_back ( theBinVector[numberOfBin-1] + 1.0 );


  // Put values to extra hidden bin. For 'dataVector', the 'value' of the
  // extra hidden bin is assumed to have the same as the one at 'numberBin-1'. 
  dataVector.push_back( theDataVector[numberOfBin-1] );

  edgeMin = binVector[0];
  edgeMax = binVector[numberOfBin-1];

  lastBin = INT_MAX;
  lastEnergy = -DBL_MAX;
  lastValue = DBL_MAX;

}  


G4PhysicsFreeVector::~G4PhysicsFreeVector(){}


void G4PhysicsFreeVector::PutValue( size_t theBinNumber, G4double theBinValue, 
			            G4double theDataValue )
{
  binVector[theBinNumber]  = theBinValue;
  dataVector[theBinNumber] = theDataValue;


  if( theBinNumber == numberOfBin-1 ) {
     edgeMax = binVector[numberOfBin-1];

     // Put values to extra hidden bin. For 'binVector', the 'edgeMin' 
     // of the extra hidden bin is assumed to have the following value. 
     // For binary search, this value is completely arbitrary if it is 
     // greater than the 'edgeMin' at 'numberOfBin-1'. 
     binVector[numberOfBin] = theBinValue + 1.0;

     // Put values to extra hidden bin. For 'dataVector', the 'value' 
     // of the extra hidden bin is assumed to have the same as the one 
     // at 'numberBin-1'. 
     dataVector[numberOfBin] = theDataValue;
  }

  if( theBinNumber == 0 ) {
     edgeMin = binVector[0];
  }
}
