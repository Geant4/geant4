//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id$
//
// 
// --------------------------------------------------------------
//      GEANT 4 class implementation file
//
//  G4PhysicsLnVector.cc
//
//  27 Apr 1999 - M.G.Pia: Created from G4PhysicsLogVector
//  19 Jun 2009 - V.Ivanchenko : removed hidden bin 
//
// --------------------------------------------------------------

#include "G4PhysicsLnVector.hh"

G4PhysicsLnVector::G4PhysicsLnVector()
  : G4PhysicsVector()
{
  type = T_G4PhysicsLnVector;
}

G4PhysicsLnVector::G4PhysicsLnVector(size_t theNbin)
  : G4PhysicsVector()
{
  type = T_G4PhysicsLnVector;

  numberOfNodes = theNbin + 1;
  dataVector.reserve(numberOfNodes);
  binVector.reserve(numberOfNodes);      

  for (size_t i=0; i<numberOfNodes; i++)
  {
     binVector.push_back(0.0);
     dataVector.push_back(0.0);
  }
}  

G4PhysicsLnVector::G4PhysicsLnVector(G4double theEmin, 
                                     G4double theEmax, size_t theNbin)
  : G4PhysicsVector()
{
  type = T_G4PhysicsLnVector;

  dBin    = std::log(theEmax/theEmin)/theNbin;
  baseBin = std::log(theEmin)/dBin;
  numberOfNodes = theNbin + 1;
  dataVector.reserve(numberOfNodes);
  binVector.reserve(numberOfNodes);      

  binVector.push_back(theEmin);
  dataVector.push_back(0.0);

  for (size_t i=1; i<numberOfNodes-1; i++)
    {
      binVector.push_back(std::exp((baseBin+i)*dBin));
      dataVector.push_back(0.0);
    }
  binVector.push_back(theEmax);
  dataVector.push_back(0.0);

  edgeMin = binVector[0];
  edgeMax = binVector[numberOfNodes-1];
}

G4PhysicsLnVector::~G4PhysicsLnVector(){}

G4bool G4PhysicsLnVector::Retrieve(std::ifstream& fIn, G4bool ascii)
{
  G4bool success = G4PhysicsVector::Retrieve(fIn, ascii);
  if (success)
  {
    G4double theEmin = binVector[0];
    dBin = std::log(binVector[1]/theEmin);
    baseBin = std::log(theEmin)/dBin;
  }
  return success;
}

void G4PhysicsLnVector::ScaleVector(G4double factorE, G4double factorV)
{
  G4PhysicsVector::ScaleVector(factorE, factorV);
  G4double theEmin = binVector[0];
  dBin = std::log(binVector[1]/theEmin);
  baseBin = std::log(theEmin)/dBin;
}

size_t G4PhysicsLnVector::FindBinLocation(G4double theEnergy) const 
{
 
  // For G4PhysicsLnVector, FindBinLocation is implemented using
  // a simple arithmetic calculation.
  //
  // Because this is a virtual function, it is accessed through a
  // pointer to the G4PhyiscsVector object for most usages. In this
  // case, 'inline' will not be invoked. However, there is a possibility 
  // that the user access to the G4PhysicsLnVector object directly and 
  // not through pointers or references. In this case, the 'inline' will
  // be invoked. (See R.B.Murray, "C++ Strategies and Tactics", Chap.6.6)

  return size_t( std::log(theEnergy)/dBin - baseBin );
}
