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
// G4PhysicsLogVector class implementation
//
// Authors:
// - 02 Dec. 1995, G.Cosmo: Structure created based on object model
// - 03 Mar. 1996, K.Amako: Implemented the 1st version
// Revisions:
// - 11 Nov. 2000, H.Kurashige : Use STL vector for dataVector and binVector
// - 19 Jun. 2009, V.Ivanchenko : removed hidden bin
// --------------------------------------------------------------------

#include "G4PhysicsLogVector.hh"
#include "G4Exp.hh"

// --------------------------------------------------------------------
G4PhysicsLogVector::G4PhysicsLogVector()
  : G4PhysicsVector()
{
  type = T_G4PhysicsLogVector;
}

// --------------------------------------------------------------------
G4PhysicsLogVector::G4PhysicsLogVector(G4double theEmin, G4double theEmax,
                                       std::size_t theNbin)
  : G4PhysicsVector()
{
  numberOfNodes = theNbin + 1;
  if(theNbin < 2 || theEmin == theEmax)
  {
    G4ExceptionDescription ed;
    ed << "G4PhysicsLogVector with wrong parameters: theNbin= " << theNbin
       << " theEmin= " << theEmin << " theEmax= " << theEmax;
    G4Exception("G4PhysicsLogVector::G4PhysicsLogVector()", "glob03",
                FatalException, ed, "theNbins should be > 2");
  }
  if(numberOfNodes < 3)
  {
    numberOfNodes = 3;
  }
  type = T_G4PhysicsLogVector;

  invdBin = 1. / (G4Log(theEmax / theEmin) / (G4double) (numberOfNodes - 1));
  baseBin = G4Log(theEmin) * invdBin;

  dataVector.reserve(numberOfNodes);
  binVector.reserve(numberOfNodes);

  binVector.push_back(theEmin);
  dataVector.push_back(0.0);

  for(std::size_t i = 1; i < numberOfNodes - 1; ++i)
  {
    binVector.push_back(G4Exp((baseBin + i) / invdBin));
    dataVector.push_back(0.0);
  }
  binVector.push_back(theEmax);
  dataVector.push_back(0.0);

  edgeMin = binVector[0];
  edgeMax = binVector[numberOfNodes - 1];
}

// --------------------------------------------------------------------
G4PhysicsLogVector::~G4PhysicsLogVector() {}

// --------------------------------------------------------------------
G4bool G4PhysicsLogVector::Retrieve(std::ifstream& fIn, G4bool ascii)
{
  G4bool success = G4PhysicsVector::Retrieve(fIn, ascii);
  if(success)
  {
    invdBin = 1. / G4Log(binVector[1] / edgeMin);
    baseBin = G4Log(edgeMin) * invdBin;
  }
  return success;
}

// --------------------------------------------------------------------
void G4PhysicsLogVector::ScaleVector(G4double factorE, G4double factorV)
{
  G4PhysicsVector::ScaleVector(factorE, factorV);
  invdBin = 1. / G4Log(binVector[1] / edgeMin);
  baseBin = G4Log(edgeMin) * invdBin;
}
