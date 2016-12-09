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
// $Id: G4PhysicsLinearVector.cc 98864 2016-08-15 11:53:26Z gcosmo $
//
// 
//--------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//  G4PhysicsLinearVector.cc
//
//  15 Feb 1996 - K.Amako : 1st version
//  19 Jun 2009 - V.Ivanchenko : removed hidden bin 
//  02 Oct 2013 - V.Ivanchenko : removed FindBinLocation
//
//--------------------------------------------------------------------

#include "G4PhysicsLinearVector.hh"

G4PhysicsLinearVector::G4PhysicsLinearVector()
  : G4PhysicsVector()
{
  type = T_G4PhysicsLinearVector;
}

G4PhysicsLinearVector::G4PhysicsLinearVector(G4double theEmin, 
                                             G4double theEmax, size_t theNbin)
  : G4PhysicsVector()
{
  type = T_G4PhysicsLinearVector;

  dBin    = (theEmax-theEmin)/(G4double)theNbin;
  baseBin = theEmin/dBin;

  numberOfNodes = theNbin + 1;
  dataVector.reserve(numberOfNodes);
  binVector.reserve(numberOfNodes);      

  binVector.push_back(theEmin);
  dataVector.push_back(0.0);

  for (size_t i=1; i<numberOfNodes-1; ++i)
    {
      binVector.push_back( theEmin + i*dBin );
      dataVector.push_back(0.0);
    }
  binVector.push_back(theEmax);
  dataVector.push_back(0.0);

  edgeMin = binVector[0];
  edgeMax = binVector[numberOfNodes-1];

}  

G4PhysicsLinearVector::~G4PhysicsLinearVector()
{}

G4bool G4PhysicsLinearVector::Retrieve(std::ifstream& fIn, G4bool ascii)
{
  G4bool success = G4PhysicsVector::Retrieve(fIn, ascii);
  if (success)
  {
    dBin = binVector[1]-edgeMin;
    baseBin = edgeMin/dBin;
  }
  return success;
}

void G4PhysicsLinearVector::ScaleVector(G4double factorE, G4double factorV)
{
  G4PhysicsVector::ScaleVector(factorE, factorV);
  dBin = binVector[1]-edgeMin;
  baseBin = edgeMin/dBin;
}

