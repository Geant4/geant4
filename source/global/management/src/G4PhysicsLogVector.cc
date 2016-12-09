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
// $Id: G4PhysicsLogVector.cc 98864 2016-08-15 11:53:26Z gcosmo $
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
//    9  Mar. 2001, H.Kurashige : added PhysicsVector type and Retrieve
//    05 Sep. 2008, V.Ivanchenko : added protections for zero-length vector
//    19 Jun. 2009, V.Ivanchenko : removed hidden bin 
//    02 Oct. 2013  V.Ivanchenko : Remove FindBinLocation method, use G4Log
//
// --------------------------------------------------------------

#include "G4PhysicsLogVector.hh"
#include "G4Exp.hh"

G4PhysicsLogVector::G4PhysicsLogVector()
  : G4PhysicsVector()
{ 
  type = T_G4PhysicsLogVector;
}

G4PhysicsLogVector::G4PhysicsLogVector(G4double theEmin, 
                                       G4double theEmax, size_t theNbin)
  : G4PhysicsVector()
{
  type = T_G4PhysicsLogVector;

  dBin     =  G4Log(theEmax/theEmin)/(G4double)theNbin;
  baseBin  =  G4Log(theEmin)/dBin;

  numberOfNodes = theNbin + 1;
  dataVector.reserve(numberOfNodes);
  binVector.reserve(numberOfNodes);

  binVector.push_back(theEmin);
  dataVector.push_back(0.0);

  for (size_t i=1; i<numberOfNodes-1; ++i)
    {
      binVector.push_back(G4Exp((baseBin+i)*dBin));
      dataVector.push_back(0.0);
    }
  binVector.push_back(theEmax);
  dataVector.push_back(0.0);

  edgeMin = binVector[0];
  edgeMax = binVector[numberOfNodes-1];
}  

G4PhysicsLogVector::~G4PhysicsLogVector()
{}

G4bool G4PhysicsLogVector::Retrieve(std::ifstream& fIn, G4bool ascii)
{
  G4bool success = G4PhysicsVector::Retrieve(fIn, ascii);
  if (success)
  {
    dBin = G4Log(binVector[1]/edgeMin);
    baseBin = G4Log(edgeMin)/dBin;
  }
  return success;
}

void G4PhysicsLogVector::ScaleVector(G4double factorE, G4double factorV)
{
  G4PhysicsVector::ScaleVector(factorE, factorV);
  dBin = G4Log(binVector[1]/edgeMin);
  baseBin = G4Log(edgeMin)/dBin;
}
