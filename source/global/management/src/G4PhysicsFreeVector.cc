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
// $Id: G4PhysicsFreeVector.cc 98864 2016-08-15 11:53:26Z gcosmo $
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
//    19 Jun. 2009, V.Ivanchenko : removed hidden bin 
//
//--------------------------------------------------------------------

#include "G4PhysicsFreeVector.hh"

G4PhysicsFreeVector::G4PhysicsFreeVector() 
  : G4PhysicsVector()
{
  type = T_G4PhysicsFreeVector;
}

G4PhysicsFreeVector::G4PhysicsFreeVector(size_t length)
  : G4PhysicsVector()
{
  type = T_G4PhysicsFreeVector;
  numberOfNodes = length;

  dataVector.reserve(numberOfNodes);
  binVector.reserve(numberOfNodes);

  for (size_t i=0; i<numberOfNodes; ++i)
  {
     binVector.push_back(0.0);
     dataVector.push_back(0.0);
  }
}  

G4PhysicsFreeVector::G4PhysicsFreeVector(const G4DataVector& theBinVector, 
                                         const G4DataVector& theDataVector)
  : G4PhysicsVector()
{
  type = T_G4PhysicsFreeVector;
  numberOfNodes = theBinVector.size();

  dataVector.reserve(numberOfNodes);
  binVector.reserve(numberOfNodes);

  for (size_t i=0; i<numberOfNodes; ++i)
  {
    binVector.push_back(theBinVector[i]);
    dataVector.push_back(theDataVector[i]);
  }
  if(numberOfNodes > 0) 
  {
    edgeMin = binVector[0];
    edgeMax = binVector[numberOfNodes-1];
  }
}  

G4PhysicsFreeVector::~G4PhysicsFreeVector()
{}

