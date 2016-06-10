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
// $Id: G4LPhysicsFreeVector.cc 74256 2013-10-02 14:24:02Z gcosmo $
//
// 
// --------------------------------------------------------------------
// Class G4LPhysicsFreeVector
// Derived from base class G4PhysicsVector
// This is a free vector for Low Energy Physics cross section data
//
// F.W. Jones, TRIUMF, 04-JUN-96
// 
// Modified:
//    19 Jun. 2009, V.Ivanchenko : removed hidden bin 
//    19 Jun. 2009, V.Ivanchenko : removed FindBinLocation 
//
// --------------------------------------------------------------------

#include "G4LPhysicsFreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

// --------------------------------------------------------------------

G4LPhysicsFreeVector::G4LPhysicsFreeVector()
   : G4PhysicsVector()
{
   type = T_G4LPhysicsFreeVector;
}

// --------------------------------------------------------------------

G4LPhysicsFreeVector::G4LPhysicsFreeVector(size_t nbin,
                                           G4double binmin,
                                           G4double binmax)
  : G4PhysicsVector()
{
   type = T_G4LPhysicsFreeVector;
   
   edgeMin = binmin;
   edgeMax = binmax;
   numberOfNodes = nbin;
   binVector.reserve(numberOfNodes);
   dataVector.reserve(numberOfNodes);
   for (size_t i=0; i<numberOfNodes; i++)
   {
     binVector.push_back(0.0);
     dataVector.push_back(0.0);
   }
}  

// --------------------------------------------------------------------

G4LPhysicsFreeVector::~G4LPhysicsFreeVector()
{
}

// --------------------------------------------------------------------

void G4LPhysicsFreeVector::DumpValues()
{
   for (size_t i = 0; i < numberOfNodes; i++)
   {
      G4cout << binVector[i] << "   " << dataVector[i]/millibarn << G4endl;
   }
}

// --------------------------------------------------------------------
