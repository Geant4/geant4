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
// $Id: G4LPhysicsFreeVector.cc,v 1.11 2007-05-14 10:07:35 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------------
// Class G4LPhysicsFreeVector
// Derived from base class G4PhysicsVector
// This is a free vector for Low Energy Physics cross section data
//
// F.W. Jones, TRIUMF, 04-JUN-96
//
// 27-MAR-97 FWJ: first version for Alpha release
// 11-NOV-00 H.Kurashige : use STL vector for dataVector and binVector
// 13-MAY-07 V.Ivanchenko : initialise vectors of length=nbin+1 (one bin
//                          is hidden for users) - fixed bug #942 
//

#include "G4LPhysicsFreeVector.hh"

#include <stdio.h>

G4LPhysicsFreeVector::G4LPhysicsFreeVector()
   : verboseLevel(0)
{
   type = T_G4LPhysicsFreeVector;

   edgeMin = 0.0;
   edgeMax = 0.0;
   numberOfBin = 0;
}

G4LPhysicsFreeVector::G4LPhysicsFreeVector(size_t nbin,
                                           G4double binmin,
                                           G4double binmax)
   : verboseLevel(0)
{
   type = T_G4LPhysicsFreeVector;

   edgeMin = binmin;
   edgeMax = binmax;
   numberOfBin = nbin;
   lastEnergy = 0.;
   lastValue = 0.;
   lastBin = 0;
   binVector.reserve(nbin+1);
   dataVector.reserve(nbin+1);
   for (size_t i=0; i<=numberOfBin; i++)
   {
     binVector.push_back(0.0);
     dataVector.push_back(0.0);
   }
}  

G4LPhysicsFreeVector::~G4LPhysicsFreeVector()
{
}

void G4LPhysicsFreeVector::DumpValues()
{
   for (size_t i = 0; i < numberOfBin; i++)
   {
      printf(" %12.4f   %7.1f\n", binVector[i], dataVector[i]/millibarn);
   }
}
