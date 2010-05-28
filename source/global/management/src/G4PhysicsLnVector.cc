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
// $Id: G4PhysicsLnVector.cc,v 1.21 2010-05-28 05:13:43 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
  : G4PhysicsVector(), dBin(0.), baseBin(0.)
{
  type = T_G4PhysicsLnVector;
}

G4PhysicsLnVector::G4PhysicsLnVector(size_t theNbin)
  : G4PhysicsVector(), dBin(0.), baseBin(0.)
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
  : G4PhysicsVector(),
    dBin(std::log(theEmax/theEmin)/theNbin),
    baseBin(std::log(theEmin)/dBin)
{
  type = T_G4PhysicsLnVector;

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

G4PhysicsLnVector::G4PhysicsLnVector(const G4PhysicsLnVector& right)
  : G4PhysicsVector(right)
{
  dBin = right.dBin;
  baseBin = right.baseBin;
}

G4PhysicsLnVector& 
G4PhysicsLnVector::operator=(const G4PhysicsLnVector& right)
{
  // Check assignment to self
  //
  if(this == &right) { return *this; }

  DeleteData();
  CopyData(right);

  dBin    = right.dBin;
  baseBin = right.baseBin;
  return *this;
}
