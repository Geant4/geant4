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
// $Id: G4PhysicsLinearVector.cc,v 1.18 2010-05-28 05:13:43 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//--------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//  G4PhysicsLinearVector.cc
//
//  15 Feb 1996 - K.Amako : 1st version
//  19 Jun 2009 - V.Ivanchenko : removed hidden bin 
//
//--------------------------------------------------------------------

#include "G4PhysicsLinearVector.hh"

G4PhysicsLinearVector::G4PhysicsLinearVector()
  : G4PhysicsVector(), dBin(0.), baseBin(0.)
{
  type = T_G4PhysicsLinearVector;
}

G4PhysicsLinearVector::G4PhysicsLinearVector(size_t theNbin)
  : G4PhysicsVector(), dBin(0.), baseBin(0.)
{
  type = T_G4PhysicsLinearVector;

  numberOfNodes = theNbin + 1;
  dataVector.reserve(numberOfNodes);
  binVector.reserve(numberOfNodes);      

  for (size_t i=0; i<numberOfNodes; i++)
  {
     binVector.push_back(0.0);
     dataVector.push_back(0.0);
  }
}  

G4PhysicsLinearVector::G4PhysicsLinearVector(G4double theEmin, 
                                             G4double theEmax, size_t theNbin)
  : G4PhysicsVector(),
    dBin((theEmax-theEmin)/theNbin),
    baseBin(theEmin/dBin)
{
  type = T_G4PhysicsLinearVector;

  numberOfNodes = theNbin + 1;
  dataVector.reserve(numberOfNodes);
  binVector.reserve(numberOfNodes);      

  binVector.push_back(theEmin);
  dataVector.push_back(0.0);

  for (size_t i=1; i<numberOfNodes-1; i++)
    {
      binVector.push_back( theEmin + i*dBin );
      dataVector.push_back(0.0);
    }
  binVector.push_back(theEmax);
  dataVector.push_back(0.0);

  edgeMin = binVector[0];
  edgeMax = binVector[numberOfNodes-1];

}  

G4PhysicsLinearVector::~G4PhysicsLinearVector(){}

G4bool G4PhysicsLinearVector::Retrieve(std::ifstream& fIn, G4bool ascii)
{
  G4bool success = G4PhysicsVector::Retrieve(fIn, ascii);
  if (success)
  {
    G4double theEmin = binVector[0];
    dBin = binVector[1]-theEmin;
    baseBin = theEmin/dBin;
  }
  return success;
}

G4PhysicsLinearVector::G4PhysicsLinearVector(const G4PhysicsLinearVector& right)
  : G4PhysicsVector(right)
{
  dBin = right.dBin;
  baseBin = right.baseBin;
}

G4PhysicsLinearVector& 
G4PhysicsLinearVector::operator=(const G4PhysicsLinearVector& right)
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
