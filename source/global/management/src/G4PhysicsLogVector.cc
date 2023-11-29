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
// --------------------------------------------------------------------

#include "G4PhysicsLogVector.hh"
#include "G4Exp.hh"
#include "G4Log.hh"

// --------------------------------------------------------------------
G4PhysicsLogVector::G4PhysicsLogVector(G4bool spline)
  : G4PhysicsVector(spline)
{
  type = T_G4PhysicsLogVector;
}

// --------------------------------------------------------------------
G4PhysicsLogVector::G4PhysicsLogVector(G4double Emin, G4double Emax,
                                       std::size_t Nbin, G4bool spline)
  : G4PhysicsVector(spline)
{
  numberOfNodes = Nbin + 1;
  if(Nbin < 2 || Emin >= Emax || Emin <= 0.0)
  {
    G4ExceptionDescription ed;
    ed << "G4PhysicsLogVector with wrong parameters: theNbin= " << Nbin
       << " Emin= " << Emin << " Emax= " << Emax;
    G4Exception("G4PhysicsLogVector::G4PhysicsLogVector()", "glob03",
                FatalException, ed, "Nbins should be > 1 and Emax > Emin > 0");
  }
  if(numberOfNodes < 3)
  {
    numberOfNodes = 3;
  }
  type = T_G4PhysicsLogVector;

  binVector.resize(numberOfNodes);
  dataVector.resize(numberOfNodes, 0.0);
  binVector[0] = Emin;
  binVector[numberOfNodes - 1] = Emax;
  Initialise();

  for(std::size_t i = 1; i <= idxmax; ++i)
  {
    binVector[i] = edgeMin*G4Exp(i / invdBin);
  }
}

// --------------------------------------------------------------------
void G4PhysicsLogVector::Initialise()
{
  idxmax  = numberOfNodes - 2;
  edgeMin = binVector[0];
  edgeMax = binVector[numberOfNodes - 1];
  invdBin = (idxmax + 1) / G4Log(edgeMax/edgeMin);
  logemin = G4Log(edgeMin);
}

// --------------------------------------------------------------------
