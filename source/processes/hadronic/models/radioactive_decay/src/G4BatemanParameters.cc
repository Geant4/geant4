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
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File:   G4BatemanParameters.cc                                            //
//  Author: D.H. Wright (SLAC)                                                //
//  Date:   15 December 2015                                                  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "G4BatemanParameters.hh"
// #include "G4ParticleDefinition.hh"
// #include "G4ParticleTable.hh"
// #include "G4DecayTable.hh"
// #include "G4DecayProducts.hh"


G4BatemanParameters::G4BatemanParameters()
 : Z(0), A(0), E(0.0), generation(0)
{}


G4BatemanParameters::G4BatemanParameters(const G4BatemanParameters& right)
{
  Z = right.Z;
  A = right.A;
  E = right.E;
  generation = right.generation;
  Acoeffs = right.Acoeffs;
  taus = right.taus;
}

G4BatemanParameters& G4BatemanParameters::operator=(const G4BatemanParameters& right)
{
  if (this != &right) { 
    Z = right.Z;
    A = right.A;
    E = right.E;
    generation = right.generation;
    Acoeffs = right.Acoeffs;
    taus = right.taus;
  }
  return *this;
}


G4BatemanParameters::~G4BatemanParameters()
{} 


void 
G4BatemanParameters::SetParameters(G4int aZ, G4int anA, G4double anE, G4int aG, 
                                   std::vector<G4double> theCoeffs, 
                                   std::vector<G4double> theTaus)
{
  Z = aZ;
  A = anA;
  E = anE;
  generation = aG;
  Acoeffs = theCoeffs; 
  taus = theTaus;
}


void G4BatemanParameters::DumpInfo()
{
  G4cout << " Z: " << Z << "  A: " << A << "  E: " << E << " Generation: "
         << generation << G4endl;

  G4cout << " A coefficients: ";
  for (G4int i = 0; i < G4int(Acoeffs.size()); i++) G4cout << Acoeffs[i];
  G4cout << G4endl;

  G4cout << " Mean lifes (tau): ";
  for (G4int i = 0; i < G4int(taus.size()); i++) G4cout << taus[i];
  G4cout << G4endl;
}

