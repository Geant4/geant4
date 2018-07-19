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
// $Id: G4AlphaEvaporationProbability.cc 99962 2016-10-12 14:09:03Z gcosmo $
//
// J.M. Quesada (August2008). Based on:
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//
// Modified:
// 03-09-2008 J.M. Quesada for external choice of inverse cross section option
// 17-11-2010 V.Ivanchenko integer Z and A

#include "G4AlphaEvaporationProbability.hh"
#include "G4SystemOfUnits.hh"

G4AlphaEvaporationProbability::G4AlphaEvaporationProbability() :
    G4EvaporationProbability(4,2,1.0,&theCoulombBarrier) 
{}

G4AlphaEvaporationProbability::~G4AlphaEvaporationProbability()
{}

G4double G4AlphaEvaporationProbability::CalcAlphaParam(const G4Fragment& fragment)
{ 
  // Data comes from 
  // Dostrovsky, Fraenkel and Friedlander
  // Physical Review, vol 116, num. 3 1959
  // 
  // const G4int size = 5;
  // G4double Zlist[5] = { 10.0, 20.0, 30.0, 50.0, 70.0};
  //	G4double Calpha[5] = { 0.10, 0.10, 0.10, 0.08, 0.06};

  G4int aZ = fragment.GetZ_asInt()-GetZ();
  G4double C;
	
  if (aZ <= 30) 
    {
      C = 0.10;
    }
  else if (aZ <= 50)
    {
      C = 0.1 - (aZ-30)*0.001;
    }
  else if (aZ < 70)
    {
      C = 0.08 - (aZ-50)*0.001;
    }
  else 
    {
      C = 0.06;
    }
  return 1.0 + C;
}
	
G4double G4AlphaEvaporationProbability::CalcBetaParam(const G4Fragment &) 
{ 
  return 0.0; 
}
