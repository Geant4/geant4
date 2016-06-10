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
// $Id: G4DeuteronEvaporationProbability.cc 89518 2015-04-15 14:43:30Z gcosmo $
//
// J.M. Quesada (August2008). Based on:
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//
// Modified:
// 03-09-2008 J.M. Quesada for external choice of inverse cross section option
// 17-11-2010 V.Ivanchenko integer Z and A


#include "G4DeuteronEvaporationProbability.hh"

G4DeuteronEvaporationProbability::G4DeuteronEvaporationProbability() :
    G4EvaporationProbability(2,1,3.0,&theCoulombBarrier)
{}

G4DeuteronEvaporationProbability::~G4DeuteronEvaporationProbability()
{}

G4double 
G4DeuteronEvaporationProbability::CalcAlphaParam(const G4Fragment& fragment) 
{
  // Data comes from 
  // Dostrovsky, Fraenkel and Friedlander
  // Physical Review, vol 116, num. 3 1959
  // 
  // const G4int size = 5;
  // G4double Zlist[5] = { 10.0, 20.0, 30.0, 50.0, 70.0};
  // G4double Cp[5] = { 0.50, 0.28, 0.20, 0.15, 0.10};
  // C for deuteron is equal to C for protons divided by 2

  G4int aZ = fragment.GetZ_asInt()-GetZ();
  G4double C;
	
  if (aZ >= 70) {
    C = 0.10;
  } else {
    C = ((((0.15417e-06*aZ) - 0.29875e-04)*aZ + 0.21071e-02)*aZ 
	 - 0.66612e-01)*aZ + 0.98375;
  }
  return 1.0 + C*0.5;
}

G4double G4DeuteronEvaporationProbability::CalcBetaParam(const G4Fragment & ) 
{
  return 0.0;
}

