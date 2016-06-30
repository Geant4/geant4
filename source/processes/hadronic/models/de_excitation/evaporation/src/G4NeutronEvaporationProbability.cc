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
// $Id: G4NeutronEvaporationProbability.cc 96634 2016-04-27 09:31:49Z gcosmo $
//
// J.M. Quesada (August2008). Based on:
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//
// Modified:
// 03-09-2008 J.M. Quesada for external choice of inverse cross section option
// 17-11-2010 V.Ivanchenko integer Z and A

#include "G4NeutronEvaporationProbability.hh"
#include "G4SystemOfUnits.hh"

G4NeutronEvaporationProbability::G4NeutronEvaporationProbability() :
    G4EvaporationProbability(1,0,2.0,&theCoulombBarrier)
{}

G4NeutronEvaporationProbability::~G4NeutronEvaporationProbability()
{}

G4double G4NeutronEvaporationProbability::CalcAlphaParam(const G4Fragment& fragment)
{ 
  return 0.76+2.2/fG4pow->Z13(fragment.GetA_asInt() - 1);
}
	
G4double G4NeutronEvaporationProbability::CalcBetaParam(const G4Fragment& fragment) 
{ 
  return (2.12/fG4pow->Z23(fragment.GetA_asInt() - 1) - 0.05)*CLHEP::MeV/
    CalcAlphaParam(fragment); 
}

