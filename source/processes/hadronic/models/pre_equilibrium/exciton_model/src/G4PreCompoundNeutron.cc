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
// $Id: G4PreCompoundNeutron.cc 90337 2015-05-26 08:34:27Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4PreCompoundNeutron
//
// Author:         V.Lara
//
// Modified:  
// 21.08.2008 J. M. Quesada add choice of options  
// 10.02.2009 J. M. Quesada set default opt3
// 20.08.2010 V.Ivanchenko added G4Pow and G4PreCompoundParameters pointers
//                         use int Z and A and cleanup
// 

#include "G4PreCompoundNeutron.hh"
#include "G4SystemOfUnits.hh"
#include "G4Neutron.hh"

G4PreCompoundNeutron::G4PreCompoundNeutron()
  : G4PreCompoundNucleon(G4Neutron::Neutron(), &theNeutronCoulombBarrier)
{}

G4PreCompoundNeutron::~G4PreCompoundNeutron()
{}

G4double G4PreCompoundNeutron::GetRj(G4int nParticles, G4int nCharged) const
{
  G4double rj = 0.0;
  if(nParticles > 0) { 
    rj = static_cast<G4double>(nParticles - nCharged)/
      static_cast<G4double>(nParticles);
  }
  return rj;
}

G4double G4PreCompoundNeutron::GetAlpha() const
{
  return 0.76+2.2/theResA13;
}

G4double G4PreCompoundNeutron::GetBeta() const
{
  return (2.12/(theResA13*theResA13)-0.05)*MeV/GetAlpha();
}

