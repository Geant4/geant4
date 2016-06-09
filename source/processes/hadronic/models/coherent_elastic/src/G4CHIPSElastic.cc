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
// $Id: G4CHIPSElastic.cc,v 1.3 2009/10/08 18:56:57 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-03 $
//
//---------------------------------------------------------------------
//
// Geant4 Class : G4CHIPSElastic
//
// Author : V.Ivanchenko 29 June 2009 
//  
// Modified:
//
//---------------------------------------------------------------------
// CHIPS model of hadron elastic scattering
//

#include "G4CHIPSElastic.hh"
#include "G4VQCrossSection.hh"
#include "G4ParticleDefinition.hh"
#include "G4QElasticCrossSection.hh"

G4VQCrossSection* G4CHIPSElastic::xsManager = 0;

G4CHIPSElastic::G4CHIPSElastic() : G4VHadronElastic("hElasticCHIPS")
{
  if(!xsManager) {xsManager = G4QElasticCrossSection::GetPointer();}
}

G4CHIPSElastic::~G4CHIPSElastic()
{}

G4double 
G4CHIPSElastic::SampleInvariantT(const G4ParticleDefinition* p, 
				 G4double plab, G4int Z, G4int A)
{
  G4int N = A - Z;
  if(Z == 1 && N == 2) N = 1;
  else if(Z == 2 && N == 1) N = 2;
  G4int projPDG = p->GetPDGEncoding();
  G4double cs = xsManager->GetCrossSection(false,plab,Z,N,projPDG);
  G4double t = 0.0;
  if(cs > 0.0) t = xsManager->GetExchangeT(Z,N,projPDG);
  else         t = G4VHadronElastic::SampleInvariantT(p, plab, Z, A);
  return t;
}

