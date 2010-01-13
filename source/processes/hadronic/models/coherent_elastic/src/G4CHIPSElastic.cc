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
// $Id: G4CHIPSElastic.cc,v 1.4 2010-01-13 15:42:06 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------
//
// Geant4 Class : G4CHIPSElastic
//
// Author : V.Ivanchenko 29 June 2009 
//  
// Modified:
// 13.01.10: M.Kosov: Use G4Q(Pr/Neut)ElasticCS instead of G4QElasticCS
//
//---------------------------------------------------------------------
// CHIPS model of hadron elastic scattering
//

#include "G4CHIPSElastic.hh"
#include "G4VQCrossSection.hh"
#include "G4ParticleDefinition.hh"
#include "G4QProtonElasticCrossSection.hh"
#include "G4QNeutronElasticCrossSection.hh"

G4VQCrossSection* G4CHIPSElastic::pxsManager = 0;
G4VQCrossSection* G4CHIPSElastic::nxsManager = 0;

G4CHIPSElastic::G4CHIPSElastic() : G4VHadronElastic("hElasticCHIPS")
{
  if(!pxsManager)
  {
    pxsManager = G4QProtonElasticCrossSection::GetPointer();
    nxsManager = G4QNeutronElasticCrossSection::GetPointer();
  }
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
  G4double cs = 0.;
  if     (projPDG==2212) cs = pxsManager->GetCrossSection(false,plab,Z,N,projPDG);
  else if(projPDG==2112) cs = nxsManager->GetCrossSection(false,plab,Z,N,projPDG);
  G4double t = 0.0;
  if(cs > 0.0)
  {
    if     (projPDG==2212) t = pxsManager->GetExchangeT(Z,N,projPDG);
    else if(projPDG==2112) t = nxsManager->GetExchangeT(Z,N,projPDG);
  }
  else         t = G4VHadronElastic::SampleInvariantT(p, plab, Z, A);
  return t;
}

