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
// $Id: G4CHIPSElasticXS.cc,v 1.2 2010-09-24 13:56:00 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:    G4CHIPSElasticXS
//
// Author  Ivantchenko, Geant4, 3-Aug-09
//
// Modifications:
//

#include "G4CHIPSElasticXS.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4Element.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4VQCrossSection.hh"
#include "G4QProtonElasticCrossSection.hh"
#include "G4QNeutronElasticCrossSection.hh"

G4CHIPSElasticXS::G4CHIPSElasticXS() 
  :  G4VCrossSectionDataSet("CHIPSElasticXS"),
     theProton(G4Proton::Proton()), 
     theNeutron(G4Neutron::Neutron()),
     thEnergy(19*CLHEP::MeV),
     isInitialized(false)
{
  //  verboseLevel = 0;
  pCManager   = G4QProtonElasticCrossSection::GetPointer();
  nCManager   = G4QNeutronElasticCrossSection::GetPointer();
}

G4CHIPSElasticXS::~G4CHIPSElasticXS()
{}

G4bool 
G4CHIPSElasticXS::IsApplicable(const G4DynamicParticle* dyn, 
			       const G4Element* elm)
{
  return (elm->GetZ() < 2.5 && dyn->GetKineticEnergy() > thEnergy); 
}

G4bool 
G4CHIPSElasticXS::IsZAApplicable(const G4DynamicParticle* dyn,
				 G4double ZZ, G4double /*AA*/)
{
  return (ZZ < 2.5 && dyn->GetKineticEnergy() > thEnergy);
}

G4bool 
G4CHIPSElasticXS::IsIsoApplicable(const G4DynamicParticle* dyn, 
				  G4int Z, G4int /*N*/)
{
  return (Z <= 2 && dyn->GetKineticEnergy() > thEnergy);
}

G4double 
G4CHIPSElasticXS::GetCrossSection(const G4DynamicParticle* aParticle,
				  const G4Element* elm,
				  G4double)
{
  G4double xs = 0.0;
  G4int Z = G4int(elm->GetZ());
  G4IsotopeVector* isv = elm->GetIsotopeVector();
  G4int ni = 0;
  if(isv) { ni = isv->size(); }

  if(ni <= 1) {
    G4int A = G4int(elm->GetN()+0.5);
    xs = GetZandACrossSection(aParticle, Z, A);
  } else {
    G4double* ab = elm->GetRelativeAbundanceVector();
    for(G4int j=0; j<ni; ++j) {
      G4int A = (*isv)[j]->GetN();
      xs += ab[j]*GetZandACrossSection(aParticle, Z, A);
    }
  }

  if(verboseLevel > 1) {
    G4cout  << "G4CHIPSElasticXS::GetCrossSection for "
	    << theParticle->GetParticleName()
	    << " on " << elm->GetName()
	    << " ekin(MeV)= " << aParticle->GetKineticEnergy()/CLHEP::MeV 
	    << ",  XSel(bn)= " << xs/CLHEP::barn << G4endl;
  }
  return xs;
}

G4double 
G4CHIPSElasticXS::GetIsoCrossSection(const G4DynamicParticle* p, 
				     const G4Isotope* iso,
				     G4double)
{
  return GetZandACrossSection(p, iso->GetZ(), iso->GetN());
}

G4double 
G4CHIPSElasticXS::GetIsoZACrossSection(const G4DynamicParticle* p, 
				       G4double ZZ,
				       G4double AA, 
				       G4double)
{
  return GetZandACrossSection(p, G4int(ZZ), G4int(AA));
}

G4double 
G4CHIPSElasticXS::GetZandACrossSection(const G4DynamicParticle* dyn, 
				       G4int Z, G4int A, G4double)
{
  G4double momentum = dyn->GetTotalMomentum();

  // only proton, deuteron and He4 x-sections
  G4int N = A - Z;
  if(Z == 1) {
    if(N > 1) { N = 1; }
  } else if(Z == 2) { N = 2; }
  
  G4double x = 0.0;
  if(theParticle == theProton) {
    x = pCManager->GetCrossSection(false,momentum,Z,N,pPDG);
  } else {
    x = nCManager->GetCrossSection(false,momentum,Z,N,pPDG);
  }
  return x;
}

void 
G4CHIPSElasticXS::BuildPhysicsTable(const G4ParticleDefinition& p)
{
  if(isInitialized) { return; }
  if(verboseLevel > 0){
    G4cout << "G4CHIPSElasticXS::BuildPhysicsTable for " 
	   << p.GetParticleName() 
	   << "  Elow(MeV)= " << thEnergy/MeV 
	   << G4endl;
  }
  isInitialized = true;
  theParticle = &p;
  if(theParticle != theProton && theParticle != theNeutron) {
    G4cout << "G4CHIPSElasticXS::BuildPhysicsTable ERROR for " 
	   << p.GetParticleName() 
	   << G4endl;
    G4Exception("G4CHIPSElasticXS", "", FatalException,"Not applicable"); 
  }
  pPDG = theParticle->GetPDGEncoding();
}

void 
G4CHIPSElasticXS::DumpPhysicsTable(const G4ParticleDefinition&)
{}

