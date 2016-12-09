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
// $Id: G4VCrossSectionDataSet.cc 98738 2016-08-09 12:53:25Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:    G4VCrossSectionDataSet
//
// Author  F.W. Jones, TRIUMF, 20-JAN-97
//
// Modifications:
// 23.01.2009 V.Ivanchenko move constructor and destructor to source
// 12.08.2011 G.Folger, V.Ivanchenko, T.Koi, D.Wright redesign the class
//

#include "G4VCrossSectionDataSet.hh"
#include "G4SystemOfUnits.hh"
#include "G4CrossSectionDataSetRegistry.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4NistManager.hh"
#include "G4HadronicException.hh"
#include "G4HadTmpUtil.hh"
#include "Randomize.hh"

G4VCrossSectionDataSet::G4VCrossSectionDataSet(const G4String& nam) :
  verboseLevel(0),minKinEnergy(0.0),maxKinEnergy(100*TeV),name(nam) 
{
  registry = G4CrossSectionDataSetRegistry::Instance();
  registry->Register(this);
}

G4VCrossSectionDataSet::~G4VCrossSectionDataSet()
{
  registry->DeRegister(this);
}

G4bool 
G4VCrossSectionDataSet::IsElementApplicable(const G4DynamicParticle*, 
					    G4int,
					    const G4Material*)
{
  return false;
}

G4bool 
G4VCrossSectionDataSet::IsIsoApplicable(const G4DynamicParticle*, 
                                        G4int, G4int,
                                        const G4Element*,  
                                        const G4Material*)
{
  return false;
}

G4double 
G4VCrossSectionDataSet::ComputeCrossSection(const G4DynamicParticle* part, 
					    const G4Element* elm,
					    const G4Material* mat)
{
  G4int Z = G4lrint(elm->GetZ());

  if (IsElementApplicable(part, Z, mat)) { 
    return GetElementCrossSection(part, Z, mat);
  }

  // isotope-wise cross section making sum over available
  // isotope cross sections, which may be incomplete, so
  // the result is corrected 
  G4int nIso = elm->GetNumberOfIsotopes();    
  G4double fact = 0.0;
  G4double xsec = 0.0;
  G4Isotope* iso = 0;

  if (0 < nIso) { 

    // user-defined isotope abundances        
    G4IsotopeVector* isoVector = elm->GetIsotopeVector();
    G4double* abundVector = elm->GetRelativeAbundanceVector();

    for (G4int j = 0; j<nIso; ++j) {
      iso = (*isoVector)[j];
      G4int A = iso->GetN();
      if(abundVector[j] > 0.0 && IsIsoApplicable(part, Z, A, elm, mat)) {
        fact += abundVector[j];
	xsec += abundVector[j]*GetIsoCrossSection(part, Z, A, iso, elm, mat);
      }
    }

  } else {

    // natural isotope abundances
    G4NistManager* nist = G4NistManager::Instance();
    G4int n0 = nist->GetNistFirstIsotopeN(Z);
    G4int nn = nist->GetNumberOfNistIsotopes(Z);
    for (G4int A = n0; A < n0+nn; ++A) {
      G4double abund = nist->GetIsotopeAbundance(Z, A);
      if(abund > 0.0 && IsIsoApplicable(part, Z, A, elm, mat)) {
        fact += abund;
	xsec += abund*GetIsoCrossSection(part, Z, A, iso, elm, mat);
      }
    }
  }
  if(fact > 0.0) { xsec /= fact; }
  return xsec;
}

G4double 
G4VCrossSectionDataSet::GetElementCrossSection(const G4DynamicParticle* dynPart,
					       G4int Z,
					       const G4Material* mat)
{
  G4cout << "G4VCrossSectionDataSet::GetCrossSection per element ERROR: "
	 << " there is no cross section for "
	 << dynPart->GetDefinition()->GetParticleName()
	 << "  E(MeV)= "  << dynPart->GetKineticEnergy()/MeV;
  if(mat) { G4cout << "  inside " << mat->GetName(); }
  G4cout << " for Z= " << Z << G4endl;
  throw G4HadronicException(__FILE__, __LINE__, 
        "G4VCrossSectionDataSet::GetElementCrossSection is absent");
}

G4double 
G4VCrossSectionDataSet::GetIsoCrossSection(const G4DynamicParticle* dynPart,
					   G4int Z, G4int A,
					   const G4Isotope*,
					   const G4Element* elm,
					   const G4Material* mat)
{
  G4cout << "G4VCrossSectionDataSet::GetCrossSection per isotope ERROR: "
	 << " there is no cross section for "
	 << dynPart->GetDefinition()->GetParticleName()
	 << "  E(MeV)= "  << dynPart->GetKineticEnergy()/MeV;
  if(mat) { G4cout << "  inside " << mat->GetName(); }
  if(elm) { G4cout << " for " << elm->GetName(); }
  G4cout << "  Z= " << Z << " A= " << A << G4endl;
  throw G4HadronicException(__FILE__, __LINE__,
        "G4VCrossSectionDataSet::GetIsoCrossSection is absent");
}

G4Isotope* 
G4VCrossSectionDataSet::SelectIsotope(const G4Element* anElement, G4double)
{
  G4int nIso = anElement->GetNumberOfIsotopes();
  G4IsotopeVector* isoVector = anElement->GetIsotopeVector();
  G4Isotope* iso = (*isoVector)[0];

  // more than 1 isotope
  if(1 < nIso) {
    G4double* abundVector = anElement->GetRelativeAbundanceVector();
    G4double sum = 0.0;
    G4double q = G4UniformRand();
    for (G4int j = 0; j<nIso; ++j) {
      sum += abundVector[j];
      if(q <= sum) {
	iso = (*isoVector)[j];
	break;
      }
    }
  }
  return iso;
}

void G4VCrossSectionDataSet::BuildPhysicsTable(const G4ParticleDefinition&)
{}

void G4VCrossSectionDataSet::DumpPhysicsTable(const G4ParticleDefinition&)
{}

void G4VCrossSectionDataSet::CrossSectionDescription(std::ostream& outFile) const
{
  outFile << "The description for this cross section data set has not been written yet.\n";
}
