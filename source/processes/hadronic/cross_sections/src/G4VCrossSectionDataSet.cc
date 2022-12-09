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
// -------------------------------------------------------------------
//
// GEANT4 Class file
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
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4NistManager.hh"
#include "Randomize.hh"
#include "G4HadronicParameters.hh"

G4VCrossSectionDataSet::G4VCrossSectionDataSet(const G4String& nam) :
  verboseLevel(0),name(nam),minKinEnergy(0.0),
  maxKinEnergy(G4HadronicParameters::Instance()->GetMaxEnergy()),
  isForAllAtomsAndEnergies(false) 
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
  G4int Z = elm->GetZasInt();

  if (IsElementApplicable(part, Z, mat)) { 
    return GetElementCrossSection(part, Z, mat);
  }

  // isotope-wise cross section making sum over available
  // isotope cross sections, which may be incomplete, so
  // the result is corrected 
  std::size_t nIso = elm->GetNumberOfIsotopes();    
  G4double fact = 0.0;
  G4double xsec = 0.0;

  // user-defined isotope abundances        
  const G4IsotopeVector* isoVector = elm->GetIsotopeVector();
  const G4double* abundVector = elm->GetRelativeAbundanceVector();
  for (std::size_t j=0; j<nIso; ++j) {
    const G4Isotope* iso = (*isoVector)[j];
    G4int A = iso->GetN();
    if(abundVector[j] > 0.0 && IsIsoApplicable(part, Z, A, elm, mat)) {
      fact += abundVector[j];
      xsec += abundVector[j]*GetIsoCrossSection(part, Z, A, iso, elm, mat);
    }
  }
  return (fact > 0.0) ? xsec/fact : 0.0; 
}

G4double 
G4VCrossSectionDataSet::ComputeCrossSectionPerElement(
                        G4double kinEnergy, G4double loge,
                        const G4ParticleDefinition* pd, 
                        const G4Element* elm, const G4Material* mat)
{
  G4int Z = elm->GetZasInt();
  std::size_t nIso = elm->GetNumberOfIsotopes();    
  G4double xsec = 0.0;
  const G4IsotopeVector* isoVector = elm->GetIsotopeVector();
  const G4double* abundVector = elm->GetRelativeAbundanceVector();
  for (std::size_t j=0; j<nIso; ++j) {
    const G4Isotope* iso = (*isoVector)[j];
    G4int A = iso->GetN();
    xsec += abundVector[j]*
      ComputeIsoCrossSection(kinEnergy, loge, pd, Z, A, iso, elm, mat);
  }
  return xsec;
}

G4double 
G4VCrossSectionDataSet::GetElementCrossSection(const G4DynamicParticle* dynPart,
					       G4int Z,
					       const G4Material* mat)
{
  G4ExceptionDescription ed;
  ed << "GetElementCrossSection is not implemented in <" << name << ">\n"
     << "Particle: " << dynPart->GetDefinition()->GetParticleName()
     << "  Ekin(MeV)= "  << dynPart->GetKineticEnergy()/MeV;
  if(nullptr != mat) { ed << "  material: " << mat->GetName(); }
  ed << " target Z= " << Z << G4endl;
  G4Exception("G4VCrossSectionDataSet::GetElementCrossSection", "had001", 
              FatalException, ed);
  return 0.0;
}

G4double 
G4VCrossSectionDataSet::GetIsoCrossSection(const G4DynamicParticle* dynPart,
					   G4int Z, G4int A,
					   const G4Isotope*,
					   const G4Element* elm,
					   const G4Material* mat)
{
  G4ExceptionDescription ed;
  ed << "GetIsoCrossSection is not implemented in <" << name << ">\n"
     << "Particle: " << dynPart->GetDefinition()->GetParticleName()
     << "  Ekin(MeV)= "  << dynPart->GetKineticEnergy()/MeV;
  if(nullptr != mat) { ed << "  material: " << mat->GetName(); }
  if(nullptr != elm) { ed << " element: " << elm->GetName(); }
  ed << " target Z= " << Z << " A= " << A << G4endl;
  G4Exception("G4VCrossSectionDataSet::GetIsoCrossSection", "had001", 
              FatalException, ed);
  return 0.0;
}

G4double 
G4VCrossSectionDataSet::ComputeIsoCrossSection(G4double kinEnergy, G4double,
			 		       const G4ParticleDefinition* pd,
					       G4int Z, G4int A,
					       const G4Isotope*,
					       const G4Element* elm,
					       const G4Material* mat)
{
  G4ExceptionDescription ed;
  ed << "GetIsoCrossSection is not implemented in <" << name << ">\n"
     << "Particle: " << pd->GetParticleName()
     << "  Ekin(MeV)= "  << kinEnergy/CLHEP::MeV;
  if(nullptr != mat) { ed << "  material: " << mat->GetName(); }
  if(nullptr != elm) { ed << " element: " << elm->GetName(); }
  ed << " target Z= " << Z << " A= " << A << G4endl;
  G4Exception("G4VCrossSectionDataSet::GetIsoCrossSection", "had001", 
              FatalException, ed);
  return 0.0;
}

const G4Isotope* 
G4VCrossSectionDataSet::SelectIsotope(const G4Element* anElement, 
                                      G4double, G4double)
{
  G4int nIso = (G4int)anElement->GetNumberOfIsotopes();
  const G4Isotope* iso = anElement->GetIsotope(0);

  // more than 1 isotope
  if(1 < nIso) {
    const G4double* abundVector = anElement->GetRelativeAbundanceVector();
    G4double sum = 0.0;
    G4double q = G4UniformRand();
    for (G4int j=0; j<nIso; ++j) {
      sum += abundVector[j];
      if(q <= sum) {
	iso = anElement->GetIsotope(j);
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
