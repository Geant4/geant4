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
// author: Vladimir.Grichine@cern.ch
//
// Implements data from: Barashenkov V.S., Nucleon-Nucleus Cross Section,
// Preprint JINR P2-89-770, p. 12, Dubna 1989 (scanned version from KEK)
// Based on G. Folger version of G4PiNuclearCrossSection class
//
// Modified: 16.08.2018 V.Ivanchenko - major update
//

#include "G4NucleonNuclearCrossSection.hh"
#include "G4SystemOfUnits.hh"
#include "G4DynamicParticle.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4ComponentBarNucleonNucleusXsc.hh"

// factory
#include "G4CrossSectionFactory.hh"
//
G4_DECLARE_XS_FACTORY(G4NucleonNuclearCrossSection);


using namespace std;

///////////////////////////////////////////////////////////////////////////////

G4NucleonNuclearCrossSection::G4NucleonNuclearCrossSection()
 : G4VCrossSectionDataSet(Default_Name()),
   fTotalXsc(0.0), fInelasticXsc(0.0), fElasticXsc(0.0)
{
  theNeutron = G4Neutron::Neutron();
  theProton  = G4Proton::Proton();
  fBarash = new G4ComponentBarNucleonNucleusXsc();  
}

///////////////////////////////////////////////////////////////////////////////
//

G4NucleonNuclearCrossSection::~G4NucleonNuclearCrossSection()
{}

////////////////////////////////////////////////////////////////////////////
//

G4bool G4NucleonNuclearCrossSection::IsElementApplicable(
       const G4DynamicParticle*, G4int Z, const G4Material*)
{
  return (Z > 1);
}

////////////////////////////////////////////////////////////////////////////
//

G4double G4NucleonNuclearCrossSection::GetElementCrossSection(
         const G4DynamicParticle* dp, G4int Z, const G4Material*)
{
  ComputeCrossSections(dp->GetDefinition(), dp->GetKineticEnergy(), Z);
  return fInelasticXsc;
}

////////////////////////////////////////////////////////////////////////////
//

void G4NucleonNuclearCrossSection::ComputeCrossSections(
            const G4ParticleDefinition* pd, 
            G4double kinEnergy, G4int Z)
{
  fBarash->ComputeCrossSections(pd, kinEnergy, Z);
  fTotalXsc = fBarash->GetTotalXsc();
  fInelasticXsc = fBarash->GetInelasticXsc();
  fElasticXsc = fBarash->GetElasticXsc();
}


////////////////////////////////////////////////////////////////////////////
//

void
G4NucleonNuclearCrossSection::CrossSectionDescription(std::ostream& outFile) const
{
  outFile << "G4NucleonNuclearCrossSection is a variant of the Barashenkov\n"
          << "cross section parameterization to be used of protons and\n"
          << "nucleons on targets heavier than hydrogen.  It is intended for\n"
          << "use as a cross section component and is currently used by\n"
          << "G4BGGNucleonInelasticXS.  It is valid for incident energies up\n"
          << "to 1 TeV.\n"; 
}

