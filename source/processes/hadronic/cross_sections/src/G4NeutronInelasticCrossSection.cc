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
// By JPW, working, but to be cleaned up. @@@
// 22 Dec 2006 - DHW added isotope dependence
// G.Folger, 25-Nov-2009: extend to 100TeV, using a constant above 20GeV
// 19 Aug 2011, V.Ivanchenko move to new design and make x-section per element
//

#include "G4NeutronInelasticCrossSection.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DynamicParticle.hh"
#include "G4HadTmpUtil.hh"
#include "G4NistManager.hh"
#include "G4Pow.hh"

G4NeutronInelasticCrossSection::G4NeutronInelasticCrossSection()
  : G4VCrossSectionDataSet("Wellisch-Laidlaw"), 
    minEnergy(19.9*MeV), maxEnergy(19.9*GeV)
{}

G4NeutronInelasticCrossSection::~G4NeutronInelasticCrossSection()
{}

void
G4NeutronInelasticCrossSection::CrossSectionDescription(std::ostream& outFile) const
{
  outFile << "G4NeutronInelasticCrossSection calculates the inelastic neutron\n"
          << "scattering cross section for nuclei using the Wellisch-Laidlaw\n"
          << "parameterization between 19.9 MeV and 19.9 GeV.  Above 19.9 GeV\n"
          << "the cross section is assumed to be constant.\n";
}

G4bool 
G4NeutronInelasticCrossSection::IsElementApplicable(
   const G4DynamicParticle* part, G4int Z, const G4Material*)
{
  G4double e = part->GetKineticEnergy();
  return (1 < Z && e > minEnergy);    
}

G4double G4NeutronInelasticCrossSection::
GetElementCrossSection(const G4DynamicParticle* aPart, 
		       G4int Z, const G4Material*)
{
  G4int A = G4int(G4NistManager::Instance()->GetAtomicMassAmu(Z));
  return GetCrossSection(aPart->GetKineticEnergy(), Z, A);
}

G4double 
G4NeutronInelasticCrossSection::GetCrossSection(G4double anEnergy, 
						G4int Z, G4int A)
{
  if(anEnergy > maxEnergy) { anEnergy = maxEnergy; }
  G4double cross_section = 0.0;
  if(anEnergy < keV) { return cross_section; }
  
  G4Pow* g4pow = G4Pow::GetInstance();
  G4double A13 = g4pow->Z13(A);
  
  G4double elog = std::log10(anEnergy/MeV);
  G4int nOfNeutrons = A - Z;
  G4double atomicNumber = G4double(A);
  static const G4double p1=1.3773;
  G4double p2 = 1. + 10./atomicNumber   - 0.0006*atomicNumber;
  G4double p3 = 0.6+ 13./atomicNumber   - 0.0005*atomicNumber;
  G4double p4 = 7.2449 - 0.018242*atomicNumber;
  G4double p5 = 1.64 - 1.8/atomicNumber - 0.0005*atomicNumber;
  G4double p6 = 1. + 200./atomicNumber + 0.02*atomicNumber;
  G4double p7 = (atomicNumber-70.)*(atomicNumber-200.)/11000.;

  G4double logN  = g4pow->logZ(nOfNeutrons);
  G4double part1 = pi*p1*p1*logN;
  G4double part2 = 1.+ A13 - p2*(1.-1./A13);

  G4double firstexp = -p4*(elog-p5);
  G4double first    = 1. + G4Exp(firstexp);
  G4double corr     = 1. + p3*(1.-1./first); 

  G4double secondexp= -p6*(elog-p7);
  G4double secondv   = 1.+G4Exp(secondexp);
  G4double corr2    = 1./secondv;

  G4double xsec = corr*corr2*part1*part2*10.*millibarn;
  if(xsec < 0.0) { xsec = 0.0; }
  return xsec;
}
