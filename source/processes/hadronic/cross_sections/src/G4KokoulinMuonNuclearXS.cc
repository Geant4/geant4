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
//
// Author:      D.H. Wright
// Date:        1 February 2011
//
// Modified:
//
// 19 Aug 2011, V.Ivanchenko move to new design and make x-section per element

// Description: use Kokoulin's parameterized calculation of virtual 
//              photon production cross section and conversion to
//              real photons.

#include "G4KokoulinMuonNuclearXS.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsVector.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4NucleiProperties.hh"
#include "G4NistManager.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

// factory
#include "G4CrossSectionFactory.hh"
//
G4_DECLARE_XS_FACTORY(G4KokoulinMuonNuclearXS);

G4PhysicsVector* G4KokoulinMuonNuclearXS::theCrossSection[] = {0};

G4KokoulinMuonNuclearXS::G4KokoulinMuonNuclearXS()
  :G4VCrossSectionDataSet(Default_Name()), 
  LowestKineticEnergy(1*GeV), HighestKineticEnergy(1*PeV),
  TotBin(60), CutFixed(0.2*GeV), isInitialized(false), isMaster(false)
{}

G4KokoulinMuonNuclearXS::~G4KokoulinMuonNuclearXS()
{
  if (isMaster) {
    for(G4int i=0; i<MAXZMUN; ++i) {
      delete theCrossSection[i];
      theCrossSection[i] = 0;
    }
  }
}


void
G4KokoulinMuonNuclearXS::CrossSectionDescription(std::ostream& outFile) const
{
    outFile << "G4KokoulinMuonNuclearXS provides the total inelastic\n"
    << "cross section for mu- and mu+ interactions with nuclei.\n"
    << "R. Kokoulin's approximation of the Borog and Petrukhin double\n"
    << "differential cross section at high energy and low Q**2 is integrated\n"
    << "over the muon energy loss to get the total cross section as a\n"
    << "function of muon kinetic energy\n" ;
}


G4bool 
G4KokoulinMuonNuclearXS::IsElementApplicable(const G4DynamicParticle*, 
					     G4int, const G4Material*)
{
  return true;
}

void 
G4KokoulinMuonNuclearXS::BuildPhysicsTable(const G4ParticleDefinition&)
{
  if(!isInitialized) { 
    isInitialized = true; 
    for(G4int i=0; i<MAXZMUN; ++i) {
      if(theCrossSection[i]) { return; }
    }
    isMaster = true; 
  }
  if(isMaster) { BuildCrossSectionTable(); }
}

void G4KokoulinMuonNuclearXS::BuildCrossSectionTable()
{
  G4double energy, A, Value;
  G4int Z;

  std::size_t nEl = G4Element::GetNumberOfElements(); 
  const G4ElementTable* theElementTable = G4Element::GetElementTable();
  G4NistManager* nistManager = G4NistManager::Instance();

  for (std::size_t j = 0; j < nEl; ++j) {
    Z = G4lrint((*theElementTable)[j]->GetZ());

    //AR-24Apr2018 Switch to treat transuranic elements as uranium  
    const G4bool isHeavyElementAllowed = true; if ( isHeavyElementAllowed && Z>92 ) Z=92;

    A  = nistManager->GetAtomicMassAmu(Z);
    if(Z < MAXZMUN && !theCrossSection[Z]) {
      theCrossSection[Z] = new G4PhysicsLogVector(LowestKineticEnergy,
						  HighestKineticEnergy, 
						  TotBin);
      for (G4int i = 0; i <= TotBin; ++i) {
	energy = theCrossSection[Z]->Energy(i);
	Value = ComputeMicroscopicCrossSection(energy, A);
	theCrossSection[Z]->PutValue(i,Value);
      }
    }
  }
}

G4double G4KokoulinMuonNuclearXS::
ComputeMicroscopicCrossSection(G4double KineticEnergy, G4double A)
{
  // Calculate cross section (differential in muon incident kinetic energy) by 
  // integrating the double differential cross section over the energy loss

  static const G4double xgi[] = 
    {0.0199,0.1017,0.2372,0.4083,0.5917,0.7628,0.8983,0.9801};
  static const G4double wgi[] = 
    {0.0506,0.1112,0.1569,0.1813,0.1813,0.1569,0.1112,0.0506};
  static const G4double ak1 = 6.9;
  static const G4double ak2 = 1.0;

  G4double Mass = G4MuonMinus::MuonMinus()->GetPDGMass();

  G4double CrossSection = 0.0;
  if (KineticEnergy <= CutFixed) return CrossSection; 

  G4double epmin = CutFixed;
  G4double epmax = KineticEnergy + Mass - 0.5*proton_mass_c2;
  if (epmax <= epmin) return CrossSection; // NaN bug correction

  G4double aaa = G4Log(epmin);
  G4double bbb = G4Log(epmax);
  G4int kkk = std::max(1,G4int((bbb-aaa)/ak1 +ak2));
  G4double hhh = (bbb-aaa)/kkk ;
  G4double epln;
  G4double ep;
  G4double x;

  for (G4int l = 0; l < kkk; ++l) {
    x = aaa + hhh*l;
    for (G4int ll = 0; ll < 8; ++ll) {
      epln=x+xgi[ll]*hhh;
      ep = G4Exp(epln);
      CrossSection += 
	ep*wgi[ll]*ComputeDDMicroscopicCrossSection(KineticEnergy, 0, A, ep);
    }
  }

  CrossSection *= hhh ;
  if (CrossSection < 0.) { CrossSection = 0.; }
  return CrossSection;
}

G4double G4KokoulinMuonNuclearXS::
ComputeDDMicroscopicCrossSection(G4double KineticEnergy, G4double,
                                 G4double A, G4double epsilon)
{
  // Calculate the double-differential microscopic cross section (in muon
  // incident kinetic energy and energy loss) using the cross section formula
  // of R.P. Kokoulin (18/01/98)

  static const G4double alam2 = 0.400*GeV*GeV;
  static const G4double alam  = 0.632456*GeV;
  static const G4double coeffn = fine_structure_const/pi;   

  G4double ParticleMass = G4MuonMinus::MuonMinus()->GetPDGMass();
  G4double TotalEnergy = KineticEnergy + ParticleMass;

  G4double DCrossSection = 0.;

  if ((epsilon >= TotalEnergy - 0.5*proton_mass_c2) ||
      (epsilon <= CutFixed) ) { return DCrossSection; }

  G4double ep = epsilon/GeV;
  G4double aeff = 0.22*A+0.78*G4Exp(0.89*G4Log(A));       //shadowing 
  G4double sigph = (49.2+11.1*G4Log(ep)+151.8/std::sqrt(ep))*microbarn; 
  
  G4double v = epsilon/TotalEnergy;
  G4double v1 = 1.-v;
  G4double v2 = v*v;
  G4double mass2 = ParticleMass*ParticleMass;

  G4double up = TotalEnergy*TotalEnergy*v1/mass2*(1.+mass2*v2/(alam2*v1));
  G4double down = 1.+epsilon/alam*(1.+alam/(2.*proton_mass_c2)+epsilon/alam);

  DCrossSection = coeffn*aeff*sigph/epsilon*
                  (-v1+(v1+0.5*v2*(1.+2.*mass2/alam2))*G4Log(up/down));

  if (DCrossSection < 0.) { DCrossSection = 0.; }
  return DCrossSection;
}

G4double G4KokoulinMuonNuclearXS::
GetElementCrossSection(const G4DynamicParticle* aPart,
		       G4int Z, const G4Material*)
{
  //AR-24Apr2018 Switch to treat transuranic elements as uranium  
  const G4bool isHeavyElementAllowed = true; if ( isHeavyElementAllowed && Z>92 ) Z=92;

  return theCrossSection[Z]->Value(aPart->GetKineticEnergy());
}

