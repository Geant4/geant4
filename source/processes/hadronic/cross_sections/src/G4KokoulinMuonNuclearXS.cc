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
// $Id: $
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
#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh"


G4KokoulinMuonNuclearXS::G4KokoulinMuonNuclearXS()
 :G4VCrossSectionDataSet("KokoulinMuonNuclearXS"), theCrossSectionTable(0),
  LowestKineticEnergy(1*GeV), HighestKineticEnergy(1*PeV),
  TotBin(60), CutFixed(0.2*GeV)
{
  BuildCrossSectionTable();
}

G4KokoulinMuonNuclearXS::~G4KokoulinMuonNuclearXS()
{
  if (theCrossSectionTable) {
    theCrossSectionTable->clearAndDestroy();
    delete theCrossSectionTable;
  }
}

G4bool 
G4KokoulinMuonNuclearXS::IsElementApplicable(const G4DynamicParticle*, 
					     G4int, const G4Material*)
{
  return true;
}


void G4KokoulinMuonNuclearXS::BuildCrossSectionTable()
{
  G4double lowEdgeEnergy;
  G4PhysicsLogVector* ptrVector;

  G4int nEl = G4Element::GetNumberOfElements(); 
  const G4ElementTable* theElementTable = G4Element::GetElementTable();

  theCrossSectionTable = new G4PhysicsTable(nEl);

  G4double AtomicNumber;
  G4double AtomicWeight;
  G4double Value;

  for (G4int j = 0; j < nEl; j++) {
    ptrVector = new G4PhysicsLogVector(LowestKineticEnergy,
                                       HighestKineticEnergy, TotBin);
    AtomicNumber = (*theElementTable)[j]->GetZ();
    AtomicWeight = (*theElementTable)[j]->GetA();

    for (G4int i = 0; i <= TotBin; ++i) {
      lowEdgeEnergy = ptrVector->Energy(i);
      Value = ComputeMicroscopicCrossSection(lowEdgeEnergy,
                                             AtomicNumber, AtomicWeight);
      ptrVector->PutValue(i,Value);
    }

    theCrossSectionTable->insertAt(j, ptrVector);
  }

  // Build (Z,element) look-up table (for use in GetZandACrossSection) 
  for (G4int j = 0; j < nEl; j++) {
    G4int ZZ = G4int((*theElementTable)[j]->GetZ());
    zelMap[ZZ] = j;
  }
}

G4double G4KokoulinMuonNuclearXS::
ComputeMicroscopicCrossSection(G4double KineticEnergy,
                               G4double AtomicNumber, G4double AtomicWeight)
{
  // Calculate cross section (differential in muon incident kinetic energy) by 
  // integrating the double differential cross section over the energy loss

  const G4double xgi[] = {0.0199,0.1017,0.2372,0.4083,0.5917,0.7628,0.8983,0.9801};
  const G4double wgi[] = {0.0506,0.1112,0.1569,0.1813,0.1813,0.1569,0.1112,0.0506};
  const G4double ak1 = 6.9;
  const G4double ak2 = 1.0;

  G4double Mass = G4MuonMinus::MuonMinus()->GetPDGMass();

  G4double CrossSection = 0.0;
  if (AtomicNumber < 1.) return CrossSection;
  if (KineticEnergy <= CutFixed) return CrossSection; 

  G4double epmin = CutFixed;
  G4double epmax = KineticEnergy + Mass - 0.5*proton_mass_c2;
  if (epmax <= epmin) return CrossSection; // NaN bug correction

  G4double aaa = std::log(epmin) ;
  G4double bbb = std::log(epmax) ;
  G4int kkk = G4int((bbb-aaa)/ak1 +ak2) ;
  G4double hhh = (bbb-aaa)/kkk ;
  G4double epln;
  G4double ep;
  G4double x;

  for (G4int l = 0; l < kkk; l++) {
    x = aaa + hhh*l;
    for (G4int ll = 0; ll < 8; ll++) {
      epln=x+xgi[ll]*hhh;
      ep = std::exp(epln);
      CrossSection += ep*wgi[ll]*ComputeDDMicroscopicCrossSection(KineticEnergy,  
                                                AtomicNumber, AtomicWeight, ep);
    }
  }

  CrossSection *= hhh ;
  if (CrossSection < 0.) CrossSection = 0.;
  return CrossSection;
}

G4double G4KokoulinMuonNuclearXS::
ComputeDDMicroscopicCrossSection(G4double KineticEnergy,
                                 G4double, G4double AtomicWeight,
                                 G4double epsilon)
{
  // Calculate the double-differential microscopic cross section (in muon
  // incident kinetic energy and energy loss) using the cross section formula
  // of R.P. Kokoulin (18/01/98)

  const G4double alam2 = 0.400*GeV*GeV;
  const G4double alam  = 0.632456*GeV;
  const G4double coeffn = fine_structure_const/pi;   

  G4double ParticleMass = G4MuonMinus::MuonMinus()->GetPDGMass();
  G4double TotalEnergy = KineticEnergy + ParticleMass;

  G4double DCrossSection = 0.;

  if ((epsilon >= TotalEnergy - 0.5*proton_mass_c2) ||
      (epsilon <= CutFixed) ) return DCrossSection;

  G4double ep = epsilon/GeV;
  G4double a = AtomicWeight/(g/mole);
  G4double aeff = 0.22*a+0.78*std::exp(0.89*std::log(a));       //shadowing 
  G4double sigph = (49.2+11.1*std::log(ep)+151.8/std::sqrt(ep))*microbarn; 
  
  G4double v = epsilon/TotalEnergy;
  G4double v1 = 1.-v;
  G4double v2 = v*v;
  G4double mass2 = ParticleMass*ParticleMass;

  G4double up = TotalEnergy*TotalEnergy*v1/mass2*(1.+mass2*v2/(alam2*v1));
  G4double down = 1.+epsilon/alam*(1.+alam/(2.*proton_mass_c2)+epsilon/alam);

  DCrossSection = coeffn*aeff*sigph/epsilon*
                  (-v1+(v1+0.5*v2*(1.+2.*mass2/alam2))*std::log(up/down));

  if (DCrossSection < 0.) DCrossSection = 0.; 
  return DCrossSection;
}

G4double G4KokoulinMuonNuclearXS::
GetElementCrossSection(const G4DynamicParticle* aPart,
		       G4int ZZ, const G4Material*)
{
  G4int j = zelMap[ZZ];
  return (*theCrossSectionTable)[j]->Value(aPart->GetKineticEnergy());
}

