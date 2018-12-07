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
// Based on G4NucleonNuclearCrossSection class
//
//

#include "G4ComponentBarNucleonNucleusXsc.hh"
#include "G4SystemOfUnits.hh"
#include "G4DynamicParticle.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4Pow.hh"
#include "G4BarashenkovData.hh"
#include "G4NistManager.hh"

using namespace std;

///////////////////////////////////////////////////////////////////////////////

G4double G4ComponentBarNucleonNucleusXsc::theA[93] = {0.0};
G4double G4ComponentBarNucleonNucleusXsc::A75[93]  = {0.0};

G4ComponentBarNucleonNucleusXsc::G4ComponentBarNucleonNucleusXsc()
 : G4VComponentCrossSection("G4ComponentBarNucleonNucleusXsc"),
   fTotalXsc(0.0), fInelasticXsc(0.0), fElasticXsc(0.0)
{
  theNeutron = G4Neutron::Neutron();
  theProton  = G4Proton::Proton();
  
  // He, Be, C
  thePimData.push_back(new G4PiData(he_m_t, he_m_in, e1, 44));
  thePipData.push_back(new G4PiData(he_m_t, he_p_in, e1, 44));

  thePimData.push_back(new G4PiData(be_m_t, be_m_in, e1, 44));
  thePipData.push_back(new G4PiData(be_m_t, be_p_in, e1, 44));

  thePimData.push_back(new G4PiData(c_m_t,  c_m_in,  e1, 44));
  thePipData.push_back(new G4PiData(c_m_t,  c_p_in,  e1, 44));

  // N, O, Na
  thePimData.push_back(new G4PiData(n_m_t,  n_m_in,  e2, 44));
  thePipData.push_back(new G4PiData(n_m_t,  n_p_in,  e2, 44));

  thePimData.push_back(new G4PiData(o_m_t,  o_m_in,  e2, 44));
  thePipData.push_back(new G4PiData(o_m_t,  o_p_in,  e2, 44));

  thePimData.push_back(new G4PiData(na_m_t, na_m_in, e2, 44));
  thePipData.push_back(new G4PiData(na_m_t, na_p_in, e2, 44));

  // Al, Si, Ca
  thePimData.push_back(new G4PiData(al_m_t, al_m_in, e3, 45));
  thePipData.push_back(new G4PiData(al_m_t, al_p_in, e3, 45));

  thePimData.push_back(new G4PiData(si_m_t, si_m_in, e3, 45));
  thePipData.push_back(new G4PiData(si_m_t, si_p_in, e3, 45));

  thePimData.push_back(new G4PiData(ca_m_t, ca_m_in, e3, 45));
  thePipData.push_back(new G4PiData(ca_m_t, ca_p_in, e3, 45));

  // Fe, Cu, Mo
  thePimData.push_back(new G4PiData(fe_m_t, fe_m_in, e4, 47));
  thePipData.push_back(new G4PiData(fe_m_t, fe_p_in, e4, 47));

  thePimData.push_back(new G4PiData(cu_m_t, cu_m_in, e4, 47));
  thePipData.push_back(new G4PiData(cu_m_t, cu_p_in, e4, 47));

  thePimData.push_back(new G4PiData(mo_m_t, mo_m_in, e4, 47));
  thePipData.push_back(new G4PiData(mo_m_t, mo_p_in, e4, 47));

  // Cd, Sn, W
  thePimData.push_back(new G4PiData(cd_m_t, cd_m_in, e5, 48));
  thePipData.push_back(new G4PiData(cd_m_t, cd_p_in, e5, 48));

  thePimData.push_back(new G4PiData(sn_m_t, sn_m_in, e5, 48));
  thePipData.push_back(new G4PiData(sn_m_t, sn_p_in, e5, 48));

  thePimData.push_back(new G4PiData(w_m_t,  w_m_in,  e5, 48));
  thePipData.push_back(new G4PiData(w_m_t,  w_p_in,  e5, 48));

  // Pb, U
  thePimData.push_back(new G4PiData(pb_m_t, pb_m_in, e6, 46));
  thePipData.push_back(new G4PiData(pb_m_t, pb_p_in, e6, 46));

  thePimData.push_back(new G4PiData(u_m_t,  u_m_in,  e6, 46));
  thePipData.push_back(new G4PiData(u_m_t,  u_p_in,  e6, 46));

  theZ.push_back(2); // He
  theZ.push_back(4); // Be
  theZ.push_back(6); // C
  theZ.push_back(7); // N
  theZ.push_back(8); // O
  theZ.push_back(11); // Na
  theZ.push_back(13); // Al
  theZ.push_back(14); // Si
  theZ.push_back(20); // Ca
  theZ.push_back(26); // Fe
  theZ.push_back(29); // Cu
  theZ.push_back(42); // Mo
  theZ.push_back(48); // Cd
  theZ.push_back(50); // Sn
  theZ.push_back(74); // W
  theZ.push_back(82); // Pb
  theZ.push_back(92); // U
  
  if(A75[0] == 0.0) {
    G4NistManager* nist = G4NistManager::Instance();
    A75[0] = theA[0] = 1.0;
    for(G4int i=1; i<93; ++i) {
      theA[i] = nist->GetAtomicMassAmu(i);
      A75[i] = G4Pow::GetInstance()->A23(theA[i]); // interpolate by square ~ A^(2/3)
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
//

G4ComponentBarNucleonNucleusXsc::~G4ComponentBarNucleonNucleusXsc()
{
  for(auto xsec : thePimData) { delete xsec; }
  for(auto xsec : thePipData) { delete xsec; }
}

////////////////////////////////////////////////////////////////////

G4double G4ComponentBarNucleonNucleusXsc::GetTotalIsotopeCrossSection(
         const G4ParticleDefinition* aParticle,
	 G4double kinEnergy, G4int Z, G4int)
{
  ComputeCrossSections(aParticle, kinEnergy, Z);
  return fTotalXsc;
}

//////////////////////////////////////////////////////////////////////

G4double G4ComponentBarNucleonNucleusXsc::GetTotalElementCrossSection(
         const G4ParticleDefinition* aParticle,
	 G4double kinEnergy, G4int Z, G4double)
{
  ComputeCrossSections(aParticle, kinEnergy, Z);
  return fTotalXsc;
}

////////////////////////////////////////////////////////////////////

G4double G4ComponentBarNucleonNucleusXsc::GetInelasticIsotopeCrossSection(
         const G4ParticleDefinition* aParticle,
	 G4double kinEnergy, G4int Z, G4int)
{
  ComputeCrossSections(aParticle, kinEnergy, Z);
  return fInelasticXsc;
}

/////////////////////////////////////////////////////////////////////

G4double G4ComponentBarNucleonNucleusXsc::GetInelasticElementCrossSection(
         const G4ParticleDefinition* aParticle,
	 G4double kinEnergy, G4int Z, G4double)
{
  ComputeCrossSections(aParticle, kinEnergy, Z);
  return fInelasticXsc;
}

//////////////////////////////////////////////////////////////////

G4double G4ComponentBarNucleonNucleusXsc::GetElasticElementCrossSection(
         const G4ParticleDefinition* aParticle,
	 G4double kinEnergy, G4int Z, G4double)
{
  ComputeCrossSections(aParticle, kinEnergy, Z);
  return fElasticXsc;
}

///////////////////////////////////////////////////////////////////

G4double G4ComponentBarNucleonNucleusXsc::GetElasticIsotopeCrossSection(
         const G4ParticleDefinition* aParticle,
	 G4double kinEnergy, G4int Z, G4int)
{
  ComputeCrossSections(aParticle, kinEnergy, Z);
  return fElasticXsc;
}



////////////////////////////////////////////////////////////////////////////
//

G4bool G4ComponentBarNucleonNucleusXsc::IsElementApplicable(
       const G4DynamicParticle*, G4int Z) 
{
  return (Z > 1);
}

////////////////////////////////////////////////////////////////////////////
//
//

void G4ComponentBarNucleonNucleusXsc::ComputeCrossSections(
     const G4ParticleDefinition* aParticle, G4double kineticEnergy, G4int ZZ) 
{
  G4int Z = std::min(ZZ, 92);
  size_t it = 0;
  size_t itmax = theZ.size() - 1;
  for(; it <= itmax; ++it) { if(Z <= theZ[it]) { break; } }
  if( it > itmax ) { it = itmax; }
  G4int Z1, Z2;
  G4double x1, x2, xt1, xt2;

  std::vector<G4PiData*>* theData = (aParticle == theNeutron) ?  &thePimData : &thePipData;

   if( theZ[it] == Z )
     {
       fInelasticXsc = (*theData)[it]->ReactionXSection(kineticEnergy);
       fTotalXsc = (*theData)[it]->TotalXSection(kineticEnergy);
     }
   else
     {
       if(0 == it) { it = 1; }
       x1  = (*theData)[it-1]->ReactionXSection(kineticEnergy);
       xt1 = (*theData)[it-1]->TotalXSection(kineticEnergy);
       Z1  = theZ[it-1];
       x2  = (*theData)[it]->ReactionXSection(kineticEnergy);
       xt2 = (*theData)[it]->TotalXSection(kineticEnergy);
       Z2  = theZ[it];

       fInelasticXsc = Interpolate(Z1, Z2, Z, x1, x2);
       fTotalXsc = Interpolate(Z1, Z2, Z, xt1, xt2);
     }

   fElasticXsc = std::max(fTotalXsc - fInelasticXsc, 0.0);
}

/////////////////////////////////////////////////////////////////////////////
//

G4double G4ComponentBarNucleonNucleusXsc::
Interpolate(G4int Z1, G4int Z2, G4int Z, G4double x1, G4double x2)
{ 
  // for tabulated data, cross section scales with A^(2/3)
  G4double r1 = x1* A75[Z] / A75[Z1];
  G4double r2 = x2* A75[Z] / A75[Z2];
  G4double alp1 = (theA[Z] - theA[Z1]);
  G4double alp2 = (theA[Z2] - theA[Z]);
  G4double result = (r1*alp2 + r2*alp1)/(alp1 + alp2);
  //       G4cout << "x1/2, z1/2 z" <<x1<<" "<<x2<<" "<<Z1<<" "<<Z2<<" "<<Z<<G4endl;
  //       G4cout << "res1/2 " << r1 <<" " << r2 <<" " << result<< G4endl;
  return result;
}

void
G4ComponentBarNucleonNucleusXsc::CrossSectionDescription(std::ostream& outFile) const
{
  outFile << "G4ComponentBarNucleonNucleusXsc is a variant of the Barashenkov\n"
          << "cross section parameterization to be used of protons and\n"
          << "neutrons on targets heavier than hydrogen.  It is intended for\n"
          << "use as a cross section component and is currently used by\n"
          << "G4BGGNucleonInelasticXS.  It is valid for incident energies up\n"
          << "to 1 TeV.\n"; 
}

