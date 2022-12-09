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


#include "G4MuNeutrinoNucleusTotXsc.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4HadTmpUtil.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4ElementVector.hh"

#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"

using namespace std;
using namespace CLHEP;

G4MuNeutrinoNucleusTotXsc::G4MuNeutrinoNucleusTotXsc()
 : G4VCrossSectionDataSet("NuMuNuclTotXsc")
{
  fCofXsc = 1.e-38*cm2/GeV;

  // G4cout<<"fCofXsc = "<<fCofXsc*GeV/cm2<<" cm2/GeV"<<G4endl;

  // PDG2016: sin^2 theta Weinberg

  fSin2tW = 0.23129; // 0.2312;

  // 9 <-> 6, 5/9 or 5/6 ?

  fCofS = 5.*fSin2tW*fSin2tW/9.;
  fCofL = 1. - fSin2tW + fCofS;

  // G4cout<<"fCosL = "<<fCofL<<", fCofS = "<<fCofS<<G4endl;

  fCutEnergy = 0.; // default value

  fBiasingFactor = 1.; // default as physics
  fEmc = 0.2*GeV;
  fIndex = 50;

  fTotXsc = 0.;
  fCcTotRatio = 0.75; // from nc/cc~0.33 ratio
  fCcFactor = fNcFactor = 1.;
  fQEratio = 0.5; // mean in the 1 GeV range

  // theMuonMinus = G4MuonMinus::MuonMinus(); 
  // theMuonPlus  = G4MuonPlus::MuonPlus(); 
}

G4MuNeutrinoNucleusTotXsc::~G4MuNeutrinoNucleusTotXsc() 
{}

//////////////////////////////////////////////////////

G4bool 
G4MuNeutrinoNucleusTotXsc::IsIsoApplicable( const G4DynamicParticle* aPart, G4int, G4int, const G4Element*, const G4Material*)
{
  G4bool result  = false;
  G4String pName = aPart->GetDefinition()->GetParticleName();
  G4double tKin = aPart->GetKineticEnergy();
  
  if(      ( pName == "nu_mu"   || pName == "anti_nu_mu")  && tKin >= fEmc ) 
  {
    result = true;
  }
  return result;
}

//////////////////////////////////////

G4double G4MuNeutrinoNucleusTotXsc::GetElementCrossSection(const G4DynamicParticle* part,
					       G4int Z,    const G4Material* mat )
{
  G4int Zi(0);
  size_t i(0), j(0);
  const G4ElementVector* theElementVector = mat->GetElementVector();
  
  for ( i = 0; i < theElementVector->size(); ++i )
  {
    Zi = (*theElementVector)[i]->GetZasInt();
    if( Zi == Z ) break;
  }
  const G4Element* elm = (*theElementVector)[i];
  size_t nIso = elm->GetNumberOfIsotopes();    
  G4double fact = 0.0;
  G4double xsec = 0.0;
  const G4Isotope* iso = nullptr;       
  const G4IsotopeVector* isoVector = elm->GetIsotopeVector();
  const G4double* abundVector = elm->GetRelativeAbundanceVector();

  for (j = 0; j<nIso; ++j)
  {
    iso = (*isoVector)[j];
    G4int A = iso->GetN();
    
    if( abundVector[j] > 0.0 && IsIsoApplicable(part, Z, A, elm, mat) )
    {
      fact += abundVector[j];
      xsec += abundVector[j]*GetIsoCrossSection( part, Z, A, iso, elm, mat);
    }
  }
  if( fact > 0.0) { xsec /= fact; }
  return xsec;
}

////////////////////////////////////////////////////
//
//

G4double G4MuNeutrinoNucleusTotXsc::GetIsoCrossSection(const G4DynamicParticle* aPart, G4int Z, G4int A,  
			      const G4Isotope*, const G4Element*, const G4Material* )
{
  fCcFactor   = fNcFactor = 1.;
  fCcTotRatio = 0.25;

  G4double ccnuXsc, ccanuXsc, ncXsc, totXsc(0.);

  G4double energy  = aPart->GetTotalEnergy();
  G4String pName   = aPart->GetDefinition()->GetParticleName();

  G4int index = GetEnergyIndex(energy);

  if( index >= fIndex )
  {
    G4double pm = proton_mass_c2;
    G4double s2 = 2.*energy*pm+pm*pm;
    G4double aa = 1.;
    G4double bb = 1.085;
    G4double mw = 80.385*GeV;
    fCcFactor   = bb/(1.+ aa*s2/mw/mw);

    G4double mz = 91.1876*GeV;
    fNcFactor   =  bb/(1.+ aa*s2/mz/mz);
  }
  ccnuXsc  = GetNuMuTotCsXsc(index, energy, Z, A);
  ccnuXsc *= fCcFactor;
  ccanuXsc = GetANuMuTotCsXsc(index, energy, Z, A);
  ccanuXsc *= fCcFactor;

  if( pName == "nu_mu" )
  {
    ncXsc = fCofL*ccnuXsc + fCofS*ccanuXsc;
    ncXsc *= fNcFactor/fCcFactor;
    totXsc = ccnuXsc + ncXsc;
    if( totXsc > 0.) fCcTotRatio = ccnuXsc/totXsc;
  }
  else if( pName == "anti_nu_mu" )
  {
    ncXsc = fCofL*ccanuXsc + fCofS*ccnuXsc;
    ncXsc *= fNcFactor/fCcFactor;
    totXsc = ccanuXsc + ncXsc;
    if( totXsc > 0.) fCcTotRatio = ccanuXsc/totXsc;
  }
  else return totXsc;

  totXsc *= fCofXsc; 
  totXsc *= energy; 
  // totXsc *= A;  // incoherent sum over  all isotope nucleons

  totXsc *= fBiasingFactor; // biasing up, if set >1

  fTotXsc = totXsc;

  return totXsc;
}

/////////////////////////////////////////////////////
//
// Return index of nu/anu energy array corresponding to the neutrino energy

G4int G4MuNeutrinoNucleusTotXsc::GetEnergyIndex(G4double energy)
{
  G4int i, eIndex = 0;

  for( i = 0; i < fIndex; i++)
  {
    if( energy <= fNuMuEnergy[i]*GeV ) 
    {
      eIndex = i;
      break;
    }
  }
  if( i >= fIndex ) eIndex = i;
  // G4cout<<"eIndex = "<<eIndex<<G4endl;
  return eIndex;
}

/////////////////////////////////////////////////////
//
// nu_mu xsc for index-1, index linear over energy

G4double G4MuNeutrinoNucleusTotXsc::GetNuMuTotCsXsc(G4int index, G4double energy, G4int zz, G4int aa)
{
  G4double xsc(0.), qexsc(0.), inxsc(0.);
  G4int nn = aa - zz;
  if(nn < 1) nn = 0;

  // if( index <= 0 || energy < theMuonMinus->GetPDGMass() ) xsc = aa*fNuMuInXsc[0] + nn*fNuMuQeXsc[0];
  if( index <= 0 || energy < fEmc ) xsc = aa*fNuMuInXsc[0] + nn*fNuMuQeXsc[0];
  else if (index >= fIndex) xsc = aa*fNuMuInXsc[fIndex-1] + nn*fNuMuQeXsc[fIndex-1];
  else
  {
    G4double x1 = fNuMuEnergy[index-1]*GeV;
    G4double x2 = fNuMuEnergy[index]*GeV;
    G4double y1 = fNuMuInXsc[index-1];
    G4double y2 = fNuMuInXsc[index];
    G4double z1 = fNuMuQeXsc[index-1];
    G4double z2 = fNuMuQeXsc[index];

    if(x1 >= x2) return aa*fNuMuInXsc[index] + nn*fNuMuQeXsc[index];
    else
    {
      G4double angle = (y2-y1)/(x2-x1);
      inxsc = y1 + (energy-x1)*angle;
      angle = (z2-z1)/(x2-x1);
      qexsc = z1 + (energy-x1)*angle; 
      qexsc *= nn;   
      xsc = inxsc*aa + qexsc;

      if( xsc > 0.) fQEratio = qexsc/xsc;
    }
  }
  return xsc;
}

/////////////////////////////////////////////////////
//
// anu_mu xsc for index-1, index linear over energy

G4double G4MuNeutrinoNucleusTotXsc::GetANuMuTotCsXsc(G4int index, G4double energy, G4int zz, G4int aa)
{
  G4double xsc(0.), qexsc(0.), inxsc(0.);

  // if( index <= 0 || energy < theMuonPlus->GetPDGMass() ) xsc = aa*fANuMuInXsc[0] + zz*fANuMuQeXsc[0];
  if( index <= 0 || energy < fEmc ) xsc = aa*fANuMuInXsc[0] + zz*fANuMuQeXsc[0];
  else if (index >= fIndex) xsc = aa*fANuMuInXsc[fIndex-1] + zz*fANuMuQeXsc[fIndex-1];
  else
  {
    G4double x1 = fNuMuEnergy[index-1]*GeV;
    G4double x2 = fNuMuEnergy[index]*GeV;
    G4double y1 = fANuMuInXsc[index-1];
    G4double y2 = fANuMuInXsc[index];
    G4double z1 = fANuMuQeXsc[index-1];
    G4double z2 = fANuMuQeXsc[index];

    if( x1 >= x2 ) return aa*fANuMuInXsc[index] + zz*fANuMuQeXsc[index];
    else
    {
      G4double angle = (y2-y1)/(x2-x1);
      inxsc = y1 + (energy-x1)*angle;

      angle = (z2-z1)/(x2-x1);
      qexsc = z1 + (energy-x1)*angle;
      qexsc *= zz;    
      xsc = inxsc*aa + qexsc;

      if( xsc > 0.) fQEratio = qexsc/xsc;
    }
  }
  return xsc;
}

////////////////////////////////////////////////////////
//
// return fNuMuTotXsc[index] if the index is in the array range

G4double G4MuNeutrinoNucleusTotXsc::GetNuMuTotCsArray( G4int index)
{
  if( index >= 0 && index < fIndex) return fNuMuInXsc[index] + fNuMuQeXsc[index];
  else 
  {
    G4cout<<"Improper index of fNuMuTotXsc array"<<G4endl;
    return 0.;
  }
}

////////////////////////////////////////////////////////
//
// return fANuMuTotXsc[index] if the index is in the array range

G4double G4MuNeutrinoNucleusTotXsc::GetANuMuTotCsArray( G4int index)
{
  if( index >= 0 && index < fIndex) return fANuMuInXsc[index] + fANuMuQeXsc[index];
  else 
  {
    G4cout<<"Improper index of fANuMuTotXsc array"<<G4endl;
    return 0.;
  }
}


///////////////////////////////////////////////////////
//
// E_nu in GeV, ( Eth = 111.603 MeV, EthW = 330.994 MeV) 

const G4double G4MuNeutrinoNucleusTotXsc::fNuMuEnergy[50] = 
{
  0.12, 0.141136, 0.165996, 0.195233, 0.229621, 
  0.270066, 0.317634, 0.373581, 0.439382, 0.516773, 
  0.607795, 0.714849, 0.84076, 0.988848, 1.16302, 
  1.36787, 1.6088, 1.89217, 2.22545, 2.61743, 
  3.07845, 3.62068, 4.25841, 5.00847, 5.89065, 
  6.9282, 8.14851, 9.58376, 11.2718, 13.2572, 
  15.5922, 18.3386, 21.5687, 25.3677, 29.8359, 
  35.0911, 41.2719, 48.5413, 57.0912, 67.147, 
  78.974, 92.8842, 109.244, 128.486, 151.117, 
  177.735, 209.04, 245.86, 289.164, 340.097 };

////////////////////////////////////////////////////
//
// XS/E arrays in 10^-38cm2/GeV

const G4double G4MuNeutrinoNucleusTotXsc::fNuMuInXsc[50] = 
{
  0, 0, 0, 0, 0, 
  0, 0, 0.0166853, 0.0649693, 0.132346, 
  0.209102, 0.286795, 0.3595, 0.423961,  0.479009, 
  0.524797, 0.562165, 0.592225, 0.61612,  0.63491, 
  0.649524, 0.660751, 0.669245, 0.675546, 0.680092, 
  0.683247, 0.685307, 0.686521, 0.687093, 0.687184, 
  0.686919, 0.686384, 0.685631, 0.684689, 0.68357, 
  0.682275, 0.680806, 0.67917, 0.677376, 0.675442, 
  0.673387, 0.671229, 0.668985, 0.666665, 0.664272, 
  0.661804, 0.65925, 0.656593, 0.65381, 0.650871    }; 

const G4double G4MuNeutrinoNucleusTotXsc::fNuMuQeXsc[50] = 
{
  0.20787, 0.411055, 0.570762, 0.705379, 0.814702, 
  0.89543, 0.944299, 0.959743, 0.942906, 0.897917, 
  0.831331, 0.750948, 0.66443, 0.578191, 0.496828, 
  0.423071, 0.358103, 0.302016, 0.254241, 0.213889, 
  0.179971, 0.151527, 0.12769, 0.107706, 0.0909373, 
  0.0768491, 0.0649975, 0.0550143, 0.0465948, 0.0394861, 
  0.0334782, 0.0283964, 0.0240945, 0.0204506, 0.0173623, 
  0.0147437, 0.0125223, 0.0106374, 0.00903737, 0.00767892, 
  0.00652531, 0.00554547, 0.0047131, 0.0040059, 0.003405, 
  0.00289436, 0.00246039, 0.00209155, 0.00177804, 0.00151152  }; 



//////////////////////////////////////////////////////////////////

const G4double G4MuNeutrinoNucleusTotXsc::fANuMuInXsc[50] = 
{
  0, 0, 0, 0, 0, 
  0, 0, 0.00437363, 0.0161485, 0.0333162, 
  0.0557621, 0.0814548, 0.108838, 0.136598, 0.163526, 
  0.188908, 0.212041, 0.232727, 0.250872, 0.26631, 
  0.279467, 0.290341, 0.299177, 0.306299, 0.311864, 
  0.316108, 0.319378, 0.321892, 0.323583, 0.324909, 
  0.325841, 0.326568, 0.327111, 0.327623, 0.32798, 
  0.328412, 0.328704, 0.328988, 0.329326, 0.329559, 
  0.329791, 0.330051, 0.330327, 0.33057, 0.330834, 
  0.331115, 0.331416, 0.331678, 0.33192, 0.332124   };

////////////////////////////////////////////////////////////////// 

const G4double G4MuNeutrinoNucleusTotXsc::fANuMuQeXsc[50] = 
{
  0.0770264, 0.138754, 0.177006, 0.202417, 0.21804, 
  0.225742, 0.227151, 0.223805, 0.21709, 0.208137, 
  0.197763, 0.186496, 0.174651, 0.162429, 0.14999, 
  0.137498, 0.125127, 0.113057, 0.101455, 0.0904642, 
  0.0801914, 0.0707075, 0.0620483, 0.0542192, 0.0472011, 
  0.0409571, 0.0354377, 0.0305862, 0.0263422, 0.0226451, 
  0.0194358, 0.0166585, 0.0142613, 0.0121968, 0.0104221, 
  0.00889912, 0.00759389, 0.00647662, 0.00552119, 0.00470487, 
  0.00400791, 0.00341322, 0.00290607, 0.00247377, 0.0021054, 
  0.00179162, 0.00152441, 0.00129691, 0.00110323, 0.000938345   }; 

////////////////////

