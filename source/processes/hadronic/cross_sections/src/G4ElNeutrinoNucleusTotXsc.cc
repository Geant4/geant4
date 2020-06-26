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
// 24.04.20 V. Grichine
//
// (nu_e,anti_nu_e)-nucleus xsc



#include "G4ElNeutrinoNucleusTotXsc.hh"
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

#include "G4Electron.hh"
#include "G4Positron.hh"

using namespace std;
using namespace CLHEP;

G4ElNeutrinoNucleusTotXsc::G4ElNeutrinoNucleusTotXsc()
 : G4VCrossSectionDataSet("NuElNuclTotXsc")
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

  fIndex = 50;

  fTotXsc = 0.;
  fCcTotRatio = 0.75; // from nc/cc~0.33 ratio
  fCcFactor = fNcFactor = 1.;

  theElectron = G4Electron::Electron(); 
  thePositron  = G4Positron::Positron(); 
}

G4ElNeutrinoNucleusTotXsc::~G4ElNeutrinoNucleusTotXsc() 
{}

//////////////////////////////////////////////////////

/*
G4bool 
G4ElNeutrinoNucleusTotXsc::IsIsoApplicable( const G4DynamicParticle* aPart, G4int, G4int, const G4Element*, const G4Material*)
{
  G4bool result  = false;
  G4String pName = aPart->GetDefinition()->GetParticleName();

  if(      pName == "nu_e"   || pName == "anti_nu_e"  ) 
  {
    result = true;
  }
  return result;
}

//////////////////////////////////////

G4double G4ElNeutrinoNucleusTotXsc::GetElementCrossSection(const G4DynamicParticle* part,
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
*/

////////////////////////////////////////////////////
//
//

G4double G4ElNeutrinoNucleusTotXsc::GetIsoCrossSection(const G4DynamicParticle* aPart, G4int, G4int A,  
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
  ccnuXsc  = GetNuElTotCsXsc(index, energy);
  ccnuXsc *= fCcFactor;
  ccanuXsc = GetANuElTotCsXsc(index, energy);
  ccanuXsc *= fCcFactor;

  if( pName == "nu_e")
  {
    ncXsc = fCofL*ccnuXsc + fCofS*ccanuXsc;
    ncXsc *= fNcFactor/fCcFactor;
    totXsc = ccnuXsc + ncXsc;
    if( totXsc > 0.) fCcTotRatio = ccnuXsc/totXsc;
  }
  else if( pName == "anti_nu_e")
  {
    ncXsc = fCofL*ccanuXsc + fCofS*ccnuXsc;
    ncXsc *= fNcFactor/fCcFactor;
    totXsc = ccanuXsc + ncXsc;
    if( totXsc > 0.) fCcTotRatio = ccanuXsc/totXsc;
  }
  else return totXsc;

  totXsc *= fCofXsc; 
  totXsc *= energy; 
  totXsc *= A;  // incoherent sum over  all isotope nucleons

  totXsc *= fBiasingFactor; // biasing up, if set >1

  fTotXsc = totXsc;

  return totXsc;
}

/////////////////////////////////////////////////////
//
// Return index of nu/anu energy array corresponding to the neutrino energy

G4int G4ElNeutrinoNucleusTotXsc::GetEnergyIndex(G4double energy)
{
  G4int i, eIndex = 0;

  for( i = 0; i < fIndex; i++)
  {
    if( energy <= fNuElEnergy[i]*GeV ) 
    {
      eIndex = i;
      break;
    }
  }
  if( i >= fIndex-1 ) eIndex = fIndex-1;
  // G4cout<<"eIndex = "<<eIndex<<G4endl;
  return eIndex;
}

/////////////////////////////////////////////////////
//
// nu_e xsc for index-1, index linear over energy

G4double G4ElNeutrinoNucleusTotXsc::GetNuElTotCsXsc(G4int index, G4double energy)
{
  G4double xsc(0.);

  if( index <= 0 || energy < theElectron->GetPDGMass() ) xsc = fNuElTotXsc[0];
  else if (index >= fIndex) xsc = fNuElTotXsc[fIndex-1];
  else
  {
    G4double x1 = fNuElEnergy[index-1]*GeV;
    G4double x2 = fNuElEnergy[index]*GeV;
    G4double y1 = fNuElTotXsc[index-1];
    G4double y2 = fNuElTotXsc[index];

    if(x1 >= x2) return fNuElTotXsc[index];
    else
    {
      G4double angle = (y2-y1)/(x2-x1);
      xsc = y1 + (energy-x1)*angle;
    }
  }
  return xsc;
}

/////////////////////////////////////////////////////
//
// anu_e xsc for index-1, index linear over energy

G4double G4ElNeutrinoNucleusTotXsc::GetANuElTotCsXsc(G4int index, G4double energy)
{
  G4double xsc(0.);

  if( index <= 0 || energy < thePositron->GetPDGMass() ) xsc = fANuElTotXsc[0];
  else if (index >= fIndex) xsc = fANuElTotXsc[fIndex-1];
  else
  {
    G4double x1 = fNuElEnergy[index-1]*GeV;
    G4double x2 = fNuElEnergy[index]*GeV;
    G4double y1 = fANuElTotXsc[index-1];
    G4double y2 = fANuElTotXsc[index];

    if( x1 >= x2 ) return fANuElTotXsc[index];
    else
    {
      G4double angle = (y2-y1)/(x2-x1);
      xsc = y1 + (energy-x1)*angle;
    }
  }
  return xsc;
}

////////////////////////////////////////////////////////
//
// return fNuElTotXsc[index] if the index is in the array range

G4double G4ElNeutrinoNucleusTotXsc::GetNuElTotCsArray( G4int index)
{
  if( index >= 0 && index < fIndex) return fNuElTotXsc[index];
  else 
  {
    G4cout<<"Inproper index of fNuElTotXsc array"<<G4endl;
    return 0.;
  }
}

////////////////////////////////////////////////////////
//
// return fANuElTotXsc[index] if the index is in the array range

G4double G4ElNeutrinoNucleusTotXsc::GetANuElTotCsArray( G4int index)
{
  if( index >= 0 && index < fIndex) return fANuElTotXsc[index];
  else 
  {
    G4cout<<"Inproper index of fANuElTotXsc array"<<G4endl;
    return 0.;
  }
}


///////////////////////////////////////////////////////
//
// E_nu in GeV

const G4double G4ElNeutrinoNucleusTotXsc::fNuElEnergy[50] = 
{
  0.000561138, 0.000735091, 0.000962969, 0.00126149, 0.00165255, 
  0.00216484, 0.00283594, 0.00371508, 0.00486676, 0.00637546, 
  0.00835185, 0.0109409, 0.0143326, 0.0187757, 0.0245962, 
  0.032221, 0.0422095, 0.0552945, 0.0724358, 0.0948908, 
  0.124307, 0.162842, 0.213323, 0.279453, 0.366084, 
  0.47957, 0.628237, 0.82299, 1.07812, 1.41233, 
  1.85016, 2.42371, 3.17505, 4.15932, 5.44871, 
  7.13781, 9.35053, 12.2492, 16.0464, 21.0208, 
  27.5373, 36.0739, 47.2568, 61.9064, 81.0973, 
  106.238, 139.171, 182.314, 238.832, 312.869 
};

/////////////////////////////////////////////////////////////
//
// nu_e CC xsc_tot/E_nu, in 10^-38 cm2/GeV

const G4double G4ElNeutrinoNucleusTotXsc::fNuElTotXsc[50] = 
{
  0.0026484, 0.00609503, 0.00939421, 0.0132163, 0.0178983, 
  0.0237692, 0.0312066, 0.0406632, 0.0526867, 0.0679357, 
  0.0871913, 0.111359, 0.141458, 0.178584, 0.223838, 
  0.27822, 0.342461, 0.416865, 0.501361, 0.596739, 
  0.713623, 0.905749, 1.20718, 1.52521, 1.75286, 
  1.82072, 1.67119, 1.50074, 1.3077, 1.14923, 
  1.0577, 0.977911, 0.918526, 0.792889, 0.702282, 
  0.678615, 0.687099, 0.725167, 0.706795, 0.678045, 
  0.649791, 0.651328, 0.651934, 0.658062, 0.660659, 
  0.662534, 0.662601, 0.660261, 0.656724, 0.65212
};



/////////////////////////////////////////////////////////////
//
// anu_e CC xsc_tot/E_nu, in 10^-38 cm2/GeV

const G4double G4ElNeutrinoNucleusTotXsc::fANuElTotXsc[50] = 
{
  0.00103385, 0.00237807, 0.00366358, 0.00515192, 0.00697434, 
  0.00925859, 0.0121508, 0.0158252, 0.0204908, 0.0263959, 
  0.0338304, 0.0431234, 0.0546346, 0.068735, 0.0857738, 
  0.106025, 0.129614, 0.15643, 0.186063, 0.21784, 
  0.251065, 0.28525, 0.319171, 0.348995, 0.369448, 
  0.378165, 0.377353, 0.371224, 0.363257, 0.355433, 
  0.348618, 0.343082, 0.338825, 0.33574, 0.333684, 
  0.332504, 0.332052, 0.332187, 0.332781, 0.333716, 
  0.33489, 0.336213, 0.337608, 0.339008, 0.340362, 
  0.341606, 0.342706, 0.343628, 0.344305, 0.344675
};
