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

#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"

using namespace std;
using namespace CLHEP;

G4MuNeutrinoNucleusTotXsc::G4MuNeutrinoNucleusTotXsc()
 : G4VCrossSectionDataSet("NuElectronTotXsc")
{
  fCofXsc = 1.e-38*cm2/GeV;

  // G4cout<<"fCofXsc = "<<fCofXsc*GeV/cm2<<" cm2/GeV"<<G4endl;

  // PDG2016: sin^2 theta Weinberg

  fSin2tW = 0.23129; // 0.2312;

  // 9 <-> 6, 5/9 or 5/6 ?

  fCofS = 5.*fSin2tW*fSin2tW/9.;
  fCofL = 1. - fSin2tW + fCofS;

  G4cout<<"fCosL = "<<fCofL<<", fCofS = "<<fCofS<<G4endl;

  fCutEnergy = 0.; // default value

  fBiasingFactor = 1.; // default as physics

  fIndex = 50;

  theMuonMinus = G4MuonMinus::MuonMinus(); 
  theMuonPlus  = G4MuonPlus::MuonPlus(); 
}

G4MuNeutrinoNucleusTotXsc::~G4MuNeutrinoNucleusTotXsc() 
{}

//////////////////////////////////////////////////////

G4bool 
G4MuNeutrinoNucleusTotXsc::IsElementApplicable( const G4DynamicParticle* aPart, G4int, const G4Material*)
{
  G4bool result  = false;
  G4String pName = aPart->GetDefinition()->GetParticleName();

  if(      pName == "nu_mu"   || pName == "anti_nu_mu"  ) 
  {
    result = true;
  }
  return result;
}

////////////////////////////////////////////////////
//
//

G4double G4MuNeutrinoNucleusTotXsc::GetIsoCrossSection(const G4DynamicParticle* aPart, G4int, G4int A,  
			      const G4Isotope*, const G4Element*, const G4Material* )
{
  G4double ccnuXsc, ccanuXsc, ncXsc, totXsc(0.);

  G4double energy  = aPart->GetTotalEnergy();
  G4String pName   = aPart->GetDefinition()->GetParticleName();

  G4int index = GetEnergyIndex(energy);

  ccnuXsc  = GetNuMuTotCsXsc(index, energy);
  ccanuXsc = GetANuMuTotCsXsc(index, energy);

  if( pName == "nu_mu")
  {
    ncXsc = fCofL*ccnuXsc + fCofS*ccanuXsc;
    totXsc = ccnuXsc + ncXsc;
  }
  else if( pName == "anti_nu_mu")
  {
    ncXsc = fCofL*ccanuXsc + fCofS*ccnuXsc;
    totXsc = ccanuXsc + ncXsc;
  }
  else return totXsc;

  // totXsc -= ncXsc; // to test experimentally available cc part

  totXsc *= fCofXsc; //*energy;
  totXsc *= energy; //  + 0.5*emass;
  totXsc *= A;  // incoherent sum over  all isotope nucleons

  totXsc *= fBiasingFactor; // biasing up, if set >1


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

G4double G4MuNeutrinoNucleusTotXsc::GetNuMuTotCsXsc(G4int index, G4double energy)
{
  G4double xsc(0.);

  if( index <= 0 || energy < theMuonMinus->GetPDGMass() ) xsc = 0.;
  else if (index >= fIndex) xsc = fNuMuTotXsc[fIndex-1];
  else
  {
    G4double x1 = fNuMuEnergy[index-1]*GeV;
    G4double x2 = fNuMuEnergy[index]*GeV;
    G4double y1 = fNuMuTotXsc[index-1];
    G4double y2 = fNuMuTotXsc[index];

    if(x1 >= x2) return fNuMuTotXsc[index];
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
// anu_mu xsc for index-1, index linear over energy

G4double G4MuNeutrinoNucleusTotXsc::GetANuMuTotCsXsc(G4int index, G4double energy)
{
  G4double xsc(0.);

  if( index <= 0 || energy < theMuonPlus->GetPDGMass() ) xsc = 0.;
  else if (index >= fIndex) xsc = fANuMuTotXsc[fIndex-1];
  else
  {
    G4double x1 = fNuMuEnergy[index-1]*GeV;
    G4double x2 = fNuMuEnergy[index]*GeV;
    G4double y1 = fANuMuTotXsc[index-1];
    G4double y2 = fANuMuTotXsc[index];

    if( x1 >= x2 ) return fANuMuTotXsc[index];
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
// return fNuMuTotXsc[index] if the index is in the array range

G4double G4MuNeutrinoNucleusTotXsc::GetNuMuTotCsArray( G4int index)
{
  if( index >= 0 && index < fIndex) return fNuMuTotXsc[index];
  else 
  {
    G4cout<<"Inproper index of fNuMuTotXsc array"<<G4endl;
    return 0.;
  }
}

////////////////////////////////////////////////////////
//
// return fANuMuTotXsc[index] if the index is in the array range

G4double G4MuNeutrinoNucleusTotXsc::GetANuMuTotCsArray( G4int index)
{
  if( index >= 0 && index < fIndex) return fANuMuTotXsc[index];
  else 
  {
    G4cout<<"Inproper index of fANuMuTotXsc array"<<G4endl;
    return 0.;
  }
}


///////////////////////////////////////////////////////
//
// E_nu in GeV

const G4double G4MuNeutrinoNucleusTotXsc::fNuMuEnergy[50] = 
{
  0.112103, 0.117359, 0.123119, 0.129443, 0.136404, 
  0.144084, 0.152576, 0.161991, 0.172458, 0.184126, 
  0.197171, 0.211801, 0.228261, 0.24684, 0.267887, 
  0.291816, 0.319125, 0.350417, 0.386422, 0.428032, 
  0.47634, 0.532692, 0.598756, 0.676612, 0.768868, 
  0.878812, 1.01062, 1.16963, 1.36271, 1.59876, 
  1.88943, 2.25002, 2.70086, 3.26916, 3.99166, 
  4.91843, 6.11836, 7.6872, 9.75942, 12.5259, 
  16.2605, 21.3615, 28.4141, 38.2903, 52.3062, 
  72.4763, 101.93, 145.6, 211.39, 312.172};

/////////////////////////////////////////////////////////////
//
// nu_mu CC xsc_tot/E_nu, in 10^-38 cm2/GeV

const G4double G4MuNeutrinoNucleusTotXsc::fNuMuTotXsc[50] = 
{
  0.0716001, 0.241013, 0.337793, 0.416033, 0.484616, 
  0.546945, 0.604661, 0.658623, 0.709277, 0.756815, 
  0.801263, 0.842519, 0.880396, 0.914642, 0.944968, 
  0.971069, 0.992655, 1.00947, 1.02133, 1.02812, 
  1.02987, 1.02671, 1.01892, 1.00689, 0.991167, 
  0.972618, 0.951518, 0.928806, 0.90521, 0.881404, 
  0.857978, 0.835424, 0.814112, 0.794314, 0.776204, 
  0.759884, 0.745394, 0.732719, 0.721809, 0.712164, 
  0.704299, 0.697804, 0.692491, 0.688137, 0.68448, 
  0.681232, 0.676128, 0.674154, 0.670553, 0.666034};



/////////////////////////////////////////////////////////////
//
// anu_mu CC xsc_tot/E_nu, in 10^-38 cm2/GeV

const G4double G4MuNeutrinoNucleusTotXsc::fANuMuTotXsc[50] = 
{
0.0291812, 0.0979725, 0.136884, 0.16794, 0.194698, 
0.218468, 0.23992, 0.259241, 0.27665, 0.292251, 
0.30612, 0.318314, 0.328886, 0.337885, 0.345464, 
0.351495, 0.356131, 0.359448, 0.361531, 0.362474, 
0.362382, 0.361365, 0.359538, 0.357024, 0.353943, 
0.350422, 0.346685, 0.342662, 0.338567, 0.334514, 
0.330612, 0.326966, 0.323668, 0.320805, 0.318451, 
0.316671, 0.315514, 0.315013, 0.315187, 0.316036, 
0.317541, 0.319667, 0.322362, 0.325556, 0.329159, 
0.332577, 0.337133, 0.341214, 0.345128, 0.347657};
