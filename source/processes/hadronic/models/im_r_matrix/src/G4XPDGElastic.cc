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
// $Id: G4XPDGElastic.cc,v 1.4 2010-03-12 15:45:18 gunter Exp $ //
// -------------------------------------------------------------------
//      
// PDG  Elastic cross section 
// PDG fits, Phys.Rev. D50 (1994), 1335
//
//
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4ios.hh"
#include "G4Log.hh"
#include "G4Pow.hh"
#include "G4SystemOfUnits.hh"
#include "G4XPDGElastic.hh"
#include "G4KineticTrack.hh"
#include "G4ParticleDefinition.hh"
#include "G4DataVector.hh"

#include "G4AntiProton.hh"
#include "G4AntiNeutron.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4KaonMinus.hh"
#include "G4KaonPlus.hh"

const G4double G4XPDGElastic::_lowLimit = 5. * GeV;
const G4double G4XPDGElastic::_highLimit = DBL_MAX;

// Parameters of the PDG Elastic cross-section fit (Rev. Particle Properties, 1998)
// Columns are: lower and higher fit range, X, Y1, Y2
const G4int G4XPDGElastic::nPar = 7;
// p pi+
const G4double G4XPDGElastic::pPiPlusPDGFit[7] =  { 2.,     200.,   0.,   11.4, -0.4,  0.079,  0. };
// p pi-
const G4double G4XPDGElastic::pPiMinusPDGFit[7] = { 2.,     360.,   1.76, 11.2, -0.64, 0.043,  0. };
// p K+
const G4double G4XPDGElastic::pKPlusPDGFit[7] =   { 2.,     175.,    5.0,  8.1, -1.8,  0.16,  -1.3 }; 
// p K-
const G4double G4XPDGElastic::pKMinusPDGFit[7] =  { 2.,     175.,    7.3,  0.,   0.,   0.29,  -2.40 };
// p p
const G4double G4XPDGElastic::ppPDGFit[7] =       { 2.,    2100.,   11.9, 26.9, -1.21, 0.169, -1.85 };
// p pbar
const G4double G4XPDGElastic::ppbarPDGFit[7] =    { 5., 1730000.,   10.2, 52.7, -1.16, 0.125, -1.28 };
// n pbar
const G4double G4XPDGElastic::npbarPDGFit[7] =    { 1.1,      5.55, 36.5,  0.,   0.,   0.,    -11.9 };


G4XPDGElastic::G4XPDGElastic() 
{
  const G4ParticleDefinition * proton = G4Proton::ProtonDefinition();
  const G4ParticleDefinition * neutron = G4Neutron::NeutronDefinition();
  const G4ParticleDefinition * piPlus = G4PionPlus::PionPlusDefinition();
  const G4ParticleDefinition * piMinus = G4PionMinus::PionMinusDefinition();
  const G4ParticleDefinition * KPlus = G4KaonPlus::KaonPlusDefinition();
  const G4ParticleDefinition * KMinus = G4KaonMinus::KaonMinusDefinition();
  const G4ParticleDefinition * antiproton = G4AntiProton::AntiProtonDefinition();
  
  std::pair<const G4ParticleDefinition *,const G4ParticleDefinition *> pp(proton,proton);
  std::pair<const G4ParticleDefinition *,const G4ParticleDefinition *> pn(proton,neutron);
  std::pair<const G4ParticleDefinition *,const G4ParticleDefinition *> piPlusp(piPlus,proton);
  std::pair<const G4ParticleDefinition *,const G4ParticleDefinition *> piMinusp(piMinus,proton);
  std::pair<const G4ParticleDefinition *,const G4ParticleDefinition *> KPlusp(KPlus,proton);
  std::pair<const G4ParticleDefinition *,const G4ParticleDefinition *> KMinusp(KMinus,proton);
  std::pair<const G4ParticleDefinition *,const G4ParticleDefinition *> nn(neutron,neutron);
  std::pair<const G4ParticleDefinition *,const G4ParticleDefinition *> ppbar(proton,antiproton);
  std::pair<const G4ParticleDefinition *,const G4ParticleDefinition *> npbar(antiproton,neutron);

  std::vector<G4double> ppData;
  std::vector<G4double> pPiPlusData;
  std::vector<G4double> pPiMinusData;
  std::vector<G4double> pKPlusData;
  std::vector<G4double> pKMinusData;
  std::vector<G4double> ppbarData;
  std::vector<G4double> npbarData;

  G4int i;
  for (i=0; i<2; i++) 
    {      
      ppData.push_back(ppPDGFit[i] * GeV); 
      pPiPlusData.push_back(pPiPlusPDGFit[i] * GeV); 
      pPiMinusData.push_back(pPiMinusPDGFit[i] * GeV); 
      pKPlusData.push_back(pKPlusPDGFit[i] * GeV); 
      pKMinusData.push_back(pKMinusPDGFit[i] * GeV); 
      ppbarData.push_back(ppbarPDGFit[i] * GeV);
      npbarData.push_back(npbarPDGFit[i] * GeV);
    }

  for (i=2; i<nPar; i++) 
    {      
      ppData.push_back(ppPDGFit[i]); 
      pPiPlusData.push_back(pPiPlusPDGFit[i]); 
      pPiMinusData.push_back(pPiMinusPDGFit[i]); 
      pKPlusData.push_back(pKPlusPDGFit[i]); 
      pKMinusData.push_back(pKMinusPDGFit[i]); 
      ppbarData.push_back(ppbarPDGFit[i]);
      npbarData.push_back(npbarPDGFit[i]);
    }

  xMap[nn] = ppData;
  xMap[pp] = ppData;
  xMap[pn] = ppData;
  xMap[piPlusp] = pPiPlusData;
  xMap[piMinusp] = pPiMinusData;
  xMap[KPlusp] = pKPlusData;
  xMap[KMinusp] = pKMinusData;
  xMap[ppbar] = ppbarData;
  xMap[npbar] = npbarData;
}


G4XPDGElastic::~G4XPDGElastic()
{ }


G4bool G4XPDGElastic::operator==(const G4XPDGElastic &right) const
{
  return (this == (G4XPDGElastic *) &right);
}


G4bool G4XPDGElastic::operator!=(const G4XPDGElastic &right) const
{
  return (this != (G4XPDGElastic *) &right);
}


G4double G4XPDGElastic::CrossSection(const G4KineticTrack& trk1, const G4KineticTrack& trk2) const
{
  // Elastic Cross-section fit, 1994 Review of Particle Properties, (1994), 1

  G4double sigma = 0.;

  const G4ParticleDefinition* def1 = trk1.GetDefinition();
  const G4ParticleDefinition* def2 = trk2.GetDefinition();
  
  G4double sqrtS = (trk1.Get4Momentum() + trk2.Get4Momentum()).mag();
  G4double m_1 = def1->GetPDGMass();
  G4double m_2 = def2->GetPDGMass();
  G4double m_max = std::max(m_1,m_2);
  //  if (m1 > m) m = m1;
  
  G4double pLab = 0.;

  if (m_max > 0. && sqrtS > (m_1 + m_2)) 
    {
      pLab = std::sqrt( (sqrtS*sqrtS - (m_1+m_2)*(m_1+m_2) ) * (sqrtS*sqrtS - (m_1-m_2)*(m_1-m_2)) ) / (2*m_max);
      
      // The PDG fit formula requires p in GeV/c
      
      // Order the pair: first is the lower mass particle, second is the higher mass one
      std::pair<const G4ParticleDefinition *,const G4ParticleDefinition *> trkPair(def1,def2);
      if (def1->GetPDGMass() > def2->GetPDGMass())
	trkPair = std::pair<const G4ParticleDefinition *,const G4ParticleDefinition *>(def2,def1);
      
      std::vector<G4double> data; 
      G4double pMinFit = 0.;
      G4double pMaxFit = 0.;
      G4double aFit = 0.;
      G4double bFit = 0.;
      G4double cFit = 0.;
      G4double dFit = 0.;
      G4double nFit = 0.;

      // Debug
//      G4cout << "Map has " << xMap.size() << " elements" << G4endl;
      // Debug end
 
      if (xMap.find(trkPair) != xMap.end())
	{
	  PairDoubleMap::const_iterator iter;
	  for (iter = xMap.begin(); iter != xMap.end(); ++iter)
	    {
	      std::pair<const G4ParticleDefinition *,const G4ParticleDefinition *> thePair = (*iter).first;
	      if (thePair == trkPair)
		{
		  data = (*iter).second;
		  pMinFit = data[0];
		  pMaxFit = data[1];
		  aFit = data[2];
		  bFit = data[3];
		  cFit = data[5];
		  dFit = data[6];
		  nFit = data[4];
	      
		  if (pLab < pMinFit) return 0.0;
		  if (pLab > pMaxFit )
		    G4cout << "WARNING! G4XPDGElastic::PDGElastic " 
			   << trk1.GetDefinition()->GetParticleName() << "-" 
			   << trk2.GetDefinition()->GetParticleName() 
			   << " elastic cross section: momentum " 
			   << pLab / GeV << " GeV outside valid fit range " 
			   << pMinFit /GeV << "-" << pMaxFit / GeV
			   << G4endl;
		  
		  pLab /= GeV;
		  if (pLab > 0.) 
		    {
		      G4double logP = G4Log(pLab);
		      sigma = aFit + bFit * G4Pow::GetInstance()->powA(pLab, nFit) + cFit * logP*logP + dFit * logP;
		      sigma = sigma * millibarn;
		    }
		}
	    }
	}
      else
	{
	  G4cout << "G4XPDGElastic::CrossSection - Data not found in Map" << G4endl;
	}
    }
  
  if (sigma < 0.)
    {
      G4cout << "WARNING! G4XPDGElastic::PDGElastic "      
	     << def1->GetParticleName() << "-" << def2->GetParticleName()
	     << " elastic cross section: momentum " 
	     << pLab << " GeV, negative cross section " 
	     << sigma / millibarn << " mb set to 0" << G4endl;
      sigma = 0.;
    }
  
  return sigma;
}


G4String G4XPDGElastic::Name() const
{
  G4String name = "PDGElastic ";
  return name;
}


G4bool G4XPDGElastic::IsValid(G4double e) const
{
  G4bool answer = InLimits(e,_lowLimit,_highLimit);

  return answer;
}











