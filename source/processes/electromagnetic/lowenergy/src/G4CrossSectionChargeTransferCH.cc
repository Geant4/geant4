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
// $Id: G4CrossSectionChargeTransferCH.cc,v 1.1 2008-02-22 16:26:33 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// Contact Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// Reference: TNS Geant4-DNA paper
// Reference for implementation model: NIM. 155, pp. 145-156, 1978

// History:
// -----------
// Date         Name              Modification
// 28 Dec 2007  M.G. Pia          Created
//
// -------------------------------------------------------------------

// Class description:
// Total cross section for incident p charge transfer in hydrocarbons
// Reference: K. Janev, J. G. Wang and T. Kato, NIFS-DATA (2001)
//            http://www-cfadc.phy.ornl.gov/astro/ps/data/
//
// The fit formula is:
// ln(cross section)=\sum C_i T_i(x)
//            x=[2ln(E)-ln(E_min)-ln(E_max)]/[ln(E_max)-ln(E_min)]
//            T_1(x)=1
//            T_2(x)=x
//          T_n+2(x)=2T_n+1(x)-T_n(x)
//        
// Where cross section is given in 10^-16 cm2 and E is in keV,
// E_min=1e-4 keV and E_max=1.e3 keV.

// Further documentation available from http://www.ge.infn.it/geant4/

// -------------------------------------------------------------------


#include "G4CrossSectionChargeTransferCH.hh"
#include "G4Track.hh"
#include "G4DynamicParticle.hh"
#include "G4Proton.hh"

G4CrossSectionChargeTransferCH::G4CrossSectionChargeTransferCH()
{
  // Default energy limits (defined for protection against anomalous behaviour only)
  name = "ChargeTransferCH";
  lowEnergyLimit = 0.1 * keV;
  highEnergyLimit = 1.0 * MeV;

  coeff.push_back(-1.060468);
  coeff.push_back(-5.662572);
  coeff.push_back(-4.376450);
  coeff.push_back(-3.567226);
  coeff.push_back(-1.433069);
  coeff.push_back(-0.5789399);
  coeff.push_back(-0.3523295);  
  coeff.push_back(-0.09956988);   
  coeff.push_back(0.01532751);
  
}


G4CrossSectionChargeTransferCH::~G4CrossSectionChargeTransferCH()
{ }
 

G4double G4CrossSectionChargeTransferCH::CrossSection(const G4Track& track)
{
  G4double crossSection = DBL_MIN;
 
  const G4DynamicParticle* particle = track.GetDynamicParticle();
  G4double e = particle->GetKineticEnergy() / keV;
  G4double eMin = lowEnergyLimit / keV;
  G4double eMax = highEnergyLimit / keV;

  if (e >=  eMin && e <= eMax)
    {

      G4double x = (2.*std::log(e)-std::log(eMin)-std::log(eMax)) / (std::log(eMax)-std::log(eMin));
      
      std::vector<G4double> t;
      t.push_back(1.);
      t.push_back(x);
      G4double cross = coeff[0] + coeff[1] * x;
      
      for (G4int i=0; i<7; i++)
	{
	  G4double tNext = t[i] + 2. * t[i+1];
	  t.push_back(tNext);
	  cross = coeff[i+2] * tNext;
	}
      crossSection = std::exp(cross);

    }      
      
 

  return crossSection;
}

