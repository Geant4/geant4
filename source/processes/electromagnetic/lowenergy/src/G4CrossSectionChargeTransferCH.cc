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
// $Id: G4CrossSectionChargeTransferCH.cc,v 1.4 2009-06-10 13:32:36 mantero Exp $
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
// http://www-pub.iaea.org/MTCD/publications/PDF/APID-VOL10.pdf 
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
  lowEnergyLimit = 0.1 * eV;
  highEnergyLimit = 1.0 * MeV;

  std::vector<G4double> coeff;
  G4String materialName;

  // CH
  materialName = "CH";
  coeff.push_back(-1.060468);
  coeff.push_back(-5.662572);
  coeff.push_back(-4.376450);
  coeff.push_back(-3.567226);
  coeff.push_back(-1.433069);
  coeff.push_back(-0.5789399);
  coeff.push_back(-0.3523295);  
  coeff.push_back(-0.09956988);   
  coeff.push_back(0.01532751);
  crossMap[materialName] = coeff;
  coeff.clear();

  // CH2
  materialName = "CH2";
  coeff.push_back(1.118652E-01);  
  coeff.push_back(-5.216803E+00);  
  coeff.push_back(-3.893841E+00);  
  coeff.push_back(-2.756865E+00);  
  coeff.push_back(-1.192043E+00);  
  coeff.push_back(-5.200059E-01);  
  coeff.push_back(-1.816781E-01);   
  coeff.push_back(1.129866E-02);  
  coeff.push_back(-3.658702E-03);
  crossMap[materialName] = coeff;
  coeff.clear();

  // CH3
  coeff.push_back(6.854778E-01);  
  coeff.push_back(-5.810741E+00);  
  coeff.push_back(-3.622136E+00);  
  coeff.push_back(-2.179920E+00);  
  coeff.push_back(-7.801006E-01);  
  coeff.push_back(-2.421371E-02);   
  coeff.push_back(2.762333E-01);  
  coeff.push_back(2.288237E-01);   
  coeff.push_back(6.164797E-02);
  crossMap[materialName] = coeff;
  coeff.clear();

  // CH4
  materialName = "CH4";
  coeff.push_back(0.3901564);
  coeff.push_back(-6.426675);
  coeff.push_back(-3.706893);
  coeff.push_back(-1.999034);
  coeff.push_back(-0.5625439);
  coeff.push_back(0.2279431);
  coeff.push_back(0.3443980);  
  coeff.push_back(0.1566892);   
  coeff.push_back(-0.05398410);
  coeff.push_back(-0.1822252);  
  coeff.push_back(-0.1593352);   
  coeff.push_back(-0.08826322);
  crossMap[materialName] = coeff;
  coeff.clear();

  // C2H
  materialName = "C2H";
  coeff.push_back(-1.986507E+00);
  coeff.push_back(-5.720283E+00);
  coeff.push_back(-3.535139E+00);
  coeff.push_back(-4.230273E+00);
  coeff.push_back(-1.254502E+00);
  coeff.push_back(-2.056898E-01);
  coeff.push_back(-4.595756E-01);
  coeff.push_back(-7.842824E-02);
  coeff.push_back(3.002537E-02);
  coeff.push_back(-5.626713E-02);
  coeff.push_back(6.455583E-02);
  crossMap[materialName] = coeff;
  coeff.clear();

  // C2H2
  materialName = "C2H2";
  coeff.push_back(2.513870E-01);
  coeff.push_back(-5.812705E+00);
  coeff.push_back(-3.338185E+00);
  coeff.push_back(-3.071630E+00);
  coeff.push_back(-1.433263E+00);
  coeff.push_back(-3.583544E-01);
  coeff.push_back(-1.456216E-01);
  coeff.push_back(-7.391778E-03);
  coeff.push_back(1.151712E-02);
  crossMap[materialName] = coeff;
  coeff.clear();
  // 

}


G4CrossSectionChargeTransferCH::~G4CrossSectionChargeTransferCH()
{ }
 

G4double G4CrossSectionChargeTransferCH::CrossSection(const G4Track& track)
{
  G4double sigma = DBL_MIN;
 
  const G4DynamicParticle* particle = track.GetDynamicParticle();
  G4double e = particle->GetKineticEnergy() / keV;

  // Get the material the track is in
  G4Material* material = track.GetMaterial();
  G4String materialName;
  materialName = material->GetName();

  // Check whether cross section data are available for the current material

  std::map<G4String,std::vector<G4double>,std::less<G4String> >::const_iterator pos;
  pos = crossMap.find(materialName);
  if (pos!= crossMap.end())
    {
      std::vector<G4double> coeff = (*pos).second;
      G4double eMin = lowEnergyLimit / keV;
      G4double eMax = highEnergyLimit / keV;
      G4double scaleFactor = 1.e-16 * cm2;

      if (e >=  eMin && e <= eMax)
	{
	  G4double x = (2.*std::log(e) - std::log(eMin) - std::log(eMax)) / (std::log(eMax) - std::log(eMin));
	  
	  std::vector<G4double> t;
	  t.push_back(1.);
	  t.push_back(x);
	  G4double cross = t[0] * coeff[0] + t[1] * coeff[1];

	  G4int nCoeff = coeff.size();
	  for (G4int i=2; i<nCoeff; i++)
	    {
	      G4double tNext = 2. * t[i-1] - t[i-2];
	      t.push_back(tNext);
	      cross = cross + coeff[i] * tNext;
	    }
	  sigma = std::exp(cross) * scaleFactor;
	}
    }      

  return sigma;
}

