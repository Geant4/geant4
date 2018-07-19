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
// $Id: G4NeutrinoElectronNcModel.cc 91806 2015-08-06 12:20:45Z gcosmo $
//
// Geant4 Header : G4NeutrinoElectronNcModel
//
// Author : V.Grichine 6.4.17
//  

#include "G4NeutrinoElectronNcModel.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "Randomize.hh"
#include "G4Electron.hh"

using namespace std;
using namespace CLHEP;

G4NeutrinoElectronNcModel::G4NeutrinoElectronNcModel(const G4String& name) 
  : G4HadronElastic(name)
{
  SetMinEnergy( 0.0*GeV );
  SetMaxEnergy( 100.*TeV );
  SetLowestEnergyLimit(1.e-6*eV);  

  theElectron = G4Electron::Electron();
  // PDG2016: sin^2 theta Weinberg

  fSin2tW = 0.23129; // 0.2312;

  fCutEnergy = 0.; // default value

}


G4NeutrinoElectronNcModel::~G4NeutrinoElectronNcModel()
{}


void G4NeutrinoElectronNcModel::ModelDescription(std::ostream& outFile) const
{

    outFile << "G4NeutrinoElectronNcModel is a neutrino-electron (neutral current) elastic scattering\n"
            << "model which uses the standard model \n"
            << "transfer parameterization.  The model is fully relativistic\n";

}

/////////////////////////////////////////////////////////

G4bool G4NeutrinoElectronNcModel::IsApplicable(const G4HadProjectile & aTrack, 
  			      G4Nucleus & targetNucleus)
{
  G4bool result  = false;
  G4String pName = aTrack.GetDefinition()->GetParticleName();
  G4double minEnergy = 0., energy = aTrack.GetTotalEnergy();

  if( fCutEnergy > 0. ) // min detected recoil electron energy
  {
    minEnergy = 0.5*(fCutEnergy+sqrt(fCutEnergy*(fCutEnergy+2.*electron_mass_c2)));
  }
  if( ( pName == "nu_e"   || pName == "anti_nu_e"   || 
        pName == "nu_mu"  || pName == "anti_nu_nu"  || 
        pName == "nu_tau" || pName == "anti_nu_tau"   ) &&
        energy > minEnergy                                 )
  {
    result = true;
  }
  G4int Z = targetNucleus.GetZ_asInt();
        Z *= 1;

  return result;
}

////////////////////////////////////////////////
//
//

G4HadFinalState* G4NeutrinoElectronNcModel::ApplyYourself(
		 const G4HadProjectile& aTrack, G4Nucleus& targetNucleus)
{
  theParticleChange.Clear();

  const G4HadProjectile* aParticle = &aTrack;
  G4double nuTkin = aParticle->GetKineticEnergy();

  if( nuTkin <= LowestEnergyLimit() ) 
  {
    theParticleChange.SetEnergyChange(nuTkin);
    theParticleChange.SetMomentumChange(aTrack.Get4Momentum().vect().unit());
    return &theParticleChange;
  }
  // sample and make final state in lab frame

  G4double eTkin = SampleElectronTkin( aParticle );

  if( eTkin > fCutEnergy )
  {
    G4double ePlab = sqrt( eTkin*(eTkin + 2.*electron_mass_c2) );

    G4double cost2  = eTkin*(nuTkin + electron_mass_c2)*(nuTkin + electron_mass_c2);
             cost2 /= nuTkin*nuTkin*(eTkin + 2.*electron_mass_c2);

    if( cost2 > 1. ) cost2 = 1.;
    if( cost2 < 0. ) cost2 = 0.;

    G4double cost = sqrt(cost2);
    G4double sint = std::sqrt( (1.0 - cost)*(1.0 + cost) );
    G4double phi  = G4UniformRand()*CLHEP::twopi;

    G4ThreeVector eP( sint*std::cos(phi), sint*std::sin(phi), cost );
    eP *= ePlab;
    G4LorentzVector lvt2( eP, eTkin + electron_mass_c2 );
    G4DynamicParticle * aSec = new G4DynamicParticle( theElectron, lvt2 );
    theParticleChange.AddSecondary( aSec );

    G4LorentzVector lvp1 = aParticle->Get4Momentum();
    G4LorentzVector lvt1(0.,0.,0.,electron_mass_c2);
    G4LorentzVector lvsum = lvp1+lvt1;

    G4LorentzVector lvp2 = lvsum-lvt2;
    G4double nuTkin2 = lvp2.e()-aParticle->GetDefinition()->GetPDGMass();
    theParticleChange.SetEnergyChange(nuTkin2);
    theParticleChange.SetMomentumChange(lvp2.vect().unit());
  }
  else if( eTkin > 0.0 ) 
  {
    theParticleChange.SetLocalEnergyDeposit( eTkin );
    nuTkin -= eTkin;

    if( nuTkin > 0. )
    {
      theParticleChange.SetEnergyChange( nuTkin );
      theParticleChange.SetMomentumChange( aTrack.Get4Momentum().vect().unit() );
    }
  }
  else 
  {
    theParticleChange.SetEnergyChange( nuTkin );
    theParticleChange.SetMomentumChange( aTrack.Get4Momentum().vect().unit() );
  }
  G4int Z = targetNucleus.GetZ_asInt();
        Z *= 1;
 
  return &theParticleChange;
}

//////////////////////////////////////////////////////
//
// sample recoil electron energy in lab frame

G4double G4NeutrinoElectronNcModel::SampleElectronTkin(const G4HadProjectile* aParticle)
{
  G4double result = 0., xi, cofL, cofR, cofL2, cofR2, cofLR;

  G4double energy = aParticle->GetTotalEnergy();
  if( energy == 0.) return result; // vmg: < th?? as in xsc 

  G4String pName  = aParticle->GetDefinition()->GetParticleName();

  if( pName == "nu_e")
  {
    cofL = 0.5 + fSin2tW;
    cofR = fSin2tW;
  }
  else if( pName == "anti_nu_e")
  {
    cofL = fSin2tW;
    cofR = 0.5 + fSin2tW;
  }
  else if( pName == "nu_mu")
  {
    cofL = -0.5 + fSin2tW;
    cofR = fSin2tW;
  }
  else if( pName == "anti_nu_mu")
  {
    cofL = fSin2tW;
    cofR = -0.5 + fSin2tW;
  }
  else if( pName == "nu_tau") // vmg: nu_tau as nu_mu ???
  {
    cofL = -0.5 + fSin2tW;
    cofR = fSin2tW;
  }
  else if( pName == "anti_nu_tau")
  {
    cofL = fSin2tW;
    cofR = -0.5 + fSin2tW;
  }
  else
  {
    return result;
  }
  xi = 0.5*electron_mass_c2/energy;

  cofL2 = cofL*cofL;
  cofR2 = cofR*cofR;
  cofLR = cofL*cofR;

  // cofs of Tkin/Enu 3rd equation

  G4double a = cofR2/3.;
  G4double b = -(cofR2+cofLR*xi);
  G4double c = cofL2+cofR2;

  G4double xMax  = 1./(1. + xi);
  G4double xMax2 = xMax*xMax;
  G4double xMax3 = xMax*xMax2;

  G4double d  = -( a*xMax3 + b*xMax2 + c*xMax );
           d *= G4UniformRand();

  // G4cout<<a<<"   "<<b<<"   "<<c<<"   "<<d<<G4endl<<G4endl;

  // cofs of the incomplete 3rd equation

  G4double p  = c/a;
           p -= b*b/a/a/3.;
  G4double q  = d/a;
           q -= b*c/a/a/3.;
           q += 2*b*b*b/a/a/a/27.;


  // cofs for the incomplete colutions

  G4double D  = p*p*p/3./3./3.;
           D += q*q/2./2.;

	   // G4cout<<"D = "<<D<<G4endl;
	   // D = -D;
	   // G4complex A1 = G4complex(- q/2., std::sqrt(-D) );
	   // G4complex A  = std::pow(A1,1./3.);

	   // G4complex B1 = G4complex(- q/2., -std::sqrt(-D) );
	   // G4complex B  = std::pow(B1,1./3.);

  G4double A1 = - q/2. + std::sqrt(D);
  G4double A = std::pow(A1,1./3.);

  G4double B1 = - q/2. - std::sqrt(D);
  G4double B = std::pow(-B1,1./3.);
           B = -B;

  // roots of the incomplete 3rd equation

  G4complex y1 =  A + B;
  // G4complex y2 = -0.5*(A + B) + 0.5*std::sqrt(3.)*(A - B)*G4complex(0.,1.);
  // G4complex y3 = -0.5*(A + B) - 0.5*std::sqrt(3.)*(A - B)*G4complex(0.,1.);
 
  G4complex x1 = y1 - b/a/3.;
  // G4complex x2 = y2 - b/a/3.;
  // G4complex x3 = y3 - b/a/3.;

  // G4cout<<"re_x1 = "<<real(x1)<<"; re_x2 = "<<real(x2)<<"; re_x3 = "<<real(x3)<<G4endl;
  // G4cout<<"im_x1 = "<<imag(x1)<<"; im_x2 = "<<imag(x2)<<"; im_x3 = "<<imag(x3)<<G4endl<<G4endl;

  result = real(x1)*energy;

  return result;
}

//
//
///////////////////////////
