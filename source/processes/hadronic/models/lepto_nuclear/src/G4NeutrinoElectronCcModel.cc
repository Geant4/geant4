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
// $Id: G4NeutrinoElectronCcModel.cc 91806 2015-08-06 12:20:45Z gcosmo $
//
// Geant4 Header : G4NeutrinoElectronCcModel
//
// Author : V.Grichine 26.4.17
//  

#include "G4NeutrinoElectronCcModel.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "Randomize.hh"
#include "G4NeutrinoE.hh"
#include "G4AntiNeutrinoE.hh"
#include "G4MuonMinus.hh"
#include "G4TauMinus.hh"

using namespace std;
using namespace CLHEP;

G4NeutrinoElectronCcModel::G4NeutrinoElectronCcModel(const G4String& name) 
  : G4HadronicInteraction(name)
{
  SetMinEnergy( 0.0*GeV );
  SetMaxEnergy( 100.*TeV );
  SetMinEnergy(1.e-6*eV);  

  theNeutrinoE = G4NeutrinoE::NeutrinoE();
  theAntiNeutrinoE = G4AntiNeutrinoE::AntiNeutrinoE();
  theMuonMinus = G4MuonMinus::MuonMinus();
  theTauMinus  = G4TauMinus::TauMinus();

  // PDG2016: sin^2 theta Weinberg

  fSin2tW = 0.23129; // 0.2312;

  fCutEnergy = 0.; // default value

}


G4NeutrinoElectronCcModel::~G4NeutrinoElectronCcModel()
{}


void G4NeutrinoElectronCcModel::ModelDescription(std::ostream& outFile) const
{

    outFile << "G4NeutrinoElectronCcModel is a neutrino-electron (neutral current) elastic scattering\n"
            << "model which uses the standard model \n"
            << "transfer parameterization.  The model is fully relativistic\n";

}

/////////////////////////////////////////////////////////

G4bool G4NeutrinoElectronCcModel::IsApplicable(const G4HadProjectile & aPart, 
					       G4Nucleus & targetNucleus)
{
  G4bool result  = false;
  G4String pName = aPart.GetDefinition()->GetParticleName();
  G4double minEnergy = 0., energy = aPart.GetTotalEnergy();
  G4double fmass, emass = electron_mass_c2;

  if(      pName == "nu_mu"   || pName == "anti_nu_mu"  ) fmass = theMuonMinus->GetPDGMass(); 
  else if( pName == "nu_tau"  || pName == "anti_nu_tau" ) fmass = theTauMinus->GetPDGMass(); 
  else fmass = emass;

  minEnergy = (fmass-emass)*(fmass+emass)/emass;
  
  if( ( pName == "nu_mu"  || pName == "anti_nu_mu"  || 
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

G4HadFinalState* G4NeutrinoElectronCcModel::ApplyYourself(
		 const G4HadProjectile& aTrack, G4Nucleus& targetNucleus)
{
  theParticleChange.Clear();

  const G4HadProjectile* aParticle = &aTrack;
  G4double energy = aParticle->GetTotalEnergy();

  if( energy <= GetMinEnergy() ) 
  {
    theParticleChange.SetEnergyChange(energy);
    theParticleChange.SetMomentumChange(aTrack.Get4Momentum().vect().unit());
    return &theParticleChange;
  }

  G4double emass = electron_mass_c2;
  G4double sTot = 2.*energy*emass + emass*emass;
 
  G4LorentzVector lvp1 = aParticle->Get4Momentum();
  G4LorentzVector lvt1(0.,0.,0.,electron_mass_c2);
  G4LorentzVector lvsum = lvp1+lvt1;
  G4ThreeVector bst = lvsum.boostVector();

  // sample and make final state in CMS frame

  G4double cost = SampleCosCMS( aParticle );
  G4double sint = std::sqrt( (1.0 - cost)*(1.0 + cost) );
  G4double phi  = G4UniformRand()*CLHEP::twopi;

  G4ThreeVector eP( sint*std::cos(phi), sint*std::sin(phi), cost );

  G4String pName  = aParticle->GetDefinition()->GetParticleName();

  G4double massf = 0.0;
  if(      pName == "nu_mu" || pName == "anti_nu_mu"  ) massf = theMuonMinus->GetPDGMass();
  else if( pName == "nu_tau" || pName == "anti_nu_tau") massf = theTauMinus->GetPDGMass();

  G4double massf2 = massf*massf;

  G4double epf = 0.5*(sTot - massf2)/sqrt(sTot);
  // G4double etf = epf*(sTot + massf2)/(sTot - massf2);

  eP *= epf;
  G4LorentzVector lvp2( eP, epf );
  lvp2.boost(bst); // back to lab frame

  G4LorentzVector lvt2 = lvsum - lvp2; // ?

  G4DynamicParticle* aNu = nullptr; 
  G4DynamicParticle* aLept = nullptr; 

  if(  pName == "nu_mu" || pName == "nu_tau")               
  {
    aNu = new G4DynamicParticle( theNeutrinoE, lvp2 );
  }
  else if( pName == "anti_nu_mu" || pName == "anti_nu_tau") 
  {
    aNu = new G4DynamicParticle( theAntiNeutrinoE, lvp2 );
  }
  if(  pName == "nu_mu" || pName == "anti_nu_mu")       
  {
    aLept = new G4DynamicParticle( theMuonMinus, lvt2 );
  }
  else if( pName == "nu_tau" || pName == "anti_nu_tau") 
  {
    aLept = new G4DynamicParticle( theTauMinus, lvt2 );
  }

  if ( aNu ) theParticleChange.AddSecondary( aNu );
  if ( aLept ) theParticleChange.AddSecondary( aLept );

  G4int Z = targetNucleus.GetZ_asInt();
        Z *= 1;
 
  return &theParticleChange;
}

//////////////////////////////////////////////////////
//
// sample recoil electron energy in lab frame

G4double G4NeutrinoElectronCcModel::SampleCosCMS(const G4HadProjectile* aParticle)
{
  G4double result = 0., cofL, cofR, cofLR, massf2, sTot, emass = electron_mass_c2, emass2;

  G4double energy = aParticle->GetTotalEnergy();
  
  if( energy == 0.) return result; // vmg: < th?? as in xsc 

  G4String pName  = aParticle->GetDefinition()->GetParticleName();

  if( pName == "nu_mu" || pName == "nu_tau")
  {
    return 2.*G4UniformRand()-1.; // uniform scattering cos in CMS
  }
  else if( pName == "anti_nu_mu" || pName == "anti_nu_tau")
  {
    emass2 = emass*emass;
    sTot = 2.*energy*emass + emass2;

    cofL = (sTot-emass2)/(sTot+emass2);

    if(pName == "anti_nu_mu") massf2 = theMuonMinus->GetPDGMass()*theMuonMinus->GetPDGMass();
    else                      massf2 = theTauMinus->GetPDGMass()*theTauMinus->GetPDGMass();

    cofR = (sTot-massf2)/(sTot+massf2);

    cofLR = cofL*cofR/3.;

    // cofs of cos 3rd equation

    G4double a = cofLR;
    G4double b = 0.5*(cofR+cofL);
    G4double c = 1.;

    G4double d  = -G4UniformRand()*2.*(1.+ cofLR);
             d += c - b + a;

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
	     if(D < 0.) D = -D;
	   // G4complex A1 = G4complex(- q/2., std::sqrt(-D) );
	   // G4complex A  = std::pow(A1,1./3.);

	   // G4complex B1 = G4complex(- q/2., -std::sqrt(-D) );
	   // G4complex B  = std::pow(B1,1./3.);

	     G4double A, B;

    G4double A1 = - q/2. + std::sqrt(D);
    if (A1 < 0.) A1 = -A1;
    A = std::pow(A1,1./3.);
    if (A1 < 0.) A = -A;

    G4double B1 = - q/2. - std::sqrt(D);
    // G4double B = std::pow(-B1,1./3.);
    if(B1 < 0.) B1 = -B1;
    B = std::pow(B1,1./3.);
    if(B1 < 0.)   B = -B;
    // G4cout<<"A1 = "<<A1<<"; A = "<<A<<"; B1 = "<<B1<<"; B = "<<B<<G4endl;
    // roots of the incomplete 3rd equation

    G4complex y1 =  A + B;
    // G4complex y2 = -0.5*(A + B) + 0.5*std::sqrt(3.)*(A - B)*G4complex(0.,1.);
    // G4complex y3 = -0.5*(A + B) - 0.5*std::sqrt(3.)*(A - B)*G4complex(0.,1.);
 
    G4complex x1 = y1 - b/a/3.;
    // G4complex x2 = y2 - b/a/3.;
    // G4complex x3 = y3 - b/a/3.;
    // G4cout<<"re_x1 = "<<real(x1)<<" + i*"<<imag(x1)<<G4endl;
    // G4cout<<"re_x1 = "<<real(x1)<<"; re_x2 = "<<real(x2)<<"; re_x3 = "<<real(x3)<<G4endl;
    // G4cout<<"im_x1 = "<<imag(x1)<<"; im_x2 = "<<imag(x2)<<"; im_x3 = "<<imag(x3)<<G4endl<<G4endl;

    result = real(x1);
  }
  else 
  {
    return result;
  }
  return result;
}

//
//
///////////////////////////
