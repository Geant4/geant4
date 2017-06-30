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


#include "G4NeutrinoElectronCcXsc.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4HadTmpUtil.hh"
#include "G4NistManager.hh"

#include "G4MuonMinus.hh"
#include "G4TauMinus.hh"

using namespace std;
using namespace CLHEP;

G4NeutrinoElectronCcXsc::G4NeutrinoElectronCcXsc()
 : G4VCrossSectionDataSet("NuElectronCcXsc")
{
  // PDG2016: Gf=1.1663787(6)e-5*(hc)^3/GeV^2
  // fCofXsc  = Gf*Gf*MeC2*2/pi

  fCofXsc  = 1.36044e-22;
  fCofXsc *= hbarc*hbarc*electron_mass_c2;
  fCofXsc /= halfpi;

  // G4cout<<"fCofXsc = "<<fCofXsc*GeV/cm2<<" cm2/GeV"<<G4endl;

  // G4cout<<"hbarc = "<<hbarc/MeV/fermi<<" MeV*fermi"<<G4endl;

  // PDG2016: sin^2 theta Weinberg

  fSin2tW = 0.23129; // 0.2312;

  fCutEnergy = 0.; // default value

  fBiasingFactor = 1.; // default as physics

  theMuonMinus = G4MuonMinus::MuonMinus(); 
  theTauMinus  = G4TauMinus::TauMinus(); 
}

G4NeutrinoElectronCcXsc::~G4NeutrinoElectronCcXsc() 
{}

//////////////////////////////////////////////////////

G4bool 
G4NeutrinoElectronCcXsc::IsElementApplicable( const G4DynamicParticle* aPart, G4int, const G4Material*)
{
  G4bool result  = false;
  G4String pName = aPart->GetDefinition()->GetParticleName();
  G4double minEnergy = 0., energy = aPart->GetTotalEnergy();
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
  return result;
}

////////////////////////////////////////////////////

G4double G4NeutrinoElectronCcXsc::
GetElementCrossSection(const G4DynamicParticle* aPart, G4int ZZ,  
		       const G4Material*) 
{
  G4double result = 0., totS, fmass, fmass2, emass=electron_mass_c2, emass2;

  G4double energy = aPart->GetTotalEnergy();
  G4String pName   = aPart->GetDefinition()->GetParticleName();

  emass2 = emass*emass;
  totS   = 2.*energy*emass + emass2;

  if( pName == "nu_mu")
  {
    fmass  = theMuonMinus->GetPDGMass();
    fmass2 = fmass*fmass;
    result = (1. - fmass2/totS)*(1. - fmass2/totS);
  }
  else if( pName == "anti_nu_mu")
  {
    fmass  = theMuonMinus->GetPDGMass();
    fmass2 = fmass*fmass;

    result  = (1.+ emass2/totS)*(1.+ fmass2/totS);
    result += (1.- emass2/totS)*(1.- fmass2/totS)/3.;
    result *= 0.25*(1. - fmass2/totS)*(1. - fmass2/totS);
  }
  else if( pName == "nu_tau") 
  {
    fmass  = theTauMinus->GetPDGMass();
    fmass2 = fmass*fmass;
    result = (1. - fmass2/totS)*(1. - fmass2/totS);
  }
  else if( pName == "anti_nu_tau")
  {
    fmass  = theTauMinus->GetPDGMass();
    fmass2 = fmass*fmass;

    result  = (1.+ emass2/totS)*(1.+ fmass2/totS);
    result += (1.- emass2/totS)*(1.- fmass2/totS)/3.;
    result *= 0.25*(1. - fmass2/totS)*(1. - fmass2/totS);
  }
  else
  {
    return result;
  }
  // if( energy <= electron_mass_c2 ) return result;

  result *= fCofXsc; //*energy;
  result *= energy + 0.5*emass;
  result *= ZZ;  // incoherent sum over  all element electrons

  result *= fBiasingFactor; // biasing up, if set >1

  return result;
}

