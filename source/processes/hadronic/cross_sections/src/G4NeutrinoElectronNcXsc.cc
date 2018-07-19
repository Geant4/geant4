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


#include "G4NeutrinoElectronNcXsc.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4HadTmpUtil.hh"
#include "G4Proton.hh"
#include "G4NistManager.hh"

using namespace std;
using namespace CLHEP;

G4NeutrinoElectronNcXsc::G4NeutrinoElectronNcXsc()
 : G4VCrossSectionDataSet("NuElectronNcXsc")
{
  // PDG2016: Gf=1.1663787(6)e-5*(hc)^3/GeV^2
  // fCofXsc  = Gf*Gf*MeC2*2/pi

  fCofXsc  = 1.36044e-22;
  fCofXsc *= hbarc*hbarc*electron_mass_c2;
  fCofXsc /= halfpi;

  // G4cout<<"hbarc = "<<hbarc/MeV/fermi<<" MeV*fermi"<<G4endl;

  // PDG2016: sin^2 theta Weinberg

  fSin2tW = 0.23129; // 0.2312;

  fCutEnergy = 0.; // default value
  fBiasingFactor = 1.; 
}

G4NeutrinoElectronNcXsc::~G4NeutrinoElectronNcXsc() 
{}

//////////////////////////////////////////////////////

G4bool 
G4NeutrinoElectronNcXsc::IsElementApplicable( const G4DynamicParticle* aPart, G4int, const G4Material*)
{
  G4bool result  = false;
  G4String pName = aPart->GetDefinition()->GetParticleName();
  G4double minEnergy = 0., energy = aPart->GetTotalEnergy();
  // Z *= 1;
  if( fCutEnergy > 0. ) // min detected recoil electron energy
  {
    minEnergy = 0.5*(fCutEnergy+sqrt(fCutEnergy*(fCutEnergy+2.*electron_mass_c2)));
  }
  if( ( pName == "nu_e"   || pName == "anti_nu_e"   || 
        pName == "nu_mu"  || pName == "anti_nu_mu"  || 
        pName == "nu_tau" || pName == "anti_nu_tau"   ) &&
        energy > minEnergy                                 )
  {
    result = true;
  }
  return result;
}

////////////////////////////////////////////////////

G4double G4NeutrinoElectronNcXsc::
GetElementCrossSection(const G4DynamicParticle* aPart, G4int ZZ,  
		       const G4Material*) 
{
  G4double result = 0., cofL, cofR, cofL2, cofR2, cofLR;

  G4double energy = aPart->GetTotalEnergy();
  G4String pName   = aPart->GetDefinition()->GetParticleName();

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
  // if( energy <= electron_mass_c2 ) return result;

  cofL2 = cofL*cofL;
  cofR2 = cofR*cofR;
  cofLR = cofL*cofR;

  if( fCutEnergy > 0. )
  {
    G4double tM  = 2.*energy*energy/(electron_mass_c2 + 2.*energy);
    G4double tM2 = tM*tM;
    G4double tM3 = tM*tM2;
    G4double tC  = fCutEnergy;
    G4double tC2 = tC*tC;
    G4double tC3 = tC*tC2;

    result  = (cofL2+cofR2)*(tM-tC);
    result -= (cofR2+cofLR*0.5*electron_mass_c2/energy)*(tM2-tC2)/energy;
    result += cofR2*(tM3-tC3)/energy/energy/3.;
  }
  else
  {
    G4double rtM  = 2.*energy/(electron_mass_c2 + 2.*energy);
    G4double rtM2 = rtM*rtM;
    G4double rtM3 = rtM*rtM2;

    result  = (cofL2+cofR2)*rtM*energy;
    result -= (cofR2*energy+cofLR*0.5*electron_mass_c2)*rtM2;
    result += cofR2*rtM3*energy/3.;
  }
  // result = cofL*cofL + cofR*cofR/3.;
  // G4cout<<"cofL2 + cofR2/3. = "<<result<<G4endl;
  // result -= 0.5*cofL*cofR*electron_mass_c2/energy;

  result *= fCofXsc; //*energy;

  result *= ZZ;  // incoherent sum over  all element electrons

  result *= fBiasingFactor;

  return result;
}

