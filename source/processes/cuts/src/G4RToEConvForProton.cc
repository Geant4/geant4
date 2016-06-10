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
// $Id: G4RToEConvForProton.cc 70745 2013-06-05 10:54:00Z gcosmo $
//
//
// --------------------------------------------------------------
//      GEANT 4 class implementation file/  History:
//    5 Oct. 2002, H.Kuirashige : Structure created based on object model
// --------------------------------------------------------------

#include "G4RToEConvForProton.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4PhysicsLogVector.hh"

#include "G4ios.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

G4RToEConvForProton::G4RToEConvForProton() 
  : G4VRangeToEnergyConverter(),
    Mass(0.0),
    Z(-1.),  
    tau0(0.0), taul(0.0), taum(0.0),
    ionpot(0.0),
    ca(0.0), cba(0.0), cc(0.0)
{    
  theParticle =  G4ParticleTable::GetParticleTable()->FindParticle("proton");
  if (theParticle ==0) {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << " G4RToEConvForProton::G4RToEConvForProton() ";
      G4cout << " proton is not defined !!" << G4endl;
    }
#endif
  } else {
    Mass = theParticle->GetPDGMass();
  } 
}

G4RToEConvForProton::~G4RToEConvForProton()
{ 
}


G4double G4RToEConvForProton::Convert(G4double rangeCut, const G4Material* )
{
  // Simple formula
  //   range = Ekin/(100*keV)*(1*mm);
  return (rangeCut/(1.0*mm)) * (100.0*keV); 
}


// **********************************************************************
// ************************* ComputeLoss ********************************
// **********************************************************************
G4double G4RToEConvForProton::ComputeLoss(G4double AtomicNumber,
                                                G4double KineticEnergy) 
{
  //  calculate dE/dx
  const G4double  z2Particle = 1.0;

  if( std::fabs(AtomicNumber-Z)>0.1 ){
    // recalculate constants
    Z = AtomicNumber;
    G4double Z13 = std::exp(std::log(Z)/3.);
    tau0 = 0.1*Z13*MeV/proton_mass_c2;
    taum = 0.035*Z13*MeV/proton_mass_c2;
    taul = 2.*MeV/proton_mass_c2;
    ionpot = 1.6e-5*MeV*std::exp(0.9*std::log(Z));
    cc = (taul+1.)*(taul+1.)*std::log(2.*electron_mass_c2*taul*(taul+2.)/ionpot)/(taul*(taul+2.))-1.;
    cc = 2.*twopi_mc2_rcl2*Z*cc*std::sqrt(taul);
    ca = cc/((1.-0.5*std::sqrt(tau0/taum))*tau0);
    cba = -0.5/std::sqrt(taum);
  }

  G4double tau = KineticEnergy/Mass;
  G4double dEdx;
  if ( tau <= tau0 ) {
    dEdx = ca*(std::sqrt(tau)+cba*tau);
  } else {
    if( tau <= taul ) {
      dEdx = cc/std::sqrt(tau);
    } else {
      dEdx = (tau+1.)*(tau+1.)*
             std::log(2.*electron_mass_c2*tau*(tau+2.)/ionpot)/(tau*(tau+2.))-1.;
      dEdx = 2.*twopi_mc2_rcl2*Z*dEdx;
    }
  }
  return dEdx*z2Particle ;
}


// **********************************************************************
// ************************* Reset       ********************************
// **********************************************************************
void G4RToEConvForProton::Reset()
{
  // do nothing because loss tables and range vectors are not used 
  return;
}
