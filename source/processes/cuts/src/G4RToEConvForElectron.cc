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
// G4RToEConvForElectron class implementation
//
// Author: H.Kurashige, 05 October 2002 - First implementation
// --------------------------------------------------------------------

#include "G4RToEConvForElectron.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Pow.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

// --------------------------------------------------------------------
G4RToEConvForElectron::G4RToEConvForElectron() 
  : G4VRangeToEnergyConverter()
{    
  theParticle = G4ParticleTable::GetParticleTable()->FindParticle("e-");
  if (theParticle == nullptr )
  {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0)
    {
      G4cout << "G4RToEConvForElectron::G4RToEConvForElectron() - ";
      G4cout << "Electron is not defined !!" << G4endl;
    }
#endif
  }
  else
  {
    fPDG = theParticle->GetPDGEncoding();
  }
}

// --------------------------------------------------------------------
G4RToEConvForElectron::~G4RToEConvForElectron() 
{}

// --------------------------------------------------------------------
G4double G4RToEConvForElectron::ComputeValue(const G4int Z,
                                             const G4double kinEnergy) 
{
  const G4double cbr1=0.02, cbr2=-5.7e-5, cbr3=1., cbr4=0.072;
  const G4double Tlow=10.*CLHEP::keV, Thigh=1.*CLHEP::GeV;
  const G4double taul = Tlow/CLHEP::electron_mass_c2;
  const G4double log05 = G4Log(0.5);
  const G4double taul12 = std::sqrt(taul);
  const G4double bremfactor = 0.1;

  const G4double Zlog = G4Pow::GetInstance()->logZ(Z);
  const G4double ionpot = 
    1.6e-5*CLHEP::MeV*G4Exp(0.9*Zlog)/CLHEP::electron_mass_c2;
  const G4double ionpotlog = G4Log(ionpot);

  const G4double tau = kinEnergy/CLHEP::electron_mass_c2;
  G4double dEdx = 0.0;

  if(tau<taul)
  {
    G4double t1 = taul+1.;
    G4double t2 = taul+2.;
    G4double tsq = taul*taul;
    G4double beta2 = taul*t2/(t1*t1);
    G4double f = 1.-beta2+G4Log(tsq/2.)
               +(0.5+0.25*tsq+(1.+2.*taul)*log05)/(t1*t1);
    dEdx = Z*(G4Log(2.*taul+4.)-2.*ionpotlog+f)/beta2;
    dEdx *= taul12/std::sqrt(tau);
  }
  else
  {
    G4double t1 = tau+1.;
    G4double t2 = tau+2.;
    G4double tsq = tau*tau;
    G4double beta2 = tau*t2/(t1*t1);
    G4double f = 1.-beta2+G4Log(tsq/2.)
                   +(0.5+0.25*tsq+(1.+2.*tau)*log05)/(t1*t1);
    dEdx = Z*(G4Log(2.*tau+4.)-2.*ionpotlog+f)/beta2;

    // loss from bremsstrahlung follows
    G4double cbrem = (cbr1+cbr2*Z)*(cbr3+cbr4*G4Log(kinEnergy/Thigh));
    dEdx += Z*(Z+1)*cbrem*bremfactor*tau/beta2;
  }

  return dEdx*CLHEP::twopi_mc2_rcl2;
}

// --------------------------------------------------------------------
