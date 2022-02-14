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
// G4RToEConvForGamma class implementation
//
// Author: H.Kurashige, 05 October 2002 - First implementation
// --------------------------------------------------------------------

#include "G4RToEConvForGamma.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4Log.hh"
#include "G4Exp.hh"
#include "G4Pow.hh"

// --------------------------------------------------------------------
G4RToEConvForGamma::G4RToEConvForGamma()  
  : G4VRangeToEnergyConverter()
{    
  theParticle = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
  if (theParticle == nullptr)
  {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0)
    {
      G4cout << " G4RToEConvForGamma::G4RToEConvForGamma() - ";
      G4cout << "Gamma is not defined !!" << G4endl;
    }
#endif
  } 
  else
  {
    fPDG = theParticle->GetPDGEncoding();
  }
}

// --------------------------------------------------------------------
G4RToEConvForGamma::~G4RToEConvForGamma()  
{}

// --------------------------------------------------------------------
G4double G4RToEConvForGamma::ComputeValue(const G4int Z, 
                                          const G4double energy) 
{
  // Compute the "absorption" cross-section of the photon "absorption".
  // Cross-section means here the sum of the cross-sections of the
  // pair production, Compton scattering and photoelectric processes

  const G4double t1keV = 1.*CLHEP::keV;
  const G4double t200keV = 200.*CLHEP::keV;
  const G4double t100MeV = 100.*CLHEP::MeV;

  G4double Zsquare = Z*Z;
  G4double Zlog = G4Pow::GetInstance()->logZ(Z);
  G4double Zlogsquare = Zlog*Zlog;

  G4double tmin = (0.552+218.5/Z+557.17/Zsquare)*CLHEP::MeV;
  G4double tlow = 0.2*G4Exp(-7.355/std::sqrt(Z))*CLHEP::MeV;

  G4double smin = (0.01239+0.005585*Zlog-0.000923*Zlogsquare)*G4Exp(1.5*Zlog);
  G4double s200keV = (0.2651-0.1501*Zlog+0.02283*Zlogsquare)*Zsquare;

  G4double cminlog = G4Log(tmin/t200keV); 
  G4double cmin = G4Log(s200keV/smin)/(cminlog*cminlog);

  G4double slowlog = G4Log(t200keV/tlow);
  G4double slow = s200keV * G4Exp(0.042*Z*slowlog*slowlog);
  G4double logtlow = G4Log(tlow/t1keV);
  G4double clow = G4Log(300.*Zsquare/slow)/logtlow;
  G4double chigh = (7.55e-5 - 0.0542e-5*Z)*Zsquare*Z/G4Log(t100MeV/tmin);

  // Calculate the cross-section (using an approximate empirical formula)
  G4double xs;
  if ( energy < tlow )
  {
    xs = (energy < t1keV) ? slow*G4Exp(clow*logtlow) :
      slow*G4Exp(clow*G4Log(tlow/energy));
  }
  else if ( energy < t200keV )
  {
    G4double x = G4Log(t200keV/energy);
    xs = s200keV * G4Exp(0.042*Z*x*x);
  }
  else if( energy<tmin )
  {
    const G4double x = G4Log(tmin/energy);
    xs = smin * G4Exp(cmin*x*x);
  }
  else
  {
    xs = smin + chigh*G4Log(energy/tmin);
  }
  return xs * CLHEP::barn;
}

// --------------------------------------------------------------------
