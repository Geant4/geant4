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
// $Id: G4HETCChargedFragment.cc 96527 2016-04-20 08:51:00Z gcosmo $
//
// by V. Lara
//
// Modified:
// 23.08.2010 V.Ivanchenko general cleanup, move constructor and destructor 
//            the source, use G4Pow
//

#include "G4HETCChargedFragment.hh"
#include "G4PhysicalConstants.hh"
#include "G4VCoulombBarrier.hh"

G4HETCChargedFragment::G4HETCChargedFragment(
  const G4ParticleDefinition* pd, G4VCoulombBarrier * aCoulombBarrier)
  : G4HETCFragment(pd, aCoulombBarrier)
{}

G4HETCChargedFragment::~G4HETCChargedFragment()
{}

G4double G4HETCChargedFragment::
SampleKineticEnergy(const G4Fragment & aFragment)
{
  G4int Pb = aFragment.GetNumberOfParticles();
  G4int H = aFragment.GetNumberOfHoles();

  G4double g0 = (6.0/pi2)*theFragA*theParameters->GetLevelDensity();

  G4double Ab = std::max(0.0,G4double(Pb*Pb+H*H+Pb-3*H)/(4.0*g0));
  G4double Emax = theMaxKinEnergy - Ab;

  G4double x = BetaRand(Pb + H, 2);
  
  return Emax - (Emax-theCoulombBarrier)*x;
}
