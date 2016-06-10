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
// $Id: G4HETCNeutron.cc 90337 2015-05-26 08:34:27Z gcosmo $
//
// by V. Lara
//
// Modified:
// 23.08.2010 V.Ivanchenko general cleanup, move constructor and destructor 
//            the source, use G4Pow

#include "G4HETCNeutron.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Neutron.hh"

G4HETCNeutron::G4HETCNeutron() 
  : G4HETCFragment(G4Neutron::Neutron(), &theNeutronCoulombBarrier)
{}

G4HETCNeutron::~G4HETCNeutron() 
{}

G4double G4HETCNeutron::GetAlpha() const
{
  return 0.76+2.2/theResA13;
}
  
G4double G4HETCNeutron::GetBeta() const
{
  return (2.12/(theResA13*theResA13)-0.05)*MeV/GetAlpha();
}

G4double G4HETCNeutron::GetSpinFactor() const
{
  // (2s+1)
  return 2.0;
}
  
G4double G4HETCNeutron::K(const G4Fragment & aFragment)
{
  // Number of protons in emitted fragment
  G4int Pa = theZ;
  // Number of neutrons in emitted fragment 
  G4int Na = theA - Pa;

  G4double r = G4double(theResZ)/G4double(theResA);
  
  G4int P = aFragment.GetNumberOfParticles();
  G4int H = aFragment.GetNumberOfHoles();
  
  G4double result = 0.0;
  if (P > 0)
    {
      result = (H + Na/(1.0-r))/P;
    }
  
  return std::max(0.0,result);
}

G4double G4HETCNeutron::SampleKineticEnergy(const G4Fragment & aFragment)
{
  G4int H = aFragment.GetNumberOfHoles();
  G4int Pb = aFragment.GetNumberOfParticles();
  G4int Nb = Pb + H;
  G4double g0 = (6.0/pi2)*theFragA*theParameters->GetLevelDensity();
  
  G4double Ab = std::max(0.0,G4double(Pb*Pb+H*H+Pb-3*H)/(4.0*g0));
  G4double Emax = theMaxKinEnergy - Ab;
  
  G4double cut = GetBeta() / (GetBeta()+Emax/G4double(Nb+1));
  G4double x(0.0);
  if (G4UniformRand() <= cut)
    {
      x = BetaRand(Nb,1);
    }
  else 
    {
      x = BetaRand(Nb,2);
    }

  return Emax * (1.0 - x);
}
