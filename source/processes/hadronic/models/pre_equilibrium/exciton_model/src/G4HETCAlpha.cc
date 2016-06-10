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
// $Id: G4HETCAlpha.cc 90337 2015-05-26 08:34:27Z gcosmo $
//
// by V. Lara
//
// Modified:
// 23.08.2010 V.Ivanchenko general cleanup, move constructor and destructor 
//            the source, use G4Pow

#include "G4HETCAlpha.hh"
#include "G4Alpha.hh"

G4HETCAlpha::G4HETCAlpha() 
  : G4HETCChargedFragment(G4Alpha::Alpha(), &theAlphaCoulombBarrier)
{}

G4HETCAlpha::~G4HETCAlpha() 
{}

G4double G4HETCAlpha::GetAlpha() const
{
  G4double C = 0.0;
  if (theFragZ <= 30) 
    {
      C = 0.10;
    } 
  else if (theFragZ <= 50) 
    {
      C = 0.1 + -((theFragZ-50.)/20.)*0.02;
    }
  else if (theFragZ < 70) 
    {
      C = 0.08 + -((theFragZ-70.)/20.)*0.02;
    }
  else 
    {
      C = 0.06;
    }
  return 1.0+C;
}
  
G4double G4HETCAlpha::GetBeta() const
{
  return theCoulombBarrier;
}
  
G4double G4HETCAlpha::GetSpinFactor() const
{
  return 1.0;
}

G4double G4HETCAlpha::K(const G4Fragment & aFragment)
{
  // Number of protons in emitted fragment
  G4int Pa = theZ;
  // Number of neutrons in emitted fragment 
  G4int Na = theA - Pa;

  G4double r = G4double(theResZ)/G4double(theResA);

  G4int P = aFragment.GetNumberOfParticles();
  G4int H = aFragment.GetNumberOfHoles();

  G4double result = 0.0;
  if (P > 3)
    {
      result = 3.0/(P*(P-1.0)*(P-2.0)*(P-3.0)) * 
	(H*(H-1.0)*(H-2.0)*(H-3.0)*r*r*(r-1.0)*(r-1.0) +
	 2.0*H*(H-1.0)*(H-2.0)*(Pa*r*(1.0-r)*(1.0-r)+Na*r*r*(1.0-r)) +
	 H*(H-1.0)*(Pa*(Pa-1.0)*(1.0-r)*(1.0-r)+4.0*Na*Pa*r*(1.0-r)+Na*(Na-1.0)*r*r) +
	 2*H*(Pa*Na*(Na-1.0)*r+Pa*(Pa-1.0)*Na*(1.0-r)) +
	 Pa*(Pa-1.0)*Na*(Na-1.0));

      result /= 6.0*r*r*(1. - r) *(1. - r);
    }
  return std::max(0.0,result);
}
