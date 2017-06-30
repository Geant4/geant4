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
// $Id: G4Generator2BS.cc 104410 2017-05-30 07:17:09Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4Generator2BS
//
// Author:        Andreia Trindade (andreia@lip.pt)
//                Pedro Rodrigues  (psilva@lip.pt)
//                Luis Peralta     (luis@lip.pt)
//
// Creation date: 2 June 2003
//
// Modifications: 
// 02 Jun 2003               First implementation acording with new design
// 05 Nov 2003  MGP          Fixed std namespace
// 17 Nov 2003  MGP          Fixed compilation problem on Windows                  
// 12 Oct 2010  V.Ivanchenko Moved RejectionFunction inline, use G4Pow to speadup
// 09 May 2011  L.Pandola    Initialize private members, to avoid Coverity warning
//
// Class Description: 
//
// Concrete base class for Bremsstrahlung Angular Distribution Generation 
// 2BS Distribution
//
// Class Description: End 
//
// -------------------------------------------------------------------
//

#include "G4Generator2BS.hh"
#include "Randomize.hh"   
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Pow.hh"   

G4Generator2BS::G4Generator2BS(const G4String&)
  : G4VEmAngularDistribution("AngularGen2BS"),fz(1),ratio(1),
    ratio1(1),ratio2(1),delta(0)
{
  g4pow = G4Pow::GetInstance();
  nwarn = 0;
}

G4Generator2BS::~G4Generator2BS() 
{}

G4ThreeVector& G4Generator2BS::SampleDirection(const G4DynamicParticle* dp,
					       G4double final_energy,
					       G4int Z,
					       const G4Material*)
{

  // Adapted from "Improved bremsstrahlung photon angular sampling in the EGS4 code system"
  // by Alex F. Bielajew, Rahde Mohan anc Chen-Shou Chui, PIRS-0203
  // Ionizing Radiation Standards
  // Institute for National Measurement Standards 
  // National Research Council of Canada
  // Departement of Medical Physics, Memorial Sloan-Kettering Cancer Center, New York

  G4double energy = dp->GetTotalEnergy();
  ratio = final_energy/energy;
  ratio1 = (1 + ratio)*(1 + ratio);
  ratio2 = 1 + ratio*ratio;

  G4double gamma = energy/electron_mass_c2;
  G4double beta  = std::sqrt((gamma - 1)*(gamma + 1))/gamma;

  // VI speadup
  fz = 0.00008116224*g4pow->Z13(Z)*g4pow->Z13(Z+1);

  // majoranta
  G4double ymax = 2*beta*(1 + beta)*gamma*gamma;
  G4double gMax = RejectionFunction(0.0);
  gMax = std::max(gMax,RejectionFunction(ymax));

  G4double y, gfun;

  do{
    G4double q = G4UniformRand();
    y = q*ymax/(1 + ymax*(1 - q));
    gfun = RejectionFunction(y);

    // violation point
    if(gfun > gMax && nwarn >= 20) {
      ++nwarn;
      G4cout << "### WARNING in G4Generator2BS: Etot(MeV)= " << energy/MeV 
	     << "  Egamma(MeV)" << (energy - final_energy)/MeV
	     << " gMax= " << gMax << "  < " << gfun
	     << "  results are not reliable!" 
	     << G4endl;
      if(20 == nwarn) { 
	G4cout << "   WARNING in G4Generator2BS is closed" << G4endl; 
      }
    }

  } while(G4UniformRand()*gMax > gfun || y > ymax);

  
  G4double cost = 1 - 2*y/ymax;
  G4double sint = std::sqrt((1 - cost)*(1 + cost));
  G4double phi  = twopi*G4UniformRand(); 

  fLocalDirection.set(sint*std::cos(phi), sint*std::sin(phi),cost);
  fLocalDirection.rotateUz(dp->GetMomentumDirection());

  return fLocalDirection;
}

void G4Generator2BS::PrintGeneratorInformation() const
{
  G4cout << "\n" << G4endl;
  G4cout << "Bremsstrahlung Angular Generator is 2BS Generator "
	 << "from 2BS Koch & Motz distribution (Rev Mod Phys 31(4), 920 (1959))" << G4endl;
  G4cout << "Sampling algorithm adapted from PIRS-0203" << G4endl;
  G4cout << "\n" << G4endl;
} 

