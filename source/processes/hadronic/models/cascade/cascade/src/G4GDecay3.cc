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
// $Id$
//
// File:    G4GDecay3.cc
// Author:  Dennis Wright (SLAC)
// Date:    19 April 2013
//
// Description: three-body phase space momentum generator based on 
//              GDECA3 of Geant3
//
// 20130620  Address Coverity #51433, initializing all data members
// 20141201  Fix error message text to show correct class name
// 20150608  M. Kelsey -- Label all while loops as terminating.

#include "G4GDecay3.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"

G4GDecay3::G4GDecay3(const G4double& pMass, const G4double& dMass0,
                     const G4double& dMass1, const G4double& dMass2)
 : loopMax(100), parentMass(pMass), mDaughter0(dMass0), mDaughter1(dMass1),
   mDaughter2(dMass2), pDaughter0(0.), pDaughter1(0.), pDaughter2(0.) {;}


G4bool G4GDecay3::CalculateMomentumMagnitudes()
{
  G4int looper = 0;
  G4bool status;

  G4double rndm;
  G4double rndm1;
  G4double rndm2;

  G4double momentummax;
  G4double momentumsum;
  G4double energy;

  G4double availableE = parentMass - mDaughter0 - mDaughter1 - mDaughter2;
  do {				/* Loop checking 08.06.2015 MHK */
    rndm1 = G4UniformRand();
    rndm2 = G4UniformRand();
    if (rndm2 > rndm1) {
      // keep randoms in descending order 
      rndm = rndm1;
      rndm1 = rndm2;
      rndm2 = rndm;
    }
    momentummax = 0.0;
    momentumsum = 0.0;

    // daughter 0
    energy = rndm2*availableE;
    pDaughter0 = std::sqrt(energy*energy + 2.0*energy*mDaughter0);
    if (pDaughter0 > momentummax) momentummax = pDaughter0;
    momentumsum += pDaughter0;

    // daughter 1
    energy = (1.-rndm1)*availableE;
    pDaughter1 = std::sqrt(energy*energy + 2.0*energy*mDaughter1);
    if (pDaughter1 > momentummax) momentummax = pDaughter1;
    momentumsum += pDaughter1;

    // daughter 2
    energy = (rndm1-rndm2)*availableE;
    pDaughter2 = std::sqrt(energy*energy + 2.0*energy*mDaughter2);
    if (pDaughter2 > momentummax) momentummax = pDaughter2;
    momentumsum += pDaughter2;
    looper++;
    status = looper < loopMax;
  } while ((momentummax > momentumsum - momentummax) && status);

  return status;
}


std::vector<G4ThreeVector> G4GDecay3::GetThreeBodyMomenta()
{

  std::vector<G4ThreeVector> pVect;

  if (CalculateMomentumMagnitudes() ) {
    
    // Calculate directions
    G4double costheta = 2.*G4UniformRand()-1.;
    G4double sintheta = std::sqrt((1.0-costheta)*(1.0+costheta));
    G4double phi = twopi*G4UniformRand();
    G4double sinphi = std::sin(phi);
    G4double cosphi = std::cos(phi);
    G4ThreeVector direction0(sintheta*cosphi, sintheta*sinphi, costheta);

    G4double costhetan = (pDaughter1*pDaughter1 - pDaughter2*pDaughter2
                          - pDaughter0*pDaughter0)/(2.0*pDaughter2*pDaughter0);
    G4double sinthetan = std::sqrt((1.0-costhetan)*(1.0+costhetan));
    G4double phin = twopi*G4UniformRand();
    G4double sinphin = std::sin(phin);
    G4double cosphin = std::cos(phin);
    G4ThreeVector direction2;
    direction2.setX(sinthetan*cosphin*costheta*cosphi -
                    sinthetan*sinphin*sinphi + costhetan*sintheta*cosphi);
    direction2.setY(sinthetan*cosphin*costheta*sinphi +
                    sinthetan*sinphin*cosphi + costhetan*sintheta*sinphi);
    direction2.setZ(-sinthetan*cosphin*sintheta + costhetan*costheta);

    // Return momentum vectors
    pVect.push_back(pDaughter0*direction0);
    pVect.push_back(-direction0*pDaughter0 - direction2*pDaughter2);
    pVect.push_back(pDaughter2*direction2);

  } else {
    G4cerr << "G4GDecay3::GetThreeBodyMomenta: " << loopMax
           << " or more loops in momentum magnitude calculation " << G4endl;
  }

  return pVect;
}

