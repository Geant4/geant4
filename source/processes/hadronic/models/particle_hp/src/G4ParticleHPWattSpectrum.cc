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
// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//

#include "G4ParticleHPWattSpectrum.hh"
#include "G4SystemOfUnits.hh"

  G4double G4ParticleHPWattSpectrum::Sample(G4double anEnergy) 
  {
    G4double a = theApar.GetY(anEnergy)*eV;
    G4double b = theBpar.GetY(anEnergy)/eV;
    G4double result=0.;
    G4double random, cut, max;
    max = std::sinh(std::sqrt(b*15.*a));

    G4int icounter=0;
    G4int icounter_max=1024;
    do
    {
      icounter++;
      if ( icounter > icounter_max ) {
	 G4cout << "Loop-counter exceeded the threshold value at " << __LINE__ << "th line of " << __FILE__ << "." << G4endl;
         break;
      }
      random = G4UniformRand();
      result = -a*G4Log(random);
      cut = G4UniformRand();
    }
    while(cut>std::sinh(std::sqrt(b*result))/max); // Loop checking, 11.05.2015, T. Koi
    return result;
  }
