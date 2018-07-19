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
// $Id: G4HadronBuilder.hh 102048 2016-12-19 09:02:38Z gcosmo $
//
// -----------------------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: 
//             Gunter Folger, August/September 2001
//               Create class; 
// -----------------------------------------------------------------------------
//

#ifndef G4HadronBuilder_h
#define G4HadronBuilder_h 1

#include "globals.hh"
#include <vector>
#include "G4ParticleDefinition.hh"
#include "G4ParticleDefinition.hh"

class G4HadronBuilder
{
  public:

     G4ParticleDefinition * Build(G4ParticleDefinition * black, G4ParticleDefinition * white);
     G4ParticleDefinition * BuildLowSpin(G4ParticleDefinition * black, G4ParticleDefinition * white);
     G4ParticleDefinition * BuildHighSpin(G4ParticleDefinition * black, G4ParticleDefinition * white);

     //  ctor
     G4HadronBuilder(G4double mesonMix, G4double barionMix,
		     std::vector<double> scalarMesonMix,
		     std::vector<double> vectorMesonMix); 

  private:

     G4HadronBuilder(); // no default ctor

     enum Spin { SpinZero=1, SpinHalf=2, SpinOne=3, SpinThreeHalf=4 };

     G4ParticleDefinition * Meson(G4ParticleDefinition * black, G4ParticleDefinition * white, Spin spin);

     G4ParticleDefinition * Barion(G4ParticleDefinition * black, G4ParticleDefinition * white, Spin spin);
     
     G4double mesonSpinMix;
     G4double barionSpinMix;
     std::vector<double> scalarMesonMixings;
     std::vector<double> vectorMesonMixings;
     
};

#endif

