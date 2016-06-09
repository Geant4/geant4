//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4HadronBuilder.hh,v 1.1 2003/10/07 11:25:40 hpw Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
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
