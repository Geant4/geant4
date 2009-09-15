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
// J. M. Quesada (July 2009) based on G4ProtonCoulombBarrier
// Coded strictly according to Furihata's GEM paper 
//

#ifndef G4ProtonGEMCoulombBarrier_h
#define G4ProtonGEMCoulombBarrier_h 1

#include "G4GEMCoulombBarrier.hh"
#include "globals.hh"

class G4ProtonGEMCoulombBarrier : public G4GEMCoulombBarrier
{
public:
  G4ProtonGEMCoulombBarrier() : G4GEMCoulombBarrier(1,1) {}
  ~G4ProtonGEMCoulombBarrier() {}
  
private:
  G4ProtonGEMCoulombBarrier(const G4ProtonGEMCoulombBarrier & right);
  
  const G4ProtonGEMCoulombBarrier & operator=(const G4ProtonGEMCoulombBarrier & right);
  G4bool operator==(const G4ProtonGEMCoulombBarrier & right) const;
  G4bool operator!=(const G4ProtonGEMCoulombBarrier & right) const;
  
private:
  G4double BarrierPenetrationFactor(const G4double aZ) const {
    // Data comes from 
    // Dostrovsky, Fraenkel and Friedlander
    // Physical Review, vol 116, num. 3 1959
    // (JMQ 190709: according to notes added on proof)
    //dataK = {{20, 0.51}, {30, 0.60}, {40, 0.66}, {50, 0.68}};
    //
    G4double K = 1.0;	
    if (aZ >= 50){
      K=0.68;     
    } else if (aZ <= 20) {
      K=0.51; 
    } else K=0.28445+0.0115956*aZ+0.000026329*aZ*aZ-2.18583*1e-6*aZ*aZ*aZ+3.7083*1e-9*aZ*aZ*aZ*aZ;	
    return K;
  }
};

#endif
