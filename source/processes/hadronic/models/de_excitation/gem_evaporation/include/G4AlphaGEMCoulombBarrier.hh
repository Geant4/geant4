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
// J. M. Quesada (July 2009) based on G4AlphaCoulombBarrier
// Coded strictly according to Furihata's GEM paper 
//


#ifndef G4AlphaGEMCoulombBarrier_h
#define G4AlphaGEMCoulombBarrier_h 1

#include "G4GEMCoulombBarrier.hh"
#include "globals.hh"

class G4AlphaGEMCoulombBarrier : public G4GEMCoulombBarrier
{
public:
  G4AlphaGEMCoulombBarrier() : G4GEMCoulombBarrier(4,2) {};
  ~G4AlphaGEMCoulombBarrier() {};
  
private:
  G4AlphaGEMCoulombBarrier(const G4AlphaGEMCoulombBarrier & right);
  
  const G4AlphaGEMCoulombBarrier & operator=(const G4AlphaGEMCoulombBarrier & right);
  G4bool operator==(const G4AlphaGEMCoulombBarrier & right) const;
  G4bool operator!=(const G4AlphaGEMCoulombBarrier & right) const;
  
public:
  G4double BarrierPenetrationFactor(G4double aZ) const
  {
    // Data comes from 
    // Dostrovsky, Fraenkel and Friedlander
    // Physical Review, vol 116, num. 3 1959
    // (JMQ 190709: according to notes added on proof)
    // dataKa = {{20, 0.81}, {30, 0.85}, {40, 0.89}, {50, 0.93}};
    //
    G4double K = 1.0;	
    if (aZ >= 50){
      K=0.93;     
    } else if (aZ <= 20) {
      K=0.81; 
    } else K=0.729802+0.00402544*aZ-1.17276*1e-6*aZ*aZ+2.31248*1e-8*aZ*aZ*aZ-1.65177*1e-10*aZ*aZ*aZ*aZ;	
    return K;
  }
};

#endif
