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
#ifndef G4HadSecondary_hh
#define G4HadSecondary_hh

#include "G4DynamicParticle.hh"

class G4HadSecondary
{
  public:
    G4HadSecondary(G4DynamicParticle * aT, G4double aWeight = 1);
    G4DynamicParticle * GetParticle() {return theP;}
    G4double GetWeight() {return theWeight;}
    void SetWeight(G4double aW){theWeight= aW;}
    void SetTime(G4double aT) {theTime = aT;}
    G4double GetTime() {return theTime;}
    
    
  private:
  
   G4HadSecondary(){};
   
   G4DynamicParticle * theP; 
   G4double theWeight;
   G4double theTime;
};

#endif
