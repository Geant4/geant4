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
#ifndef VANAPlot_h
#define VANAPlot_h

#include <iostream>
#include "HadronNuclear/batch_monitor/Particle.h"
#include "globals.hh"

class VANAPlot
{
  public:
  
    virtual G4bool Insert(G4int aPdg, G4double anX, G4double anEntryToAccumulate) = 0;
    virtual void DumpInfo(ostream &, G4String aPreFix) = 0;
    virtual void SetNevents(G4int aNumber) = 0;
    virtual G4bool Filter(ANAParticle * aParticle) = 0;
    
  private:
};

#endif
