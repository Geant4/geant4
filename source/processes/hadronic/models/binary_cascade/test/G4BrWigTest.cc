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
#include "globals.hh"
#include "G4ios.hh"

#include "CLHEP/Hist/TupleManager.h"
#include "CLHEP/Hist/HBookFile.h"
#include "CLHEP/Hist/Histogram.h"
#include "CLHEP/Hist/Tuple.h"
#include "Randomize.hh"

#include <fstream>
#include <iomanip>
#include <assert.h>
#include <iostream>
#include <G4KineticTrack.hh>

#include "G4PionMinus.hh"
#include "G4VShortLivedParticle.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4ParticleTable.hh"
#include "G4ShortLivedTable.hh"

int main()
{
  // MGP ---- HBOOK initialization
  HepTupleManager* hbookManager;
  hbookManager = new HBookFile("brwig.hbook", 58);
  assert (hbookManager != 0);  

  // MGP ---- Book histograms
 
  HepHistogram* hbrwig;
  HepHistogram* hbmass;
  hbrwig = hbookManager->histogram("Breit Wigner", 1000,0.,1.);
  assert (hbrwig != 0);
  hbmass = hbookManager->histogram("random mass", 1000,0.,2500.);
  assert (hbmass != 0);

  G4ParticleDefinition* pionMinus = G4PionMinus::PionMinusDefinition();
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  G4ShortLivedConstructor ShortLived;
  ShortLived.ConstructParticle();

  G4ParticleDefinition* deltaPlus = particleTable->FindParticle("delta+");
  G4cout << deltaPlus->GetParticleName() << " created, type is " 
	 << deltaPlus->GetParticleType() << G4endl;

  G4ParticleDefinition* definition1;
  definition1 = deltaPlus;
  G4double mass1 = definition1->GetPDGMass();
  G4double width1 = definition1->GetPDGWidth();

  int nevt, DEBUG, nwid;
  cout << "Enter Debug 0/1" << G4endl;
  G4cin >> DEBUG;
  cout << " Enter number of events " << G4endl;
  G4cin >> nevt;
  cout << " Enter number of width " << G4endl;
  G4cin >> nwid;

  G4KineticTrack wig;
  G4double max=wig.BrWig(width1, mass1, mass1);
  int i=0;
  while(i<nevt){
  // rand mass-nwid*width:mass+nwid*width
   G4double mass=(mass1+nwid*width1)-(G4UniformRand()*2*nwid*width1);
   G4double check=max*G4UniformRand();
   G4double wigt=wig.BrWig(width1, mass1, mass);
   // fill histo
   hbrwig->accumulate(wigt);
   if(check < wigt){
     hbmass->accumulate(mass);
     i++;
   }
   if (DEBUG == 1)cout << " w/m/xm/brwig " << width1 << "  " << mass << "  " << mass1 << "  " << wigt << G4endl;
  } 
  hbookManager->write();
  return EXIT_SUCCESS;
}










