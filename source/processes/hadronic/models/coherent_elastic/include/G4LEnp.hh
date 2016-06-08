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
 // G4 Low energy model: n-p scattering
 // F.W. Jones, L.G. Greeniaus, H.P. Wellisch
 //  
 // For further comments see G4LEnpData.hh and G4LEnp.cc
 //

#ifndef G4LEnp_h
#define G4LEnp_h 1
 
#include "globals.hh"
#include "Randomize.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"
#include "G4ElementTable.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsVector.hh"
#include "G4LPhysicsFreeVector.hh"
#include "G4Gamma.hh"
#include "G4Step.hh"
#include "G4TrackStatus.hh"
#include "G4HadronicInteraction.hh"



class G4LEnp : public G4HadronicInteraction
{
 private:

   enum { NENERGY=21, NANGLE=180 };

 public:

   G4LEnp();

   ~G4LEnp();
 
   G4VParticleChange* ApplyYourself(const G4Track& aTrack,
                                    G4Nucleus& targetNucleus);


 private:

 // The following arrays are declared static to allow the use of initializers.
 // They are initialized in G4LEnpData.hh

   static G4float sig[NENERGY][NANGLE];
   static G4float pcm[NENERGY], elab[NENERGY], 
     dsigmax[NENERGY], sigtot[NENERGY];

};

#endif
