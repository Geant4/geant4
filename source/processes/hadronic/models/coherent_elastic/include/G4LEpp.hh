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
 // G4 Low energy model: n-n or p-p scattering
 // F.W. Jones, L.G. Greeniaus, H.P. Wellisch
 //  
 // For further comments see G4LEppData.hh and G4LEpp.cc
 //

#ifndef G4LEpp_h
#define G4LEpp_h 1
 
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


class G4LEpp : public G4HadronicInteraction
{
 private:

   enum { NENERGY=22, NANGLE=180 };

 public:

   G4LEpp();

   ~G4LEpp();
 
   G4VParticleChange* ApplyYourself(const G4Track& aTrack,
                                    G4Nucleus& targetNucleus);

   void SetCoulombEffects(G4int State);
  
 private:

  G4float * sig[NANGLE];
  G4float * elab;

 // The following arrays are declared static to allow the use of initializers.
 // They are initialized in G4LEppData.hh

 // Coulomb effects suppressed:
   static G4float Sig[NENERGY][NANGLE];
   static G4float Pcm[NENERGY], Elab[NENERGY], 
     dSigmax[NENERGY], Sigtot[NENERGY];

 // Coulomb effects not suppressed:
   static G4float SigCoul[NENERGY][NANGLE];
   static G4float PcmCoul[NENERGY], ElabCoul[NENERGY], 
     dSigmaxCoul[NENERGY], SigtotCoul[NENERGY];


};

#endif
