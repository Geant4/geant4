// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VShortLivedParticle.hh,v 1.2 1999-04-13 08:18:29 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file 
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on object model of
//      Hisaya Kurashige, 27 June 1998
// ----------------------------------------------------------------

#ifndef G4VShortLivedParticle_h
#define G4VShortLivedParticle_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4ParticleDefinition.hh"

class G4VShortLivedParticle : public G4ParticleDefinition
{
  //  A virtual class for short lived particles. 
  //  ShortLivedParticle particles will not be tracked by the TrackingManager
  //  So, G4VShortLivedParticle is not derived from G4ParticleWIthCuts 
  private:

   const G4VShortLivedParticle & operator=(const G4VShortLivedParticle &right);

  public:

   G4VShortLivedParticle(const G4String&  aName,  
               G4double         mass,     
               G4double         width,
               G4double         charge,   
               G4int            iSpin,
               G4int            iParity,
               G4int            iConjugation,
               G4int            iIsospin,   
               G4int            iIsospinZ, 
               G4int            gParity,
               const G4String&  pType,
               G4int            lepton,
               G4int            baryon,
               G4int            encoding,
               G4bool           stable,
               G4double         lifetime,
               G4DecayTable     *decaytable);

   virtual ~G4VShortLivedParticle() {};

   G4int operator==(const G4VShortLivedParticle &right) const;
   G4int operator!=(const G4VShortLivedParticle &right) const;

  public:
      // These methods concerning cut values are not supported for shortlives
      virtual void              ResetCuts();
      virtual void              SetCuts(G4double );
      virtual void              ReCalcCuts();
      virtual G4double      	GetLengthCuts() const;
      virtual G4double*	        GetEnergyCuts() const;
      virtual G4double      	GetEnergyThreshold(const G4Material* ) const;

};

#endif




