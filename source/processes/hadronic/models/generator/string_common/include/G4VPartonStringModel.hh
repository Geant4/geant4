// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VPartonStringModel.hh,v 1.2 1998/10/09 13:33:29 maxim Exp $
// GEANT4 tag $Name: geant4-00 $
//
#ifndef G4VPartonStringModel_h
#define G4VPartonStringModel_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4VPartonStringModel ----------------
//             by Gunter Folger, May 1998.
//      abstract class for all Parton String Models
// ------------------------------------------------------------

// Modified at 8-Oct-1998 by Maxim Komogorov. Method EnergyAndMomentumCorrector was added.

#include "G4StringModel.hh"
#include "G4VParticipants.hh"
#include "G4ReactionProductVector.hh"
#include "G4ExcitedString.hh"
#include "G4ExcitedStringVector.hh"
#include "G4VStringFragmentation.hh"
#include "G4V3DNucleus.hh"

class G4VPartonStringModel : public G4VHighEnergyGenerator

{
  public:
      G4VPartonStringModel();
      virtual ~G4VPartonStringModel();

  private:
      G4VPartonStringModel(const G4VPartonStringModel &right);
      const G4VPartonStringModel & operator=(const G4VPartonStringModel &right);
      int operator==(const G4VPartonStringModel &right) const;
      int operator!=(const G4VPartonStringModel &right) const;

  public:
      void SetFragmentationModel(G4VStringFragmentation * aModel);
      G4KineticTrackVector * Scatter(const G4Nucleus &theNucleus, const G4DynamicParticle &thePrimary);
      virtual G4V3DNucleus * GetWoundedNucleus() const = 0;

  protected:
        
      virtual void Init(const G4Nucleus &theNucleus, const G4DynamicParticle &thePrimary) = 0;
      virtual G4ExcitedStringVector * GetStrings() = 0;
      void SetThisPointer(G4VPartonStringModel * aPointer);

      G4bool EnergyAndMomentumCorrector(G4KineticTrackVector* Output, G4LorentzVector& TotalCollisionMomentum);   

  private:
  
      G4VStringFragmentation * stringFragmentationModel;
      G4VPartonStringModel * theThis;

};

inline void G4VPartonStringModel::SetFragmentationModel(G4VStringFragmentation * aModel)
{
  stringFragmentationModel = aModel;
}

inline
void G4VPartonStringModel::SetThisPointer(G4VPartonStringModel * aPointer)
{
	theThis=aPointer;
}
#endif


