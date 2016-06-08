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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4VPartonStringModel.hh,v 1.5 2001/08/01 17:08:59 hpw Exp $
// GEANT4 tag $Name: geant4-04-01-patch-01 $
//
#ifndef G4VPartonStringModel_h
#define G4VPartonStringModel_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
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
#include "G4KineticTrackVector.hh"

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


