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
// $Id: G4VPartonStringModel.hh 100828 2016-11-02 15:25:59Z gcosmo $
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
    G4VPartonStringModel(const G4String& modelName = "Parton String Model");
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
    virtual void ModelDescription(std::ostream& outFile) const;
    virtual G4V3DNucleus * GetProjectileNucleus() const;

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

inline void G4VPartonStringModel::SetThisPointer(G4VPartonStringModel * aPointer)
{
  theThis=aPointer;
}

#endif

