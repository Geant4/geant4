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
// $Id: G4ITReactionChange.hh 64057 2012-10-30 15:04:49Z gcosmo $
//
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr) 
//
// WARNING : This class is released as a prototype.
// It might strongly evolve or even disapear in the next releases.
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#ifndef G4ITReactionChange_H
#define G4ITReactionChange_H

#include "globals.hh"
#include "G4ParticleChange.hh"
#include "G4TrackFastVector.hh"

/** Similar to G4ParticleChange, but deal with two tracks
  * rather than one.
  */

class G4ITReactionChange
{
public:
    /** Default constructor */
    G4ITReactionChange();
    /** Default destructor */
    virtual ~G4ITReactionChange();

    //_____________________________________________________
    // To be used by reaction processes
    void Initialize(const G4Track&,
                    const G4Track&,
                    G4VParticleChange* particleChangeA = 0,
                    G4VParticleChange* particleChangeB = 0
                    ) ;
    void AddSecondary(G4Track* aSecondary);
    inline void KillParents(G4bool);

    // If both parents are not killed therefore
    // we can used the particle change
    // you will have to give the pointers
    // of your particleChange while intializing
    // but it might tell you that energy is not
    // conserved for instance
    G4VParticleChange* GetParticleChange(const G4Track*);

    //_____________________________________________________
    // Not to be used in reaction processes
    void UpdateStepInfo(G4Step*, G4Step*);
    G4Track* GetSecondary(G4int) const;
    G4TrackFastVector* GetfSecondary() ;

    G4int GetNumberOfSecondaries() const;
    G4bool WereParentsKilled() const;

    const G4Track* GetTrackA();
    const G4Track* GetTrackB();

protected:
    /** Copy constructor
     *  \param other Object to copy from
     */
    G4ITReactionChange(const G4ITReactionChange& other);
    /** Assignment operator
     *  \param other Object to assign from
     *  \return A reference to this
     */
    G4ITReactionChange& operator=(const G4ITReactionChange& other);
    // equal/unequal operator
    G4bool operator==(const G4ITReactionChange &right) const;
    G4bool operator!=(const G4ITReactionChange &right) const;
    // "equal" means that the objects have the same pointer.
protected:
    std::map<const G4Track*, G4VParticleChange*> fParticleChange;
    G4TrackFastVector* fSecondaries ;
    G4int fNumberOfSecondaries;
    G4bool fKillParents ;
    G4bool fParticleChangeIsSet;
};

inline G4Track* G4ITReactionChange::GetSecondary(G4int anIndex) const
{
    return (*fSecondaries)[anIndex];
}

inline G4int G4ITReactionChange::GetNumberOfSecondaries() const
{
    return fNumberOfSecondaries;
}

inline void G4ITReactionChange::KillParents(G4bool kill)
{
    fKillParents = kill;
}

inline G4bool G4ITReactionChange::WereParentsKilled() const
{
    return fKillParents ;
}

inline G4TrackFastVector* G4ITReactionChange::GetfSecondary()
{
    return fSecondaries;
}

#endif // G4ITReactionChange_H
