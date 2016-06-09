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
// $Id: G4DNAMolecularDecay.hh 65022 2012-11-12 16:43:12Z gcosmo $
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

#ifndef G4MOLECULARDECAYPROCESS_HH
#define G4MOLECULARDECAYPROCESS_HH

#include "G4VITRestProcess.hh"
#include <map>

class G4ParticleChange;
class G4VMolecularDecayDisplacer;
class G4MoleculeDefinition;

/**
  * G4DNAMolecularDecay should be called only for molecules.
  * It will dissociate the molecules using the decay associated to
  * this molecule and if a displacement scheme has been registered,
  * it will place the products to the expected position.
  */

class G4DNAMolecularDecay: public G4VITRestProcess
{
public:
    G4DNAMolecularDecay(const G4String& processName = "DNAMolecularDecay",
                            G4ProcessType type = fDecay);

    virtual ~G4DNAMolecularDecay();

    G4IT_ADD_CLONE(G4VITProcess, G4DNAMolecularDecay)

    virtual G4bool IsApplicable(const G4ParticleDefinition&);

    inline G4double AtRestGetPhysicalInteractionLength(
        const G4Track& track,
        G4ForceCondition* condition
        );
    inline G4VParticleChange* AtRestDoIt(
        const G4Track& track,
        const G4Step& step
        );

    inline void SetVerbose(G4int);

    //__________________________________________________________________
    void SetDecayDisplacer(const G4ParticleDefinition*, G4VMolecularDecayDisplacer*);
    G4VMolecularDecayDisplacer* GetDecayDisplacer(const G4ParticleDefinition*);

protected:
    //__________________________________________________________________
    // Make the decay
    virtual G4VParticleChange* DecayIt(const G4Track& ,const G4Step&);
    virtual G4double GetMeanLifeTime(const G4Track&,G4ForceCondition*);

private:
    G4DNAMolecularDecay();
    G4DNAMolecularDecay(const G4DNAMolecularDecay &right);
    G4DNAMolecularDecay & operator=(const G4DNAMolecularDecay &right);

private:
    G4bool fDecayAtFixedTime ;

    typedef std::map<const G4ParticleDefinition*, G4VMolecularDecayDisplacer*>  DecayDisplacementMap;
    DecayDisplacementMap fDecayDisplacementMap;

    G4int fVerbose;
};

inline void G4DNAMolecularDecay::SetVerbose(G4int verbose)
{
    fVerbose = verbose ;
}

inline G4double G4DNAMolecularDecay::AtRestGetPhysicalInteractionLength(
    const G4Track& track,
    G4ForceCondition* condition)
{
    if(fDecayAtFixedTime)
    {
        return GetMeanLifeTime(track, condition);
    }

    return G4VITRestProcess::AtRestGetPhysicalInteractionLength(track, condition);
}

inline G4VParticleChange* G4DNAMolecularDecay::AtRestDoIt(
    const G4Track& track,
    const G4Step& step
    )
{
    ClearNumberOfInteractionLengthLeft();
    ClearInteractionTimeLeft();
    return DecayIt(track, step);
}


#endif /* G4MOLECULARDECAYPROCESS_HH */
