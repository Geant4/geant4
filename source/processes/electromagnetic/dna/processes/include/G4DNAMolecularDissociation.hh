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
// $Id: G4DNAMolecularDissociation.hh 93936 2015-11-04 09:37:59Z gcosmo $
//
// Author: Mathieu Karamitros, kara@cenbg.in2p3.fr

// The code is developed in the framework of the ESA AO7146
//
// We would be very happy hearing from you, send us your feedback! :)
//
// In order for Geant4-DNA to be maintained and still open-source,
// article citations are crucial. 
// If you use Geant4-DNA chemistry and you publish papers about your software, 
// in addition to the general paper on Geant4-DNA:
//
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// we would be very happy if you could please also cite the following
// reference papers on chemistry:
//
// J. Comput. Phys. 274 (2014) 841-882
// Prog. Nucl. Sci. Tec. 2 (2011) 503-508 


#ifndef G4MOLECULARDECAYPROCESS_HH
#define G4MOLECULARDECAYPROCESS_HH

#include "G4VITRestProcess.hh"
#include "G4VITRestDiscreteProcess.hh"
#include <map>

class G4ParticleChange;
class G4VMolecularDecayDisplacer;
class G4MoleculeDefinition;

/**
  * G4DNAMolecularDissociation should be called only for molecules.
  * It will dissociate the molecules using the decay associated to
  * this molecule and if a displacement scheme has been registered,
  * it will place the products to the expected position.
  */

class G4DNAMolecularDissociation: public G4VITRestDiscreteProcess //G4VITRestProcess
{
public:
    G4DNAMolecularDissociation(const G4String& processName = "DNAMolecularDecay",
                            G4ProcessType type = fDecay);

    virtual ~G4DNAMolecularDissociation();

    G4IT_ADD_CLONE(G4VITProcess, G4DNAMolecularDissociation)

    virtual G4bool IsApplicable(const G4ParticleDefinition&);


    virtual G4double
    PostStepGetPhysicalInteractionLength(const G4Track& track,
                                         G4double previousStepSize,
                                         G4ForceCondition* condition);

    inline G4double AtRestGetPhysicalInteractionLength(
        const G4Track& track,
        G4ForceCondition* condition
        );
    inline G4VParticleChange* AtRestDoIt(
        const G4Track& track,
        const G4Step& step
        );

    inline G4VParticleChange* PostStepDoIt(
        const G4Track& track,
        const G4Step& step
        )
    {
      return AtRestDoIt(track, step);
    }

    inline void SetVerbose(G4int);

    //__________________________________________________________________
    void SetDecayDisplacer(const G4ParticleDefinition*, G4VMolecularDecayDisplacer*);
    G4VMolecularDecayDisplacer* GetDecayDisplacer(const G4ParticleDefinition*);

protected:
    //__________________________________________________________________
    // Make the decay
    virtual G4VParticleChange* DecayIt(const G4Track& ,const G4Step&);
    virtual G4double GetMeanLifeTime(const G4Track&,G4ForceCondition*);
    virtual G4double GetMeanFreePath(const G4Track& ,
                                     G4double ,
                                     G4ForceCondition* )
    {
      return 0;
    }

private:
    G4DNAMolecularDissociation();
    G4DNAMolecularDissociation(const G4DNAMolecularDissociation &right);
    G4DNAMolecularDissociation & operator=(const G4DNAMolecularDissociation &right);

private:
    G4bool fDecayAtFixedTime ;

    typedef std::map<const G4ParticleDefinition*, G4VMolecularDecayDisplacer*>  DecayDisplacementMap;
    DecayDisplacementMap fDecayDisplacementMap;

    G4int fVerbose;
};

inline void G4DNAMolecularDissociation::SetVerbose(G4int verbose)
{
    fVerbose = verbose ;
}

inline G4double G4DNAMolecularDissociation::AtRestGetPhysicalInteractionLength(
    const G4Track& track,
    G4ForceCondition* condition)
{
    if(fDecayAtFixedTime)
    {
        return GetMeanLifeTime(track, condition);
    }

    return G4VITRestDiscreteProcess::AtRestGetPhysicalInteractionLength(track, condition);
//    return G4VITRestProcess::AtRestGetPhysicalInteractionLength(track, condition);
}

inline G4VParticleChange* G4DNAMolecularDissociation::AtRestDoIt(
    const G4Track& track,
    const G4Step& step
    )
{
    ClearNumberOfInteractionLengthLeft();
    ClearInteractionTimeLeft();
    return DecayIt(track, step);
}


#endif /* G4MOLECULARDECAYPROCESS_HH */
