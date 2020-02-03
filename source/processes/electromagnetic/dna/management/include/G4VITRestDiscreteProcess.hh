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
/// \brief Identical to G4VRestDiscreteProcess with dependency from G4VITProcess
//
// WARNING : This class is released as a prototype.
// It might strongly evolve or even disapear in the next releases.
//
// Author: Mathieu Karamitros

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

#pragma  once
#include "G4VITProcess.hh"

/** Identical to G4VRestDiscreteProcess with dependency on G4VITProcess */

class G4VITRestDiscreteProcess : public G4VITProcess
{
    //  Abstract class which defines the public behavior of
    //  rest + discrete physics interactions.
public:
    G4VITRestDiscreteProcess() = delete;
    G4VITRestDiscreteProcess(const G4String&, G4ProcessType aType = fNotDefined);
    G4VITRestDiscreteProcess(const G4VITRestDiscreteProcess&) = delete;
    G4VITRestDiscreteProcess& operator=(const G4VITRestDiscreteProcess& right) = delete;
    ~G4VITRestDiscreteProcess() override;

public:
    // with description
    G4double
    PostStepGetPhysicalInteractionLength(const G4Track& track,
                                         G4double previousStepSize,
                                         G4ForceCondition* condition) override;

    G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&) override;

    G4double AtRestGetPhysicalInteractionLength(const G4Track&,
                                                G4ForceCondition*) override;

    G4VParticleChange* AtRestDoIt(const G4Track&, const G4Step&) override;

    //  no operation in  AlongStepDoIt
    G4double AlongStepGetPhysicalInteractionLength(const G4Track&,
                                                   G4double,
                                                   G4double,
                                                   G4double&,
                                                   G4GPILSelection*) override
    {
        return -1.0;
    }

    //  no operation in  AlongStepDoIt
    G4VParticleChange* AlongStepDoIt(const G4Track&, const G4Step&) override
    {
        return nullptr;
    }

protected:
    // with description
    virtual G4double GetMeanFreePath(const G4Track& aTrack,
                                     G4double previousStepSize,
                                     G4ForceCondition* condition) = 0;
    //  Calculates from the macroscopic cross section a mean
    //  free path, the value is returned in units of distance.

    virtual G4double GetMeanLifeTime(const G4Track& aTrack,
                                     G4ForceCondition* condition) = 0;
    //  Calculates the mean life-time (i.e. for decays) of the
    //  particle at rest due to the occurrence of the given process,
    //  or converts the probability of interaction (i.e. for
    //  annihilation) into the life-time of the particle for the
    //  occurrence of the given process.
};
