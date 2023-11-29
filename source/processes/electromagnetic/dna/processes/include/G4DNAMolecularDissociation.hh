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


#pragma once

#include "G4VITRestDiscreteProcess.hh"
#include <map>

class G4ParticleChange;
class G4VMolecularDissociationDisplacer;
class G4MoleculeDefinition;
class G4VUserBrownianAction;
/**
  * G4DNAMolecularDissociation should be called only for molecules.
  * It will dissociate the molecules using the decay associated to
  * this molecule and if a displacement scheme has been registered,
  * it will place the products to the expected position.
  */

class G4DNAMolecularDissociation : public G4VITRestDiscreteProcess
{
public:
    G4DNAMolecularDissociation(const G4String& processName = "DNAMolecularDecay",
                               G4ProcessType type = fDecay);
    G4DNAMolecularDissociation() = delete;
    G4DNAMolecularDissociation(const G4DNAMolecularDissociation& right) = delete;
    G4DNAMolecularDissociation& operator=(const G4DNAMolecularDissociation& right) = delete;
    ~G4DNAMolecularDissociation() override;

    using Species = const G4MoleculeDefinition;
    using Displacer = G4VMolecularDissociationDisplacer;

    G4bool IsApplicable(const G4ParticleDefinition&) override;

    G4double PostStepGetPhysicalInteractionLength(const G4Track& track,
                                                  G4double previousStepSize,
                                                  G4ForceCondition* condition) override;

    G4double AtRestGetPhysicalInteractionLength(const G4Track& track,
                                                G4ForceCondition* condition) override;

    G4VParticleChange* AtRestDoIt(const G4Track& track, const G4Step& step) override;

    G4VParticleChange* PostStepDoIt(const G4Track& track, const G4Step& step) override;

    void SetVerbose(G4int);

    inline void SetUserBrownianAction(G4VUserBrownianAction* pBrownianAction)
    {
      fpBrownianAction = pBrownianAction;
    }

    void SetDisplacer(Species*, Displacer*);
    Displacer* GetDisplacer(Species*);

protected:
    virtual G4VParticleChange* DecayIt(const G4Track&, const G4Step&);

    G4double GetMeanLifeTime(const G4Track&, G4ForceCondition*) override;

    G4double GetMeanFreePath(const G4Track&,
                             G4double,
                             G4ForceCondition*) override;

private:
    using DisplacementMap = std::map<Species*, std::unique_ptr<Displacer>>;
    G4VUserBrownianAction* fpBrownianAction = nullptr;
    G4bool fDecayAtFixedTime;
    DisplacementMap fDisplacementMap;
    G4int fVerbose;
};

