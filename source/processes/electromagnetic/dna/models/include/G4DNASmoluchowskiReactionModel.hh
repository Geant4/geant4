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
// $Id: G4DNASmoluchowskiReactionModel.hh 65695 2012-11-27 11:39:12Z gunter $
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

#ifndef G4DiffControlledReactionRadius_
#define G4DiffControlledReactionRadius_

#include "AddClone_def.hh"
#include "G4VDNAReactionModel.hh"
#include <vector>

class G4DNAMolecularReactionData;

/**
  * G4DNASmoluchowskiReactionModel should be used
  * for very fast reactions (high reaction rate) : the reactions between
  * reactants occuring at encounter.
  * When the time step is constrained this model
  * uses brownian bridge : "Absorbing (Smoluchowski) boundary condition"
  * Reference : On the simulation of diffusion processes close to boundaries,
  * N. J. B. Green, Molecular Physics, 65: 6, 1399 — 1408(1988)
  */

class G4DNASmoluchowskiReactionModel : public G4VDNAReactionModel
{
public :
    G4DNASmoluchowskiReactionModel();
    virtual ~G4DNASmoluchowskiReactionModel();

    G4DNASmoluchowskiReactionModel(const G4DNASmoluchowskiReactionModel&);

    G4IT_ADD_CLONE(G4VDNAReactionModel, G4DNASmoluchowskiReactionModel)

    virtual void Initialise(const G4Molecule*, const G4Track&) ;
    virtual void InitialiseToPrint(const G4Molecule*) ;
    virtual G4double GetReactionRadius(const G4Molecule*, const G4Molecule*);
    virtual G4double GetReactionRadius(const G4int);

    virtual G4bool FindReaction(const G4Track&,
                                const G4Track&,
                                const G4double /*reactionRadius*/,
                                G4double& /*separationDistance*/,
                                const G4bool /*alongStepInteraction*/) ;

private :
    const std::vector<const G4DNAMolecularReactionData*>* fReactionData ;
    G4DNASmoluchowskiReactionModel& operator=(const G4DNASmoluchowskiReactionModel&);
};

#endif
