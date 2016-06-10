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
// $Id: G4DNASmoluchowskiReactionModel.hh 91584 2015-07-27 13:01:48Z gcosmo $
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

    virtual void Initialise(G4MolecularConfiguration*, const G4Track&) ;
    virtual void InitialiseToPrint(G4MolecularConfiguration*) ;
    virtual G4double GetReactionRadius(G4MolecularConfiguration*,
                                       G4MolecularConfiguration*);
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
