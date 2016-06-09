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
// $Id: G4DNAMolecularDecayDisplacer.hh 64057 2012-10-30 15:04:49Z gcosmo $
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

#ifndef G4DNAMolecularDecayDisplacer_h
#define G4DNAMolecularDecayDisplacer_h 1

#include "G4VMolecularDecayDisplacer.hh"
#include "G4MoleculeID.hh"

class G4DNAMolecularDecayDisplacer : public G4VMolecularDecayDisplacer
{
public :
    G4DNAMolecularDecayDisplacer();
    virtual ~G4DNAMolecularDecayDisplacer() ;

    virtual std::vector<G4ThreeVector> GetProductsDisplacement(const G4MolecularDecayChannel*) const;
    virtual G4ThreeVector GetMotherMoleculeDisplacement(const G4MolecularDecayChannel*) const;
    G4ThreeVector radialDistributionOfElectron() const;

#if defined G4EM_ALLOC_EXPORT
    G4DLLEXPORT static const DisplacementType Ionisation_DissociationDecay;
    G4DLLEXPORT static const DisplacementType A1B1_DissociationDecay;
    G4DLLEXPORT static const DisplacementType B1A1_DissociationDecay;
    G4DLLEXPORT static const DisplacementType AutoIonisation;
    G4DLLEXPORT static const DisplacementType DissociativeAttachment;
#else
    G4DLLIMPORT static const DisplacementType Ionisation_DissociationDecay;
    G4DLLIMPORT static const DisplacementType A1B1_DissociationDecay;
    G4DLLIMPORT static const DisplacementType B1A1_DissociationDecay;
    G4DLLIMPORT static const DisplacementType AutoIonisation;
    G4DLLEXPORT static const DisplacementType DissociativeAttachment;
#endif

private :
    G4ThreeVector radialDistributionOfProducts(G4double) const;
};
#endif

