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

#ifndef G4VMolecularDecayDisplacer_h
#define G4VMolecularDecayDisplacer_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include <vector>

class G4Molecule;
class G4MolecularDecayChannel;

typedef int DisplacementType;

class G4VMolecularDecayDisplacer
{
public :
    virtual std::vector<G4ThreeVector> GetProductsDisplacement(const G4MolecularDecayChannel*) const = 0;
    virtual G4ThreeVector GetMotherMoleculeDisplacement(const G4MolecularDecayChannel*) const = 0;
    inline void SetVerbose(G4int);
    virtual ~G4VMolecularDecayDisplacer();

#if defined G4EM_ALLOC_EXPORT
    G4DLLEXPORT static const DisplacementType NoDisplacement;
#else
    G4DLLIMPORT static const DisplacementType NoDisplacement;
#endif

protected :
    G4VMolecularDecayDisplacer();
    G4int fVerbose ;
    static DisplacementType AddDisplacement();
    static DisplacementType Last;
};

void G4VMolecularDecayDisplacer :: SetVerbose(G4int verbose)
{
    fVerbose = verbose ;
}
#endif


