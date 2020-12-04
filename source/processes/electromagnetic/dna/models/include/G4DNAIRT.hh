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
/*
 * G4DNAIRT.hh
 *
 *  Created on: Jul 23, 2019
 *      Author: W. G. Shin
 *              J. Ramos-Mendez and B. Faddegon
*/

#ifndef G4DNAIRT_HH_
#define G4DNAIRT_HH_


#include "globals.hh"
#include "G4ThreeVector.hh"

#include "G4DNAMolecularReaction.hh"
#include "G4DNAMolecularReactionTable.hh"
#include "G4MoleculeTable.hh"

#include "G4VDNAReactionModel.hh"
#include "G4VITReactionProcess.hh"

#include "G4ITReactionTable.hh"
#include "G4ITTrackHolder.hh"
#include "G4ITReaction.hh"

#include "G4Molecule.hh"
#include "G4VITReactionProcess.hh"
#include "G4ParticleChange.hh"

#include "AddClone_def.hh"
#include <vector>
#include <map>

class G4DNAMolecularReactionTable;
class G4VDNAReactionModel;
class G4ErrorFunction;

class G4DNAIRT : public G4VITReactionProcess
{
public:

    G4DNAIRT();
    explicit G4DNAIRT(G4VDNAReactionModel*);
    ~G4DNAIRT() override;
    G4DNAIRT(const G4DNAIRT& other) = delete;
    G4DNAIRT& operator=(const G4DNAIRT& other) = delete;

    G4bool TestReactibility(const G4Track&,
                            const G4Track&,
                            double ,
                            bool ) override;
    std::vector<std::unique_ptr<G4ITReactionChange>> FindReaction(G4ITReactionSet*, const double, const double, const bool) override;
    std::unique_ptr<G4ITReactionChange> MakeReaction(const G4Track&, const G4Track&) override;

    void SetReactionModel(G4VDNAReactionModel*);

    void Initialize() override;
    void SpaceBinning();
    void IRTSampling();
    void Sampling(G4Track*);

    G4double GetIndependentReactionTime(const G4MolecularConfiguration*, const G4MolecularConfiguration*, G4double);
    G4int FindBin(G4int, G4double, G4double, G4double);
    G4double SamplePDC(G4double , G4double );

protected:
    const G4DNAMolecularReactionTable*& fMolReactionTable;
    G4VDNAReactionModel* fpReactionModel;

private:
    G4ITTrackHolder* fTrackHolder;
    G4ITReactionSet* fReactionSet;
    G4ErrorFunction* erfc;

    std::map<G4int,std::map<G4int,std::map<G4int,std::vector<G4Track*>>>> spaceBinned;

    G4double fRCutOff;
    G4double timeMin;
    G4double timeMax;

    G4double fXMin, fYMin, fZMin;
    G4double fXMax, fYMax, fZMax;
    G4int fNx, fNy, fNz;
    G4int xiniIndex, yiniIndex, ziniIndex;
    G4int xendIndex, yendIndex, zendIndex;

};

#endif /* G4DNAIRT_HH_ */
