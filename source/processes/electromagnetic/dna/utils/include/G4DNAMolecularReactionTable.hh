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
// $Id: G4DNAMolecularReactionTable.hh 64057 2012-10-30 15:04:49Z gcosmo $
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

#ifndef G4MolecularReactionTable_h
#define G4MolecularReactionTable_h 1

#include "G4ITReactionTable.hh"
#include "G4Molecule.hh"
#include <vector>
#include <map>
#include "G4ReferenceCast.hh"

class G4VDNAReactionModel ;
class G4DNAMolecularReactionTable;

/**
  * G4DNAMolecularReactionData contains the information
  * relative to a given reaction (eg : °OH + °OH -> H2O2)
  */

class G4DNAMolecularReactionData
{
public :
    G4DNAMolecularReactionData(G4double reactionRate,
                            const G4Molecule* reactive1,
                            const G4Molecule* reactive2);
    ~G4DNAMolecularReactionData();

    const G4Molecule* GetReactive1() const { return fReactive1.get(); }
    const G4Molecule* GetReactive2() const { return fReactive2.get(); }

    G4double GetReactionRate() const {return fReactionRate;}
    G4double GetReducedReactionRadius() const {return fReducedReactionRadius;}

    //_____________________________________________________
    void SetReactive1(const G4Molecule* reactive) ;
    void SetReactive2(const G4Molecule* reactive) ;
    void SetReactive(const G4Molecule* reactive1, const G4Molecule* reactive2);
    void AddProduct(const G4Molecule* molecule);

    G4int GetNbProducts() const
    {
        if(fProducts) return fProducts->size();
        return 0;
    }

    const G4Molecule* GetProduct(G4int i) const
    {
        if(fProducts) return (*fProducts)[i].get();
        return 0;
    }

protected :
    G4DNAMolecularReactionData();
    G4MoleculeHandle fReactive1;
    G4MoleculeHandle fReactive2;
    G4double fReactionRate;
    G4double fReducedReactionRadius;

    std::vector<G4MoleculeHandle>* fProducts;
};

struct compMoleculeP
{
    bool operator () (const G4Molecule* molecule1, const G4Molecule* molecule2) const
    {
        return *molecule1<*molecule2;
    }
};

/**
  * G4DNAMolecularReactionTable sorts out the G4DNAMolecularReactionData
  * for bimolecular reaction
  */

class G4DNAMolecularReactionTable : public G4ITReactionTable
{
protected:
    G4DNAMolecularReactionTable();
    static G4DNAMolecularReactionTable* fInstance;

public :
    static G4DNAMolecularReactionTable* GetReactionTable();
    static void DeleteInstance();
    virtual ~G4DNAMolecularReactionTable();

    /**
    * Define a reaction :
    * First argument : reaction rate
    * Second argument : reactant 1
    * Third argument : reactant 2
    * Fourth argument : a std std::vector holding the molecular products
    * if this last argument is NULL then it will be interpreted as
    * a reaction giving no products
    */
    void SetReaction(G4double observedReactionRate,
                     const G4Molecule* reactive1, const G4Molecule* reactive2);

    void SetReaction(G4DNAMolecularReactionData*);

    const G4DNAMolecularReactionData* GetReactionData(const G4Molecule*, const G4Molecule*) const;

    //_________________________________________________________________
    /**
     * Given a molecule's type, it returns with which a reaction is allowed
     */
    const std::vector<const G4Molecule*>* CanReactWith(const G4Molecule* aMolecule) const ;
    const std::map<const G4Molecule*, const G4DNAMolecularReactionData*, compMoleculeP>* GetReativesNData(const G4Molecule* aMolecule) const ;
    const std::vector<const G4DNAMolecularReactionData*>* GetReactionData(const G4Molecule*) const ;

    //_________________________________________________________________
    void PrintTable(G4VDNAReactionModel* = 0);

protected :
    const G4MoleculeHandleManager* fMoleculeHandleManager;
    G4bool fVerbose;

    //_________________________________________________
    typedef std::map<const G4Molecule*,
                        std::map<const G4Molecule*,
                                 const G4DNAMolecularReactionData*,
                                 compMoleculeP>,
                    compMoleculeP > ReactionDataMap;
    typedef std::map<const G4Molecule*,std::vector<const G4Molecule*>,compMoleculeP> ReactivesMV ;
    typedef std::map<const G4Molecule*,std::vector<const G4DNAMolecularReactionData*>,compMoleculeP> ReactionDataMV ;

    ReactionDataMap fReactionData;
    ReactivesMV     fReactivesMV;
    ReactionDataMV  fReactionDataMV;
};
#endif /*G4MolecularReactionTable_HH*/

