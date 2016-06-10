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
// $Id: G4DNAMolecularReactionTable.hh 90769 2015-06-09 10:33:41Z gcosmo $
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
  
  G4DNAMolecularReactionData(G4double reactionRate,
                             const G4String& reactive1,
                             const G4String& reactive2);
  ~G4DNAMolecularReactionData();
  
  const G4Molecule* GetReactive1() const { return fReactive1; }
  const G4Molecule* GetReactive2() const { return fReactive2; }
  
  G4double GetReactionRate() const {return fReactionRate;}
  G4double GetReducedReactionRadius() const {return fReducedReactionRadius;}
  
  //_____________________________________________________
  void SetReactive1(const G4Molecule* reactive) ;
  void SetReactive2(const G4Molecule* reactive) ;
  void SetReactive(const G4Molecule* reactive1, const G4Molecule* reactive2);
  void AddProduct(const G4Molecule* molecule);
  
  void SetReactive1(const G4String& reactive) ;
  void SetReactive2(const G4String& reactive) ;
  void SetReactive(const G4String& reactive1, const G4String& reactive2);
  void AddProduct(const G4String& molecule);
  
  G4int GetNbProducts() const
  {
    if(fProducts) return fProducts->size();
    return 0;
  }
  
  const G4Molecule* GetProduct(G4int i) const
  {
    if(fProducts) return (*fProducts)[i];
    return 0;
  }
  
  protected :
  G4DNAMolecularReactionData();
  const G4Molecule* fReactive1;
  const G4Molecule* fReactive2;
  G4double fReactionRate;
  G4double fReducedReactionRadius;
  
  std::vector<const G4Molecule*>* fProducts;
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
  //    static G4ThreadLocal G4DNAMolecularReactionTable* fInstance;
  
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
  
  const G4DNAMolecularReactionData* GetReactionData(const G4Molecule*,
                                                    const G4Molecule*) const;
  
  //_________________________________________________________________
  /**
   * Given a molecule's type, it returns with which a reaction is allowed
   */
  const std::vector<const G4Molecule*>* CanReactWith(const G4Molecule*) const ;
  const std::map<const G4Molecule*,
                 const G4DNAMolecularReactionData*,
                 compMoleculeP>* GetReativesNData(const G4Molecule*) const ;
  const std::vector<const G4DNAMolecularReactionData*>*
                                      GetReactionData(const G4Molecule*) const ;
  
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
  typedef std::map<const G4Molecule*,
                   std::vector<const G4Molecule*>,
                   compMoleculeP> ReactivesMV ;
  typedef std::map<const G4Molecule*,
                   std::vector<const G4DNAMolecularReactionData*>,
                   compMoleculeP> ReactionDataMV ;
  
  ReactionDataMap fReactionData;
  ReactivesMV     fReactivesMV;
  ReactionDataMV  fReactionDataMV;
};
#endif /*G4MolecularReactionTable_HH*/

