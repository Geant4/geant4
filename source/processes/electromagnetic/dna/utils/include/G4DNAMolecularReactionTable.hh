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
// $Id: G4DNAMolecularReactionTable.hh 100802 2016-11-02 14:55:27Z gcosmo $
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


#ifndef G4MolecularReactionTable_h
#define G4MolecularReactionTable_h 1

#include "G4ITReactionTable.hh"
#include "G4MolecularConfiguration.hh"
#include <vector>
#include <map>
#include <functional>
#include "G4ReferenceCast.hh"

class G4VDNAReactionModel ;
class G4DNAMolecularReactionTable;
class G4ReactionTableMessenger;

/**
 * G4DNAMolecularReactionData contains the information
 * relative to a given reaction (eg : °OH + °OH -> H2O2)
 */

class G4DNAMolecularReactionData
{
  public :
  //----------------------------------------------------------------------------

  G4DNAMolecularReactionData(G4double reactionRate,
                             G4MolecularConfiguration* reactive1,
                             G4MolecularConfiguration* reactive2);
  
  G4DNAMolecularReactionData(G4double reactionRate,
                             const G4String& reactive1,
                             const G4String& reactive2);
  ~G4DNAMolecularReactionData();
  

  //----------------------------------------------------------------------------
  inline int GetReactionID() const { return fReactionID; }
  inline void SetReactionID(int ID) { fReactionID = ID; }

  //----------------------------------------------------------------------------
  inline std::pair<G4MolecularConfiguration*, G4MolecularConfiguration*>
  GetReactants()
  {
    return std::make_pair(fReactant1, fReactant2);
  }

  inline G4MolecularConfiguration* GetReactant1() const
  {
    return fReactant1;
  }
  inline G4MolecularConfiguration* GetReactant2() const
  {
    return fReactant2;
  }

  inline void SetObservedReactionRateConstant(G4double rate)
  {
    fObservedReactionRate = rate;
  }

  inline G4double GetObservedReactionRateConstant() const
  {
    return fObservedReactionRate;
  }

  inline G4double GetEffectiveReactionRadius() const
  {
    return fEffectiveReactionRadius;
  }
  
  inline void SetEffectiveReactionRadius(G4double radius)
  {
    fEffectiveReactionRadius = radius;
  }

  //_____________________________________________________

  void SetReactant1(G4MolecularConfiguration* reactive) ;
  void SetReactant2(G4MolecularConfiguration* reactive) ;
  
  void SetReactants(G4MolecularConfiguration* reactive1,
                   G4MolecularConfiguration* reactive2);

  void AddProduct(G4MolecularConfiguration* molecule);

  void SetReactant1(const G4String& reactive) ;
  void SetReactant2(const G4String& reactive) ;
  void SetReactants(const G4String& reactive1, const G4String& reactive2);
  void AddProduct(const G4String& molecule);
  
  inline G4int GetNbProducts() const
  {
    if(fProducts) return fProducts->size();
    return 0;
  }
  
  inline G4MolecularConfiguration* GetProduct(G4int i) const
  {
    if(fProducts) return (*fProducts)[i];
    return 0;
  }
  
  inline const std::vector<G4MolecularConfiguration*>* GetProducts() const
  {
    return fProducts;
  }

  inline void RemoveProducts()
  {
    if(fProducts)
    {
      fProducts->clear();
      delete fProducts;
    }
  }

  //----------------------------------------------------------------------------
  // Temperature scaling
  typedef std::function<double(double)> RateParam;

  static double PolynomialParam(double temp_K, std::vector<double> P);
  static double ArrehniusParam(double temp_K, std::vector<double> P);
  static double ScaledParameterization(double temp_K,
                                       double temp_init,
                                       double rateCste_init);

  void SetPolynomialParameterization(const std::vector<double>& P);

  void SetArrehniusParameterization(double A0, double E_R);
  void SetScaledParameterization(double temperature_K,
      double rateCste);

  void ScaleForNewTemperature(double temp_K);

  protected :
  G4DNAMolecularReactionData();
  G4MolecularConfiguration* fReactant1;
  G4MolecularConfiguration* fReactant2;
  G4double fObservedReactionRate;
  G4double fEffectiveReactionRadius;
  
  std::vector<G4MolecularConfiguration*>* fProducts;
  // G4DNAReactionType fReactionType;
  RateParam fRateParam;
  int fReactionID;
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
  static G4DNAMolecularReactionTable* Instance();
  static void DeleteInstance();
  virtual ~G4DNAMolecularReactionTable();
  
  /**
   * Define a reaction :
   * First argument : reaction rate
   * Second argument : reactant 1
   * Third argument : reactant 2
   * Fourth argument : a std::vector holding the molecular products
   * if this last argument is NULL then it will be interpreted as
   * a reaction giving no products
   */
  void SetReaction(G4double observedReactionRate,
                   G4MolecularConfiguration* reactive1,
                   G4MolecularConfiguration* reactive2);
  
  void SetReaction(G4DNAMolecularReactionData*);
  
  const G4DNAMolecularReactionData* GetReactionData(G4MolecularConfiguration*,
                                                    G4MolecularConfiguration*) const;
  
  const G4DNAMolecularReactionData* GetReactionData(const G4String&,
                                                    const G4String&) const;

  const G4DNAMolecularReactionData* GetReaction(int reactionID) const;

  size_t GetNReactions() const
  { return fVectorOfReactionData.size(); }

  //_________________________________________________________________
  /**
   * Given a molecule's type, it returns with which a reaction is allowed
   */
  const std::vector<G4MolecularConfiguration*>*
  CanReactWith(G4MolecularConfiguration*) const ;

  const std::map<G4MolecularConfiguration*, const G4DNAMolecularReactionData*>*
  GetReativesNData(G4MolecularConfiguration*) const;

  const std::vector<const G4DNAMolecularReactionData*>*
  GetReactionData(G4MolecularConfiguration*) const;
  
  inline const std::map<G4MolecularConfiguration*,
        std::map<G4MolecularConfiguration*,
            const G4DNAMolecularReactionData*> >&
  GetAllReactionData()
  {
    return fReactionData;
  }

  inline const std::vector<const G4DNAMolecularReactionData*>&
  GetVectorOfReactionData()
  {
    return fVectorOfReactionData;
  }

  void ScaleReactionRateForNewTemperature(double temp_K);

  //_________________________________________________________________
  void PrintTable(G4VDNAReactionModel* = 0);
  
  protected :
  G4bool fVerbose;
  
  //_________________________________________________
  typedef std::map<G4MolecularConfiguration*,
      std::map<G4MolecularConfiguration*,
          const G4DNAMolecularReactionData*> > ReactionDataMap;
  typedef std::map<G4MolecularConfiguration*,
      std::vector<G4MolecularConfiguration*> > ReactivesMV;
  typedef std::map<G4MolecularConfiguration*,
      std::vector<const G4DNAMolecularReactionData*> > ReactionDataMV;
  
  ReactionDataMap fReactionData;
  ReactivesMV     fReactantsMV;
  ReactionDataMV  fReactionDataMV;
  std::vector<const G4DNAMolecularReactionData*> fVectorOfReactionData;
  G4ReactionTableMessenger* fpMessenger;
};
#endif /*G4MolecularReactionTable_HH*/

