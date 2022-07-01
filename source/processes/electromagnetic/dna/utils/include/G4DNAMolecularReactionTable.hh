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

#include "G4ITReactionTable.hh"
#include "G4MolecularConfiguration.hh"
#include "G4ReferenceCast.hh"
#include "G4VDNAMolecularGeometry.hh"
#include <vector>
#include <map>
#include <functional>
#include <memory>

class G4VDNAReactionModel;
class G4DNAMolecularReactionTable;
class G4ReactionTableMessenger;

/**
 * G4DNAMolecularReactionData contains the information
 * relative to a given reaction (eg : °OH + °OH -> H2O2)
 */
class G4DNAMolecularReactionData
{
public:
    //----------------------------------------------------------------------------

    G4DNAMolecularReactionData(G4double reactionRate,
                               const G4MolecularConfiguration* reactive1,
                               const G4MolecularConfiguration* reactive2);

    G4DNAMolecularReactionData(G4double reactionRate,
                               const G4String& reactive1,
                               const G4String& reactive2);
    ~G4DNAMolecularReactionData();

    using Reactant = const G4MolecularConfiguration;
    using ReactantPair = std::pair<Reactant*, Reactant*>;
    using ReactionProducts = std::vector<Reactant*>;

    int GetReactionID() const;
    void SetReactionID(int ID);

    ReactantPair GetReactants();

    Reactant* GetReactant1() const;
    Reactant* GetReactant2() const;

    void SetObservedReactionRateConstant(G4double rate);
    G4double GetObservedReactionRateConstant() const;
    G4double GetActivationRateConstant() const;
    G4double GetDiffusionRateConstant() const;

    void SetReactionRadius(G4double radius);
    G4double GetReactionRadius() const;

    void SetEffectiveReactionRadius(G4double radius);
    G4double GetEffectiveReactionRadius() const;
    G4double GetOnsagerRadius() const;

    void SetProbability(G4double prob);
    G4double GetProbability() const;

    void SetReactionType(G4int type);
    G4int GetReactionType() const;

    void SetReactant1(Reactant* reactive);
    void SetReactant2(Reactant* reactive);

    void SetReactants(Reactant* reactive1,
                      Reactant* reactive2);

    void AddProduct(Reactant* molecule);

    void SetReactant1(const G4String& reactive);
    void SetReactant2(const G4String& reactive);
    void SetReactants(const G4String& reactive1, const G4String& reactive2);
    void AddProduct(const G4String& molecule);

    G4int GetNbProducts() const;
    Reactant* GetProduct(G4int i) const;

    const ReactionProducts* GetProducts() const;
    void RemoveProducts();

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

    void ComputeEffectiveRadius();

protected:
    G4DNAMolecularReactionData();
    Reactant* fpReactant1;
    Reactant* fpReactant2;

    G4double fObservedReactionRate;
    G4double fActivationRate;
    G4double fDiffusionRate;

    G4double fOnsagerRadius;

    G4double fReactionRadius;
    G4double fEffectiveReactionRadius;

    G4double fProbability;
    G4int fType;

    ReactionProducts fProducts;
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
    static G4DNAMolecularReactionTable* fpInstance;

public:
    static G4DNAMolecularReactionTable* GetReactionTable();
    static G4DNAMolecularReactionTable* Instance();
    static void DeleteInstance();
    virtual ~G4DNAMolecularReactionTable();

    using Reactant = const G4MolecularConfiguration;
    using Data = const G4DNAMolecularReactionData;
    using ReactantList = std::vector<Reactant*>;
    using DataList = std::vector<Data*>;
    using SpecificDataList = std::map<Reactant*, Data*>;

    using ReactionDataMap = std::map<Reactant*, SpecificDataList>;
    using ReactivesMV = std::map<Reactant*, ReactantList>;
    using ReactionDataMV = std::map<Reactant*, DataList>;

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
                     Reactant* reactive1,
                     Reactant* reactive2);

    void SetReaction(G4DNAMolecularReactionData*);

    void SetGeometry(G4VDNAMolecularGeometry* geometry){fGeometry = geometry;};
    G4VDNAMolecularGeometry* GetGeometry() const;

    Data* GetReactionData(Reactant*, Reactant*) const;

    Data* GetReactionData(const G4String&, const G4String&) const;

    Data* GetReaction(int reactionID) const;

    size_t GetNReactions() const;

    //_________________________________________________________________
    /**
     * Given a molecule's type, it returns with which a reaction is allowed
     */
    const ReactantList* CanReactWith(Reactant*) const;

    const SpecificDataList* GetReativesNData(const G4MolecularConfiguration*) const;

    const DataList* GetReactionData(const G4MolecularConfiguration*) const;

    const ReactionDataMap& GetAllReactionData();

    DataList GetVectorOfReactionData();

    void ScaleReactionRateForNewTemperature(double temp_K);

    //_________________________________________________________________
    void PrintTable(G4VDNAReactionModel* = 0);

    void Reset();

protected:
    G4bool fVerbose;

    G4VDNAMolecularGeometry* fGeometry;
    ReactionDataMap fReactionData;
    ReactivesMV     fReactantsMV;
    ReactionDataMV  fReactionDataMV;
    std::vector<std::unique_ptr<Data>> fVectorOfReactionData;
    std::unique_ptr<G4ReactionTableMessenger> fpMessenger;
};
