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
// HadrontherapyRBE.hh;
//

#ifndef HadrontherapyRBE_H
#define HadrontherapyRBE_H 1

#include "globals.hh"
#include <vector>
#include <valarray>
#include <map>
#include "G4Pow.hh"

class G4GenericMessenger;

/**
 * @brief Main class of the RBE calculation.
 *
 * The calculation has to be explicitly enabled
 * (use macro command)
 *
 * Available macro commands:
 *
 *  - /rbe/calculation 0/1 : enable or disable RBE calculation
 *  - /rbe/verbose 0/1/2 : level of screen output detail
 *  - /rbe/loadLemTable [path] : read a CSV file with alphas, betas, ...
 *  - /rbe/cellLine [name] : select one of the cell lines from the data file
 *  - /rbe/doseScale [number] : factor to make the survival/RBE calculation correct
 *  - /rbe/accumulate 0/1 : enable or disable data summing over multiple runs
 *  - /rbe/reset : clear accumulated data back to 0.
 */
class HadrontherapyRBE
{
public:
    virtual ~HadrontherapyRBE();

    // Make a new & get instance pointer
    static HadrontherapyRBE* CreateInstance(G4int nX, G4int nY, G4int nZ, G4double massOfVoxel);
    static HadrontherapyRBE* GetInstance();

    // If this is false (set with macro), nothing happens
    G4bool IsCalculationEnabled() const { return fCalculationEnabled; }

    // If this is true, dose (and alpha/beta parameters) are accumulated over multiple runs
    G4bool IsAccumulationEnabled() const { return fAccumulate; }

    // Initialization of data from a CSV file
    void LoadLEMTable(G4String path);

    // Select the cell and update the pointer
    void SetCellLine(G4String name);

    // Calculate alpha and beta for single deposition, {0,0} if not applicable
    std::tuple<G4double, G4double> GetHitAlphaAndBeta(G4double E, G4int Z);

    // Parameter setting
    void SetDoseScale(G4double scale);
    void SetCalculationEnabled(G4bool enabled);
    void SetAccumulationEnabled(G4bool accumulate);

    // Verbosity for output
    void SetVerboseLevel(G4int level) { fVerboseLevel = level; }
    G4int GetVerboseLevel() const { return fVerboseLevel; }
    
    // Alias for matrix type
    using array_type = std::valarray<G4double>;

    // Calculation
    void ComputeAlphaAndBeta();
    void ComputeRBE();

    // Update the class with accumulated data
    // (To be used from HadrontherapyRBEAccumulable)
    void SetAlphaNumerator(const array_type alpha);
    void SetBetaNumerator(const array_type beta);
    void SetEnergyDeposit(const array_type eDep);
    void SetDenominator(const array_type denom);
    
    // Accumulation variants necessary for multi-run sumation
    void AddAlphaNumerator(const array_type alpha);
    void AddBetaNumerator(const array_type beta);
    void AddEnergyDeposit(const array_type eDep);
    void AddDenominator(const array_type denom);
    
    // Clear accumulated data
    void Reset();

    // Output to text files (called at the end of run)
    void StoreAlphaAndBeta();
    void StoreRBE();

    // Information about voxels
    size_t GetNumberOfVoxelsAlongX() const { return fNumberOfVoxelsAlongX; }
    size_t GetNumberOfVoxelsAlongY() const { return fNumberOfVoxelsAlongY; }
    size_t GetNumberOfVoxelsAlongZ() const { return fNumberOfVoxelsAlongZ; }

    // Some basic output to the screen
    void PrintParameters();

protected:
    inline G4int Index(G4int i, G4int j, G4int k) {return (i * fNumberOfVoxelsAlongY + j) * fNumberOfVoxelsAlongZ + k;}

    // Interpolation
    // G4int GetRowVecEnergy();
    // G4bool NearLookup(G4double E, G4double DE);
    // G4bool LinearLookup(G4double E, G4double DE, G4int Z);
    // void interpolation_onE(G4int k,G4int m, G4int indexE, G4double E, G4int Z);
    // G4bool interpolation_onLET1_onLET2_onE(G4int k,G4int m, G4int indexE, G4double E, G4double LET);
    // void InitDynamicVec(std::vector<G4double> &vecEnergy, G4int matrix_energy);

    // Messenger initialization
    void CreateMessenger();

private:
    HadrontherapyRBE(G4int numberOfVoxelX, G4int numberOfVoxelY, G4int numberOfVoxelZ, G4double massOfVoxel);

    G4GenericMessenger* fMessenger;
    G4Pow* g4pow = G4Pow::GetInstance();

    static HadrontherapyRBE* instance;
    G4int fVerboseLevel { 1 };

    // Parameters for calculation
    G4double fAlphaX { 0.0 };
    G4double fBetaX { 0.0 };
    G4double fDoseCut { 0.0 };
    G4double fDoseScale { 1.0 };

    // Output paths (TODO: Change to analysis tools)
    G4String fAlphaBetaPath { "AlphaAndBeta.out" };
    G4String fRBEPath { "RBE.out" };

    // Voxelization
    G4int fNumberOfVoxelsAlongX, fNumberOfVoxelsAlongY, fNumberOfVoxelsAlongZ;
    G4int fNumberOfVoxels;
    G4double fMassOfVoxel;

    G4double* x; // TODO: Currently not used (that much)

    G4bool fCalculationEnabled { false };
    G4bool fAccumulate { false };

    // Matrices to be set when accumulated
    array_type fAlpha;
    array_type fBeta;
    array_type fDose; // Note: this is calculated from energyDeposit, massOfVoxel and doseScale
    
    array_type fAlphaNumerator;
    array_type fBetaNumerator;
    array_type fDenominator;

    // Matrices of calculated values
    array_type fLnS;
    array_type fSurvival;
    array_type fDoseX;
    array_type fRBE;

    // Available tables and associated values.
    using vector_type = std::map<G4int, std::vector<G4double>>;
    std::map<G4String, vector_type> fTablesEnergy;
    // std::map<G4String, vector_type> fTablesLet;
    std::map<G4String, vector_type> fTablesAlpha;
    std::map<G4String, vector_type> fTablesBeta;
    std::map<G4String, G4double> fTablesAlphaX;
    std::map<G4String, G4double> fTablesBetaX;
    std::map<G4String, G4double> fTablesDoseCut;

    // Selected tables and associated values.
    // (changed when the cell line is set)
    G4String fActiveCellLine = "";
    vector_type* fActiveTableEnergy { nullptr };
    // vector_type* fActiveTableLet { nullptr };
    vector_type* fActiveTableAlpha { nullptr };
    vector_type* fActiveTableBeta { nullptr };
    std::map<G4int, G4double> fMaxEnergies;
    std::map<G4int, G4double> fMinEnergies;
    G4int fMinZ;
    G4int fMaxZ;
};
#endif

