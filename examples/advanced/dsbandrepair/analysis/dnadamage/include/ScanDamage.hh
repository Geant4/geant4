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
//
/// \file ScanDamage.hh
/// \brief Definition of the ScanDamage class

#ifndef ScanDamage_h
#define ScanDamage_h

#include <map>
#include <vector>
#include <string>
#include <tuple>
#include <set>
#include "Damage.hh"
#include <filesystem>
namespace fs = std::filesystem;
class TFile;
using ullint = unsigned long long int;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

struct VoxelData
{
    VoxelData(int chromo, int domain, ullint firstBpCN)
    {
        fChromosome = chromo;
        fDomain = domain;
        fFirstBpCopyNum = firstBpCN;
    }

    ~VoxelData() {}

    int fChromosome{0};
    int fDomain{0};
    ullint fFirstBpCopyNum{0};
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

typedef std::vector<std::vector<ullint> > Table;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

class Trier {
public:
    bool operator()(const std::vector<ullint>& a, const std::vector<ullint>& b)
    {
        bool bb = false;
        if(a[0] < b[0]) bb = true;
        return bb;
    }
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

class ScanDamage
{
public:
    ScanDamage();
    ~ScanDamage() = default;  
    std::map<unsigned int,std::map<unsigned int,std::vector<Damage> > > ExtractDamage();
    void SetThresholdEnergy(double e) {fThresholdEnergy = e;}
    void SetProbabilityForIndirectSBSelection(double p) {fProbabilityForIndirectSB = p;}
    double GetThresholdEnergy() {return fThresholdEnergy;}
    double GetProbabilityForIndirectSBSelection() {return fProbabilityForIndirectSB;}
    std::map<int, Table> GetMergedSBData() {return fMergedTables;}
    double GetEdepSumInNucleus() {return fEdepSumInNucleus;} //eV
    double GetTotalNbBpPlacedInGeo() {return fTotalNbBpPlacedInGeo;}
    double GetTotalNbHistonePlacedInGeo() {return fTotalNbHistonePlacedInGeo;}
    double GetNucleusVolume() {return fNucleusVolume;} 
    double GetNucleusMassDensity() {return fNucleusMassDensity;}
    double GetNucleusMass() {return fNucleusMass;}
    std::map<int,ullint> GetChromosomeBpSizesMap() {return fChromosomeBpMap;}
    void SkipScanningIndirectDamage() {fSkipScanningIndirectDamage = true;}
    bool SkippedScanningIndirectDamage() {return fSkipScanningIndirectDamage;}
private:
    void ScanDamageFromPhys();
    void ScanDamageFromChem();
    void RetrieveVoxelBp();
    void FillVoxelData();
    void AnaPhysRootFile(const std::string fileName);
    void AnaChemRootFile(fs::directory_entry entry);
    void AnaPhysRootTree1(TFile*);
    void AnaPhysRootTree2(TFile*);
    void SortPhysTableWithSelection();
    void SortChemTableWithSelection();
    void ReadCellandVoxelDefFilePaths();
    void MergeDamageFromPhysChem();
    std::tuple<unsigned int, unsigned int> GetEventNberAndVoxelNberFromChemRoot(const std::string fileNam);
    double fThresholdEnergy{17.5};//eV
    std::string fCellDefFilePath{""};
    std::set<std::string> fVoxelDefFilesList;
    std::map<std::string, int> fBpPerVoxel;
    std::vector<VoxelData> fVoxels;
    std::map<int, Table> fphysTables, fphysSlectedTables, fchemTables, fchemSlectedTables,fMergedTables;
    std::map<unsigned int,std::map<unsigned int,std::vector<Damage> > > fDamage;
    double fEdepSumInNucleus{0}; //eV
    int corruptedFiles = 0;
    double fTotalNbBpPlacedInGeo{0};
    double fTotalNbHistonePlacedInGeo{0};
    double fNucleusVolume{0};
    double fNucleusMassDensity{0};
    double fNucleusMass{0};
    double fProbabilityForIndirectSB{0.4};
    std::map<int,ullint> fChromosomeBpMap; //Store number of Bp in each Chomosomes;
    bool fSkipScanningIndirectDamage{false};
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

#endif


