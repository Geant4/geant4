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
/// \file InformationKeeper.cc
/// \brief Implementation of the InformationKeeper class

#include "InformationKeeper.hh"

#ifdef USE_MPI
#include "G4MPImanager.hh"
#endif

#include <fstream>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

InformationKeeper* InformationKeeper::fInstance = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

InformationKeeper* InformationKeeper::Instance()
{
    if (fInstance == nullptr) {
        static InformationKeeper parParser;
        fInstance = &parParser;
    }
    return fInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

InformationKeeper::InformationKeeper() 
{
    fPhysOutFolderName = "phys_output";
    fChemInputFolderName = "chem_input";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void InformationKeeper::WritePhysGeo()
{
    //Create a file containing geopath and Physlist in "FS build directory" for chemStage and anlysis
    std::string fname = "imp.info";
    std::ofstream fout(fname);
    fout<<"====> Auto_generated file. Do not delete me !!!\n";
    fout<<"====> The file conveys some information for chem_geo and Analysis modules!!!\n";
    fout<<"_geocellpath "<<fCellDefFilePath<<"\n";
    for (const auto & entry : fVoxelDefFilesList) {
        fout<<"_geovolxelpath "<<std::string(entry)<<"\n";
    }
    fout<<"_physList "<<fPhysDNAName<<"\n";
    fout<<"_numberOfBasepairs "<<fTotalNbBpPlacedInGeo<<"\n";
    fout<<"_numberOfHistones "<<fTotalNbHistonePlacedInGeo<<"\n";
    fout<<"_nucleusVolume "<<fNucleusVolume<<"\n"; // m3
    fout<<"_nucleusMassDensity "<<fNucleusMassDensity<<"\n"; // kg/m3
    fout<<"_nucleusMass "<<fNucleusVolume*fNucleusMassDensity; //kg
    fout.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......