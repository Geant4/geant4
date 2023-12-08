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
/// \file ChemGeoImport.cc
/// \brief Implementation of the ChemGeoImport class

#include "ChemGeoImport.hh"
#include "G4Filesystem.hh"

ChemGeoImport::ChemGeoImport()
{
    GetVoxelDefFilePathList();
    fpGun = new UserMoleculeGun();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ChemGeoImport::~ChemGeoImport()
{
    if(fpGun)
        delete fpGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChemGeoImport::InsertMoleculeInWorld()
{
    // The idea is to add all the molecules specified in the input files

    if(fIsParsed)
    {
        // Create the molecules

        // Loop on all the parsed molecules
        for(G4int i=0, ie=fMolecules.size(); i<ie; i++)
        {
            // Retrieve general molecule informations
            //
            G4String name = fMolecules[i].fName;

            G4ThreeVector moleculePosition = fMolecules[i].fPosition;
            G4int copyNum = fMolecules[i].fCopyNumber;

            G4int strand = fMolecules[i].fStrand;
            ChemMolecule Bmolecule(name, copyNum, moleculePosition, strand, -1, -1, -1);

            if(name=="phosphate1") 
            {
                name="Phosphate";
                Bmolecule.fName=name;
            }
            else if(name=="phosphate2") 
            {
                name="Phosphate";
                Bmolecule.fName=name;
            }
            else if(name=="deoxyribose1") 
            {
                name="Deoxyribose";
                Bmolecule.fName=name;
            }
            else if(name=="deoxyribose2") 
            {
                name="Deoxyribose";
                Bmolecule.fName=name;
            }
            else if(name=="base_adenine")
            {
                name="Adenine";
                Bmolecule.fName="base";
            }
            else if(name=="base_guanine") 
            {
                name="Guanine";
                Bmolecule.fName="base";
            }
            else if(name=="base_thymine")
            {
                name="Thymine";
                Bmolecule.fName="base";
            }
            else if(name=="base_cytosine") 
            {
                name="Cytosine";
                Bmolecule.fName="base";
            }
            else if(name=="histone") 
            {
                name="Histone";
                Bmolecule.fName=name;
            }
            else if(name=="solvatedElectron") 
            {
                name=G4Electron_aq::Definition()->GetName();
            }
            else if(name=="water") 
            {
                name=G4H2O::Definition()->GetName();
            }
            else
            {
                G4String msg = 
                "The name "+ name+" is not specified in the listed chemical molecules";
                G4Exception("ChemGeoImport::BuildGeometry", "", FatalException, msg);
            }
            
            

            // Check if the molecule is on the "remove list"
            G4bool toBeRemoved = IsMoleculeInTheRemoveTable(Bmolecule);
            if(!toBeRemoved)
            {
                // Molecule is not in the "remove list" and we can add it to the simulation
                // Check the molecule to be added is not a water molecule (special case)
                if(name != G4H2O::Definition()->GetName() )
                    fpGun->AddMolecule(name, moleculePosition, 1.e-12*s, copyNum, strand);
                else // Water molecule case
                    fpGun->AddWaterMolecule(moleculePosition, fMolecules.at(i).fTrackId, 
                        ElectronicModification(fMolecules.at(i).fState), 
                        fMolecules.at(i).fElectronicLevel);
            }

        }
        G4DNAChemistryManager::Instance()->SetGun(fpGun);
    }
    else
    {
        G4String msg = 
        "ChemGeoImport::BuildGeometry: The parse method needs to be called first.";
        G4Exception("ChemGeoImport::BuildGeometry", "", FatalException, msg);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChemGeoImport::Reset()
{
    // Clear the containers
    if(fpGun){
        delete fpGun;
        fpGun = new UserMoleculeGun();
    }
    fMolecules.clear();
    fIsParsed = false;
    fFactor = 1.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChemGeoImport::ParseFiles(const G4String& chemInputFile)
{
    G4fs::path aP{std::string(chemInputFile)};
    if (G4fs::exists(aP)) {
        Reset();
        ParseChemInputFile(chemInputFile);
        auto geoPathFileName = GetVoxelDefFilePath(fGeoNameFromChemInput);

        ParseGeoFile(geoPathFileName);
        fIsParsed = true;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChemGeoImport::ParseChemInputFile(const G4String& fileName)
{
    // Setup the input stream
    std::ifstream file;
    file.open(fileName.c_str() );

    if(!file.good() )
    {
        // Geant4 exception
        G4String msg = fileName+" could not be opened";
        G4Exception("ChemGeoImport::ParseChemInputFile", "", FatalException, msg);
    }

    // Define the line string variable
    G4String line;

    // Read the file line per line
    while(std::getline(file, line) )
    {
        // Check the line to determine if it is empty
        if(line.empty() )
            continue; // skip the line if it is empty

        // Data string stream
        std::istringstream issLine(line);

        // String to determine the first letter/word
        G4String firstItem;

        // Put the first letter/word within the string
        issLine >> firstItem;

        // Check first letter to determine if the line is data or comment
        if(firstItem=="#")
            continue; // skip the line if it is comment

        else if(firstItem=="_input")
        {
            G4int type(-1), state(-1), electronicLevel(-1), parentTrackId(-1);
            G4double x, y, z;
            issLine >> type >> state >> electronicLevel;
            issLine >> x >> y >> z;
            issLine >> parentTrackId;

            x *= fFactor*nm;
            y *= fFactor*nm;
            z *= fFactor*nm;

            G4String name;
            if(type==1)
                name="water";
            else if(type==2)
                name="solvatedElectron";
            else
            {
                G4ExceptionDescription description;
                description <<  "The type " << type <<" is not recognized";
                G4Exception("ChemGeoImport::ParseFile", "Fatal", FatalException, description, "");
            }

            ChemMolecule molecule(name, -1, G4ThreeVector(x,y,z), -1, 
                state, electronicLevel, parentTrackId);

            fMolecules.push_back(molecule);
        }

        else if(firstItem=="_remove")
        {
            G4String name;
            issLine >> name;

            G4int copyNumber;
            issLine >> copyNumber;

            G4int strand;
            issLine >> strand;

            fToBeRemovedMol.push_back(ChemMolecule(name,copyNumber,G4ThreeVector(),strand,-1,-1,-1));
        }

        else if(firstItem=="_eventNum")
        {
            // Nothing
        }

        else if(firstItem=="_voxelType")
        {
            issLine >> fGeoNameFromChemInput;
        }

        else if(firstItem=="_voxelCopyNumber")
        {
            // Nothing
        }

        else if(firstItem=="_Version")
        {
            // Nothing
        }

        else
        {
            // Geant4 exception
            G4String msg = 
            firstItem+" is not defined in the parser. Check the input file: "+fileName+".";
            G4Exception("ChemGeoImport::ParseChemInputFile", "Geo_WrongParse",FatalException, msg);
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChemGeoImport::ParseGeoFile(const G4String& fileName)
{
    // Setup the input stream
    std::ifstream file(fileName.c_str());

    // Check if the file was correctly opened
    if(!file.is_open() )
    {
        // Geant4 exception
        G4String msg = fileName+" could not be opened";
        G4Exception("ChemGeoImport::ParseGeoFile", "", FatalException, msg);
    }

    // Define the line string variable
    G4String line;

    // Read the file line per line
    while(std::getline(file, line) )
    {
        // Check the line to determine if it is empty
        if(line.empty() )
            continue; // skip the line if it is empty

        // Data string stream
        std::istringstream issLine(line);

        // String to determine the first letter/word
        G4String firstItem;

        // Put the first letter/word within the string
        issLine >> firstItem;

        // Check first letter to determine if the line is data or comment
        if(firstItem=="#")
            continue; // skip the line if it is comment

        // Use the file
        else if(firstItem=="_Name")
        {
            G4String name;
            issLine >> name;
        }
        else if(firstItem=="_Size")
        {
            G4double size;
            issLine >> size;
            size *= fFactor*nm;

            fSize = size;
        }
        else if(firstItem=="_Number")
        {
            // Nothing
        }
        else if(firstItem=="_Radius")
        {
            // Nothing
        }
        else if(firstItem=="_Version")
        {
            // Nothing
        }
        else if(firstItem=="_pl")
        {
            G4String name;
            issLine >> name;

            G4String material;
            issLine >> material;

            G4int strand;
            issLine >> strand;

            G4int copyNumber;
            issLine >> copyNumber;

            G4double x;
            issLine >> x;
            x *= fFactor*nm;

            G4double y;
            issLine >> y;
            y *= fFactor*nm;

            G4double z;
            issLine >> z;
            z *= fFactor*nm;

            ChemMolecule molecule(name, copyNumber, G4ThreeVector(x, y, z), strand, -1, -1, -1);

            fMolecules.push_back(molecule);
        }

        else
        {
            // Geant4 exception
            G4String msg = 
            firstItem+" is not defined in the parser. Check the input file: "+fileName+".";
            G4Exception("ChemGeoImport::ParseGeoFile", "Geo_WrongParse", FatalException, msg);
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool ChemGeoImport::IsMoleculeInTheRemoveTable(const ChemMolecule& molecule)
{
    if(std::find(fToBeRemovedMol.begin(),fToBeRemovedMol.end(),molecule) != fToBeRemovedMol.end())
        return true;
    else
        return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String ChemGeoImport::GetVoxelDefFilePath(G4String bareName)
{
    G4String strRes = "";
    for (auto const &entry : fVoxelDefFilesList) {
        G4fs::path voxelP{std::string(entry)};
        if (voxelP.stem().string() == bareName) {
            strRes = entry;
        }
    }
    return strRes;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChemGeoImport::GetVoxelDefFilePathList()
{
    G4fs::path thisP = G4fs::current_path();
    G4bool doesWantedFileExist = false;
    for (const auto &entry : G4fs::directory_iterator(thisP)){
        if (entry.path().filename() == "imp.info") {
            std::ifstream file(entry.path().c_str());
            if(!file.good() ){
                G4String msg = 
                "File imp.info is broken. Check its content or try to rerun the PhysicalStage?";
                G4Exception("ChemGeoImport::GetVoxelDefFilePathList()", "", FatalException, msg);
            }
            doesWantedFileExist = true;
            G4String line;
            while(std::getline(file, line) ){
                std::istringstream iss(line);
                G4String flag;
                G4String voxelDefFile;
                iss >> flag;
                if ( flag == "_geovolxelpath") {
                    iss >> voxelDefFile;
                    fVoxelDefFilesList.insert(voxelDefFile);
                }
            }
            file.close();
        }
    }

    if (!doesWantedFileExist) {
        G4String msg = "File imp.info does not exist. Did you run the Physical Stage?";
        G4Exception("ChemGeoImport::GetVoxelDefFilePathList()", "", FatalException, msg);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......