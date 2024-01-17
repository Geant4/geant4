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
/// \file ParametersParser.cc
/// \brief Implementation of the ParametersParser class

#include "ParametersParser.hh"

#include <fstream>
#include <iostream>
#include <sstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

ParametersParser* ParametersParser::fInstance = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

ParametersParser* ParametersParser::Instance()
{
    if (fInstance == nullptr) {
        static ParametersParser parParser;
        fInstance = &parParser;
    }
    return fInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

ParametersParser::ParametersParser()
{
    fSDDfileName = "SDDformat_"+fOutputName;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ParametersParser::LoadParameters(const std::string &fileName)
{
    std::ifstream file;
    file.open(fileName.c_str());
    if (!file.is_open()) {
        std::cout<<"ParametersParser::LoadParameters Error in openning file!!!"<<std::endl;
    } else {
        std::string line;
        while(std::getline(file, line))
        {
            std::istringstream iss(line);
            std::string flag;
            iss >> flag;
            std::string tvalue;
            iss >> tvalue;
            if (flag == "/ana/thresholdFordirectSBSelection") fThresholdE = (tvalue);
            if (flag == "/ana/probForIndirectSBSelection") fProbabilityForIndirectSB = (tvalue);
            if (flag == "/ana/BpForDSB") BpForDSB = std::stoi(tvalue);
            if (flag == "/ana/TLK/lambda1") TLKLambda1 = (tvalue);
            if (flag == "/ana/TLK/lambda2") TLKLambda2 = (tvalue);
            if (flag == "/ana/TLK/beta1") TLKBeta1 = (tvalue);
            if (flag == "/ana/TLK/beta2") TLKBeta2 = (tvalue);
            if (flag == "/ana/TLK/eta") TLKEta = (tvalue);
            if (flag == "/ana/TLK/doseMax") TLKdoseMax = (tvalue);
            if (flag == "/ana/TLK/deltaDose") TLKdeltaDose = (tvalue);
            if (flag == "/ana/LEMIV/loopLength") LEMIVLoopLength = tvalue;
            if (flag == "/ana/LEMIV/Ni") LEMIVNi = tvalue;
            if (flag == "/ana/LEMIV/Nc") LEMIVNc = tvalue;
            if (flag == "/ana/LEMIV/NDSB") LEMIVNDSB = tvalue;
            if (flag == "/ana/LEMIV/Funrej") LEMIVFunrej = tvalue;
            if (flag == "/ana/LEMIV/Tfast") LEMIVTfast = tvalue;
            if (flag == "/ana/LEMIV/Tslow") LEMIVTslow = tvalue;
            if (flag == "/ana/LEMIV/timeMax") LEMIVtimeMax = tvalue;
            if (flag == "/ana/LEMIV/deltaTime") LEMIVdeltaTime = tvalue;
            if (flag == "/ana/BELOV/Nirrep") BELOVNirrep = tvalue;
            if (flag == "/ana/BELOV/Dz") BELOVDz = tvalue;

            if (flag == "/ana/TLK/used") useTLK = tvalue;
            if (flag == "/ana/LEMIV/used") useLEMIV = tvalue;
            if (flag == "/ana/BELOV/used") useBELOV = tvalue;

            if (flag == "/ana/ouputName") fOutputName = tvalue;
            if (flag == "/ana/cellNucleusName") fCellNucleusName = tvalue;
            if (flag == "/ana/loadDamagesFromSDD") {
                fSDDfileName = tvalue;
                fLoadDamagesFromSDD = true;
            }
            if (flag == "/ana/unitOfNormalization") fUnitOfNormalization = std::stoi(tvalue);

            if (flag == "/ana/skipIndirectDamages") fSkipScanningIndirectDamage = true;

            if (flag == "/gps/particle") fParticleName = (tvalue);
            if (flag == "/gps/energy") {
                fParticleEnergy = std::stof(tvalue);
                iss >> tvalue;
                fEnergyUnit = tvalue;
            }
            if (flag == "/run/beamOn") fNumberOfParticles = std::stoi(tvalue);
            if (flag == "/scheduler/endTime") {
                fEndTimeForChemReactions = tvalue;
                iss >> tvalue;
                fEndTimeForChemReactions += tvalue;
            }
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

bool ParametersParser::UseTLK()
{
    bool used = false;
    if (useTLK == "true" || useTLK == "TRUE" || useTLK == "True" || useTLK == "1" || 
        useTLK == "yes" || useTLK == "YES" || useTLK == "Yes") used = true;
    return used;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

bool ParametersParser::UseLEMIV()
{
    bool used = false;
    if (useLEMIV == "true" || useLEMIV == "TRUE" || useLEMIV== "True" || useLEMIV == "1" || 
        useLEMIV == "yes" || useLEMIV == "YES"|| useLEMIV == "Yes") used = true;
    return used;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

bool ParametersParser::UseBelov()
{
    bool used = false;
    if (useBELOV == "true" || useBELOV == "TRUE" || useBELOV == "True" || useBELOV == "1" || 
        useBELOV == "yes" || useBELOV == "YES"|| useBELOV == "Yes") used = true;
    return used;
}