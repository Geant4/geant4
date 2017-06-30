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
// Authors: S. Meylan and C. Villagrasa (IRSN, France)
// Models come from
// M. Bug et al, Rad. Phys and Chem. 130, 459-479 (2017)

#include "G4DNAPTBIonisationStructure.hh"
#include "G4SystemOfUnits.hh"

G4DNAPTBIonisationStructure::G4DNAPTBIonisationStructure()
{
    energyConstant["G4_WATER"].push_back(10.79*eV);
    energyConstant["G4_WATER"].push_back(13.39*eV);
    energyConstant["G4_WATER"].push_back(16.05*eV);
    energyConstant["G4_WATER"].push_back(32.30*eV);
    energyConstant["G4_WATER"].push_back(539.0*eV);

    energyConstant["THF"].push_back(9.74*eV);
    energyConstant["THF"].push_back(12.31*eV);
    energyConstant["THF"].push_back(12.99*eV);
    energyConstant["THF"].push_back(13.57*eV);
    energyConstant["THF"].push_back(13.60*eV);
    energyConstant["THF"].push_back(15.11*eV);
    energyConstant["THF"].push_back(15.97*eV);
    energyConstant["THF"].push_back(16.28*eV);
    energyConstant["THF"].push_back(18.19*eV);
    energyConstant["THF"].push_back(18.69*eV);
    energyConstant["THF"].push_back(22.14*eV);
    energyConstant["THF"].push_back(22.25*eV);
    energyConstant["THF"].push_back(27.21*eV);
    energyConstant["THF"].push_back(28.97*eV);
    energyConstant["THF"].push_back(36.97*eV);
    energyConstant["THF"].push_back(305.07*eV);
    energyConstant["THF"].push_back(305.08*eV);
    energyConstant["THF"].push_back(306.17*eV);
    energyConstant["THF"].push_back(306.17*eV);
    energyConstant["THF"].push_back(557.94*eV);

    energyConstant["PY"].push_back(9.73*eV);
    energyConstant["PY"].push_back(10.96*eV);
    energyConstant["PY"].push_back(11.54*eV);
    energyConstant["PY"].push_back(12.58*eV);
    energyConstant["PY"].push_back(15.96*eV);
    energyConstant["PY"].push_back(16.27*eV);
    energyConstant["PY"].push_back(16.53*eV);
    energyConstant["PY"].push_back(17.98*eV);
    energyConstant["PY"].push_back(19.37*eV);
    energyConstant["PY"].push_back(20.52*eV);
    energyConstant["PY"].push_back(24.55*eV);
    energyConstant["PY"].push_back(24.64*eV);
    energyConstant["PY"].push_back(29.75*eV);
    energyConstant["PY"].push_back(33.02*eV);
    energyConstant["PY"].push_back(36.57*eV);
    energyConstant["PY"].push_back(305.92*eV);
    energyConstant["PY"].push_back(307.09*eV);
    energyConstant["PY"].push_back(307.09*eV);
    energyConstant["PY"].push_back(307.52*eV);
    energyConstant["PY"].push_back(423.44*eV);
    energyConstant["PY"].push_back(423.44*eV);

    energyConstant["PU"].push_back(9.58*eV);
    energyConstant["PU"].push_back(10.57*eV);
    energyConstant["PU"].push_back(10.97*eV);
    energyConstant["PU"].push_back(12.22*eV);
    energyConstant["PU"].push_back(12.92*eV);
    energyConstant["PU"].push_back(13.44*eV);
    energyConstant["PU"].push_back(15.05*eV);
    energyConstant["PU"].push_back(16.56*eV);
    energyConstant["PU"].push_back(17.18*eV);
    energyConstant["PU"].push_back(17.88*eV);
    energyConstant["PU"].push_back(17.90*eV);
    energyConstant["PU"].push_back(19.11*eV);
    energyConstant["PU"].push_back(20.09*eV);
    energyConstant["PU"].push_back(21.70*eV);
    energyConstant["PU"].push_back(23.52*eV);
    energyConstant["PU"].push_back(24.35*eV);
    energyConstant["PU"].push_back(25.41*eV);
    energyConstant["PU"].push_back(29.34*eV);
    energyConstant["PU"].push_back(32.44*eV);
    energyConstant["PU"].push_back(33.67*eV);
    energyConstant["PU"].push_back(36.26*eV);
    energyConstant["PU"].push_back(38.22*eV);
    energyConstant["PU"].push_back(306.53*eV);
    energyConstant["PU"].push_back(307.19*eV);
    energyConstant["PU"].push_back(307.64*eV);
    energyConstant["PU"].push_back(308.14*eV);
    energyConstant["PU"].push_back(308.17*eV);
    energyConstant["PU"].push_back(423.31*eV);
    energyConstant["PU"].push_back(423.43*eV);
    energyConstant["PU"].push_back(423.64*eV);
    energyConstant["PU"].push_back(423.98*eV);

    energyConstant["TMP"].push_back(10.81*eV);
    energyConstant["TMP"].push_back(10.81*eV);
    energyConstant["TMP"].push_back(12.90*eV);
    energyConstant["TMP"].push_back(13.32*eV);
    energyConstant["TMP"].push_back(13.32*eV);
    energyConstant["TMP"].push_back(13.59*eV);
    energyConstant["TMP"].push_back(14.33*eV);
    energyConstant["TMP"].push_back(14.33*eV);
    energyConstant["TMP"].push_back(15.90*eV);
    energyConstant["TMP"].push_back(17.09*eV);
    energyConstant["TMP"].push_back(17.09*eV);
    energyConstant["TMP"].push_back(17.13*eV);
    energyConstant["TMP"].push_back(17.85*eV);
    energyConstant["TMP"].push_back(17.85*eV);
    energyConstant["TMP"].push_back(18.44*eV);
    energyConstant["TMP"].push_back(19.37*eV);
    energyConstant["TMP"].push_back(19.37*eV);
    energyConstant["TMP"].push_back(21.40*eV);
    energyConstant["TMP"].push_back(26.20*eV);
    energyConstant["TMP"].push_back(26.20*eV);
    energyConstant["TMP"].push_back(27.43*eV);
    energyConstant["TMP"].push_back(35.23*eV);
    energyConstant["TMP"].push_back(37.67*eV);
    energyConstant["TMP"].push_back(37.67*eV);
    energyConstant["TMP"].push_back(39.64*eV);
    energyConstant["TMP"].push_back(152.42*eV);
    energyConstant["TMP"].push_back(152.42*eV);
    energyConstant["TMP"].push_back(152.44*eV);
    energyConstant["TMP"].push_back(209.59*eV);
    energyConstant["TMP"].push_back(306.92*eV);
    energyConstant["TMP"].push_back(306.92*eV);
    energyConstant["TMP"].push_back(306.92*eV);
    energyConstant["TMP"].push_back(557.34*eV);
    energyConstant["TMP"].push_back(559.40*eV);
    energyConstant["TMP"].push_back(559.40*eV);
    energyConstant["TMP"].push_back(559.41*eV);
    energyConstant["TMP"].push_back(2178.05*eV);

    std::map<G4String, std::vector<G4double> >::iterator it;
    for(it=energyConstant.begin();it!=energyConstant.end();it++)
    {
        nLevels[it->first] = (it->second).size();
    }
}


G4DNAPTBIonisationStructure::~G4DNAPTBIonisationStructure()
{ }


G4double G4DNAPTBIonisationStructure::IonisationEnergy(G4int level, const G4String& materialName)
{
    G4String matNameModif = ReplaceMaterial(materialName);

    // check if the material exist in the map
    if(energyConstant.find(matNameModif)==energyConstant.end())
    {
        std::ostringstream oss;
        oss << "Material name was not found in energyConstantMap. Problematic material is: "<<matNameModif;
        G4Exception("G4DNAPTBIonisationStructure::IonisationEnergy","em0002",
                    FatalException, oss.str().c_str());
    }

    G4double ionisation = 0.;

    if (level >=0 && level < nLevels[matNameModif]) ionisation = energyConstant[matNameModif][level];

    return ionisation;
}

G4int G4DNAPTBIonisationStructure::NumberOfLevels(const G4String& materialName)
{
    G4String matNameModif = ReplaceMaterial(materialName);

    // check if the material exist in the map
    if(nLevels.find(matNameModif)==nLevels.end())
    {
        std::ostringstream oss;
        oss << "Material name was not found in energyConstantMap. Problematic material is: "<<matNameModif;
        G4Exception("G4DNAPTBIonisationStructure::NumberOfLevels","em0002",
                    FatalException, oss.str().c_str());
    }

    return nLevels[matNameModif];
}

G4String G4DNAPTBIonisationStructure::ReplaceMaterial(const G4String& materialName)
{
    G4String materialNameModified (materialName);

    if(materialName=="backbone_THF") materialNameModified  = "THF";
    else if(materialName=="backbone_TMP") materialNameModified  = "TMP";
    else if(materialName=="adenine_PU") materialNameModified  = "PU";
    else if(materialName=="guanine_PU") materialNameModified  = "PU";
    else if(materialName=="thymine_PY") materialNameModified  = "PY";
    else if(materialName=="cytosine_PY") materialNameModified  = "PY";

    return materialNameModified;
}
