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
/// \file ParametersParser.hh
/// \brief Definition of the ParametersParser class

#ifndef ParametersParser_h
#define ParametersParser_h 1
#include <string>
class ParametersParser
{
public:
    static ParametersParser* Instance();
    ~ParametersParser() = default;
    void LoadParameters(const std::string &fileName);
    std::string GetTLKLambda1() {return TLKLambda1;}
    std::string GetTLKLambda2() {return TLKLambda2;}
    std::string GetTLKBeta1() {return TLKBeta1;}
    std::string GetTLKBeta2() {return TLKBeta2;}
    std::string GetTLKEta() {return TLKEta;}
    std::string GetTLKdoseMax() {return TLKdoseMax;}
    std::string GetTLKdeltaDose() {return TLKdeltaDose;}
    std::string GetEMIVLoopLength() {return LEMIVLoopLength;}
    std::string GetEMIVNi() {return LEMIVNi;}
    std::string GetEMIVNc() {return LEMIVNc;}
    std::string GetEMIVNDSB() {return LEMIVNDSB;}
    std::string GetEMIVFunrej() {return LEMIVFunrej;}
    std::string GetEMIVTFast() {return LEMIVTfast;}
    std::string GetEMIVTSlow() {return LEMIVTslow;}
    std::string GetLEMtimeMax() {return LEMIVtimeMax;}
    std::string GetLEMdeltaTime() {return LEMIVdeltaTime;}
    std::string GetBELOVNirrep() {return BELOVNirrep;}
    std::string GetBELOVDz() {return BELOVDz;}
    std::string GetThresholdE() {return fThresholdE;}
    std::string GetProbabilityForIndirectSB() {return fProbabilityForIndirectSB;}
    std::string GetParticleName() {return fParticleName;}
    float GetParticleEnergy() {return fParticleEnergy;}
    std::string GetEnergyUnit() {return fEnergyUnit;}
    std::string GetEndTimeForChemReactions() {return fEndTimeForChemReactions;}
    int GetNumberOfParticles() {return fNumberOfParticles;}
    int GetBpForDSB() {return BpForDSB;}
    bool UseTLK();
    bool UseLEMIV();
    bool UseBelov();
    bool WannaLoadDamagesFromSDD() {return fLoadDamagesFromSDD;}
    std::string GetOutputName() {return fOutputName;};
    std::string GetSDDFileName() {return fSDDfileName;};
    std::string GetCellNucleusName() {return fCellNucleusName;};
    int GetUnitTypeOfNormalization() {return fUnitOfNormalization;}
    bool WannaSkipScanningIndirectDamage() {return fSkipScanningIndirectDamage;}
private:
    explicit ParametersParser();
    static ParametersParser* fInstance;
    std::string fOutputName{"Output.dat"};
    std::string fSDDfileName{""};
    std::string fCellNucleusName{"Undefined"};
    std::string fThresholdE{""};
    std::string fProbabilityForIndirectSB{""};
    // num of bp to consider a DSB for the MakeCluster function default value is 10
    int BpForDSB{0}; 
    // TLK parameters
    std::string useTLK{"true"};
	// simple DSB repair probability (h-1)	
	std::string TLKLambda1{""};
	// complex DSB repair probability (h-1)	
	std::string TLKLambda2{""};
	// simple DSB misrepair probability (h-1)
	std::string TLKBeta1{""};
	// complex DSB misrepair probability (h-1)
	std::string TLKBeta2{""};
	// binary misrepair probability (h-1)
	std::string TLKEta{""};
    std::string TLKdoseMax{""}, TLKdeltaDose{""};

    // LEMIV parameters
    std::string useLEMIV{"true"};
    std::string LEMIVLoopLength{""}; // length of the loop in Mbp, 
	//isolated DSB yield in Gy-1
	std::string LEMIVNi{""};
	//clustered DSB yield in Gy-1
	std::string LEMIVNc{""};
	// DSB yield in Gy-1
	std::string LEMIVNDSB{""};
	std::string LEMIVFunrej{""};
	// constant time in h-1
	std::string LEMIVTfast{""};
	std::string LEMIVTslow{""};
    std::string LEMIVtimeMax{""}, LEMIVdeltaTime{""};

    // BELOV parameters
    std::string useBELOV{"false"};
    std::string BELOVNirrep{""};
    std::string BELOVDz{""};

    // source info
    std::string fParticleName{""};
    float fParticleEnergy{0};
    std::string fEnergyUnit{""};
    int fNumberOfParticles{0};

    std::string fEndTimeForChemReactions{""};

    bool fLoadDamagesFromSDD{false};
    int fUnitOfNormalization{1}; // unit type for normization: 2 : [Gy-1]; 1: [Gy-1 * Gbp-1]
    bool fSkipScanningIndirectDamage{false};
};

#endif