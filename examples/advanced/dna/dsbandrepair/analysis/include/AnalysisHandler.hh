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
/// \file AnalysisHandler.hh
/// \brief Definition of the AnalysisHandler class

#ifndef AnalysisHandler_h
#define AnalysisHandler_h 1

#include <memory>
#include "ScanDamage.hh"
#include "TLKModel.hh"
#include "LEMIVModel.hh"
#include "BelovModel.hh"

class AnalysisHandler
{
    public:
    AnalysisHandler(/* args */);
    ~AnalysisHandler() = default;
    void SetThresholdEnergy(double e);
    void GetAllDamageAndScanSB();
    void GiveMeSBs();
    void ApplyDNAModel(const std::string dnamodel);
    void SetBpForDSB(unsigned int pVal);
    void SetParametersForTLKModel(double pLambda1 = 3.0,double pLambda2=0.03, 
                                    double pBeta1=0.01, double pBeta2=0.06,double pEta=0.002);
    void SetParametersForLEMIVModel(double pLoopLength=2e6,double pFunrej=0,
                                    double pTfast=0.24,double pTslow=2.81);
    void CreateSDD(std::string filename);
private:
    ///void GetDoseFromEdep();
    std::unique_ptr<ScanDamage> fScanDamage;
    std::unique_ptr<TLKModel> fTLKModel;
    std::unique_ptr<LEMIVModel> fLEMIVModel;
    std::unique_ptr<BelovModel> fBelovModel;
    std::vector<Damage> fAllDamage;
    
    std::pair<float,float> fNsDSBandError = {0,0};
    std::pair<float,float> fNcDSBandError = {0,0};
    std::pair<float,float> fNDSBandError = {0,0};
    std::pair<float,float> fNDSBdirandError = {0,0}; // DSB has contribution from at least one direct  damage
    std::pair<float,float> fNDSBIndandError = {0,0}; // DSB has contribution from at least one indirect  damage
    std::pair<float,float> fNDSBdirIandError = {0,0}; // DSB has contribution from both direct and indirect damage
    std::pair<float,float> fNSSBandError = {0,0};
    std::pair<float,float> fNSBandError = {0,0};
    std::pair<float,float> fNdirSBandError = {0,0};
    std::pair<float,float> fNindirSBandError = {0,0};

    bool fIsSBScanned = false;
    // num of bp to consider a DSB for the MakeCluster function
	// default value is 10
	unsigned int fBpForDSB{10};
    // store dose deposited in nucleus cell
	double fDose{0};
    //TLK: Compute a SF Curve
    double pTLKDoseMax{0}, pTLKDeltaDose{0};
    //LEMIV: compute fraction of unrejoined DSB up to pLEMIVTimeMax (h) and  pLEMIVDeltaTime steps
    double pLEMIVtimeMax{0}, pLEMIVdeltaTime{0};

    double fNBp{0};// number of base pairs
    double fEdepInNucleus{0}; // eV
    double fNucleusVolume{0};

    std::map<int,unsigned long long int> fChromosomeBpMap; //Store number of Bp in each Chomosomes;
};

#endif
