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
/// \file AnalysisHandler.cc
/// \brief Implementation of the AnalysisHandler class

#include "AnalysisHandler.hh"
#include "ScanDamage.hh"
#include "DamageClassifier.hh"
#include "ParametersParser.hh"
#include "SDDData.hh"

#include <cmath>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

AnalysisHandler::AnalysisHandler(): pTLKDoseMax(6.0), pTLKDeltaDose(0.25), fNBp (-1)
{
    fScanDamage = std::make_unique<ScanDamage>();
    fTLKModel = std::make_unique<TLKModel>();
    fLEMIVModel = std::make_unique<LEMIVModel>();
    fBelovModel = std::make_unique<BelovModel>();
    fBpForDSB = 10;
   
    if (ParametersParser::Instance()->GetThresholdE() != "") {
        auto e = std::stod(ParametersParser::Instance()->GetThresholdE());
        if (e > 0) fScanDamage->SetThresholdEnergy(e);
    }
    if (ParametersParser::Instance()->GetProbabilityForIndirectSB() != "") {
        auto p = std::stod(ParametersParser::Instance()->GetProbabilityForIndirectSB());
        if (p > 0 && p <= 100) fScanDamage->SetProbabilityForIndirectSBSelection(p/100.);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void AnalysisHandler::SetThresholdEnergy(double e)
{
    fScanDamage->SetThresholdEnergy(e);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void AnalysisHandler::GetAllDamageAndScanSB()
{
    std::map<unsigned int,std::map<unsigned int,std::vector<Damage> > > dmMap; 
    if (ParametersParser::Instance()->WannaLoadDamagesFromSDD()) {
        SDDData sdddata(ParametersParser::Instance()->GetSDDFileName());
        dmMap = sdddata.GetAllDamage();
        fDose = sdddata.GetDose();
        fChromosomeBpMap = sdddata.GetChromosomeBpSizesMap(fNBp);
    } else {
        if (ParametersParser::Instance()->WannaSkipScanningIndirectDamage()) {
            fScanDamage->SkipScanningIndirectDamage();
        }
        dmMap = fScanDamage->ExtractDamage();
        fEdepInNucleus = fScanDamage->GetEdepSumInNucleus();//eV
        double nuclesumass = fScanDamage->GetNucleusMass(); // kg
        double eVtoJ = 1.60E-19;
        fDose = fEdepInNucleus*eVtoJ/nuclesumass;
        fNBp = fScanDamage->GetTotalNbBpPlacedInGeo();
        fChromosomeBpMap = fScanDamage->GetChromosomeBpSizesMap();
    }
    DamageClassifier damClass;
    std::map<int,int> ndsbMap, ncdsbMap, nssbMap, nsbMap; 
    std::map<int,int> ndirsbMap, ndsbdirMap, ndsbdirIMap, ndsbInMap;
    for (const auto& [chromo,evtDm] : dmMap) {
        for (const auto& [evt, dmV] : evtDm) {
            fAllDamage.insert(fAllDamage.end(),dmV.begin(),dmV.end()); // to write SDD file and for LEM-IV
            std::vector<Damage> tmpV{dmV};
            auto classifiedDamage = damClass.MakeCluster(tmpV,fBpForDSB,false);
			
            for (auto dm : dmV) {
                if (dm.GetDamageType() == Damage::Damage::fBackbone ) {
                    if ( (nsbMap.find(evt) == nsbMap.end()) ) {
                        nsbMap.insert({evt,1});
                    } else {
                        nsbMap[evt] ++;
                    }
                    if (dm.GetCause() == Damage::Damage::fDirect) {
                        if ( (ndirsbMap.find(evt) == ndirsbMap.end()) ) {
                            ndirsbMap.insert({evt,1});
                        } else {
                            ndirsbMap[evt] ++;
                        }
                    }
                }
            }
			
            int NumDSBForThisCluster = damClass.GetNumDSB(classifiedDamage);
            if ( (ndsbMap.find(evt) == ndsbMap.end()) ) {
                if (NumDSBForThisCluster > 0) ndsbMap.insert({evt,NumDSBForThisCluster});
            } else {
                ndsbMap[evt] += NumDSBForThisCluster;
            }

            int NumcDSBForThisCluster = damClass.GetNumComplexDSB(classifiedDamage);
            if ( (ncdsbMap.find(evt) == ncdsbMap.end()) ) {
                if (NumcDSBForThisCluster > 0) ncdsbMap.insert({evt,NumcDSBForThisCluster});
            } else {
                ncdsbMap[evt] += NumcDSBForThisCluster;
            }

            int NumSSBForThisCluster = damClass.GetNumSSB(classifiedDamage);
            if ( (nssbMap.find(evt) == nssbMap.end()) ) {
                if (NumSSBForThisCluster > 0) nssbMap.insert({evt,NumSSBForThisCluster});
            } else {
                nssbMap[evt] += NumSSBForThisCluster;
            }

            int NumDSBdirForThisCluster = damClass.GetNumDSBwithDirectDamage(classifiedDamage);
            if ( (ndsbdirMap.find(evt) == ndsbdirMap.end()) ) {
                if (NumDSBdirForThisCluster > 0) ndsbdirMap.insert({evt,NumDSBdirForThisCluster});
            } else {
                ndsbdirMap[evt] += NumDSBdirForThisCluster;
            }

            int NumDSBInForThisCluster = damClass.GetNumDSBwithIndirectDamage(classifiedDamage);
            if ( (ndsbInMap.find(evt) == ndsbInMap.end()) ) {
                if (NumDSBInForThisCluster > 0) ndsbInMap.insert({evt,NumDSBInForThisCluster});
            } else {
                ndsbInMap[evt] += NumDSBInForThisCluster;
            }

            int NumDSBdirInForThisCluster = damClass.GetNumDSBwithBothDirectIndirectDamage(classifiedDamage);
            if ( (ndsbdirIMap.find(evt) == ndsbdirIMap.end()) ) {
                if (NumDSBdirInForThisCluster > 0) ndsbdirIMap.insert({evt,NumDSBdirInForThisCluster});
            } else {
                ndsbdirIMap[evt] += NumDSBdirInForThisCluster;
            }
        }
    }

    // DSB and its error:
    float xxtotal=0, rms, xtotal = 0;
    if (ndsbMap.size() >0) {
        for (auto const& [evt,numdsb] : ndsbMap) {
        xtotal += (float)numdsb;
        xxtotal += float(numdsb*numdsb);
        }
        if (ndsbMap.size() == 1) {
            // try to estimate error using poisson distribution
            rms = std::sqrt(xtotal);
        } else rms = std::sqrt(std::fabs(xxtotal - xtotal*xtotal)/float(ndsbMap.size()));
        fNDSBandError.first = xtotal;
        fNDSBandError.second = rms;
    }
   
    
    // cDSB and its error:
    if (ncdsbMap.size() > 0) {
        xxtotal=0, rms = 0, xtotal = 0;
        for (auto const& [evt,numcdsb] : ncdsbMap) {
            xtotal += (float)numcdsb;
            xxtotal += float(numcdsb*numcdsb);
        }
        if (ncdsbMap.size() == 1) {
            // try to estimate error using poisson distribution
            rms = std::sqrt(xtotal);
        } else rms = std::sqrt(std::fabs(xxtotal - xtotal*xtotal)/float(ncdsbMap.size()));
        fNcDSBandError.first = xtotal;
        fNcDSBandError.second = rms;
    }
    // sDSB and its error, using error propagation method:
    if (fNDSBandError.first > 0) {
        fNsDSBandError.first = fNDSBandError.first -  fNcDSBandError.first;
        if (fNcDSBandError.first > 0) {
            fNsDSBandError.second = (fNsDSBandError.first)*std::sqrt(
                (fNDSBandError.second/fNDSBandError.first)*(fNDSBandError.second/fNDSBandError.first) + 
                (fNcDSBandError.second/fNcDSBandError.first)*(fNcDSBandError.second/fNcDSBandError.first));
        }
        else fNsDSBandError.second = (fNsDSBandError.first)*(fNDSBandError.second/fNDSBandError.first);
    }

    // DSBdir and its error:
    if (ndsbdirMap.size() > 0) {
        xxtotal=0, rms = 0, xtotal = 0;
        for (auto const& [evt,numdsbdir] : ndsbdirMap) {
            xtotal += (float)numdsbdir;
            xxtotal += float(numdsbdir*numdsbdir);
        }
        if (ndsbdirMap.size() == 1) {
            // try to estimate error using poisson distribution
            rms = std::sqrt(xtotal);
        } else rms = std::sqrt(std::fabs(xxtotal - xtotal*xtotal)/float(ndsbdirMap.size()));
        fNDSBdirandError.first = xtotal;
        fNDSBdirandError.second = rms;
    }

    // DSBIn and its error:
    if (ndsbInMap.size() > 0) {
        xxtotal=0, rms = 0, xtotal = 0;
        for (auto const& [evt,numdsbIn] : ndsbInMap) {
            xtotal += (float)numdsbIn;
            xxtotal += float(numdsbIn*numdsbIn);
        }
        if (ndsbInMap.size() == 1) {
            // try to estimate error using poisson distribution
            rms = std::sqrt(xtotal);
        } else rms = std::sqrt(std::fabs(xxtotal - xtotal*xtotal)/float(ndsbInMap.size()));
        fNDSBIndandError.first = xtotal;
        fNDSBIndandError.second = rms;
    }

    // DSBdirIn and its error:
    if (ndsbdirIMap.size() > 0) {
        xxtotal=0, rms = 0, xtotal = 0;
        for (auto const& [evt,numdsbdirIn] : ndsbdirIMap) {
            xtotal += (float)numdsbdirIn;
            xxtotal += float(numdsbdirIn*numdsbdirIn);
        }
        if (ndsbdirIMap.size() == 1) {
            // try to estimate error using poisson distribution
            rms = std::sqrt(xtotal);
        } else rms = std::sqrt(std::fabs(xxtotal - xtotal*xtotal)/float(ndsbdirIMap.size()));
        fNDSBdirIandError.first = xtotal;
        fNDSBdirIandError.second = rms;
    }


    // SSB and its error:
    if (nssbMap.size() > 0) {
        xxtotal=0, rms = 0, xtotal = 0;
        for (auto const& [evt,numssb] : nssbMap) {
            xtotal += (float)numssb;
            xxtotal += float(numssb*numssb);
        }
        if (nssbMap.size() == 1) {
            // try to estimate error using poisson distribution
            rms = std::sqrt(xtotal);
        } else rms = std::sqrt(std::fabs(xxtotal - xtotal*xtotal)/float(nssbMap.size()));
        fNSSBandError.first = xtotal;
        fNSSBandError.second = rms;
    }

    // SB and its error:
    if (nsbMap.size() > 0) {
        xxtotal=0, rms = 0, xtotal = 0;
        for (auto const& [evt,numsb] : nsbMap) {
            xtotal += (float)numsb;
            xxtotal += float(numsb*numsb);
        }
        if (nsbMap.size() == 1) {
            // try to estimate error using poisson distribution
            rms = std::sqrt(xtotal);
        } else rms = std::sqrt(std::fabs(xxtotal - xtotal*xtotal)/float(nsbMap.size()));
        fNSBandError.first = xtotal;
        fNSBandError.second = rms;
    }

    // direct SB and its error:
    if (ndirsbMap.size() > 0) {
        xxtotal=0, rms = 0, xtotal = 0;
        for (auto const& [evt,numdirsb] : ndirsbMap) {
            xtotal += (float)numdirsb;
            xxtotal += float(numdirsb*numdirsb);
        }
        if (ndirsbMap.size() == 1) {
            // try to estimate error using poisson distribution
            rms = std::sqrt(xtotal);
        } else rms = std::sqrt(std::fabs(xxtotal - xtotal*xtotal)/float(ndirsbMap.size()));
        fNdirSBandError.first = xtotal;
        fNdirSBandError.second = rms;
    }

    // indirect SB its error, using error propagation method:
    if (fNSBandError.first > 0) {
        fNindirSBandError.first = fNSBandError.first -  fNdirSBandError.first;
        if (fNdirSBandError.first > 0) {
            fNindirSBandError.second = (fNindirSBandError.first)*std::sqrt(
                (fNSBandError.second/fNSBandError.first)*(fNSBandError.second/fNSBandError.first) + 
                (fNdirSBandError.second/fNdirSBandError.first)*(fNdirSBandError.second/fNdirSBandError.first));
        }
        else fNindirSBandError.second = (fNindirSBandError.first)*(fNSBandError.second/fNSBandError.first);
    }

    // clear Maps
    ndsbMap.clear();
    ncdsbMap.clear();
    nssbMap.clear();
    nsbMap.clear();
    ndirsbMap.clear();
    ndsbdirMap.clear();
    ndsbdirIMap.clear();
    ndsbInMap.clear();
    dmMap.clear();

    fIsSBScanned = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void AnalysisHandler::GiveMeSBs()
{
    if (!fIsSBScanned) GetAllDamageAndScanSB();
    std::string outName = ParametersParser::Instance()->GetOutputName();
    std::fstream file;
	file.open(outName.c_str(), std::ios_base::out);
    double norm = 1.0;
    std::string normunit = "";
    if (ParametersParser::Instance()->GetUnitTypeOfNormalization() == 2) {
        norm = 1.0/fDose;
        normunit = "[SB/Gy]";
    } else {
        double BbToGb = 1e-9; // convert Bb to Gb
        norm = 1.0/(fDose*fNBp*BbToGb);
        normunit = "[SB/Gy/Gbp]";
    }
    file <<"#=========================== Strand Breaks ============================#\n";
    file <<"Name of Cell Nucleus: "<<ParametersParser::Instance()->GetCellNucleusName()<<"\n";
    if (ParametersParser::Instance()->WannaLoadDamagesFromSDD()) {
        std::string outstrtmp = "No info from SDD file!!!\n";
        file <<"Volume of Cell Nucleus: "<<outstrtmp;
        file <<"Mass Density of Cell Nucleus: "<<outstrtmp;
        file <<"Mass of Cell Nucleus: "<<outstrtmp;
        file <<"Energy deposited in Cell Nucleus: "<<outstrtmp;
        file <<"Dose delivered in Cell Nucleus: "<<fDose <<" (Gy)\n";
        file <<"Minimum Distance between two clusters: "<<fBpForDSB <<" (bp)\n";
        file <<"Number of basepairs in Cell Nucleus: "<<fNBp <<" (bp)\n";
        file <<"Threshold Energy for direct damage selection: "<<outstrtmp;
        file <<"Propability for indirect damage selection: "<<outstrtmp;
    } else {
        file <<"Volume of Cell Nucleus: "<<fScanDamage->GetNucleusVolume()<<" (m3)\n";
        file <<"Mass Density of Cell Nucleus: "<<fScanDamage->GetNucleusMassDensity()<<" (kg/m3)\n";
        file <<"Mass of Cell Nucleus: "<<fScanDamage->GetNucleusMass()<<" (kg)\n";
        file <<"Energy deposited in Cell Nucleus: "<<fEdepInNucleus <<" (eV)\n";
        file <<"Dose delivered in Cell Nucleus: "<<fDose <<" (Gy)\n";
        file <<"Minimum Distance between two clusters: "<<fBpForDSB <<" (bp)\n";
        file <<"Number of basepairs in Cell Nucleus: "<<fNBp <<" (bp)\n";
        file <<"Threshold Energy for direct damage selection: "<<fScanDamage->GetThresholdEnergy() <<" (eV)\n";
        if (fScanDamage->SkippedScanningIndirectDamage()) file <<"Propability for indirect SB selection: "
                                                                <<" Skipped the indirect analysis\n";
        else file <<"Propability for indirect damage selection: "
                <<fScanDamage->GetProbabilityForIndirectSBSelection()*100.<<" (%)\n";
    }
    
    file <<"#======================================================================#\n";
    file << "\n";
    file <<"#Un-normalized results:\n";
    file << "TotalSB  [SB]   \t" << fNSBandError.first <<"\t+/-\t"<<fNSBandError.second<< "\n";
	file << "DirSB    [SB]   \t" << fNdirSBandError.first <<"\t+/-\t"<<fNdirSBandError.second<< "\n";
	file << "IndirSB  [SB]   \t" << fNindirSBandError.first <<"\t+/-\t"<<fNindirSBandError.second<< "\n";
    file << "SSB      [SB]   \t" << fNSSBandError.first <<"\t+/-\t"<<fNSSBandError.second<< "\n";
	file << "DSB      [SB]   \t" << fNDSBandError.first <<"\t+/-\t"<<fNDSBandError.second<< "\n";
	file << "cDSB     [SB]   \t" << fNcDSBandError.first <<"\t+/-\t"<<fNcDSBandError.second<< "\n";
    file << "sDSB     [SB]   \t" << fNsDSBandError.first <<"\t+/-\t"<<fNsDSBandError.second<< "\n";
    file << "DSBdir   [SB]   \t" << fNDSBdirandError.first <<"\t+/-\t"<<fNDSBdirandError.second<< "\n";
    file << "DSBind   [SB]   \t" << fNDSBIndandError.first <<"\t+/-\t"<<fNDSBIndandError.second<< "\n";
    file << "DSBdirIn [SB]   \t" << fNDSBdirIandError.first <<"\t+/-\t"<<fNDSBdirIandError.second<< "\n";
    file << "\n";
    file <<"#Normalized results:\n";
    file << "TotalSB  " + normunit +"    \t" << fNSBandError.first * norm <<"\t+/-\t"<<fNSBandError.second * norm<< "\n";
	file << "DirSB    " + normunit +"    \t" << fNdirSBandError.first * norm<<"\t+/-\t"<<fNdirSBandError.second * norm<< "\n";
	file << "IndirSB  " + normunit +"    \t" << fNindirSBandError.first * norm<<"\t+/-\t"<<fNindirSBandError.second * norm<< "\n";
    file << "SSB      " + normunit +"    \t" << fNSSBandError.first * norm<<"\t+/-\t"<<fNSSBandError.second * norm<< "\n";
	file << "DSB      " + normunit +"    \t" << fNDSBandError.first * norm<<"\t+/-\t"<<fNDSBandError.second * norm<< "\n";
	file << "cDSB     " + normunit +"    \t" << fNcDSBandError.first * norm<<"\t+/-\t"<<fNcDSBandError.second * norm<< "\n";
    file << "sDSB     " + normunit +"    \t" << fNsDSBandError.first * norm<<"\t+/-\t"<<fNsDSBandError.second * norm<< "\n";
    file << "DSBdir   " + normunit +"    \t" << fNDSBdirandError.first * norm<<"\t+/-\t"<<fNDSBdirandError.second * norm<< "\n";
    file << "DSBind   " + normunit +"    \t" << fNDSBIndandError.first * norm<<"\t+/-\t"<<fNDSBIndandError.second * norm<< "\n";
    file << "DSBdirIn " + normunit +"    \t" << fNDSBdirIandError.first * norm<<"\t+/-\t"<<fNDSBdirIandError.second * norm<< "\n";
    file <<"#======================================================================#\n";
    file << "where: \n";
    file << "-----> TotalSB: Total strand-breaks\n";
    file << "-----> DirSB: Direct strand-breaks\n";
    file << "-----> IndirSB: Indirect strand-breaks\n";
    file << "-----> SSB: Single strand-breaks\n";
    file << "-----> DSB: Double strand-breaks\n";
    file << "-----> cDSB: Complex DSB\n";
    file << "-----> sDSB: Simple DSB\n";
    file << "-----> DSBdir: DSB that contains at least one direct SB\n";
    file << "-----> DSBdind: DSB that contains at least one indirect SB\n";
    file << "-----> DSBdirIn: DSB that contains at both direct and indirect SB\n";
    file <<"#============================== End ===================================#\n";
	file.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void AnalysisHandler::ApplyDNAModel(std::string dnaModel)
{
    if (!fIsSBScanned) GetAllDamageAndScanSB();

    if (dnaModel == "TLK") {
        std::cout << "Invoking TLK Model" << std::endl;
        //fTLKModel->SetDose(fDose);
        SetParametersForTLKModel();
        //fTLKModel->ComputeAndSetDamageInput(fAllDamage);
        fTLKModel->SetSingleDSBYield(fNsDSBandError.first/fDose);
        fTLKModel->SetComplexDSBYield(fNcDSBandError.first/fDose);
        if (ParametersParser::Instance()->GetTLKdoseMax() != "") {
            auto val = std::stod(ParametersParser::Instance()->GetTLKdoseMax());
            if (val != pTLKDoseMax) pTLKDoseMax = val;
        }
        if (ParametersParser::Instance()->GetTLKdeltaDose() != "") {
            auto val = std::stod(ParametersParser::Instance()->GetTLKdeltaDose());
            if (val != pTLKDeltaDose) pTLKDeltaDose = val;
        }
        fTLKModel->CalculateRepair(pTLKDoseMax,pTLKDeltaDose);
        std::string outname = "TLK_"+ParametersParser::Instance()->GetOutputName();
        fTLKModel->WriteOutput(outname);
    }

    if (dnaModel == "LEMIV") {
        std::cout << "Invoking LEMIV Model" << std::endl;
        fLEMIVModel->SetChromosomeBpSizesMap(fChromosomeBpMap);
        fLEMIVModel->SetDose(fDose);
        SetParametersForLEMIVModel();
        fLEMIVModel->ComputeAndSetDamageInput(fAllDamage);
        if (ParametersParser::Instance()->GetLEMtimeMax() != "") {
            auto val = std::stod(ParametersParser::Instance()->GetLEMtimeMax());
            if (val != pLEMIVtimeMax) pLEMIVtimeMax = val;
        }
        if (ParametersParser::Instance()->GetLEMdeltaTime() != "") {
            auto val = std::stod(ParametersParser::Instance()->GetLEMdeltaTime());
            if (val != pLEMIVdeltaTime) pLEMIVdeltaTime = val;
        }
        fLEMIVModel->CalculateRepair(pLEMIVtimeMax,pLEMIVdeltaTime);
        std::string outname = "LEMIV_"+ParametersParser::Instance()->GetOutputName();
        fLEMIVModel->WriteOutput(outname);
    }

    if (dnaModel == "BELOV") {
        std::cout << "Invoking Belov's Model" << std::endl;
        fBelovModel->SetDSBandComDSBandDose(fNDSBandError.first,fNcDSBandError.first,fDose);
        if (ParametersParser::Instance()->GetBELOVNirrep() != "") {
            auto Nirrep = std::stod(ParametersParser::Instance()->GetBELOVNirrep());
            fBelovModel->SetNirrep(Nirrep);
        }
        double Dz = 1.0; 
        if (ParametersParser::Instance()->GetBELOVDz() != "") {
            Dz = std::stod(ParametersParser::Instance()->GetBELOVDz());
        }
        fBelovModel->CalculateRepair(Dz);
        std::string outname = "BELOV_"+ParametersParser::Instance()->GetOutputName();
        fBelovModel->WriteOutput(outname);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void AnalysisHandler::CreateSDD(std::string filename)
{
    std::string str_tmp;
    std::fstream ofile;
	ofile.open(filename.c_str(), std::ios_base::out);
    if (ofile.is_open()) {
        // Create header
        ofile
            <<"SDD version, SDDv1.0;\n"
            <<"Software, dsbandrepair;\n"
            <<"Author contact, Le Tuan Anh - anh.letuan@irsn.fr, , ;\n"
            <<"Simulation Details, DNA damages from direct and indirect effects;\n"
            <<"Source, ;\n"
            <<"Source type, ;\n"
            <<"Incident particles, "<<0<<";\n"
            <<"Mean particle energy ("<<ParametersParser::Instance()->GetEnergyUnit()<<"), "
            <<ParametersParser::Instance()->GetParticleEnergy()<<";\n"
            <<"Energy distribution, , ;\n"
            <<"Particle fraction, 0;\n"
            <<"Dose or fluence, 1, "<<fDose<<";\n"
            <<"Dose rate, 0;\n"
            <<"Irradiation target, ;\n"
            <<"Volumes, 0;\n";
        ofile<<"Chromosome sizes, "<<fChromosomeBpMap.size();
        for (auto const& [chroID, nBps] :fChromosomeBpMap) {
            float nMBps = nBps*1E-6;// convert from Bp to MBp
            ofile<<", "<<nMBps;
        }
        ofile<<";\n";
        ofile<<"DNA Density, 0;\n"
            <<"Cell Cycle Phase, 0;\n"
            <<"DNA Structure, 0;\n"
            <<"In vitro / in vivo, ;\n"
            <<"Proliferation status, ;\n"
            <<"Microenvironment, 0, 0;\n"
            <<"Damage definition, 0;\n"
            <<"Time, 0;\n"
            <<"Damage and primary count, "+std::to_string(fAllDamage.size())+", 0;\n"
            <<"Data entries, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0;\n"
            <<"Data field explaination, Field 1: [1]-eventID, Field 3: [0]-Chromatin "
            <<"type [1]-ChromosomeID [3]-strand, Field 4:chrom position (copynb), Field 5: "
            <<"Cause (direct: [0]=0)  (indirect: [0]=1), Field 6: Damage types (Base:[0]>0) (Backbone: [1]>0);\n"
            <<"\n"
            <<"***EndOfHeader***;\n"
            <<"\n";
        
        // Data Section
        int prevEvt = -1;
        for (auto &damage : fAllDamage) {
            //Field 1 Calassification
            int newEvtFlag = 0; // = 2 if new event;
            if (prevEvt != damage.GetEvt()) {
                newEvtFlag = 2;
                prevEvt = damage.GetEvt();
            }
            ofile<<newEvtFlag<<", "<<damage.GetEvt()<<"; ";
            //Field 3 Chromosome IDs
            ofile<<damage.GetDamageChromatin()<<", "<<damage.GetChromo()<<", "<<0<<", "<<damage.GetStrand()<<"; ";
            //Field 4, Chromosome position 
            ofile<<damage.GetCopyNb()<<"; ";
            //Field 5, Cause: Unknown = -1, Direct = 0, Indirect = 1
            ofile<<damage.GetCause()<<", "<<0<<", "<<0<<"; ";
            //Field 6, Damage types:
            int firstval = 0, secval = 0;
            if (damage.GetDamageType() == Damage::DamageType::fBase) firstval = 1;
            if (damage.GetDamageType() == Damage::DamageType::fBackbone) secval = 1;
            ofile<<firstval<<", "<<secval<<", "<<0<<"; ";
            ofile<<"\n";
        }
    }
    ofile.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void AnalysisHandler::SetBpForDSB(unsigned int pVal)
{
    if (pVal == fBpForDSB) return;
    fBpForDSB = pVal;
    fTLKModel->SetBpForDSB(fBpForDSB);
    fLEMIVModel->SetBpForDSB(fBpForDSB);
    fBelovModel->SetBpForDSB(fBpForDSB);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void AnalysisHandler::SetParametersForTLKModel(double pLambda1,double pLambda2, double pBeta1, double pBeta2,double pEta)
{
    if (ParametersParser::Instance()->GetTLKLambda1() != "") 
        pLambda1 = std::stod(ParametersParser::Instance()->GetTLKLambda1());
    if (ParametersParser::Instance()->GetTLKLambda2() != "") 
        pLambda2 = std::stod(ParametersParser::Instance()->GetTLKLambda2());
    if (ParametersParser::Instance()->GetTLKBeta1() != "") 
        pBeta1 = std::stod(ParametersParser::Instance()->GetTLKBeta1());
    if (ParametersParser::Instance()->GetTLKBeta2() != "") 
        pBeta2 = std::stod(ParametersParser::Instance()->GetTLKBeta2());
    if (ParametersParser::Instance()->GetTLKEta() != "") 
        pEta = std::stod(ParametersParser::Instance()->GetTLKEta());
    if (pBeta1 != fTLKModel->GetBeta1()) fTLKModel->SetBeta1(pBeta1);
    if (pBeta2 != fTLKModel->GetBeta2()) fTLKModel->SetBeta2(pBeta2);
    if (pLambda1 != fTLKModel->GetLambda1()) fTLKModel->SetLambda1(pLambda1);
    if (pLambda2 != fTLKModel->GetLambda2()) fTLKModel->SetLambda2(pLambda2);
    if (pEta != fTLKModel->GetEta()) fTLKModel->SetEta(pEta);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void AnalysisHandler::SetParametersForLEMIVModel(double pLoopLength,double pFunrej,double pTfast,double pTslow)
{
    if (ParametersParser::Instance()->GetEMIVLoopLength() != "") 
        pLoopLength = std::stod(ParametersParser::Instance()->GetEMIVLoopLength());
    if (ParametersParser::Instance()->GetEMIVFunrej() != "") 
        pFunrej = std::stod(ParametersParser::Instance()->GetEMIVFunrej());
    if (ParametersParser::Instance()->GetEMIVTFast() != "") 
        pTfast = std::stod(ParametersParser::Instance()->GetEMIVTFast());
    if (ParametersParser::Instance()->GetEMIVTSlow() != "") 
        pTslow = std::stod(ParametersParser::Instance()->GetEMIVTSlow());
    
    if (pLoopLength != fLEMIVModel->GetLoopLength()) fLEMIVModel->SetLoopLength(pLoopLength);
    if (pFunrej != fLEMIVModel->GetFunrej()) fLEMIVModel->SetFunrej(pFunrej);
    if (pTfast != fLEMIVModel->GetTfast()) fLEMIVModel->SetTfast(pTfast);
    if (pTslow != fLEMIVModel->GetTslow()) fLEMIVModel->SetTslow(pTslow);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....