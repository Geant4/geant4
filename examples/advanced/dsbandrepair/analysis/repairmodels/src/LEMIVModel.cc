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
/// \file LEMIVModel.cc
/// \brief Implementation of the LEMIVModel class

#include "LEMIVModel.hh"

#include "ClassifiedDamage.hh"
#include "Damage.hh"
#include "DamageClassifier.hh"

#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

LEMIVModel::LEMIVModel(double pLoopLength,double pNi, double pNc, 
double pNDSB,double pFunrej,double pTfast,double pTslow):
fLoopLength(pLoopLength),
fNi(pNi),
fNc(pNc),
fNDSB(pNDSB),
fFunrej(pFunrej),
fTfast(pTfast),
fTslow(pTslow)
{
	fBpForDSB = 10;
	fDefaultsChromosomeSizes={250, 250, 242, 242, 198, 198, 190, 190, 182, 
					182, 171, 171, 159, 159, 145, 145, 138, 138, 134, 
					134, 135, 135, 133, 133, 114, 114, 107, 107, 102,
					102, 90, 90, 83, 83, 80, 80, 59, 59, 64, 64, 47, 47, 51, 51, 156, 57}; // in Mbp
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

double LEMIVModel::ComputeUnrej(double pTime)
{
	double lambdac = (fNDSB-fNi)/fNc;
	double Ffast = fNi/fNDSB;
	double Fslow = fNc*lambdac/fNDSB;

	return Ffast*std::exp(-std::log(2)*pTime/fTfast)+
	(Fslow-fFunrej)*std::exp(-std::log(2)*pTime/fTslow)+
	fFunrej;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LEMIVModel::ComputeAndSetDamageInput(std::vector<Damage> vecDamage)
{

	double nidsb=0;
	double ncdsb=0;
	double ndsb=0;

	DamageClassifier damclass = DamageClassifier();

	auto sortedDamage = damclass.SortDamageByChromo(vecDamage);
	// Check chromosome sizes:
	if (fChromosomeBpMap.size() == 0) {// using default values
		std::cout<<"=====> LEMIV Calculation will Using Default Choromosome sizes"<<std::endl;
		for (int ii=0;ii<fDefaultsChromosomeSizes.size();ii++) {
			unsigned long long int nBp = fDefaultsChromosomeSizes[ii]*1E6; // convert MBp tp Bp
			fChromosomeBpMap.insert({ii,nBp});
		}
	}
	// Loop on each chromosome
	int i=0;
	for(auto it=sortedDamage.begin();it!=sortedDamage.end();it++)
	{
		i++;
		// Damage are now sorted by event, push all the damage in the same vector
		std::vector<Damage> chromoDamage;
		for(auto itt=it->second.begin();itt!=it->second.end();itt++)
		{
			std::move(itt->second.begin(), itt->second.end(), std::back_inserter(chromoDamage));
		}
		
		// sort the list of damage by ascending bp
		std::sort(chromoDamage.begin(), chromoDamage.end(),
			[](const Damage& a, const Damage& b) {
				return a.GetCopyNb() < b.GetCopyNb();
			});

		// for each loop inside the chromosome
		auto chromID = it->first;
		if (fChromosomeBpMap.find(chromID) == fChromosomeBpMap.end()) {
			std::cerr<<"**** Fatal Error *****"<<std::endl;
			std::cerr<<"LEMIVModel::ComputeAndSetDamageInput: Cannot find size info for chrom ID "
					<<chromID<<std::endl;
			std::cerr<<"*************** *****"<<std::endl;
			exit(EXIT_FAILURE);
		}
		for(auto startLoop=0;startLoop<fChromosomeBpMap[chromID];startLoop+=fLoopLength)
		{

			int n = GetDSBPerLoop(chromoDamage,startLoop);

			ndsb+=n;

			if(n==1)
				nidsb+=1.0;
			if(n>=2)
				ncdsb+=1.0;
		} 
	}

	// Set in the model input parameters
	SetNumDSB(ndsb/fDose);
	SetNumDomainIsolated(nidsb/fDose);
	SetNumDomainClustered(ncdsb/fDose);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int LEMIVModel::GetDSBPerLoop(std::vector<Damage> vecDamage,unsigned int startLoop)
{
	// Start to fill a vector with damage having bp between startLopp and startLoop+2Mbp
	std::vector<Damage> loopDamage;

	for(int i=0;i<vecDamage.size();i++)
	{
		if((vecDamage[i].GetCopyNb()>startLoop)&&(vecDamage[i].GetCopyNb()<startLoop+fLoopLength))
		{
			loopDamage.push_back(vecDamage[i]);
		}
	}

	// Make cluster
	DamageClassifier dam;
	auto classifiedDamage = dam.MakeCluster(loopDamage,fBpForDSB,false);

	// Return the number of DSB in this loop
	return dam.GetNumDSB(classifiedDamage);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LEMIVModel::CalculateRepair(double pTMax, double pDeltaT)
{
	fUCurve.clear();

	for(double time=0.;time<=pTMax;time+=pDeltaT)
	{
		fUCurve.push_back(std::make_pair(time,ComputeUnrej(time)));
	}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LEMIVModel::WriteOutput(std::string pFileName)
{
	std::fstream file;
	file.open(pFileName.c_str(), std::ios_base::out);
	//Header part
	file <<"#============================================= LEMIV  MODEL =============================================#\n";
	file << "                                  LEMIV Model, CalculateRepair with:\n";
	file << "#Number of DSBs: " << fNDSB * fDose    << " DSBs.\n";
	file << "#Number of domains with clustered DSB, Nc         = " << fNc * fDose    << " domains.\n";
	file << "#Number of domains with isolated DSB, Ni    = " << fNi * fDose<< " domains.\n";
	file << "#Funrej     = " << fFunrej << "\n";
	file << "#Tfast 	 = " << fTfast  << " h-1              " << "#Tslow = " << fTslow << " h-1\n";
	file << "#LoopLength (length of domain) = " << fLoopLength<<" bp \n";
	file <<"#========================================================================================================#\n";
	file << "Time (h)\tU\n";
	//End Header part
	for(int i=0;i<fUCurve.size();i++)
	{
		file << fUCurve[i].first << "\t" << fUCurve[i].second << "\n";
	}

	file.close();
}