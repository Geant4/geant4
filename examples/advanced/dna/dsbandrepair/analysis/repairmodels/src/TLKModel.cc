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
/// \file TLKModel.cc
/// \brief Implementation of the TLKModel class

#include "TLKModel.hh"

#include "ClassifiedDamage.hh"
#include "Damage.hh"
#include "DamageClassifier.hh"
#include "ODESolver.hh"
#include <cmath>
#include <iostream>
#include <fstream>
#include <functional>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

TLKModel::TLKModel(double pLambda1,double pLambda2, double pBeta1, double pBeta2,double pEta):
fLambda1(pLambda1),
fLambda2(pLambda2),
fBeta1(pBeta1),
fBeta2(pBeta2),
fEta(pEta)
{
	fSingleDSBYield = 0;
	fComplexDSBYield = 0;
	fBpForDSB = 10;
	fStartTime = 0.0;
	fStopTime = 480.0;
	fStepTime = 0.048;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector<double> TLKModel::TLK_odes_system(double t,std::vector<double> y){

	std::vector<double> dxdt;
	double dxdt1 = -fLambda1*y[0]-fEta*y[0]*(y[0]+y[1]);
	double dxdt2 = -fLambda2*y[1]-fEta*y[1]*(y[0]+y[1]);
	double dxdt3 = fBeta1*fLambda1*y[0]+fBeta2*fLambda2*y[1]+0.25*fEta*(y[0]+y[1])*(y[0]+y[1]);
	dxdt.push_back(dxdt1);
	dxdt.push_back(dxdt2);
	dxdt.push_back(dxdt3);
	return dxdt;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TLKModel::ComputeAndSetDamageInput(std::vector<Damage> vecDamage)
{
	DamageClassifier damClass;
	auto classifiedDamage = damClass.MakeCluster(vecDamage,fBpForDSB,false);

	fComplexDSBYield = damClass.GetNumComplexDSB(classifiedDamage);
	fSingleDSBYield = damClass.GetNumDSB(classifiedDamage)-fComplexDSBYield;


	fComplexDSBYield = fComplexDSBYield/fDose;
	fSingleDSBYield = fSingleDSBYield/fDose;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double TLKModel::ComputeSF(double pDose)
{
	std::vector<double> y( 3 );

	y[0] = fSingleDSBYield*pDose;
	y[1] = fComplexDSBYield*pDose;
	y[2] = 0;
	std::function<std::vector<double>(double,std::vector<double>)> 
	func = [this] (double t,std::vector<double> y) -> std::vector<double> {
		return TLK_odes_system(t,y);
	};
	ODESolver odeSolver;
	odeSolver.RungeKutta4(func,y,fStartTime,fStopTime,fStepTime);
	return std::exp(-y[2]);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TLKModel::CalculateRepair(double pDoseMax, double pDeltaDose)
{
	fSFCurve.clear();

	for(double dose=0.;dose<=pDoseMax;dose+=pDeltaDose)
	{
		fSFCurve.push_back(std::make_pair(dose,ComputeSF(dose)));
	}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TLKModel::WriteOutput(std::string pFileName)
{
	std::fstream file;
	file.open(pFileName.c_str(), std::ios_base::out);
	//Header part
	file <<"#============================================= TLK  MODEL =============================================#\n";
	file << "                                  TLK Model, CalculateRepair with:\n";
	file << "#Single DSB  = " << fSingleDSBYield   << " (DSB/Gy)                     " << "#Complex DSB    = " << fComplexDSBYield << " (DSB/Gy)\n";
	file << "#Lambda1  = " << fLambda1     << "                     " << "#Lambda2    = " << fLambda2 << "\n";
	file << "#Beta1    = " << fBeta1       << "                  " << "#Beta2      = " << fBeta2 << "\n";
	file << "#Eta	  = " << fEta         << "\n";
	file <<"#========================================================================================================#\n";
	
	file << "Dose (Gy)\tSF\n";
	//End Header part
	for(int i=0;i<fSFCurve.size();i++)
	{
		file << fSFCurve[i].first << "\t" << fSFCurve[i].second << "\n";
	}

	file.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......