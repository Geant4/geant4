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
/// \file TLKModel.hh
/// \brief Definition of the TLKModel class

#ifndef TLKMODEL_HH
#define TLKMODEL_HH

//#include <boost/numeric/odeint.hpp>
#include <vector>
#include <string>
class Damage;
class TLKModel
{
public:
	
  /// \brief constructor
	TLKModel(double pLambda1 =-1,double pLambda2=-1, double pBeta1=-1, double pBeta2=-1,double pEta=-1);
  /// \brief destructor
	~TLKModel() = default;

	double GetLambda1(){return fLambda1;};
	void SetLambda1(double pVal){fLambda1 = pVal;};
	double GetLambda2(){return fLambda2;};
	void SetLambda2(double pVal){fLambda2 = pVal;};
	double GetBeta1(){return fBeta1;};
	void SetBeta1(double pVal){fBeta1 = pVal;};
	double GetBeta2(){return fBeta2;};
	void SetBeta2(double pVal){fBeta2 = pVal;};
	double GetEta(){return fEta;};
	void SetEta(double pVal){fEta = pVal;};

	double GetSingleDSBYield(){return fSingleDSBYield;};
	void SetSingleDSBYield(double pVal){fSingleDSBYield = pVal;};

	double GetComplexDSBYield(){return fComplexDSBYield;};
	void SetComplexDSBYield(double pVal){fComplexDSBYield = pVal;};

	unsigned int GetBpForDSB(){return fBpForDSB;};
	void SetBpForDSB(unsigned int pVal){fBpForDSB = pVal;};

	double GetStartTime(){return fStartTime;};
	void SetStartTime(double pVal){fStartTime = pVal;};
	double GetStopTime(){return fStopTime;};
	void SetStopTime(double pVal){fStopTime = pVal;};
	double GetStepTime(){return fStepTime;};
	void SetStepTime(double pVal){fStepTime = pVal;};

	// Compute damage inputs required by TLK model
	void ComputeAndSetDamageInput(std::vector<Damage>);

	// Compute SF for a given absorbed dose expressed in Gy
	double ComputeSF(double pDose);

	// Compute a SF Curve
	void CalculateRepair(double pDoseMax, double pDeltaDose);

	// Write output, dose expressed in Gy
	void WriteOutput(std::string pFileName);
	
	void SetDose(double d) {fDose = d;}
private:
	std::vector<double> TLK_odes_system(double t,std::vector<double> y);
	// TLK parameters
	// simple DSB repair probability (h-1)	
	double fLambda1{0};
	// complex DSB repair probability (h-1)	
	double fLambda2{0};
	// simple DSB misrepair probability (h-1)
	double fBeta1{0};
	// complex DSB misrepair probability (h-1)
	double fBeta2{0};
	// binary misrepair probability (h-1)
	double fEta{0};

	// DSB Yields in Gy-1
	double fSingleDSBYield{0};
	double fComplexDSBYield{0};

	// num of bp to consider a DSB for the MakeCluster function
	// default value is 10
	unsigned int fBpForDSB{10};

	// Time for integration
	double fStartTime{0};
	double fStopTime{0};
	double fStepTime{0};

	// SF curve
  	std::vector<std::pair<double,double>> fSFCurve;
	// store dose deposited in nucleus cell
	double fDose{0};
};

#endif