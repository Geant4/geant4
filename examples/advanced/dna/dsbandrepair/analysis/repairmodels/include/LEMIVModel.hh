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
/// \file LEMIVModel.hh
/// \brief Definition of the LEMIVModel class

#ifndef LEMIVMODEL_HH
#define LEMIVMODEL_HH

#include <string>
#include <vector>
#include <map>

class Damage;

class LEMIVModel
{
public:
	
  /// \brief constructor
	// pTfast and pTslow in h-1
	LEMIVModel(double pLoopLength = 2E6,double pNi = 0, double pNc = 0, 
		double pNDSB = 0,double pFunrej =0,double pTfast =-1,double pTslow =-1);
  /// \brief destructor
	~LEMIVModel() = default;

	double ComputeUnrej(double pTime);

	double GetLoopLength() {return fLoopLength;};
	void SetLoopLength(double pVal){fLoopLength=pVal;};
	double GetNumDSB() {return fNDSB;};
	void SetNumDSB(double pVal){fNDSB=pVal;};
	double GetNumDomainIsolated(){return fNi;};
	void SetNumDomainIsolated(double pVal){fNi=pVal;};
	double GetNumDomainClustered(){return fNc;};
	void SetNumDomainClustered(double pVal){fNc=pVal;};

	double GetFunrej(){return fFunrej;};
	void SetFunrej(double pVal){fFunrej=pVal;};

	// Tfast in h-1
	double GetTfast(){return fTfast;};
	void SetTfast(double pVal){fTfast=pVal;};

	// Tslow in h-1
	double GetTslow(){return fTslow;};
	void SetTslow(double pVal){fTslow=pVal;};

	// Computes and sets input damage parameters of LEMIV
	void ComputeAndSetDamageInput(std::vector<Damage>);

	// Write U=f(t) curve
	// pTMax and pDeltaT in h
	void CalculateRepair(double pTMax, double pDeltaT);

	// Write output
	void WriteOutput(std::string pFileName);

	unsigned int GetBpForDSB(){return fBpForDSB;};
	void SetBpForDSB(unsigned int pVal){fBpForDSB = pVal;};
	void SetDose(double d) {fDose = d;}

	void SetChromosomeBpSizesMap(std::map<int,unsigned long long int> chroSizes) {fChromosomeBpMap = chroSizes;}
private:

	double fLoopLength; // length of the loop in bp, default one is 2 Mbp

	//isolated DSB yield in Gy-1
	double fNi{0};
	//clustered DSB yield in Gy-1
	double fNc{0};
	// DSB yield in Gy-1
	double fNDSB{0};

	double fFunrej{0};

	// constant time in h-1
	double fTfast{0};
	double fTslow{0};

	// Compute the number of DSB for a given loop startint at bp pStartLop
	int GetDSBPerLoop(std::vector<Damage> pVecDamage,unsigned int pStartLoop);

	// U=f(t) curve, time in h
  	std::vector<std::pair<double,double>> fUCurve;
	// num of bp to consider a DSB for the MakeCluster function
	// default value is 10
	unsigned int fBpForDSB{0};
	// store dose deposited in nucleus cell
	double fDose{0};
	std::vector<double> fDefaultsChromosomeSizes;// Chomosome defaults sizes
	std::map<int,unsigned long long int> fChromosomeBpMap; //Store number of Bp in each Chomosomes;
};

#endif