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
// Authors: O. Belov and M. Batmunkh
// January 2017
// last edit: L.T. Anh (2023)
/// \file BelovModel.cc
/// \brief Implementation of the BelovModel class

#include "BelovModel.hh"
#include "DamageClassifier.hh"
#include "ODESolver.hh"
#include <iostream>
#include <fstream>
#include <functional>
#include <limits>
#include <cmath>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <time.h>
#include <ctime>
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BelovModel::BelovModel():
fDz(0.),
falpha(0.),
fNirrep(0.),
fTime(0.)
{
  for(int i=0;i<5;i++){
    frepairsim[i].clear();
  }
  fdnarepair.clear();
  fBpForDSB = 10;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BelovModel::Initialize()
{
  for(int i=0;i<5;i++){
    frepairsim[i].clear();
  }
  fdnarepair.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool BelovModel::CalculateRepair(double Dz)
{
	// Recalling model parameters
	fDz = Dz;

	std::cout << "Belov Model, CalculateRepair with:" << std::endl;
	std::cout << "   - Dz = " << fDz << " Gy" << std::endl;
	std::cout << "   - alpha = " << falpha << " Gy-1" << std::endl;
	std::cout << "   - Nirrep = " << fNirrep << std::endl;

	if(falpha==0)
	{
	std::cout << " falpha=0:\n"
				<<"- please use ComputeAndSetDamageInput function before calculating repair" 
				<< std::endl
				<<"- if above checked, then the reason might be: nDSBYiels is zero !!!" 
				<< std::endl;
	return false;
	}

	if(fNirrep==0)
	{
	std::cout << " fNirrep=0:\n"
				<<"- please use ComputeAndSetDamageInput function before calculating repair" 
				<< std::endl
				<<"- if above checked, then the reason might be: nComplexDSBYiels is zero !!!" 
				<< std::endl;
	return false;
	}

	// INITIAL CONDITIONS

	int NbEquat = 29;   // Total number of model equations
	std::vector<double> Y(NbEquat,0);

	//---- Initial conditions for NHEJ -----
	
	Y[0] = falpha; 
	Y[1] = Y[2] = Y[3] = Y[4] = Y[5] = Y[6] = Y[7] = Y[8] = Y[9] = 0.; 

	//---- Initial conditions for HR -------
	Y[10] = Y[11] = Y[12] = Y[13] = Y[14] = Y[15] = Y[16] = Y[17] 
	= Y[18] = Y[19] = 0.; 

	//---- Initial conditions for SSA -----
	Y[20] = Y[21] = Y[22] = Y[23] = Y[24] = 0.; 

	//---- Initial conditions for Alt-NHEJ (MMEJ) -----
	Y[25] = Y[26] = Y[27] = Y[28] = 0.; 
	
	// Integration parameters
	double t0 = 0.0;   // Starting time point (dimensionless)
	double t1 = 45.3;  // Final time point (dimensionless)
	double dt = 2.e-6; // Intergration time step (dimensionless)

	double K8 = 0.552; // [h-1], scaling variable

	std::function<std::vector<double>(double,std::vector<double>)> 
	func = [this] (double t,std::vector<double> y) -> std::vector<double> {
		return Belov_odes_system(t,y);
	};

	std::vector<std::vector<double>> Y_vec;
	std::vector<double> times;
	double epsilon = 0.1;
	ODESolver odeSolver;
	odeSolver.SetNstepsForObserver(16000);
	odeSolver.Embedded_RungeKutta_Fehlberg(func,Y,t0,t1,dt,epsilon,&times,&Y_vec);
	//odeSolver.RungeKutta4(func,Y,t0,t1,dt,&times,&Y_vec);
	size_t steps = Y_vec.size();

	// Output options for different repair stages
	size_t Nfoci=5;
	std::string FociName;
	double maxY= -1e-9;

	for (size_t ifoci=0;ifoci<Nfoci;ifoci++)
	{
	maxY = -1e-9;

	if(ifoci==0)FociName = std::string("Ku"     ); 
	if(ifoci==1)FociName = std::string("DNAPKcs"); 
	if(ifoci==2)FociName = std::string("RPA"    ); 
	if(ifoci==3)FociName = std::string("Rad51"  ); 
	if(ifoci==4)FociName = std::string("gH2AX"  ); 
	for( size_t istep=0; istep<steps; istep+=1 )
	{
		double val = 0;
		if(FociName=="Ku"     ) val = Y_vec[istep][1];
		if(FociName=="DNAPKcs") {
			val = Y_vec[istep][3] +Y_vec[istep][4] +Y_vec[istep][5]+Y_vec[istep][6]+Y_vec[istep][7];
		}
		if(FociName=="RPA"    ) val = Y_vec[istep][14]+Y_vec[istep][15]+Y_vec[istep][20];
		if(FociName=="Rad51"  ) val = Y_vec[istep][15]+Y_vec[istep][16]+Y_vec[istep][17];
		if(FociName=="gH2AX"  ) val = Y_vec[istep][9];
		if (maxY < val ) maxY = val;
		double time = times[istep]*K8;
		frepairsim[ifoci].push_back(std::make_pair(time,val));
	}
	fdnarepair.insert(make_pair(FociName,frepairsim[ifoci]));
	}
	Y_vec.clear();Y_vec.shrink_to_fit();
	return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BelovModel::ComputeAndSetDamageInput(std::vector<Damage> vecDamage)
{

  DamageClassifier damClass;
  auto classifiedDamage = damClass.MakeCluster(vecDamage,fBpForDSB,false);

  ComplexDSBYield = damClass.GetNumComplexDSB(classifiedDamage);
  DSBYield = damClass.GetNumDSB(classifiedDamage);
  falpha = DSBYield/fDose;
  
  fNirrep = ComplexDSBYield/DSBYield;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector<double> BelovModel::Belov_odes_system(double t,std::vector<double> Y)
{
	// DSBRepairPathways
	std::vector<double> YP;
	// Concentrations of repair enzymes set to be constant
	double X1; // [Ku]
	double X2; // [DNAPKcsArt]
	double X3; // [LigIV/XRCC4/XLF]
	double X4; // [PNKP]
	double X5; // [Pol]
	double X6; // [H2AX]
	double X7; // [MRN/CtIP/ExoI/Dna2]
	double X8; // [ATM]
	double X9; // [RPA]
	double X10; // [Rad51/Rad51par/BRCA2]
	double X11; // [DNAinc]
	double X12; // [Rad52]
	double X13; // [ERCC1/XPF]
	double X14; // [LigIII]
	double X15; // [PARP1]
	double X16; // [Pol]
	double X17; // [LigI]

	X1 = X2 = X3 = X4 = X5 = X6 = X7 = X8 = 
	X9 = X10 = X11 = X12 = X13 = X14 =
	X15 = X16 = X17 = 400000.; 

	fTime = t;  // Recalling t

	//  DIMENSIONAL REACTION RATES
	//
	//------------NHEJ--------------
	double K1 = 11.052;             // M-1*h-1 
	double Kmin1 = 6.59999*1e-04;   // h-1
	double K2 = 18.8305*(1.08517-std::exp(-21.418/std::pow(fDz,1.822))); // M-1*h-1 
	double Kmin2 = 5.26*1e-01;      //h-1
	double K3 = 1.86;               // h-1
	double K4 = 1.38*1e+06;          // M-1*h-1
	double Kmin4 = 3.86*1e-04;      // h-1
	double K5 = 15.24;              // M-1*h-1
	double Kmin5 = 8.28;            // h-1
	double K6 = 18.06;              // M-1*h-1
	double Kmin6 = 1.33;            // h-1
	double K7 = 2.73*1e+05;          // M-1*h-1
	double Kmin7 = 3.2;             // h-1
	double K8 = 5.52*1e-01;         // h-1
	double K9 = 1.66*1e-01;         // h-1
	double K10 = (1.93*1e-07)/fNirrep; // M
	double K11 = 7.50*1e-02;        // h-1
	double K12 = 11.1;              // h-1
	//
	//------------HR--------------
	double P1 = 1.75*1e+03;          // M-1*h-1
	double Pmin1 = 1.33*1e-04;      // h-1
	double P2 = 0.39192;            // h-1
	double Pmin2 = 2.7605512*1e+02;  // h-1
	double P3 = 1.37*1e+04;          // M-1*h-1
	double Pmin3 = 2.34;            // h-1
	double P4 = 3.588*1e-02;       // h-1
	double P5 = 1.20*1e+05;          // M-1*h-1
	double Pmin5 = 8.82*1e-05;      // h-1
	double P6 = 1.54368*1e+06;       // M-1*h-1
	double Pmin6 = 1.55*1e-03;      // h-1
	double P7 = 1.4904;              // h-1
	double P8 = 1.20*1e+04;          // M-1*h-1
	double Pmin8 = 2.49*1e-04;      // h-1
	double P9 = 1.104;            //h-1
	double P10 = 7.20*1e-03;        // h-1
	double P11 = 6.06*1e-04;        // h-1
	double P12 = 2.76*1e-01;        // h-1
	//
	//------------SSA--------------
	double Q1 = 1.9941*1e+05;       // M-1*h-1
	double Qmin1 = 1.71*1e-04;      // h-1
	double Q2 = 4.8052*1e+04;       // M-1*h-1
	double Q3 = 6*1e+03;             // M-1*h-1
	double Qmin3 = 6.06*1e-04;      // h-1
	double Q4 = 1.62*1e-03;         // h-1
	double Q5 = 8.40*1e+04;          // M-1*h-1
	double Qmin5 = 4.75*1e-04;      // h-1
	double Q6 = 11.58;              // h-1
	//
	//-------alt-NHEJ (MMEJ)--------
	double R1 = 2.39*1e+03;          // M-1*h-1
	double Rmin1 = 12.63;           // h-1
	double R2 = 4.07*1e+04;          // M-1*h-1
	double R3 = 9.82;               // h-1
	double R4 = 1.47*1e+05;          // M-1*h-1
	double Rmin4 = 2.72;            // h-1
	double R5 = 1.65*1e-01;         //h-1
	//
	// Scalling rate XX1
	double XX1 = 9.19*1e-07; // M
	//
	// DIMENSIONLESS REACTION RATES 
	//
	//------------NHEJ--------------
	double k1 = K1*XX1/K8;
	double kmin1 = Kmin1/K8;
	double k2 = K2*XX1/K8; 
	double kmin2 = Kmin2/K8;     
	double k3 = K3/K8;               
	double k4 = K4*XX1/K8;         
	double kmin4 = Kmin4/K8;      
	double k5 = K5*XX1/K8;             
	double kmin5 = Kmin5/K8;            
	double k6 = K6*XX1/K8;              
	double kmin6 = Kmin6/K8;            
	double k7 = K7*XX1/K8;          
	double kmin7 = Kmin7/K8;           
	double k8 = K8/K8;         
	double k9 = K9/K8;         
	double k10 = K10/XX1; 
	double k11 = K11/K8;       
	double k12 = K12/K8;            
	//
	//------------HR--------------
	double p1 = P1*XX1/K8;       
	double pmin1 = Pmin1/K8;     
	double p2 = P2/K8;           
	double pmin2 = Pmin2/K8;  
	double p3 = P3*XX1/K8;          
	double pmin3 = Pmin3/K8;        
	double p4 = P4/K8;       
	double p5 = P5*XX1/K8;       
	double pmin5 = Pmin5/K8;       
	double p6 = P6*XX1/K8;        
	double pmin6 = Pmin6/K8;       
	double p7 = P7/K8;              
	double p8 = P8*XX1/K8;          
	double pmin8 = Pmin8/K8;      
	double p9 = P9/K8;          
	double p10 = P10/K8;        
	double p11 = P11/K8;        
	double p12 = P12/K8;       
	//
	//------------SSA--------------
	double q1= Q1*XX1/K8;       
	double qmin1 = Qmin1/K8;      
	double q2 = Q2*XX1/K8;       
	double q3 = Q3*XX1/K8;        
	double qmin3 = Qmin3/K8;     
	double q4 = Q4/K8;         
	double q5 = Q5*XX1/K8;          
	double qmin5 = Qmin5/K8;     
	double q6 = Q6/K8;             
	//
	//-------alt-NHEJ (MMEJ)--------
	double r1 = R1*XX1/K8;         
	double rmin1 = Rmin1/K8;        
	double r2 = R2*XX1/K8;         
	double r3 = R3/K8;            
	double r4 = R4*XX1/K8;         
	double rmin4 = Rmin4/K8;           
	double r5 = R5/K8;        
	//------------------------------------

	// SYSTEM OF DIFFERENTIAL EQUATIONS

	// ----- NHEJ ----------
	YP.push_back( fNirrep - k1*Y[0]*X1 + kmin1*Y[1] - p1*Y[0]*X1 + pmin1*Y[10]); // [DSB]

	YP.push_back( k1*Y[0]*X1 - kmin1*Y[1] - k2*Y[1]*X2 + kmin2*Y[2]); // [DBS * Ku]

	YP.push_back( k2*Y[1]*X2 - k3*Y[2] - kmin2*Y[2]); // [DSB * DNA-PK/Art]

	YP.push_back( k3*Y[2] - k4*(Y[3]*Y[3]) + kmin4*Y[4]); // [DSB * DNA-PK/ArtP]

	YP.push_back( k4*(Y[3]*Y[3]) - kmin4*Y[4] - k5*Y[4]*X3 + kmin5*Y[5]); // [Bridge]

	YP.push_back( kmin6*Y[6] + k5*Y[4]*X3 - kmin5*Y[5] - k6*Y[5]*X4);
	// [Bridge * LigIV/XRCC4/XLF]

	YP.push_back( -kmin6*Y[6] - k7*Y[6]*X5 +  kmin7*Y[7] + k6*Y[5]*X4); 
	// [Bridge * LigIV/XRCC4/XLF * PNKP]

	YP.push_back( k7*Y[6]*X5 - k8*Y[7] - kmin7*Y[7]);
	// [Bridge * LigIV/XRCC4/XLF * PNKP * Pol]

	YP.push_back( r5*Y[28] + k8*Y[7] + p12*Y[18] + p11*Y[19] + q6*Y[24]); // [dsDNA]

	YP.push_back( (k9*(Y[3] + Y[4] + Y[5] + Y[6] + Y[7])*X6)/(k10 + Y[3] + Y[4] + Y[5] 
			+ Y[6] + Y[7]) - k11*Y[8] - k12*Y[9]); // [gH2AX foci]

	// ----- HR ---------- 
	YP.push_back( p1*Y[0]*X7 - pmin1*Y[10] - p3*Y[10]*Y[11] + pmin3*Y[12]); 
	// [MRN/CtIP/ExoI/Dna2]

	YP.push_back( p2*X8 - pmin2*Y[11] - p3*Y[10]*Y[11] + p4*Y[12] + pmin3*Y[12]); 
	// [ATMP]

	YP.push_back( p3*Y[10]*Y[11] - p4*Y[12] - pmin3*Y[12]); 
	// [DSB * MRN/CtIP/ExoI/Dna2 * ATMP]

	YP.push_back( rmin1*Y[25] + p4*Y[12] - r1*X15*Y[13] - p5*Y[13]*X9 + pmin5*Y[14]); 
	// [ssDNA]

	YP.push_back( pmin6*Y[15] + p5*Y[13]*X9 - pmin5*Y[14] - p6*Y[14]*X10 - 
			q1*Y[14]*X12 + qmin1*Y[20]); // [ssDNA * RPA]

	YP.push_back( -p7*Y[15] - pmin6*Y[15] + p6*Y[14]*X10);
	// [ssDNA * RPA * Rad51/Rad51par/BRCA2]

	YP.push_back( p7*Y[15] - p8*Y[16]*X11 + pmin8*Y[17]); // [Rad51 filament]

	YP.push_back( p8*Y[16]*X11 - p9*Y[17] - pmin8*Y[17]); // [Rad51 filament * DNAinc]

	YP.push_back( p9*Y[17] - p10*Y[18] - p12*Y[18]); // [D-loop]

	YP.push_back( p10*Y[18] - p11*Y[19]); // [dHJ]

	// ----- SSA ---------- 
	YP.push_back( q1*Y[14]*X12 - qmin1*Y[20] -  q2*(Y[20]*Y[20]));
	// [ssDNA * RPA * Rad52]

	YP.push_back( q2*(Y[20]*Y[20]) - q3*Y[21]*X13 + qmin3*Y[22]); // [Flap]

	YP.push_back( q3*Y[21]*X13 - q4*Y[22] - qmin3*Y[22]); // [Flap * ERCC1/XPF]

	YP.push_back( q4*Y[22] - q5*Y[23]*X14 + qmin5*Y[24]); // [dsDNAnicks]

	YP.push_back( q5*Y[23]*X14 - q6*Y[24] - qmin5*Y[24]); // [dsDNAnicks * LigIII]

	// ----- MMEJ ---------- 
	YP.push_back( -rmin1*Y[25] - r2*Y[25]*X16 + r1*X15*Y[13]); // [ssDNA * PARP1]

	YP.push_back( r2*Y[25]*X16 - r3*Y[26]); // [ssDNA * Pol]

	YP.push_back( r3*Y[26] - r4*Y[27]*X17 + rmin4*Y[28]); // [MicroHomol]

	YP.push_back( r4*Y[27]*X17 - r5*Y[28] - rmin4*Y[28]); // [MicroHomol * LigI]

	//---------------------------------------------------
	return YP;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector<std::pair<double,double>> BelovModel::GetDNARepair(std::string NameFoci)
{
	decltype(fdnarepair)::iterator it = fdnarepair.find(NameFoci);
  if (it != fdnarepair.end()) {
    return  it->second;
  }
  else{
	std::cerr<<"There is no Foci with name: "<<NameFoci<<" !!!"<<std::endl;
	exit(0);			//exception needed
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BelovModel::WriteOutput(std::string pFileName)
{
  std::fstream file;
  file.open(pFileName.c_str(), std::ios_base::out);
  //Header part
  file <<"#===================================== BELOV  MODEL ========================================#\n";
  file << "                               Belov Model, CalculateRepair with:\n";
  file << "#DSB  = " <<  DSBYield   << " (SB)        " << "#Complex DSB= " << ComplexDSBYield << " (SB)\n";
  file << "#Dz     = " << fDz << " Gy\n";
  file << "#Nirrep = " << fNirrep << "\n";
  file <<"#===========================================================================================#\n";
  file << "Time\t";
  for(auto it=fdnarepair.begin();it!=fdnarepair.end();it++)
  {
    file << it->first << "\t";
  }
  file << "\n";
  //End header part
  int nVal = fdnarepair.begin()->second.size();

  for(int i=0;i<nVal;i++)
  {
    file << fdnarepair["DNAPKcs"][i].first << "\t";

    for(auto it=fdnarepair.begin();it!=fdnarepair.end();it++)
    {
      auto data = it->second;
      file << data[i].second << "\t";
    }
    file << "\n";

  }

  file.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BelovModel::SetDSBandComDSBandDose(double dsb,double cdsb, double d)
{
	SetDose(d);
	ComplexDSBYield = cdsb;
	DSBYield = dsb;
	falpha = DSBYield/fDose;	
	fNirrep = ComplexDSBYield/DSBYield;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
