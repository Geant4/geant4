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
/// \file BelovModel.hh
/// \brief Definition of the BelovModel class


#ifndef BelovModel_H
#define BelovModel_H 1

#include <iostream>
#include <fstream>
#include <vector>
#include <map>

class Damage;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class BelovModel 
{ 
public:  
  BelovModel();
  void Initialize();
  bool CalculateRepair(double Dz);

  void SetAlpha(double value){falpha=value;};
  void SetNirrep(double value){fNirrep=value;};
  
  virtual     
  ~BelovModel() = default;  

 //Computes and sets input damage parameters
  void ComputeAndSetDamageInput(std::vector<Damage>);

  std::vector<std::pair<double,double>> GetDNARepair(std::string NameFoci);

  void WriteOutput(std::string pFileName);
  
  unsigned int GetBpForDSB(){return fBpForDSB;};
	void SetBpForDSB(unsigned int pVal){fBpForDSB = pVal;};
  void SetDose(double d) {fDose = d;}
  void SetDSBandComDSBandDose(double dsby,double cdsby,double d);
private:
  double fDz;
  double falpha;
  double fNirrep; 
  double fTime;
  double ComplexDSBYield;
  double DSBYield;
  std::vector<double> Belov_odes_system(double t,std::vector<double> y);
  std::vector<std::pair<double,double>> frepairsim[5];
  std::map<std::string,std::vector<std::pair<double,double>>> fdnarepair;
  // num of bp to consider a DSB for the MakeCluster function
	// default value is 10
	unsigned int fBpForDSB;
  // store dose deposited in nucleus cell
	double fDose;
};
#endif
