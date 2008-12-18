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
#include "G4AdjointCSMatrix.hh"
#include <iomanip>
#include <fstream>

#include "G4AdjointInterpolator.hh"

///////////////////////////////////////////////////////
//
G4AdjointCSMatrix::G4AdjointCSMatrix(G4bool aBool){
	theLogPrimEnergyVector.clear();
	theLogCrossSectionVector.clear();
	theLogSecondEnergyMatrix.clear();
	theLogProbMatrix.clear();
	theLogProbMatrixIndex.clear();
	log0Vector.clear();
	nb_of_PrimEnergy=0;
	is_scat_proj_to_proj_case  =aBool;
}
///////////////////////////////////////////////////////
//
G4AdjointCSMatrix::~G4AdjointCSMatrix(){
	theLogPrimEnergyVector.clear();
	theLogCrossSectionVector.clear();
	theLogSecondEnergyMatrix.clear();
	theLogProbMatrix.clear();
}
///////////////////////////////////////////////////////
//
void G4AdjointCSMatrix::Clear()
{
	theLogPrimEnergyVector.clear();
	theLogCrossSectionVector.clear();
	theLogSecondEnergyMatrix.clear();
	theLogProbMatrix.clear();
	theLogProbMatrixIndex.clear();
	log0Vector.clear();
	nb_of_PrimEnergy=0;
}
///////////////////////////////////////////////////////
//
 void G4AdjointCSMatrix::AddData(G4double aLogPrimEnergy,G4double aLogCS, std::vector< G4double>* aLogSecondEnergyVector,
	 							       std::vector< G4double>* aLogProbVector,size_t n_pro_decade){
	
	G4AdjointInterpolator* theInterpolator=G4AdjointInterpolator::GetInstance();
	//Add this time we consider that the energy are given monotically
	
	theLogPrimEnergyVector.push_back(aLogPrimEnergy);
	theLogCrossSectionVector.push_back(aLogCS);
	theLogSecondEnergyMatrix.push_back(aLogSecondEnergyVector);
	//G4cout<<"Test Add Data "<<this<<'\t'<<aSecondEnergyVector->size()<<std::endl;
	//G4cout<<theSecondEnergyMatrix.size()<<std::endl;
	theLogProbMatrix.push_back(aLogProbVector);
	//G4cout<<"Test Add Data 1 "<<this<<'\t'<<aSecondEnergyVector->size()<<std::endl;
	//G4cout<<theSecondEnergyMatrix.size()<<std::endl;
	std::vector< size_t>* aLogProbVectorIndex = 0;
	dlog =0;
	if (n_pro_decade > 0 && aLogProbVector->size()>0) {
		aLogProbVectorIndex = new std::vector< size_t>();
		dlog=std::log(10.)/n_pro_decade;
		G4double log_val = int(std::min((*aLogProbVector)[0],aLogProbVector->back())/dlog)*dlog;
		log0Vector.push_back(log_val);
		while(log_val<0.) {
			aLogProbVectorIndex->push_back(theInterpolator->FindPosition(log_val,(*aLogProbVector)));
			log_val+=dlog;
		} 
	}
	else {
		log0Vector.push_back(0.);
	}
	theLogProbMatrixIndex.push_back(aLogProbVectorIndex);
	
	
	nb_of_PrimEnergy++;
	
	
}
///////////////////////////////////////////////////////
//
bool G4AdjointCSMatrix::GetData(unsigned int i, G4double& aLogPrimEnergy,G4double& aLogCS,G4double& log0, std::vector< G4double>*& aLogSecondEnergyVector,
	 							       std::vector< G4double>*& aLogProbVector, std::vector< size_t>*& aLogProbVectorIndex)
{	if (i>= nb_of_PrimEnergy) return false;
	//G4cout<<"Test Get Data "<<std::endl;
	aLogPrimEnergy = theLogPrimEnergyVector[i];
	aLogCS = theLogCrossSectionVector[i];
	aLogSecondEnergyVector = theLogSecondEnergyMatrix[i];
	//G4cout<<"Test Get Data "<<this<<'\t'<<theSecondEnergyMatrix[i]->size()<<std::endl;
	//G4cout<<"Test Get Data "<<this<<'\t'<<aSecondEnergyVector->size()<<std::endl;
	//G4cout<<"Test Get Data "<<this<<'\t'<<aSecondEnergyVector<<std::endl;
	aLogProbVector = theLogProbMatrix[i];
	aLogProbVectorIndex = theLogProbMatrixIndex[i];
	log0=log0Vector[i];
	//G4cout<<"Test Get Data 1 "<<this<<'\t'<<theProbMatrix[i]->size()<<std::endl;
	//G4cout<<"Test Get Data 1 "<<this<<'\t'<<aProbVector->size()<<std::endl;
	//G4cout<<"Test Get Data 1 "<<this<<'\t'<<aLogProbVectorIndex<<std::endl;
	return true;
	
}
///////////////////////////////////////////////////////
//
void G4AdjointCSMatrix::Write(G4String file_name)
{	std::fstream FileOutput(file_name, std::ios::out);			  
 	FileOutput<<std::setiosflags(std::ios::scientific);
 	FileOutput<<std::setprecision(6);
	FileOutput<<theLogPrimEnergyVector.size()<<std::endl;
	for (size_t i=0;i<theLogPrimEnergyVector.size();i++){
		FileOutput<<std::exp(theLogPrimEnergyVector[i])/MeV<<'\t'<<std::exp(theLogCrossSectionVector[i])<<std::endl;
		size_t j1=0;
		FileOutput<<theLogSecondEnergyMatrix[i]->size()<<std::endl;
		for (size_t j=0;j<theLogSecondEnergyMatrix[i]->size();j++){
			FileOutput<<std::exp((*theLogSecondEnergyMatrix[i])[j]);
			j1++;
			if (j1<10) FileOutput<<'\t';
			else  {
				FileOutput<<std::endl;
				j1=0;
			}	
		}
		if (j1>0) FileOutput<<std::endl;
		j1=0;
		FileOutput<<theLogProbMatrix[i]->size()<<std::endl;
		for (size_t j=0;j<theLogProbMatrix[i]->size();j++){
			FileOutput<<std::exp((*theLogProbMatrix[i])[j]);
			j1++;
			if (j1<10) FileOutput<<'\t';
			else  {
				FileOutput<<std::endl;
				j1=0;
			}	
		}
		if (j1>0) FileOutput<<std::endl;
		
		
	}
	
}
///////////////////////////////////////////////////////
//
void G4AdjointCSMatrix::Read(G4String file_name)
{	std::fstream FileOutput(file_name, std::ios::in);
	size_t n1,n2;
	
	
	theLogPrimEnergyVector.clear();
	theLogCrossSectionVector.clear();
	theLogSecondEnergyMatrix.clear();
	theLogProbMatrix.clear();
	FileOutput>>n1;
	for (size_t i=0; i<n1;i++){
		G4double E,CS;
		FileOutput>>E>>CS;
		theLogPrimEnergyVector.push_back(E);
		theLogCrossSectionVector.push_back(CS);
		FileOutput>>n2;
		theLogSecondEnergyMatrix.push_back(new std::vector<double>());
		theLogProbMatrix.push_back(new std::vector<double>());
		
		for (size_t j=0; j<n2;j++){
			G4double E1;
			FileOutput>>E1;
			theLogSecondEnergyMatrix[i]->push_back(E1);
		}
		FileOutput>>n2;
		for (size_t j=0; j<n2;j++){
			G4double prob;
			FileOutput>>prob;
			theLogProbMatrix[i]->push_back(prob);
		}
		
		
		
	}
	
				  
 	
	
}
