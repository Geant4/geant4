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
// $Id: G4AdjointCSMatrix.cc 91870 2015-08-07 15:21:40Z gcosmo $
//

#include <iomanip>
#include <fstream>

#include "G4AdjointCSMatrix.hh"
#include "G4SystemOfUnits.hh"
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
	dlog =0;
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
 void G4AdjointCSMatrix::AddData(G4double aLogPrimEnergy,G4double aLogCS, std::vector< double>* aLogSecondEnergyVector,
	 							       std::vector< double>* aLogProbVector,size_t n_pro_decade){
	
	G4AdjointInterpolator* theInterpolator=G4AdjointInterpolator::GetInstance();
	
	//At this time we consider that the energy is increasing monotically
	theLogPrimEnergyVector.push_back(aLogPrimEnergy);
	theLogCrossSectionVector.push_back(aLogCS);
	theLogSecondEnergyMatrix.push_back(aLogSecondEnergyVector);
	theLogProbMatrix.push_back(aLogProbVector);
	
	std::vector< size_t>* aLogProbVectorIndex = 0;
	dlog =0;
	
	if (n_pro_decade > 0 && aLogProbVector->size()>0) {
		aLogProbVectorIndex = new std::vector< size_t>();
		dlog=std::log(10.)/n_pro_decade;
		G4double log_val = int(std::min((*aLogProbVector)[0],aLogProbVector->back())/dlog)*dlog;
		log0Vector.push_back(log_val);
		
                // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
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
G4bool G4AdjointCSMatrix::GetData(unsigned int i, G4double& aLogPrimEnergy,G4double& aLogCS,G4double& log0, std::vector< double>*& aLogSecondEnergyVector,
	 							       std::vector< double>*& aLogProbVector, std::vector< size_t>*& aLogProbVectorIndex)
{	if (i>= nb_of_PrimEnergy) return false;
	//G4cout<<"Test Get Data "<<G4endl;
	aLogPrimEnergy = theLogPrimEnergyVector[i];
	aLogCS = theLogCrossSectionVector[i];
	aLogSecondEnergyVector = theLogSecondEnergyMatrix[i];
	aLogProbVector = theLogProbMatrix[i];
	aLogProbVectorIndex = theLogProbMatrixIndex[i];
	log0=log0Vector[i];
	return true;
	
}
///////////////////////////////////////////////////////
//
void G4AdjointCSMatrix::Write(G4String file_name)
{	std::fstream FileOutput(file_name, std::ios::out);			  
 	FileOutput<<std::setiosflags(std::ios::scientific);
 	FileOutput<<std::setprecision(6);
	FileOutput<<theLogPrimEnergyVector.size()<<G4endl;
	for (size_t i=0;i<theLogPrimEnergyVector.size();i++){
		FileOutput<<std::exp(theLogPrimEnergyVector[i])/MeV<<'\t'<<std::exp(theLogCrossSectionVector[i])<<G4endl;
		size_t j1=0;
		FileOutput<<theLogSecondEnergyMatrix[i]->size()<<G4endl;
		for (size_t j=0;j<theLogSecondEnergyMatrix[i]->size();j++){
			FileOutput<<std::exp((*theLogSecondEnergyMatrix[i])[j]);
			j1++;
			if (j1<10) FileOutput<<'\t';
			else  {
				FileOutput<<G4endl;
				j1=0;
			}	
		}
		if (j1>0) FileOutput<<G4endl;
		j1=0;
		FileOutput<<theLogProbMatrix[i]->size()<<G4endl;
		for (size_t j=0;j<theLogProbMatrix[i]->size();j++){
			FileOutput<<std::exp((*theLogProbMatrix[i])[j]);
			j1++;
			if (j1<10) FileOutput<<'\t';
			else  {
				FileOutput<<G4endl;
				j1=0;
			}	
		}
		if (j1>0) FileOutput<<G4endl;
		
		
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
		theLogSecondEnergyMatrix.push_back(new std::vector<G4double>());
		theLogProbMatrix.push_back(new std::vector<G4double>());
		
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
