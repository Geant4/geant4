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
// $Id: G4FluoDataTest.cc,v ????
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//
//      File name:     G4DataSetManagerTest
//
//      Author:        Elena Guardincerri
// 
//      Creation date: 6 August 2001
//
//      Modifications: 26 April 2002 -- AM
//
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4DataVector.hh"
#include "G4ios.hh" 
#include <fstream>
#include <iomanip>
#include "AIDA/AIDA.h"
#include "G4FluoData.hh"
#include <assert.h>


int main(int argc, char* argv[])

{
  G4cout.setf( std::ios::scientific, std::ios::floatfield );

  if (argc == 2) {
    G4int Z = atoi(argv[1]);
    G4FluoData* dataSet = new G4FluoData();
    dataSet->LoadData(Z);
    G4int vac= dataSet->NumberOfVacancies();
    AIDA::ITree* tree;
    AIDA::IAnalysisFactory* analysisFactory;
    AIDA::IHistogramFactory* cloudFactory;
    AIDA::ICloud1D* cloudFluo;

    AIDA::ITree* treeTuple;
    AIDA::ITupleFactory* tupleFactory;
    AIDA::ITuple* tupleFluo;

    analysisFactory = AIDA_createAnalysisFactory();
    AIDA::ITreeFactory* treeFactory = analysisFactory->createTreeFactory();
    G4String zString;
    std::ostringstream stream;
    stream << "fluorescence" << Z << ".xml";
    G4String fileName = stream.str();
    tree = treeFactory->create(fileName,"xml",false,true);
    cloudFactory = analysisFactory->createHistogramFactory(*tree);
    cloudFluo = cloudFactory->createCloud1D("c1");
    assert(cloudFluo);
    stream.str("");
    stream << "fluorescenceTuple" << Z << ".hbk";
    G4String fileNameTuple = stream.str();
    treeTuple = treeFactory->create(fileNameTuple,"hbook",false,true);
    tupleFactory = analysisFactory->createTupleFactory(*treeTuple);
    // Book tuple column names
    std::vector<std::string> columnNames;
    // Book tuple column types
    std::vector<std::string> columnTypes;

    columnNames.push_back("Z");
    columnNames.push_back("DataType");
    columnNames.push_back("ShellStart");
    columnNames.push_back("ShellStop");
    columnNames.push_back("ShellOrigAuger");
    columnNames.push_back("Energy");
    columnNames.push_back("Probability");
    columnNames.push_back("EnUnc");
    columnNames.push_back("Type");
    
    columnTypes.push_back("int");
    columnTypes.push_back("int");
    columnTypes.push_back("int");
    columnTypes.push_back("int");
    columnTypes.push_back("int");
    columnTypes.push_back("double");
    columnTypes.push_back("double");
    columnTypes.push_back("double");
    columnTypes.push_back("int");

    tupleFluo = tupleFactory->create("10", "Total Tuple", columnNames, columnTypes, "");
    assert(tupleFluo);

    for (G4int vacancyIndex=0; vacancyIndex < vac; vacancyIndex++) {
      G4int n = dataSet->NumberOfTransitions(vacancyIndex);
      for (G4int initIndex=0; initIndex < n; initIndex++){
	G4int id = dataSet->VacancyId(vacancyIndex);
	G4double startingShellEnergy = dataSet->StartShellEnergy(initIndex,vacancyIndex);
	G4int startingShellId = dataSet->StartShellId(initIndex,vacancyIndex);
	G4double startingShellProb = dataSet-> StartShellProb(initIndex,vacancyIndex);
	cloudFluo->fill(startingShellEnergy);

	tupleFluo->fill(0,Z);
	tupleFluo->fill(1,0);
	tupleFluo->fill(2,id);
	tupleFluo->fill(3,startingShellId);
	tupleFluo->fill(4,0);
	tupleFluo->fill(5,startingShellEnergy);
	tupleFluo->fill(6,startingShellProb);
	G4double error = 0;
	if (startingShellEnergy/eV < 100) {
	  error = startingShellEnergy*0.3;
	}
	else {error = startingShellEnergy*0.10;
	}
	tupleFluo->fill(7,error);
	tupleFluo->fill(8,0);

	tupleFluo->addRow();
      }
    }
	tree->commit(); // Write histos in file. 
	tree->close();
	treeTuple->commit(); // Write histos in file. 
	treeTuple->close();
    delete dataSet;      
  }
  else {

    G4cout << "Enter Z" << G4endl;
    G4int Z;
    G4cin >> Z;
  
    G4FluoData* dataSet = new G4FluoData();

    dataSet->LoadData(Z);
 
    G4int vac= dataSet->NumberOfVacancies();
    G4cout << "The atom of atomic number "<<Z<<" has "
	   << vac<<" vacancies "<<G4endl;
    G4cout << "Enter the index of the vacancy" << G4endl;
    G4int vacancyIndex;
    G4cin >> vacancyIndex;
  
 
    G4int n = dataSet->NumberOfTransitions(vacancyIndex);
    G4int id = dataSet->VacancyId(vacancyIndex);
    G4cout << " The shell whose index is " <<vacancyIndex  
	   << " has identity " << id << G4endl;
    G4cout <<" Electrons can reach it from "<< n <<" shells."<<G4endl;
    G4cout << "Enter the index of the starting shell" << G4endl;
    G4int initIndex;
    G4cin >>initIndex;
    G4int startingShellId = dataSet->StartShellId(initIndex,vacancyIndex);
    G4double startingShellEnergy = dataSet-> StartShellEnergy(initIndex,vacancyIndex);
    G4double startingShellProb = dataSet-> StartShellProb(initIndex,vacancyIndex);
    G4cout <<" The identity of the starting shell is "<<startingShellId<<G4endl;
    G4cout<<" The energy of the transition to the final shell is "
	  << startingShellEnergy<< " MeV "<<G4endl;
    G4cout<<" The probability of the transition to the final shell is "
	  <<startingShellProb <<G4endl;

    G4cout << "PRINT DATA"<<G4endl;

    dataSet->PrintData();
    delete dataSet;
  } 


G4cout << "END OF THE MAIN PROGRAM" << G4endl;
}













