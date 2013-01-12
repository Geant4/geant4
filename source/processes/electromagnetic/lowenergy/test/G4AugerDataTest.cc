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
// $Id: G4AugerDataTest.cc,v ????
//
// -------------------------------------------------------------------
//
//      File name:     G4AugerDataTest
//
//      Author:        Alfonso Mantero (based on work By Elena Guardincerri)
// 
//      Creation date:  18 April 2002
//
//      Modifications: 
//
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4DataVector.hh"
#include "G4ios.hh"
#include <fstream>
#include <iomanip>
#include "AIDA/AIDA.h"
#include "G4AugerData.hh"
#include <assert.h>

int main(int argc, char* argv[])
{
  G4cout.setf( std::ios::scientific, std::ios::floatfield );

  if (argc == 2) {
    G4int Z = atoi(argv[1]);
    G4AugerData* dataSet = new G4AugerData();
    //    dataSet->LoadData(Z);
    G4int vac= dataSet->NumberOfVacancies(Z);

    AIDA::ITree* tree;
    AIDA::IAnalysisFactory* analysisFactory;
    AIDA::IHistogramFactory* cloudFactory;
    AIDA::ICloud1D* cloudAuger;

    AIDA::ITree* treeTuple;
    AIDA::ITupleFactory* tupleFactory;
    AIDA::ITuple* tupleAuger;

    analysisFactory = AIDA_createAnalysisFactory();
    AIDA::ITreeFactory* treeFactory = analysisFactory->createTreeFactory();
    G4String zString;
    std::ostringstream stream;
    stream << "auger" << Z << ".xml";
    G4String fileName = stream.str();
    tree = treeFactory->create(fileName,"xml",false,true);
    cloudFactory = analysisFactory->createHistogramFactory(*tree);
    cloudAuger = cloudFactory->createCloud1D("c1");
    assert(cloudAuger);

    stream.str("");
    stream << "augerTuple" << Z << ".hbk";
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

    tupleAuger = tupleFactory->create("10", "Total Tuple", columnNames, columnTypes, "");
    assert(tupleAuger);

  for (G4int vacancyIndex = 0; vacancyIndex<= vac-1; vacancyIndex++)
    {

      G4int n = dataSet->NumberOfTransitions(Z, vacancyIndex);
      G4int id = dataSet->VacancyId(Z, vacancyIndex);
      for (G4int initIndex = 0; initIndex < n; initIndex++){
       	G4int startingShellId = dataSet->StartShellId(Z, vacancyIndex, initIndex);
       	G4int nAuger = dataSet->NumberOfAuger(Z, vacancyIndex, startingShellId);

	for (G4int augerIndex = 0; augerIndex < nAuger; augerIndex++){

	G4double startingShellEnergy = dataSet-> StartShellEnergy(Z, vacancyIndex, startingShellId, augerIndex);
	G4int augerShellId = dataSet-> AugerShellId(Z, vacancyIndex, startingShellId, augerIndex);
	G4double startingShellProb = dataSet-> StartShellProb(Z, vacancyIndex, startingShellId, augerIndex);

	cloudAuger->fill(startingShellEnergy);
	tupleAuger->fill(0,Z);
	tupleAuger->fill(1,0);
	tupleAuger->fill(2,id);
	tupleAuger->fill(3,startingShellId);
	tupleAuger->fill(4,augerShellId);
	tupleAuger->fill(5,startingShellEnergy);
	tupleAuger->fill(6,startingShellProb);
	tupleAuger->fill(7,startingShellEnergy*0.15);
	tupleAuger->fill(8,0);
	tupleAuger->addRow();

	}
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
  
  G4AugerData* dataSet = new G4AugerData();

  G4cout << "G4AugerData created" << G4endl;
  
 
  G4int vac = dataSet->NumberOfVacancies(Z);

  G4cout << "The atom of atomic number "<<Z<<" has "
	 << vac<<" vacancies "<<G4endl;
  G4cout << "Enter the index of the main vacancy" << G4endl;
  G4int a, b;
  G4int vacancyIndex;
  G4cin >> a;

  if (a == -1)
    {
      a = 0;
      b = vac-1;
    }
  else { b = a;} 

  for (vacancyIndex = a; vacancyIndex<=b; vacancyIndex++)
    {

      G4int n = dataSet->NumberOfTransitions(Z, vacancyIndex);

      G4cout << " Testing VacancyId..."<< G4endl;


      G4int id = dataSet->VacancyId(Z, vacancyIndex);
      G4cout << " The shell whose index is " <<vacancyIndex  
	     << " has identity " << id << G4endl;
      G4cout <<" Electrons can reach it from "<< n <<" shells."<<G4endl;

      G4int a1 = 0;
      G4int nMax = 0;
      if (a == b) {
	G4cout << "Enter the index of the starting shell of the electron transition" << G4endl;
	G4cin >> a1;
       	nMax = n;
	n = a1+1;

	if (a1 >= nMax) G4cout << "max Index number must be less than number of available shells" << G4endl;
      }
      for (G4int initIndex = a1; initIndex < n; initIndex++){

	G4int startingShellId = dataSet->StartShellId(Z, vacancyIndex, initIndex);
	
	G4cout << " The shell whose index is " <<initIndex  
	       << " has identity " << startingShellId << G4endl;
	
	G4int nAuger = dataSet->NumberOfAuger(Z, vacancyIndex, startingShellId);
	
	G4int a2 = 0;
	if (a == b) {	
	  G4cout <<" Being a transition electron from here, an auger electron could came from  "
		 << nAuger <<" shells."<<G4endl;	  
	  G4cout << "Enter the index of the auger electron originating  shell" << G4endl;	
	  G4cin >> a2;
	  nMax=nAuger;
	  nAuger = a2 +1;
	  
	  if (a2 >= nMax) G4cout << ("max Index number must be less than number of available shells") << G4endl;
	}	

	for (G4int augerIndex = a2; augerIndex < nAuger; augerIndex++){

	G4cout << " Testing StartShellEnergy..."<< G4endl;
	
	G4double startingShellEnergy = dataSet-> StartShellEnergy(Z, vacancyIndex, startingShellId, augerIndex);
	
	G4cout << " Testing StartShellProb..."<< G4endl;
	
	G4double startingShellProb = dataSet-> StartShellProb(Z, vacancyIndex, startingShellId, augerIndex);
	G4int augerShellId = dataSet-> AugerShellId(Z, vacancyIndex, startingShellId, augerIndex);
	G4cout <<" The identity of the starting shell is "<<augerShellId<<G4endl;
	G4cout<<" The energy of the transition to the final shell is "
	      << startingShellEnergy<< " MeV "<<G4endl;
	G4cout<<" The probability of the transition to the final shell is "
	      <<startingShellProb <<G4endl;
	
	
	G4cout <<" The identity of the auger originating shell is "<<augerShellId<<G4endl;

	}
      }
    }
  
  /*
    G4cout << "PRINT DATA"<<G4endl;
    
    dataSet->PrintData(Z);
  */
  delete dataSet;
  }
 
  G4cout << "END OF THE MAIN PROGRAM" << G4endl;
}


