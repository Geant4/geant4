// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// -------------------------------------------------------------------
//      GEANT 4 class file 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      CERN, Geneva, Switzerland
//
//      File name:     G4EpdlTables
//
//      Author:        Alessandra Forti (Alessandra.Forti@cern.ch)
// 
//      Creation date: 2 February 1999
//
//      Modifications: 
//      
// -------------------------------------------------------------------

// This class header
#include "G4EpdlTables.hh"

// Other Class Headers
#include "G4VDataFile.hh"
#include "G4DataVector.hh"
#include "G4PhysicsFreeVector.hh"
#include "CLHEP/String/Strings.h"

// C++ Headers
#include <iostream.h>
#include <fstream.h>

// Constructors
G4EpdlTables::G4EpdlTables(G4VDataFile& DFile):
  G4VTables(), 
  datfile(DFile)
{   
  theDataTable1 = 0;
  theDataTable2 = 0;
  theDataTable3 = 0;
  //  allElementList = 0;
}

// Destructor  
G4EpdlTables::~G4EpdlTables()
{

}

// Member Functions
void G4EpdlTables::FillDataTable() {

  // line counters
  G4int numTable = 0;

  // variables to flag 68 characters lines
  G4bool lineMatch = FALSE;

  // list of data vectors to be filled
  G4FirstLevel vecList; 

  G4int numBin = 100;

  if(theDataTable1){    
    
    theDataTable1->clearAndDestroy(); delete theDataTable1;
  }

  if(theDataTable2){    
    
    theDataTable2->clearAndDestroy(); delete theDataTable2;
  }

  if(theDataTable3){    
    
    theDataTable3->clearAndDestroy(); delete theDataTable3;
  }

  theDataTable1 = new G4PhysicsTable(numBin);
  theDataTable2 = new G4PhysicsTable(numBin);
  theDataTable3 = new G4PhysicsTable(numBin);

  //open input file
  datfile.OpenFile();

  // loop on the stream
  for(;;){

    datfile.Eof();

    datfile.GetLine();

    // lines counters
    G4int llength = datfile.LineLength();
    if(llength == 0) break;
  
    if(llength == 70){

      HepString AtomicNum(datfile.GetBuf());
      G4int numAtom = AtomicNum(0,3).toInt();
    }

    // search for the process flags line
    if(llength == 68 || llength == 69) {

	lineMatch = datfile.FindTheProcess();
	continue;
    }

    G4double lvl;

    if(llength < 68){

      if(lineMatch == TRUE){

	 //list of values in one line
	G4Data values; 
	
	datfile.GetDataValues(values); 
	lvl = values.length(); 

	if(!vecList.entries()){
	  
	  for(G4int k = 0; k < lvl; k++){
	  
	    vecList.insert(new G4Data);
	  }
	}

	for(G4int h = 0; h < lvl; h++){

	  vecList[h]->append(values[h]);
	}


	// Clear the temporary list
	values.clear();
      }
    }

    if(llength == 72 || llength == 73){

      // build the G4PhysicsTables

      if(lineMatch == TRUE){

	if(lvl >= 1){

	  G4PhysicsFreeVector* freevec;

	  freevec = new G4PhysicsFreeVector(*vecList[0],*vecList[1]);
	  theDataTable1->insertAt(numTable, freevec);

          if(lvl == 3){
	    
	    freevec = new G4PhysicsFreeVector(*vecList[0],*vecList[2]);
            theDataTable2->insertAt(numTable, freevec);
          }
	  
          if(lvl == 4){

	    freevec = new G4PhysicsFreeVector(*vecList[0],*vecList[2]);
            theDataTable2->insertAt(numTable, freevec);

	    freevec = new G4PhysicsFreeVector(*vecList[0],*vecList[3]);
            theDataTable3->insertAt(numTable, freevec);

	  }
	  
	}

	numTable++;
	lineMatch = FALSE;
	vecList.clear();

	if(numTable == 99){
	  
	  break;
	}
      }
    }
  }// end for(;;)

  if(theDataTable1->length() == 0){
    delete theDataTable1;
  }
  
  if(theDataTable2->length() == 0){
    delete theDataTable2;
  }
  
  if(theDataTable3->length() == 0){
    delete theDataTable3;
  }
  
} // end FillDataTable


//G4SecondLevel* G4EpdlTables::GetGlobalList(){
  
//return new G4SecondLevel((*allElementList));
//////}

G4SecondLevel* G4EpdlTables::FillTheTable(G4int numEl) {

  // line counters
  G4int numTable = 0;

  // variables to flag 68 characters lines
  G4bool lineMatch = FALSE;

  // list of data vectors to be filled
  G4FirstLevel* vecList = new G4FirstLevel();

  //  if(allElementList){ 
    
  //delete allElementList;
  //}

  G4SecondLevel* allElementList = new G4SecondLevel();

  //open input file
  datfile.OpenFile();

  // loop on the stream
  G4int subSh = 0;

  for(;;){

    datfile.Eof();

    datfile.GetLine();

    // lines counters
    G4int llength = datfile.LineLength();

    if(llength == 0) break;

    G4int numAtom;

    if(llength == 70){
      
      HepString AtomicNum(datfile.GetBuf());
      numAtom = AtomicNum(0,3).toInt();
      
    }

    // search for the process flags line
    if(llength == 68 || llength == 69) {

      if(numEl){
	
	if(numEl != numAtom){ 

	  continue;
	}

	else{

	  lineMatch = datfile.FindOneElemProc(subSh);
	}
      }
      else{

	lineMatch = datfile.FindTheProcess();
      }
      continue;
    }

    G4double lvl;

    if(llength < 68){

      if(lineMatch == TRUE){

	 //list of values in one line
	G4Data values; 
	
	datfile.GetDataValues(values); 

	lvl = values.length(); 

	if(!vecList->entries()){
	
	  for(G4int k = 0; k < lvl; k++){
		
	    vecList->insert(new G4Data);

	  }
	}

	for(G4int h = 0; h < lvl; h++){

	  (*vecList)[h]->insert(values[h]);
	}

	// Clear the temporary list
	values.clear();
      }
    }

    if(llength == 72 || llength == 73){

      // build the G4PhysicsTables
      if(lineMatch == TRUE){

	allElementList->insert(vecList);
	numTable++;
	lineMatch = FALSE;
	vecList = new G4FirstLevel();

	if(numTable == 99){
	  break;
	}
      }
    }
  }// end for(;;)
  return allElementList;
} // end FillDataTable






     

   




