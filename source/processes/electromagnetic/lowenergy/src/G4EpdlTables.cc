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

#include "CLHEP/Hist/HBookFile.h"
#include "CLHEP/Hist/TupleManager.h"
#include "CLHEP/Hist/HBookTuple.h"

// C++ Headers
#include <iostream.h>
#include <fstream.h>

//RW Headers
#include <rw/tpslist.h>

// Constructors
G4EpdlTables::G4EpdlTables(G4VDataFile& DFile):
  datfile(DFile)
{   
  theDataTable1 = 0;
  theDataTable2 = 0;
  theDataTable3 = 0;
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
  RWTPtrSlist<G4DataVector> vecList; 

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
	G4DataVector values; 
	
	datfile.GetDataValues(values); 
	lvl = values.length(); 

	if(!vecList.entries()){
	  
	  for(G4int k = 0; k < lvl; k++){
	  
	    vecList.insert(new G4DataVector);
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

} // end FillDataTable


RWTPtrSlist< RWTPtrSlist<G4DataVector> >* G4EpdlTables::GetGlobalList(){
  
  return new RWTPtrSlist< RWTPtrSlist<G4DataVector> >(allElementList);
}

void G4EpdlTables::FillTheTable(G4int numEl) {

  // line counters
  G4int numTable = 0;

  // variables to flag 68 characters lines
  G4bool lineMatch = FALSE;

  // list of data vectors to be filled
  RWTPtrSlist<G4DataVector> vecList; 

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
	G4DataVector values; 
	
	datfile.GetDataValues(values); 

	lvl = values.length(); 

	if(!vecList.entries()){

	  for(G4int k = 0; k < lvl; k++){
	
	    vecList.insert(new G4DataVector);
	    if(subSh){
	      vecList[k]->append(subSh);
	    }
	    else{

	      vecList[k]->append(numAtom);
	    }
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

	allElementList.append(new RWTPtrSlist<G4DataVector>(vecList));
	vecList.clear();
	numTable++;
	lineMatch = FALSE;
	if(numTable == 99){
	  break;
	}
      }
    }
  }// end for(;;)

} // end FillDataTable


void G4EpdlTables::MakeNtuple(const char* fname){


  G4PhysicsVector* outVec = 0;
  G4PhysicsVector* outVec2 = 0;
  G4PhysicsVector* outVec3 = 0;

  HepString name(fname);
  G4int N;

  if(theDataTable1->length()){
    
    HepTupleManager* tupleManager;
    HepString hbname = name + "1" + ".hbook";
    tupleManager = new HBookFile(hbname, 100);
    cout<<"HBOOK file name is: "<<((HBookFile*) tupleManager)->filename()<<endl;

    for(N = 0; N < 50; N++){

      HepString numElem(N);
      HepString tupleName = name + numElem;
      HepTuple* tuple = tupleManager->ntuple(tupleName);

      outVec = (*theDataTable1)(N);

      if(theDataTable2->length()){

        outVec2 = (*theDataTable2)(N);
      }

      if(theDataTable3->length()){

        outVec3 = (*theDataTable3)(N);
      }

      for(G4int i = 0; i < outVec->GetVectorLength(); i++){
	
        HepString arg = "val1" + numElem;
        tuple->column(arg, outVec->GetLowEdgeEnergy(i));

        HepString fval = "val2" + numElem;
        tuple->column(fval, (*outVec)(i));

        if(theDataTable2->length()){

          HepString fval2 = "val3" + numElem;
          tuple->column(fval2, (*outVec2)(i));

        }

        if(theDataTable3->length()){

          HepString fval3 = "val4" + numElem;
          tuple->column(fval3, (*outVec3)(i));
        }

        tuple->column("index",i);      
        tuple->dumpData();
        
      }
    }

    tupleManager->write();
    delete tupleManager;

//*****************************************************************************
 
    hbname = name + "2" +".hbook";
    tupleManager = new HBookFile(hbname, 100);
    cout<<"HBOOK file name is: "<<((HBookFile*) tupleManager)->filename()<<endl;

    for(N = 0; N < theDataTable1->length()-50; N++){

      HepString numElem(N);
      HepString tupleName = name + numElem;
      HepTuple* tuple = tupleManager->ntuple(tupleName);

      outVec = (*theDataTable1)(N+50);

      if(theDataTable2->length()){
        outVec2 = (*theDataTable2)(N+50);
      }

      if(theDataTable3->length()){
        outVec3 = (*theDataTable3)(N+50);
      }

      for(G4int i = 0; i < outVec->GetVectorLength(); i++){

        HepString arg = "val1" + numElem;
        tuple->column(arg, outVec->GetLowEdgeEnergy(i));
        HepString fval = "val2" + numElem;
        tuple->column(fval, (*outVec)(i));


        if(theDataTable2->length()){

          HepString fval2 = "val3" + numElem;
          tuple->column(fval2, (*outVec2)(i));
        }

        if(theDataTable3->length()){

          HepString fval3 = "val4" + numElem;
          tuple->column(fval3, (*outVec3)(i));
        }

        tuple->column("index",i);      
        tuple->dumpData();
      }
    }


    tupleManager->write();
    delete tupleManager;
  } // end if(DataTable1)
}

















     

   




