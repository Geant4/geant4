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
//      File name:     G4EpdlData
//
//      Author:        Alessandra Forti (Alessandra.Forti@cern.ch)
// 
//      Creation date: 29 October 1998
//
//      Modifications: 
//      
// -------------------------------------------------------------------

// This class header
#include "G4EpdlData.hh"

// Other classes headers
#include "G4PhysicsFreeVector.hh"
#include "CLHEP/Hist/HBookFile.h"
#include "CLHEP/Hist/TupleManager.h"
#include "CLHEP/Hist/HBookTuple.h"
#include "G4DataVector.hh"

// C++ headers
#include <iostream.h>
#include <fstream.h>

// Constructors
G4EpdlData::G4EpdlData(HepString& dataFile, G4double* paramVec):
  _filename(dataFile),
  _userflags(paramVec)
{   

  theDataTable1 = 0;
  theDataTable2 = 0;
  theDataTable3 = 0;
}

// Destructor  
G4EpdlData::~G4EpdlData()
{

  if (theDataTable1) {
    
    //theDataTable1->clearAndDestroy(); 
    delete theDataTable1;
  }

  if (theDataTable2) {
    
    //theDataTable1->clearAndDestroy(); 
    delete theDataTable1;
  }

  if (theDataTable3) {

    //theDataTable1->clearAndDestroy(); 
    delete theDataTable1;
  }
}

// Member functions
void G4EpdlData::FillDataTable(G4double unit1, G4double unit2) {

  // open the stream
  //  HepString dir("data/");
  //  HepString dir_file = dir+_filename;

  ifstream istr(_filename, ios::in | ios::nocreate);

  if(!istr.is_open()){
    
    cout<<"G4EpdlData error: file "<<_filename<<" not found."<<endl;
    exit(1);
  }
 
  // create a buffer of character to read the stream
  char buf[72] = {0};
  G4int numBin = 100;

  // line counters
  G4int llenght = 0, numTable = 0;

  // data vectors
  G4DataVector dv0, dv1, dv2, dv3;

  // stream pointer position to override the EOF at the end of each table
  streampos pos71 = istr.tellg();

  // variables to flag 68 characters lines
  G4bool lineMatch = FALSE;

  if(theDataTable1){    

    theDataTable1->clearAndDestroy(); delete theDataTable1;
  }
  // 1 type of data
  theDataTable1 = new G4PhysicsTable(numBin);

  if(theDataTable2){    

    theDataTable2->clearAndDestroy(); delete theDataTable2;
  }

  if(_userflags[1] == 21 || _userflags[1] == 22 || _userflags[1] == 931){ // two types of data
    theDataTable2 = new G4PhysicsTable(numBin);
  }

  if(theDataTable3){    

    theDataTable3->clearAndDestroy(); delete theDataTable3;
  }
 
  if(_userflags[1] == 932){ //three types of data
    theDataTable3 = new G4PhysicsTable(numBin);
  }

  // loop on the stream
  for(;;){
     
    if(istr.eof()) {

      istr.close();
      break;
    }

    // by-pass the EOF at the end of each table
    if(pos71){

      istr.seekg(pos71);
      pos71 = 0;
    }

    istr.getline(buf, 72);
    
    // lines counters
    llenght = strlen(buf);

    // search for the first type of flags line
    if(llenght == 70){

      G4double par[4];
      CharToInt(buf, par);
      continue;
    }
    
    // search for the second type of flags line
    if(llenght == 68) {

      G4double par[4];
      CharToInt(buf, par);

      if(par[0]  == _userflags[0] && par[1]  == _userflags[1] && 
	 par[2]  == _userflags[2] && par[3]  == _userflags[3]) {
	
	numTable++;
	lineMatch = TRUE;
	continue;
      }

      else continue;
    }

    if(llenght < 68){
      
      if(lineMatch == TRUE){

	G4double values[4] = {0}; 
	GetDataValues(buf, values);
	if(values[0]) {
	  
	  values[0] *= unit1;
	  dv0.insert(values[0]); 
	}
	else continue;
	
	if(values[1]) {
	  
	  values[1] *= unit2;
	  dv1.insert(values[1]);
	}
	else continue;
	
	if(values[2]) {
	  
	  dv2.insert(values[2]);
	}
	else continue;
	
	if(values[3]) {
	  
	  dv3.insert(values[3]);
	}
	else continue;
      }
      else continue;
    }

    if(llenght == 71){

      // build the G4PhysicsFreeVector
      if(lineMatch == TRUE){

	if(!dv0.isEmpty()){

	  if(!dv1.isEmpty()){
	    G4PhysicsFreeVector* tmpVec = 0;

	    if(_userflags[1] == 941){
	    
	      tmpVec = new G4PhysicsFreeVector(dv1,dv0);
	    }
	    else{

	      tmpVec = new G4PhysicsFreeVector(dv0,dv1);
	    }

	    theDataTable1->insertAt(numTable-1, tmpVec);
	    dv1.clear();
	  }

	  if(!dv2.isEmpty()){

	    G4PhysicsFreeVector* tmpVec = new G4PhysicsFreeVector(dv0,dv2);
	    theDataTable2->insertAt(numTable-1, tmpVec);
	    dv2.clear();	
	  }

	  if(!dv3.isEmpty()){

	    G4PhysicsFreeVector* tmpVec = new G4PhysicsFreeVector(dv0,dv3);
	    theDataTable3->insertAt(numTable-1, tmpVec);
	    dv3.clear();
	  }
	  
	  dv0.clear();
	} 

	lineMatch = FALSE;
      }

      // make the PhysicsVector append it to the table, lineMatch = FALSE}
      pos71 = istr.tellg();
      istr.close();
      istr.open(_filename);
    }
    //    if(numTable == 100) cout<<"End Of the Table: "<<numTable<<endl;
  }// end while(getline)
} // end SearchDataTable


void G4EpdlData::GetDataValues(char* istrBuf, G4double* val){
  
  HepString parts, minus;
  char* token = 0;
  G4double floatTok = 0;
  HepString dec, dec2, tot;

  token = strtok(istrBuf," ");
  if(token){

    parts = token;
    minus = parts(8,1);
    if(minus == "-" || minus == "+"){
      dec = parts(0,8); // kind of process
      dec2 = parts(8,2); // kind of process
      tot = dec + "E" + dec2;
    }
    else{
      dec = parts(0,7); // kind of process
      dec2 = parts(7,3); // kind of process
      tot = dec + "E" + dec2;
    }
    floatTok = tot.toFloat();
    val[0] = floatTok;
    //    cout<<"HEPS 1: "<<dec<<", "<<dec2<<", "<<tot<<",  ";
    //    cout<<"ARGOMENTO 1: "<<floatTok<<endl;;
    //  if(token) val[0] = strtod(token,NULL); // first data column 
  }

  token=strtok(NULL," ");
  if(token){
    
    parts = token;
    minus = parts(8,1);
    if(minus == "-" || minus == "+"){
      dec = parts(0,8); // kind of process
      dec2 = parts(8,2); // kind of process
      tot = dec + "E" + dec2;
    }
    else{
      dec = parts(0,7); // kind of process
      dec2 = parts(7,3); // kind of process
      tot = dec + "E" + dec2;
    }
    floatTok = tot.toFloat();
    val[1] = floatTok; //initial value
    //    cout<<"HEPS 2: "<<dec<<", "<<dec2<<", "<<tot<<",  ";
    //    cout<<"ARGOMENTO 2: "<<floatTok<<endl;;
  }  

  token=strtok(NULL," ");
  if(token){

    parts = token;
    minus = parts(8,1);
    if(minus == "-" || minus == "+"){
      dec = parts(0,8); // kind of process
      dec2 = parts(8,2); // kind of process
      tot = dec + "E" + dec2;
    }
    else{
      dec = parts(0,7); // kind of process
      dec2 = parts(7,3); // kind of process
      tot = dec + "E" + dec2;
    }
    floatTok = tot.toFloat();
    val[2] = floatTok; //initial value
    //    cout<<"HEPS 3: "<<dec<<", "<<dec2<<", "<<tot<<",  ";
    //  cout<<"ARGOMENTO 3: "<<floatTok<<endl;;
  }
  
  token=strtok(NULL," ");
  if(token){

    parts = token;
    minus = parts(8,1);
    if(minus == "-" || minus == "+"){
      dec = parts(0,8); // kind of process
      dec2 = parts(8,2); // kind of process
      tot = dec + "E" + dec2;
    }
    else{
      dec = parts(0,7); // kind of process
      dec2 = parts(7,3); // kind of process
      tot = dec + "E" + dec2;
    }
    floatTok = tot.toFloat();
    val[3] = floatTok; //initial value
    // cout<<"HEPS 4: "<<dec<<", "<<dec2<<", "<<tot<<",  ";
    // cout<<"ARGOMENTO 4: "<<floatTok<<endl;;
  }
}

void G4EpdlData::CharToInt(char* linebuf, G4double* _localPar){
  
  G4int llenght = strlen(linebuf);
  //  cout<<"G4ED: Lunghezza del buffer: "<<llenght<<endl;

  if(llenght == 70){

    HepString flag(linebuf);
    HepString  Z = flag(0,3);
    //    cout<<"Elemento: "<<Z<<endl;
  }
  if(llenght == 68){

    HepString flag(linebuf);

    HepString C = flag(0,2); // kind of process
    _localPar[0] = C.toInt();

    HepString I = flag(2,3); // reaction property 
    _localPar[1] = I.toInt();

    HepString S = flag(5,3); // if subshell or full information
    _localPar[2] = S.toInt();
    // cout<<"G4ED: S parameter:"<<S<<endl;
    if(S == "  0"){

      //    HepString Xi = flag(22,3); // if subshell or full information
      _localPar[3] = 0;
    }

    else{

      HepString Xi1 = flag(22,1); // shells flags
      HepString Xi2 = flag(24,1); // shells flags
      HepString Xi3 = flag(31,1); // shells flags

      if(Xi3 == "0") _localPar[3] = Xi1.toInt();
      else if(Xi3 == "1"){
	
	HepString Xi = Xi1+Xi2;
	_localPar[3] = Xi.toInt();
      }
    }
    //    cout<<_localPar[0]<<" "<<_localPar[1]<<" "<<_localPar[2]<<" "<<_localPar[3]<<endl;
  }
}


void G4EpdlData::MakeNtuple(const char* fname){


  G4PhysicsVector* outVec = 0;
  G4PhysicsVector* outVec2 = 0;
  G4PhysicsVector* outVec3 = 0;

  HepString name(fname);
  G4int N;

  if(theDataTable1){
    
    HepTupleManager* tupleManager;
    HepString hbname = name + "1" + ".hbook";
    tupleManager = new HBookFile(hbname, 100);
    cout<<"HBOOK file name is: "<<((HBookFile*) tupleManager)->filename()<<endl;

    for(N = 1; N < 51; N++){

      HepString numElem(N);
      HepString tupleName = name + numElem;
      HepTuple* tuple = tupleManager->ntuple(tupleName);

      outVec = (*theDataTable1)(N-1);
      
      if(theDataTable2){
	outVec2 = (*theDataTable2)(N-1);
      }

      if(theDataTable3){
	outVec3 = (*theDataTable3)(N-1);
      }

      for(G4int i = 0; i < outVec->GetVectorLength(); i++){

	HepString arg = "val1" + numElem;
	tuple->column(arg, outVec->GetLowEdgeEnergy(i));

	HepString fval = "val2" + numElem;
	tuple->column(fval, (*outVec)(i));

	if(theDataTable2){

	  HepString fval2 = "val3" + numElem;
	  tuple->column(fval2, (*outVec2)(i));
	}

	if(theDataTable3){

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

    for(N = 51; N < 101; N++){

      HepString numElem(N);
      HepString tupleName = name + numElem;
      HepTuple* tuple = tupleManager->ntuple(tupleName);

      outVec = (*theDataTable1)(N-1);

      if(theDataTable2){
	outVec2 = (*theDataTable2)(N-1);
      }

      if(theDataTable3){
	outVec3 = (*theDataTable3)(N-1);
      }

      for(G4int i = 0; i < outVec->GetVectorLength(); i++){

	HepString arg = "val1" + numElem;
	tuple->column(arg, outVec->GetLowEdgeEnergy(i));
	HepString fval = "val2" + numElem;
	tuple->column(fval, (*outVec)(i));

	if(theDataTable2){

	  HepString fval2 = "val3" + numElem;
	  tuple->column(fval2, (*outVec2)(i));
	}

	if(theDataTable3){

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












