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
//      File name:     G4Epdl97File
//
//      Author:        Alessandra Forti (Alessandra.Forti@cern.ch)
// 
//      Creation date: 2 February 1999
//
//      Modifications: 
//      
// -------------------------------------------------------------------

// This Class Header
#include "G4Epdl97File.hh"

//C++ Headers
#include "CLHEP/String/Strings.h"

// Constructors
G4Epdl97File::G4Epdl97File(const G4String& filename, G4int* paramVec):
  G4VDataFile(filename),
  _flags(paramVec)
{   
  SetBufferSize(74);
}

// Destructor  
G4Epdl97File::~G4Epdl97File()
{
}

G4bool G4Epdl97File::FindTheElement(G4int numZ){

  G4double llength = LineLength();
  G4bool elementFound = FALSE;
  HepString flag(GetBuf());

  if(numZ){
    if(llength == 70){
      
    }
  }
  return elementFound;
}

G4bool G4Epdl97File::FindTheProcess(){
  
  G4double llength = LineLength();
  G4bool tableFound = FALSE;
  HepString flag(GetBuf());

  if(llength == 68 || llength == 69){

    if(_flags[0] == flag(0,2).toInt()){

      if(_flags[1] == flag(2,3).toInt()){

	if(_flags[2] == flag(5,3).toInt()){
	  
	  G4int subsh;
	  G4int Xi3 = flag(31,1).toInt();
	  
	  if(Xi3 == 0){
	    
	    subsh = flag(22,1).toInt();
	  }
	  else if(Xi3 == 1){
	    
	    subsh = (flag(22,1) + flag(24,1)).toInt();
	  }
	  
	  if(_flags[3] == subsh){
	    
	    tableFound = TRUE;
	  }
	}
      }
    }
  }

  return tableFound;
}

G4bool G4Epdl97File::FindOneElemProc(G4int& subsh){
  
  G4double llength = LineLength();
  G4bool tableFound = FALSE;
  
  if(llength == 68){
    
    HepString flag(GetBuf());
    
    if(_flags[0] == flag(0,2).toInt()){
      
      if(_flags[1] == flag(2,3).toInt()){
	
	if(_flags[2] == flag(5,3).toInt()){
	  
	  G4int Xi3 = flag(31,1).toInt(); 
	  
	  if(Xi3 == 0){
	    
	    subsh = flag(22,1).toInt();
	  }
	  else if(Xi3 == 1){
	    
	    subsh = (flag(22,1) + flag(24,1)).toInt();
	  }
	  
	  tableFound = TRUE;
	}
      }
    }
  }

  return tableFound;
}

G4int* G4Epdl97File::GetTheProcFlags(){ return _flags; }

void G4Epdl97File::GetDataValues(G4Data& valList){

  char* token = 0;
  G4int i = 0;

  do{

   if(i == 0){

     token = strtok(GetBuf()," ");
   }
   else{

     token = strtok(NULL," ");
   }

   if(token) {

     valList.append(GetOneData(token));
   }

   i++;
  }while(token);
}

G4double G4Epdl97File::GetOneData(const char* token){

  HepString parts, tot;
  G4double floatTok = 0;
  
  if(token){

    parts = token;

    if(parts(8,1) == "-" || parts(8,1) == "+"){

      tot = parts(0,8) + "E" + parts(8,2);
    }
    else{

      tot = parts(0,7) + "E" + parts(7,3);
    }

    floatTok = tot.toFloat();
  }

  return floatTok;
}








