// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LowEnergyUtilities.cc,v 1.8 2001-02-05 17:45:21 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// GEANT 4 class implementation file
//
// File name:     G4LowEnergyUtilitie
//
// Author:        A.Forti 
// 
// Creation date: 2 March 1999
//
// Modifications: 16/11/2000 MG Pia    Replaced HepString with G4String
//      
// --------------------------------------------------------------

// This Class Header
#include "G4LowEnergyUtilities.hh"

// Collaborating Class Headers
#include "G4Element.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "g4std/fstream"
#include "g4std/fstream"
#include "g4std/strstream"

G4LowEnergyUtilities::G4LowEnergyUtilities()
{}

G4LowEnergyUtilities::~G4LowEnergyUtilities()
{}

G4SecondLevel* G4LowEnergyUtilities::BuildSecondLevelTables(const G4int TableInd, 
						 const G4int ParNum, 
						 const char* prename){

  G4String prenameStr(prename);
  //  HepString name, prenameStr(prename);

  char nameChar[100] = {""};
  G4std::ostrstream ost(nameChar, 100, G4std::ios::out);

  if(TableInd != 0){
    ost << prenameStr << TableInd << ".dat";
    //    HepString Znum(TableInd);
    //    name = prenameStr + Znum + ".dat";
  }
  else{
    ost << prenameStr << ".dat";
    //   name = prenameStr+ ".dat";
  }
  G4String name(nameChar);

  char* path = getenv("G4LEDATA");
  if(!path){ 
    G4String excep = "Error!!! G4LEDATA (Low Energy Electromagnetic processes data directory) environment variable not set";
    G4Exception(excep);
  }

  G4String path_string(path);
  G4String dir_file = path_string + "/" + name;
  G4std::ifstream file(dir_file);
  G4std::filebuf* lsdp = file.rdbuf();

  if(!lsdp->is_open()){
    
      G4String excep = "Error!!!! data file: " + dir_file + " NOT found";
      G4Exception(excep);
  }
  
  oneAtomTable* oneAtomPar = new oneAtomTable();
  oneShellTable* oneShellPar = new oneShellTable();
  
  for(G4int j = 0; j < ParNum; j++){ 
    
    oneShellPar->insertAt(j,new G4DataVector());
  }
  
  G4double a = 0;
  G4int k = 1, s = 0;
  
  do{
    
    file>>a;
    
    if(a == -1){
      
      if(s == 0){
	
	oneAtomPar->insert(oneShellPar);
	oneShellPar = new oneShellTable();
	
	for(G4int j = 0; j < ParNum; j++){ 
	  
	  oneShellPar->insertAt(j,new G4DataVector());
	}
      }
      
      s++;
      
      if(s == ParNum){
	
	s = 0;
      }
    }

    else if(a == -2){
      
      delete oneShellPar;
    }

    else{
      
      if(k%ParNum != 0){	
	
	(*oneShellPar)[k-1]->push_back(a);
	k++;
      }
      else if(k%ParNum == 0){
	
	(*oneShellPar)[k-1]->push_back(a);
	k = 1;
      }
    }

  }while(a != -2); //end for on file
  
  file.close();
  return oneAtomPar;
}


G4FirstLevel* G4LowEnergyUtilities::BuildFirstLevelTables(const G4int TableInd, 
						const G4int ParNum, 
						const char* prename){

  G4String prenameStr(prename);
  //  HepString name, prenameStr(prename);

  char nameChar[100] = {""};
  G4std::ostrstream ost(nameChar, 100, G4std::ios::out);

  if(TableInd != 0){
    ost << prenameStr << TableInd << ".dat";
    //    HepString Znum(TableInd);
    //    name = prenameStr + Znum + ".dat";
  }
  else{
    ost << prenameStr << ".dat";
    //   name = prenameStr+ ".dat";
  }
  G4String name(nameChar);

  char* path = getenv("G4LEDATA");
  if(!path){ 
    G4String excep = "Error!!! G4LEDATA (Low Energy Electromagnetic processes data directory) environment variable not set";
    G4Exception(excep);
  }
  
  G4String path_string(path);
  G4String dir_file = path_string + "/" + name;
  G4std::ifstream file(dir_file);
  G4std::filebuf* lsdp = file.rdbuf();

  if(!lsdp->is_open()){
    
      G4String excep = "Error!!!! data file: " + dir_file + " NOT found";
      G4Exception(excep);
  }
  
  G4FirstLevel* oneAtomPar = new G4FirstLevel();
  
  for(G4int j = 0; j < ParNum; j++){ 
    
    oneAtomPar->insertAt(j,new G4DataVector());
  }
  
  G4double a = 0;
  G4int k = 1;
  
  do{
    
    file>>a;
    
    if(a == -1 || a == -2){

    }
    else{
      
      if(k%ParNum != 0){	
	
	(*oneAtomPar)[k-1]->push_back(a);
	k++;
      }
      else if(k%ParNum == 0){
	
	(*oneAtomPar)[k-1]->push_back(a);
	k = 1;
      }
    }

  }while(a != -2); //end for on file
  
  file.close();
  return oneAtomPar;
}








