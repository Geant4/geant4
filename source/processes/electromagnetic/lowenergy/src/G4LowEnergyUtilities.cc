// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LowEnergyUtilities.cc,v 1.2 1999-09-28 13:15:45 aforti Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      ------------ G4LowEnergyUtilities physics process --------
//                   by Michel Maire, April 1996
// **************************************************************
// 12-06-96, Added SelectRandomAtom() method, by M.Maire
// 21-06-96, SetCuts implementation, M.Maire
// 17-09-96, PartialSumSigma(i)
//           split of ComputeBindingEnergy, M.Maire
// 08-01-97, crossection table + meanfreepath table, M.Maire
// 13-03-97, adapted for the new physics scheme, M.Maire
// 28-03-97, protection in BuildPhysicsTable, M.Maire
// 04-06-98, in DoIt, secondary production condition: range>min(threshold,safety)
// --------------------------------------------------------------
// This Class Header
#include "G4LowEnergyUtilities.hh"

// Collaborating Class Headers
#include "G4Element.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "CLHEP/String/Strings.h"
#include <fstream.h>

G4LowEnergyUtilities::G4LowEnergyUtilities()
{}

G4LowEnergyUtilities::~G4LowEnergyUtilities()
{}

G4SecondLevel* G4LowEnergyUtilities::BuildSecondLevelTables(const G4int TableInd, 
						 const G4int ParNum, 
						 const char* prename){

  
  HepString name, prenameStr(prename);
  if(TableInd != 0){

    HepString Znum(TableInd);
    name = prenameStr + Znum + ".dat";
  }
  else{

    name = prenameStr+ ".dat";
  }

  char* path = getenv("G4LEDATA");
  if(!path){ 
    HepString excep = "Error!!! G4LEDATA (Low Energy Electromagnetic processes data directory) environment variable not set";
    G4Exception(excep);
  }
  
  HepString path_string(path);
  HepString dir_file = path_string + "/" + name;
  ifstream file(dir_file);
  filebuf* lsdp = file.rdbuf();

  if(!lsdp->is_open()){
    
      HepString excep = "Error!!!! data file: " + dir_file + " NOT found";
      G4Exception(excep);
  }
  
  oneAtomTable* oneAtomPar = new oneAtomTable();
  oneShellTable* oneShellPar = new oneShellTable();
  
  for(G4int j = 0; j < ParNum; j++){ 
    
    oneShellPar->insertAt(j,new G4Data());
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
	  
	  oneShellPar->insertAt(j,new G4Data());
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
	
	(*oneShellPar)[k-1]->insert(a);
	k++;
      }
      else if(k%ParNum == 0){
	
	(*oneShellPar)[k-1]->insert(a);
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

  
  HepString name, prenameStr(prename);
  if(TableInd != 0){

    HepString Znum(TableInd);
    name = prenameStr + Znum + ".dat";
  }
  else{

    name = prenameStr+ ".dat";
  }

  char* path = getenv("G4LEDATA");
  if(!path){ 
    HepString excep = "Error!!! G4LEDATA (Low Energy Electromagnetic processes data directory) environment variable not set";
    G4Exception(excep);
  }
  
  HepString path_string(path);
  HepString dir_file = path_string + "/" + name;
  ifstream file(dir_file);
  filebuf* lsdp = file.rdbuf();

  if(!lsdp->is_open()){
    
      HepString excep = "Error!!!! data file: " + dir_file + " NOT found";
      G4Exception(excep);
  }
  
  G4FirstLevel* oneAtomPar = new G4FirstLevel();
  
  for(G4int j = 0; j < ParNum; j++){ 
    
    oneAtomPar->insertAt(j,new G4Data());
  }
  
  G4double a = 0;
  G4int k = 1;
  
  do{
    
    file>>a;
    
    if(a == -1 || a == -2){

    }
    else{
      
      if(k%ParNum != 0){	
	
	(*oneAtomPar)[k-1]->insert(a);
	k++;
      }
      else if(k%ParNum == 0){
	
	(*oneAtomPar)[k-1]->insert(a);
	k = 1;
      }
    }

  }while(a != -2); //end for on file
  
  file.close();
  return oneAtomPar;
}








