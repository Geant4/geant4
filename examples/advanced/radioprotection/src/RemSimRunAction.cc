//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
// ********************************************************************

//
// --------------------------------------------------------------
//                 GEANT 4 - RemSimtherapy example
// --------------------------------------------------------------
//
// Code developed by:
//  S.Guatelli
//
//
//    *******************************
//    *                             *
//    *    RemSimRunAction.cc       *
//    *                             *
//    *******************************
//
// $Id: RemSimRunAction.cc,v 1.1 2004-01-30 12:25:44 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "RemSimRunAction.hh"

#ifdef G4ANALYSIS_USE
#include "RemSimAnalysisManager.hh"
#endif

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4ios.hh"
#include "RemSimRunAction.hh"
#include <fstream>
#include <strstream>

RemSimRunAction::RemSimRunAction()
{
 //Read the input file concerning the generation of primary particles
 energies = new G4DataVector;
 data = new G4DataVector;
 ReadData(MeV,"extragalacticprotons");
}

RemSimRunAction::~RemSimRunAction()
{   
 }
void RemSimRunAction::BeginOfRunAction(const G4Run*)
{ 
#ifdef G4ANALYSIS_USE
  RemSimAnalysisManager* analysis = RemSimAnalysisManager::getInstance();
  analysis->book();
#endif  

 }

void RemSimRunAction::EndOfRunAction(const G4Run*)
{  
#ifdef G4ANALYSIS_USE      
  RemSimAnalysisManager* analysis = RemSimAnalysisManager::getInstance();
  analysis->finish();
#endif
}

void RemSimRunAction::ReadData(G4double unitE, G4String fileName)
{
  char nameChar[100] = {""};
  std::ostrstream ost(nameChar, 100, std::ios::out);
  
  ost << fileName <<".dat";
  
  G4String name(nameChar);
  
  //char* path = getenv("G4INSTALL");
  
  G4String pathString = "/afs/cern.ch/user/g/guatelli/REMSIM/N01";
  G4String dirFile = pathString + "/" + name;
  std::ifstream file(dirFile);
  std::filebuf* lsdp = file.rdbuf();
  
  if (! (lsdp->is_open()) )
    {
	  G4String excep = "RemSimRunAction - data file: " + dirFile + " not found";
	  G4Exception(excep);
    }
  G4double a = 0;
  G4int k = 1;
  
  do
    {
      file >> a;
      G4int nColumns = 2;
      // The file is organized into two columns:
      // 1st column is the energy
      // 2nd column is the corresponding value
      // The file terminates with the pattern: -1   -1
      //                                       -2   -2
      if (a == -1 || a == -2)
	{
	  
	}
      else
	{
	  if (k%nColumns != 0)
	    {	
	      G4double e = a * unitE;
	      G4cout<<e<<"energy"<<G4endl;
	      energies->push_back(e);
	      
	      k++;
	      
	    }
	  else if (k%nColumns == 0)
	    {
	      G4double value = a;
	      data->push_back(value);
	       G4cout<<a<<"flux"<<G4endl;
	      k = 1;
	    }
	}
      
    } while (a != -2); // end of file
  
  file.close();
}
G4DataVector* RemSimRunAction::GetPrimaryParticleEnergy()
{
  return energies;
}
G4DataVector* RemSimRunAction::GetPrimaryParticleEnergyDistribution()
{
  return data;
}
G4double RemSimRunAction::GetPrimaryParticleEnergyDistributionSum()
{
  G4double sum = 0;
  size_t size = data->size();
  for (size_t i = 0; i <size; i++)
    {
      sum+=(*data)[i];
    }
  return sum;
}

