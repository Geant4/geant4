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
//
//
//    *******************************
//    *                             *
//    *    RemSimRunAction.cc       *
//    *                             *
//    *******************************
//
// Code developed by: S.Guatelli, guatelli@ge.infn.it
//
// $Id: RemSimRunAction.cc,v 1.10 2004/11/23 11:43:21 guatelli Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//

#include "RemSimRunAction.hh"
#include "RemSimDetectorConstruction.hh"
#include "RemSimPrimaryGeneratorAction.hh"
#ifdef G4ANALYSIS_USE
#include "RemSimAnalysisManager.hh"
#endif
#include "RemSimRunMessenger.hh"
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
 messenger = new RemSimRunMessenger(this);
 file = "";
}

RemSimRunAction::~RemSimRunAction()
{
  delete messenger;
  delete energies;
  delete data;  
 }

void RemSimRunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun -> GetRunID() << " start." << G4endl;
#ifdef G4ANALYSIS_USE
  G4int runNb = aRun -> GetRunID();
  if (runNb == 0) 
    {
     RemSimAnalysisManager* analysis = RemSimAnalysisManager::getInstance();
     analysis -> book();
    }
#endif
}

void RemSimRunAction::EndOfRunAction(const G4Run* aRun)
{  
 G4double numberEvents = aRun -> GetNumberOfEvent();
 G4cout<< "Number of events:" << numberEvents << G4endl;
}

void RemSimRunAction::Read(G4String name)
{ 
  file = name;   
  ReadData(MeV,name);
  G4cout << name << "  is the input file!" << G4endl;
}

void RemSimRunAction::ReadData(G4double unitE, G4String fileName)
{
  char nameChar[100] = {""};
  std::ostrstream ost(nameChar, 100, std::ios::out);
 
  ost << fileName;
  
  G4String name(nameChar);
  
  std::ifstream file(fileName);
  std::filebuf* lsdp = file.rdbuf();
  
  if (! (lsdp->is_open()) )
    {
	  G4String excep = "RemSimRunAction - data file: not found";
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
	      energies->push_back(e);  
	      //              G4cout<<e<<"energy";
	      
	      k++;
	      
	    }
	  else if (k%nColumns == 0)
	    {
	      G4double value = a;
	      data->push_back(value);
	      //G4cout<<" "<<a<<"flux"<<G4endl;
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
  size_t size = data -> size();
  for (size_t i = 0; i < size; i++)
    {
      sum+=(*data)[i];
    }
  return sum;
}

G4bool RemSimRunAction::GetFile()
{
  if (file == "") return false;
  else return true;
}
