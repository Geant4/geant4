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
// $Id: G4AugerDataTest.cc,v ????
// GEANT4 tag $Name: not supported by cvs2svn $
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
#include "g4std/fstream"
#include "g4std/iomanip"

#include "G4AugerData.hh"

int main()
{
  G4cout.setf( G4std::ios::scientific, G4std::ios::floatfield );

  G4cout << "Enter Z" << G4endl;
  G4int Z;
  G4cin >> Z;
  
  G4AugerData* dataSet = new G4AugerData();

  G4cout << "G4AugerData created" << G4endl;
  
 
  G4int vac = dataSet->NumberOfVacancies(Z);

  G4cout << "The atom of atomic number "<<Z<<" has "
	 << vac<<" vacancies "<<G4endl;
  G4cout << "Enter the index of the main vacancy" << G4endl;
  G4int vacancyIndex;
  G4cin >> vacancyIndex;

   G4int n = dataSet->NumberOfTransitions(Z, vacancyIndex);

  G4cout << " Testing VacancyId..."<< G4endl;


   G4int id = dataSet->VacancyId(Z, vacancyIndex);
   G4cout << " The shell whose index is " <<vacancyIndex  
	  << " has identity " << id << G4endl;
   G4cout <<" Electrons can reach it from "<< n <<" shells."<<G4endl;

   G4cout << "Enter the index of the starting shell of the electron transition" << G4endl;
   G4int initIndex;
   G4cin >>initIndex;

   G4int startingShellId = dataSet->StartShellId(Z, vacancyIndex, initIndex);

   G4cout << " The shell whose index is " <<initIndex  
	  << " has identity " << startingShellId << G4endl;

   G4int nAuger = dataSet->NumberOfAuger(Z, vacancyIndex, startingShellId);


   G4cout <<" Being a transition electron from here, an auger electron could came from  "
	  << nAuger <<" shells."<<G4endl;


   G4cout << "Enter the index of the auger electron originating  shell" << G4endl;
   G4int augerIndex;
   G4cin >>augerIndex;

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

   /*
   G4cout << "PRINT DATA"<<G4endl;
   
   dataSet->PrintData(Z);
   */
   delete dataSet;
   
   G4cout << "END OF THE MAIN PROGRAM" << G4endl;
}


