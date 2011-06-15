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
#include "DicomRunAction.hh"
#include "DicomRun.hh"

//-- In order to obtain detector information.
#include <fstream>
#include <iomanip>
#include "G4THitsMap.hh"

#include "G4UnitsTable.hh"
//=======================================================================
// DicomRunAction
//  
//
//
//=======================================================================
// Constructor
DicomRunAction::DicomRunAction():
  FieldName(15),
  FieldValue(14)
{
  // - Prepare data member for DicomRun.
  //   vector represents a list of MultiFunctionalDetector names.
  theSDName.push_back(G4String("phantomSD"));
}

// Destructor.
DicomRunAction::~DicomRunAction()
{
  theSDName.clear();
}

//
//== 
G4Run* DicomRunAction::GenerateRun()
{
  // Generate new RUN object, which is specially
  // dedicated for MultiFunctionalDetector scheme.
  //  Detail description can be found in DicomRun.hh/cc.
  return new DicomRun(theSDName);
}

//
//==
void DicomRunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
}

//
//== 
void DicomRunAction::EndOfRunAction(const G4Run* aRun)
{
    G4cout << " ###### EndOfRunAction  " <<G4endl;
  //- DicomRun object.
  DicomRun* re02Run = (DicomRun*)aRun;
  //--- Dump all scored quantities involved in DicomRun.
  for ( G4int i = 0; i < (G4int)theSDName.size(); i++ ){
    //
    //---------------------------------------------
    // Dump accumulated quantities for this RUN.
    //  (Display only central region of x-y plane)
    //      0       ConcreteSD/DoseDeposit
    //---------------------------------------------
    G4THitsMap<G4double>* DoseDeposit = re02Run->GetHitsMap(theSDName[i]+"/DoseDeposit");

    G4cout << "=============================================================" <<G4endl;
    G4cout << " Number of event processed : "<< aRun->GetNumberOfEvent() << G4endl;
    G4cout << "=============================================================" <<G4endl;

    std::ofstream fileout;
    G4String fname = "dicom.out";
    fileout.open(fname);
    G4cout << " opened file " << fname << " for dose output" << G4endl;


    if( DoseDeposit && DoseDeposit->GetMap()->size() != 0 ) {
      std::ostream *myout = &G4cout;
      PrintHeader(myout);
      std::map<G4int,G4double*>::iterator itr = DoseDeposit->GetMap()->begin();
      for(; itr != DoseDeposit->GetMap()->end(); itr++) {
        fileout <<  itr->first
               << "     "  << *(itr->second)
               << G4endl;
	G4cout << "    " << itr->first
	       << "     " << std::setprecision(6) << *(itr->second) << " Gy"
	       << G4endl;
      }
      G4cout << "============================================="<<G4endl;
    }
    fileout.close();
    G4cout << " closed file " << fname << " for dose output" << G4endl;
  
  }
}
//
// --

void DicomRunAction::PrintHeader(std::ostream *out)
{
  std::vector<G4String> vecScoreName;
  vecScoreName.push_back("DoseDeposit");

  // head line
  //std::string vname = FillString("Volume", ' ', FieldName+1);
  //*out << vname << '|';
  std::string vname;
  *out << std::setw(10) << "Voxel" << " |";
  for (std::vector<G4String>::iterator it = vecScoreName.begin();
       it != vecScoreName.end(); it++) {
      //vname = FillString((*it),
//		       ' ', 
//		       FieldValue+1, 
//		       false);
//    *out << vname << '|';
      *out << std::setw(FieldValue) << (*it) << "  |";
  }
  *out << G4endl;  
}

std::string DicomRunAction::FillString(const std::string &name, 
				       char c, G4int n, G4bool back)
{
  std::string fname("");
  G4int k = n - name.size();
  if (k > 0) {
    if (back) {
      fname = name;
      fname += std::string(k,c);
    }
    else {
      fname = std::string(k,c);
      fname += name;
    }
  }
  else {
    fname = name;
  }
  return fname;
}
