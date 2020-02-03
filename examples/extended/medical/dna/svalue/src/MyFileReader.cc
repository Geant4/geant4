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
// Derived from 
//  https://twiki.cern.ch/twiki/bin/view/Geant4/QuickMigrationGuideForGeant4V10
// Courtesy of A. Dotti
//
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software 
// shall cite the following Geant4-DNA collaboration publications:
// Med. Phys. 37 (2010) 4692-4708
// Phys. Med. 31 (2015) 861-874
// The Geant4-DNA web site is available at http://geant4-dna.org
//
/// \file medical/dna/svalue/src/MyFileReader.cc
/// \brief Implementation of the MyFileReader class

#include "MyFileReader.hh"
#include "G4ThreeVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MyFileReader::MyFileReader()
{ 
  //******************************************************************
  //*** Specify file name containing list of incident energies (in eV)
  //******************************************************************
  G4String fileName= "spectrum.txt";
  //
  
  fInputFile.open(fileName.data()); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MyFileReader::~MyFileReader()
{ fInputFile.close(); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double MyFileReader::GetAnEvent()
{
  
  //*******************************************
  //*** Specify number of lines to read at once
  //*******************************************
  G4int numberOfLinesToRead = 100;
  //
  
  if( fEvList.size() == 0 )
  {
    for(G4int i=0;i<numberOfLinesToRead;i++)
    {
      G4double nrj;
      fInputFile >> nrj;
      fEvList.push_back(nrj);
    }
  }

  //get first element from list
  G4double ev = fEvList.front();
  
  //remove this element from list
  fEvList.pop_front();
  
  return ev;
}
