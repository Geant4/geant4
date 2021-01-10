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
/// \file eventgenerator/HepMC/HepMCEx01/src/HepMC3G4AsciiReader.cc
/// \brief Implementation of the HepMC3G4AsciiReader class
//
//

#include "HepMC3G4AsciiReader.hh"
#include "HepMC3G4AsciiReaderMessenger.hh"

#include <iostream>
#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
HepMC3G4AsciiReader::HepMC3G4AsciiReader()
  :  filename("xxx.dat"), verbose(0)
{
  asciiInput= new HepMC3::ReaderAsciiHepMC2(filename.c_str());

  messenger= new HepMC3G4AsciiReaderMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
HepMC3G4AsciiReader::~HepMC3G4AsciiReader()
{
  delete asciiInput;
  delete messenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HepMC3G4AsciiReader::Initialize()
{
  delete asciiInput;

  asciiInput= new HepMC3::ReaderAsciiHepMC2(filename.c_str());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
HepMC3::GenEvent* HepMC3G4AsciiReader::GenerateHepMCEvent()
{
  HepMC3::GenEvent* evt= new HepMC3::GenEvent();
  asciiInput-> read_event(*evt);
  if (asciiInput->failed() ) { delete evt; return nullptr;}
  

  if(verbose>0) HepMC3::Print::content(*evt);

  return evt;
}
