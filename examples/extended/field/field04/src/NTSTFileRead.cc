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
#include <iomanip>
#include <fstream>
#include "NTSTFileRead.hh"
#include "globals.hh"

NTSTFileRead::NTSTFileRead(const char* FileName, G4bool echo)
  : _LineLength(0), _Istr(0), _echo(0), stuff(0){
  _Istr = new std::ifstream(FileName);
  if (!*_Istr) {
    G4cerr << "Whoops! No such input file: " << FileName << G4endl;
  } else {
    G4cout << "Opened input file " << FileName << G4endl;
    _LineLength = sizeof _Line;
  }
  _echo=echo;
}

NTSTFileRead::~NTSTFileRead(){
  delete _Istr;
  delete stuff;
}

char* NTSTFileRead::ReadLine(){
  do {_Istr->getline(_Line, _LineLength-1);
  if (_echo) G4cout << _Line << G4endl;
  } while (_Line[0]=='#');
  return _Line;
}

std::istringstream &NTSTFileRead::StreamLine(){
  delete stuff;
  stuff = new std::istringstream(ReadLine());
  return *stuff;
}

    
