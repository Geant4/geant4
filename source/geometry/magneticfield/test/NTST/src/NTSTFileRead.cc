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
#include <iomanip.h>
#include <fstream.h>
#include "NTSTFileRead.hh"
#include "globals.hh"
NTSTFileRead::NTSTFileRead(const char* FileName, G4bool echo)
  : _LineLength(0), _Istr(0), _echo(0), stuff(0){
  _Istr = new ifstream(FileName);
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

istrstream &NTSTFileRead::StreamLine(){
  //   return istrstream(ReadLine(), _LineLength);
  delete stuff;
  stuff = new istrstream( ReadLine(), _LineLength );
  return *stuff;
}

    
