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

    
