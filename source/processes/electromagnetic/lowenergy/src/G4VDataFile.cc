// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// -------------------------------------------------------------------
//      GEANT 4 class file 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      CERN, Geneva, Switzerland
//
//      File name:     G4VDataFile
//
//      Author:        Alessandra Forti (Alessandra.Forti@cern.ch)
// 
//      Creation date: 2 February 1999
//
//      Modifications: 
//      
// -------------------------------------------------------------------

// This Class Header
#include "G4VDataFile.hh"

// Other Class Headers

// C++ Headers
#include "g4std/iostream"
#include "g4std/fstream"
#include <string.h>

// Constructors
G4VDataFile::G4VDataFile(const G4String& dataFile):
  _filename(dataFile)
{   
}

// Destructor  
G4VDataFile::~G4VDataFile()
{
  if(buf){
    delete [] buf;
  }
}
// Member Functions
void G4VDataFile::OpenFile(){
  
 // open the stream
  char* path = getenv("G4LEDATA");
  if(!path){ 

    G4Exception("G4LEDATA environment variable not set");
  }

  G4String path_string(path);
  G4String dir_file = path_string + "/" + _filename;
  #ifndef G4USE_STD_NAMESPACE
    _istr.open(dir_file.data(), ios::in | ios::nocreate);
  #else
    _istr.open(dir_file.data(), G4std::ios::in);
  #endif
  G4std::filebuf* lsdp = _istr.rdbuf();

  if(!lsdp->is_open()){

    G4String excep = "Error!!!! data file: " + dir_file + " NOT found";
    G4Exception(excep);
  }
}

void G4VDataFile::CloseFile(){

  _istr.close();
}

void G4VDataFile::Eof(){
  
  if(_istr.eof()) {
    
    _istr.close();
  
  }
}

G4std::streampos G4VDataFile::TellPos(){

  return _istr.tellg();
}

G4bool G4VDataFile::IsOpen(){

  return TRUE;//_istr.is_open();
}

void G4VDataFile::SeekPos(G4std::streampos pos){

  _istr.seekg(pos);
}

void G4VDataFile::SetBufferSize(G4int sz){

  _bufSize = sz;
  buf = new char[_bufSize+1];
  
}

void G4VDataFile::GetLine(){
  
  _istr.getline(buf, _bufSize);

  if(strlen(buf) >= _bufSize){

    G4String excep = "Error!!!! G4VDataFile::GetLine() buffer out of boundaries";
    G4Exception(excep);
  }
}

G4int G4VDataFile::LineLength(){

  return strlen(buf);

}

char* G4VDataFile::GetBuf(){

  return buf;
}








