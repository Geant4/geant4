// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4strstreambuf.hh,v 1.6 1999-12-16 08:15:57 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ---------------------------------------------------------------
// GEANT 4 class header file
//
// G4strstrambuf.hh
//
// ---------------------------------------------------------------
#ifndef G4STRSTREAM_HH
#define G4STRSTREAM_HH

#include "g4std/strstream"
#include "globals.hh"     
#include "G4coutDestination.hh"

class G4strstreambuf;
extern G4strstreambuf G4coutbuf;
extern G4strstreambuf G4cerrbuf;

class G4strstreambuf : public G4std::streambuf {
public:
  G4strstreambuf() {
    destination = NULL;
    count = 0;
    size = 4095;
    buffer = new char[size+1];
  }
  ~G4strstreambuf() {
    delete buffer;
  }
  void SetDestination(G4coutDestination * value) {
    destination = value;
  }
  int overflow(int c=EOF) {
    int result = 0;
    if(count>=size) {
      buffer[count] = '\0';
      count = 0;
      result = ReceiveString ();
    }
    buffer[count] = c;
    count++;
    if(c=='\n') {
      buffer[count] = '\0';
      count = 0;
      result = ReceiveString ();
    }
    return result;
  }
  int sync() {
    buffer[count] = '\0';
    count = 0;
    return ReceiveString ();

  }
#ifdef WIN32
  int underflow() {
    return 0;
  }
#endif
  int ReceiveString () {
    G4String stringToSend = buffer;
    int result;
    if(this == & G4coutbuf && destination != NULL) {
      result =  destination->ReceiveG4cout(stringToSend);
    } else if(this == & G4cerrbuf && destination != NULL) {
      result =  destination->ReceiveG4cerr(stringToSend);
    } else if(this == & G4coutbuf && destination == NULL) {
      G4std::cout << stringToSend << G4std::flush;
      result =0;
    } else if(this == & G4cerrbuf && destination == NULL) {
      G4std::cerr << stringToSend << G4std::flush;
      result =0;
    }
    return result;
  };

private:
  G4coutDestination * destination;
  char* buffer;
  int count,size;
};


 
#endif
