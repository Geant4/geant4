// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4strstreambuf.hh,v 1.7 2000-11-20 17:26:47 gcosmo Exp $
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

class G4strstreambuf : public G4std::streambuf
{
  public:

    G4strstreambuf();
    ~G4strstreambuf();

    inline void SetDestination(G4coutDestination * value);

    inline G4int overflow(G4int c=EOF);
    inline G4int sync();
#ifdef WIN32
    inline G4int underflow();
#endif

    inline G4int ReceiveString ();

  private:

    G4strstreambuf(const G4strstreambuf&);
    G4strstreambuf& operator=(const G4strstreambuf&);

  private:

    G4coutDestination * destination;
    char* buffer;
    G4int count,size;
};

#include "G4strstreambuf.icc"

#endif
