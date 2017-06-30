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
// $Id: G4strstreambuf.hh 103661 2017-04-20 14:57:11Z gcosmo $
//
// ====================================================================
//
//   G4strstreambuf
//
// ====================================================================
#ifndef G4_STR_STREAM_BUF_HH
#define G4_STR_STREAM_BUF_HH

#include <streambuf>

#include "globals.hh"
#include "G4coutDestination.hh"

class G4strstreambuf;

#ifdef G4MULTITHREADED

  extern G4GLOB_DLL G4ThreadLocal G4strstreambuf *G4coutbuf_p;
  extern G4GLOB_DLL G4ThreadLocal G4strstreambuf *G4cerrbuf_p;
  #define G4coutbuf (*G4coutbuf_p)
  #define G4cerrbuf (*G4cerrbuf_p)

#else  // Sequential

  extern G4GLOB_DLL G4strstreambuf G4coutbuf;
  extern G4GLOB_DLL G4strstreambuf G4cerrbuf;

#endif

class G4strstreambuf : public std::basic_streambuf<char>
{
  public:

    G4strstreambuf();
    ~G4strstreambuf();
    
    virtual G4int overflow(G4int c=EOF);
    virtual G4int sync();

#ifdef WIN32
    virtual G4int underflow();
#endif

    void SetDestination(G4coutDestination* dest);
    inline G4coutDestination* GetDestination() const;
    inline G4int ReceiveString ();
  
  private:

    char* buffer;
    G4int count, size;
    G4coutDestination* destination;

    // hidden...
    G4strstreambuf(const G4strstreambuf&);
    G4strstreambuf& operator=(const G4strstreambuf&);
};

#include "G4strstreambuf.icc"

#endif
