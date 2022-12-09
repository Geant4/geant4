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
// G4strstreambuf
//
// Class description:
//
// Buffer for cout/cerr streaming

// Authors: H.Yoshida, M.Nagamatu - November 1998
// Revisions: G.Cosmo, 1998-2013
// --------------------------------------------------------------------
#ifndef G4STRSTREAMBUF_HH
#define G4STRSTREAMBUF_HH 1

#include <streambuf>

#include "G4coutDestination.hh"
#include "globals.hh"

class G4strstreambuf;

#ifdef G4MULTITHREADED

extern G4GLOB_DLL G4strstreambuf*& _G4coutbuf_p();
extern G4GLOB_DLL G4strstreambuf*& _G4cerrbuf_p();
#  define G4coutbuf (*_G4coutbuf_p())
#  define G4cerrbuf (*_G4cerrbuf_p())

#else  // Sequential

extern G4GLOB_DLL G4strstreambuf G4coutbuf;
extern G4GLOB_DLL G4strstreambuf G4cerrbuf;

#endif

class G4strstreambuf : public std::basic_streambuf<char>
{
 public:
  G4strstreambuf();
  ~G4strstreambuf() override;

  G4int overflow(G4int c = EOF) override;
  G4int sync() override;

#ifdef WIN32
  virtual G4int underflow();
#endif

  void SetDestination(G4coutDestination* dest);
  inline G4coutDestination* GetDestination() const;
  inline G4int ReceiveString();

 private:
  char* buffer = nullptr;
  G4int count = 0, size = 0;
  G4coutDestination* destination = nullptr;

  // hidden...
  G4strstreambuf(const G4strstreambuf&);
  G4strstreambuf& operator=(const G4strstreambuf&);
};

#include "G4strstreambuf.icc"

#endif
