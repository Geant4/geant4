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
// G4ios implementation
//
// Authors: H.Yoshida, M.Nagamatu - November 1998
// --------------------------------------------------------------------

#include "G4ios.hh"
#include "G4strstreambuf.hh"
#include <iostream>

#ifdef G4MULTITHREADED

G4strstreambuf*& _G4coutbuf_p()
{
  G4ThreadLocalStatic auto* _instance = new G4strstreambuf();
  return _instance;
}

G4strstreambuf*& _G4cerrbuf_p()
{
  G4ThreadLocalStatic auto* _instance = new G4strstreambuf();
  return _instance;
}

std::ostream*& _G4cout_p()
{
  G4ThreadLocalStatic auto* _instance = new std::ostream(_G4coutbuf_p());
  return _instance;
}

std::ostream*& _G4cerr_p()
{
  G4ThreadLocalStatic auto* _instance = new std::ostream(_G4cerrbuf_p());
  return _instance;
}

#  define G4coutbuf (*_G4coutbuf_p())
#  define G4cerrbuf (*_G4cerrbuf_p())
#  define G4cout (*_G4cout_p())
#  define G4cerr (*_G4cerr_p())

void G4iosInitialization()
{
  if(_G4coutbuf_p() == nullptr)
  {
    _G4coutbuf_p() = new G4strstreambuf;
  }
  if(_G4cerrbuf_p() == nullptr)
  {
    _G4cerrbuf_p() = new G4strstreambuf;
  }
  if(_G4cout_p() == &std::cout || _G4cout_p() == nullptr)
  {
    _G4cout_p() = new std::ostream(_G4coutbuf_p());
  }
  if(_G4cerr_p() == &std::cerr || _G4cerr_p() == nullptr)
  {
    _G4cerr_p() = new std::ostream(_G4cerrbuf_p());
  }
}

void G4iosFinalization()
{
  delete _G4cout_p();
  _G4cout_p() = &std::cout;
  delete _G4cerr_p();
  _G4cerr_p() = &std::cerr;
  delete _G4coutbuf_p();
  _G4coutbuf_p() = nullptr;
  delete _G4cerrbuf_p();
  _G4cerrbuf_p() = nullptr;
}

// These two functions are guaranteed to be called at load and
// unload of the library containing this code.
namespace
{
#  ifndef WIN32
  void setupG4ioSystem() __attribute__((constructor));
  void cleanupG4ioSystem() __attribute__((destructor));
#  endif
  void setupG4ioSystem() { G4iosInitialization(); }
  void cleanupG4ioSystem() { G4iosFinalization(); }
}  // namespace

#else  // Sequential

G4strstreambuf G4coutbuf;
G4strstreambuf G4cerrbuf;
std::ostream G4cout(&G4coutbuf);
std::ostream G4cerr(&G4cerrbuf);

void G4iosInitialization() {}
void G4iosFinalization() {}

#endif
