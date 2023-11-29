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

#include "G4coutDestination.hh"

#include <iostream>

namespace
{
// Concrete streambuf redirecting output to G4coutDestination via G4cout etc
// Templated on two policy types to determine:
// - DestinationPolicy: which member member function of G4coutDestination to redirect to
// - DefaultPolicy: what to do if G4coutDestination is default (nullptr) 
template <typename DestinationPolicy, typename DefaultPolicy>
class G4strstreambuf : public std::basic_streambuf<char>
{
 public:
  G4strstreambuf()
  {
    size = 4095;
    buffer = new char[size + 1];
  }

  ~G4strstreambuf() override
  {
    // flushing buffer...
    // std::cout is used because destination object may not be alive.
    if (count != 0) {
      buffer[count] = '\0';
      std::cout << buffer;
    }
    delete[] buffer;
  }

  G4strstreambuf(const G4strstreambuf&) = delete;
  G4strstreambuf& operator=(const G4strstreambuf&) = delete;

  G4int overflow(G4int c = EOF) override
  {
    G4int result = 0;
    if (count >= size) result = sync();

    buffer[count] = (char)c;
    count++;

    return result;
  }

  G4int sync() override
  {
    buffer[count] = '\0';
    count = 0;
    return ReceiveString();
  }

#ifdef WIN32
  virtual G4int underflow() { return 0; }
#endif

  void SetDestination(G4coutDestination* dest) { destination = dest; }

  inline G4int ReceiveString()
  {
    G4String stringToSend(buffer);
    if (destination != nullptr) {
      return DestinationPolicy::PostMessage(destination, stringToSend);
    }
    return DefaultPolicy::PostMessage(stringToSend);
  }

 private:
  char* buffer = nullptr;
  G4int count = 0;
  G4int size = 0;
  G4coutDestination* destination = nullptr;
};

// Policies
struct PostToG4debug
{
  static inline G4int PostMessage(G4coutDestination* d, const G4String& s)
  {
    return d->ReceiveG4debug_(s);
  }
};

struct PostToG4cout
{
  static inline G4int PostMessage(G4coutDestination* d, const G4String& s)
  {
    return d->ReceiveG4cout_(s);
  }
};

struct PostToG4cerr
{
  static inline G4int PostMessage(G4coutDestination* d, const G4String& s)
  {
    return d->ReceiveG4cerr_(s);
  }
};

struct DefaultToCout
{
  static inline G4int PostMessage(const G4String& s)
  {
    std::cout << s << std::flush;
    return 0;
  }
};

struct DefaultToCerr
{
  static inline G4int PostMessage(const G4String& s)
  {
    std::cerr << s << std::flush;
    return 0;
  }
};

using G4debugstreambuf = G4strstreambuf<PostToG4debug, DefaultToCout>;
using G4coutstreambuf = G4strstreambuf<PostToG4cout, DefaultToCout>;
using G4cerrstreambuf = G4strstreambuf<PostToG4cerr, DefaultToCerr>;
}  // namespace

#ifdef G4MULTITHREADED
// --- StreamBuffers
G4debugstreambuf*& _G4debugbuf_p()
{
  G4ThreadLocalStatic auto* _instance = new G4debugstreambuf;
  return _instance;
}

G4coutstreambuf*& _G4coutbuf_p()
{
  G4ThreadLocalStatic auto* _instance = new G4coutstreambuf;
  return _instance;
}

G4cerrstreambuf*& _G4cerrbuf_p()
{
  G4ThreadLocalStatic auto* _instance = new G4cerrstreambuf;
  return _instance;
}

// --- Streams
std::ostream*& _G4debug_p()
{
  G4ThreadLocalStatic auto* _instance = new std::ostream(_G4debugbuf_p());
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

void G4iosInitialization()
{
  // --- Stream Buffers
  if (_G4debugbuf_p() == nullptr) {
    _G4debugbuf_p() = new G4debugstreambuf;
  }
  if (_G4coutbuf_p() == nullptr) {
    _G4coutbuf_p() = new G4coutstreambuf;
  }
  if (_G4cerrbuf_p() == nullptr) {
    _G4cerrbuf_p() = new G4cerrstreambuf;
  }

  // --- Streams
  if (_G4debug_p() == &std::cout || _G4debug_p() == nullptr) {
    _G4debug_p() = new std::ostream(_G4debugbuf_p());
  }
  if (_G4cout_p() == &std::cout || _G4cout_p() == nullptr) {
    _G4cout_p() = new std::ostream(_G4coutbuf_p());
  }
  if (_G4cerr_p() == &std::cerr || _G4cerr_p() == nullptr) {
    _G4cerr_p() = new std::ostream(_G4cerrbuf_p());
  }
}

void G4iosFinalization()
{
  // Reverse order
  // --- Streams
  delete _G4debug_p();
  _G4debug_p() = &std::cout;
  delete _G4cout_p();
  _G4cout_p() = &std::cout;
  delete _G4cerr_p();
  _G4cerr_p() = &std::cerr;

  // --- Stream Buffers
  delete _G4debugbuf_p();
  _G4debugbuf_p() = nullptr;
  delete _G4coutbuf_p();
  _G4coutbuf_p() = nullptr;
  delete _G4cerrbuf_p();
  _G4cerrbuf_p() = nullptr;
}

#  define G4debugbuf (*_G4debugbuf_p())
#  define G4coutbuf (*_G4coutbuf_p())
#  define G4cerrbuf (*_G4cerrbuf_p())

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

G4debugstreambuf G4debugbuf;
G4coutstreambuf G4coutbuf;
G4cerrstreambuf G4cerrbuf;

std::ostream G4debug(&G4debugbuf);
std::ostream G4cout(&G4coutbuf);
std::ostream G4cerr(&G4cerrbuf);

void G4iosInitialization() {}
void G4iosFinalization() {}

#endif

void G4iosSetDestination(G4coutDestination* sink)
{
  G4debugbuf.SetDestination(sink);
  G4coutbuf.SetDestination(sink);
  G4cerrbuf.SetDestination(sink);
}
