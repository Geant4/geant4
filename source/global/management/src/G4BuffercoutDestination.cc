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
// G4BuffercoutDestination class implementation
//
// Author: A.Dotti (SLAC), 14 April 2017
// --------------------------------------------------------------------

#include "G4BuffercoutDestination.hh"

#include "G4AutoLock.hh"

#include <iostream>

// Private class to implement buffering of logging via an ostringstream
class G4BuffercoutDestination::BufferImpl
{
 public:
  using FlushFn_t = std::function<void(const std::string&)>;

 public:
  explicit BufferImpl(std::size_t maxSize) : m_maxSize(maxSize) {}
  explicit BufferImpl(std::size_t maxSize, FlushFn_t&& f) : m_maxSize(maxSize), m_flushFn(f) {}

  ~BufferImpl() = default;

  // Set number of characters to hold before Flush() will be called
  // If buffer exceeds new maximum, Flush() will not be called until next call to Receive()
  void SetMaxSize(std::size_t n) { m_maxSize = n; }

  // Reset buffer without flushing
  void Reset()
  {
    m_buffer.str("");
    m_buffer.clear();
    m_currentSize = 0;
  }

  G4int Receive(const G4String& msg)
  {
    m_currentSize += msg.size();
    m_buffer << msg;

    if (m_maxSize > 0 && m_currentSize > m_maxSize) {
      return Flush();
    }
    return 0;
  }

  // Flush buffer to destination and reset it
  G4int Flush()
  {
    m_flushFn(m_buffer.str());
    Reset();
    return 0;
  }

 private:
  std::size_t m_maxSize = 0;
  std::ostringstream m_buffer;
  std::size_t m_currentSize = 0;
  FlushFn_t m_flushFn = [](auto& s) { std::cout << s << std::flush; };
};

// --------------------------------------------------------------------
G4BuffercoutDestination::G4BuffercoutDestination(std::size_t max)
  : m_maxSize(max),
    m_buffer_dbg(std::make_unique<BufferImpl>(max)),
    m_buffer_out(std::make_unique<BufferImpl>(max)),
    m_buffer_err(std::make_unique<BufferImpl>(max, [](auto& s) { std::cerr << s << std::flush; }))
{}

// --------------------------------------------------------------------
G4BuffercoutDestination::~G4BuffercoutDestination() { Finalize(); }

// --------------------------------------------------------------------
void G4BuffercoutDestination::Finalize()
{
  FlushG4cerr();
  FlushG4cout();
  FlushG4debug();
}

// --------------------------------------------------------------------
G4int G4BuffercoutDestination::ReceiveG4debug(const G4String& msg)
{
  return m_buffer_dbg->Receive(msg);
}

// --------------------------------------------------------------------
G4int G4BuffercoutDestination::ReceiveG4cout(const G4String& msg)
{
  return m_buffer_out->Receive(msg);
}

// --------------------------------------------------------------------
G4int G4BuffercoutDestination::ReceiveG4cerr(const G4String& msg)
{
  return m_buffer_err->Receive(msg);
}

// --------------------------------------------------------------------
G4int G4BuffercoutDestination::FlushG4debug()
{
  return m_buffer_dbg->Flush();
}

// --------------------------------------------------------------------
G4int G4BuffercoutDestination::FlushG4cout()
{
  return m_buffer_out->Flush();
}

// --------------------------------------------------------------------
G4int G4BuffercoutDestination::FlushG4cerr()
{
  return m_buffer_err->Flush();
}

// --------------------------------------------------------------------
void G4BuffercoutDestination::SetMaxSize(std::size_t max)
{
  m_maxSize = max;
  m_buffer_dbg->SetMaxSize(m_maxSize);
  m_buffer_out->SetMaxSize(m_maxSize);
  m_buffer_err->SetMaxSize(m_maxSize);
}
