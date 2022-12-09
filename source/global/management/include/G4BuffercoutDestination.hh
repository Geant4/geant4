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
// G4BuffercoutDestination
//
// Class Description:
//
// A cout destination that buffers received stream in a buffer and
// dumps the buffer only when full.
// Optionally allows for an infinite buffer.

//      ---------------- G4BuffercoutDestination ----------------
//
// Author: A.Dotti (SLAC), 14 April 2017
// --------------------------------------------------------------------
#ifndef G4BUFFERCOUTDESTINATION_HH
#define G4BUFFERCOUTDESTINATION_HH

#include <sstream>

#include "G4coutDestination.hh"

class G4BuffercoutDestination : public G4coutDestination
{
 public:
  explicit G4BuffercoutDestination(std::size_t maxSize = 0);
  ~G4BuffercoutDestination() override;

  G4int ReceiveG4cout(const G4String& msg) override;
  G4int ReceiveG4cerr(const G4String& msg) override;
  // Flush buffer to std output
  virtual G4int FlushG4cout();
  // Flush buffer to std error
  virtual G4int FlushG4cerr();
  // Flsuh both buffers

  virtual void Finalize();

  // Set maximum size of buffer, when buffer grows to specified size,
  // it will trigger flush. Dimension in char
  void SetMaxSize(std::size_t max) { m_maxSize = max; }
  std::size_t GetMaxSize() const { return m_maxSize; }
  std::size_t GetCurrentSizeOut() const { return m_currentSize_out; }
  std::size_t GetCurrentSizeErr() const { return m_currentSize_err; }

 protected:
  void ResetCout();
  void ResetCerr();

  std::ostringstream m_buffer_out;
  std::ostringstream m_buffer_err;
  std::size_t m_currentSize_out = 0;
  std::size_t m_currentSize_err = 0;
  std::size_t m_maxSize         = 0;
};

#endif
