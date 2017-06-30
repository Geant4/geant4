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
// $Id: G4BuffercoutDestination.hh 103582 2017-04-18 17:24:45Z adotti $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class header file 
//
// Class Description:
//
// A cout destination that buffers received stream in a buffer and
// dumps the buffer only when full.
// Optionally allows for an infinite buffer.

//      ---------------- G4BuffercoutDestination ----------------
//
// Author: A.Dotti (SLAC), April 2017
// --------------------------------------------------------------------
#ifndef G4BUFFERCOUTDESTINATION_HH_
#define G4BUFFERCOUTDESTINATION_HH_

#include <sstream>

#include "G4coutDestination.hh"

class G4BuffercoutDestination : public G4coutDestination
{
  public:

    explicit G4BuffercoutDestination(size_t maxSize = 0);
    virtual ~G4BuffercoutDestination();

    virtual G4int ReceiveG4cout(const G4String& msg) override;
    virtual G4int ReceiveG4cerr(const G4String& msg) override;
    // Flush buffer to std output
    virtual G4int FlushG4cout();
    // Flush buffer to std error
    virtual G4int FlushG4cerr();
    // Flsuh both buffers

    virtual void Finalize();

    // Set maximum size of buffer, when buffer grows to specified size,
    // it will trigger flush. Dimension in char
    void SetMaxSize(size_t max) { m_maxSize = max; }
    size_t GetMaxSize() const { return m_maxSize; }
    size_t GetCurrentSizeOut() const { return m_currentSize_out; }
    size_t GetCurrentSizeErr() const { return m_currentSize_err; }

  protected:

    void ResetCout();
    void ResetCerr();

    std::ostringstream m_buffer_out;
    std::ostringstream m_buffer_err;
    size_t m_currentSize_out;
    size_t m_currentSize_err;
    size_t m_maxSize;
};

#endif /* G4BUFFERCOUTDESTINATION_HH_ */
