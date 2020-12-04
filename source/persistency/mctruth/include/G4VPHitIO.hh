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
// G4VPHitIO
//
// Class Description:
//
// Abstract base class for storing and retrieving hit collections.

// Author: Youhei Morita, 10.08.2001
// --------------------------------------------------------------------
#ifndef G4VPHITIO_HH
#define G4VPHITIO_HH 1

#include "G4HCofThisEvent.hh"
#include "G4HCIOcatalog.hh"
#include "G4VPHitsCollectionIO.hh"

class G4VPHitIO
{
  public:

    G4VPHitIO();
      // Constructor

    virtual ~G4VPHitIO() {}
      // Destructor

    G4VPHitIO* GetVPHitIO() { return f_G4VPHitIO; }
      // Returns the pointer of the hit collection I/O manager

    virtual G4bool Store(const G4HCofThisEvent*) = 0;
      // Pure virtual method for storing hit collections of this event.
      // Each persistency package should implement a concrete method
      // of storing the hit collection of this event with this signature

    virtual G4bool Retrieve(G4HCofThisEvent*&) = 0;
      // Pure virtual method for retrieving hit collections of this event.
      // Each persistency package should implement a concrete method
      // of storing the hit collection of this event with this signature

    void SetVerboseLevel(G4int v);
      // Sets verbose level

  protected:

    void SetG4VPHitIO(G4VPHitIO* hitMan) { f_G4VPHitIO = hitMan; }
      // Registers the hit collection I/O manager

  protected:

    G4int m_verbose = 0;
    static G4ThreadLocal G4VPHitIO* f_G4VPHitIO;
    G4HCIOcatalog* f_catalog = nullptr;
};

#endif
