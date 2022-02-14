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
// G4VDCIOentry
//
// Class Description:
//
// Abstract base class for digits collection I/O manager entry.

// Author: Youhei Morita, 12.09.2001
// --------------------------------------------------------------------
#ifndef G4VDCIOENTRYT_HH
#define G4VDCIOENTRYT_HH 1

#include "G4Types.hh"
#include "G4PersistencyCenter.hh"

class G4VDCIOentry
{
  public:

    G4VDCIOentry(const G4String& n);
      // Constructor

    virtual ~G4VDCIOentry() {}
      // Destructor

    void SetVerboseLevel(G4int v) { m_verbose = v; }
      // Sets verbose level

    const G4String& GetName() { return m_name; }
      // Returns the name of the DC I/O manager entry

    virtual void CreateDCIOmanager(const G4String&, const G4String&) {}
      // virtual method for creating DC I/O manager for the detector

  protected:

    G4int m_verbose = 0;

  private:

    G4String m_name;
};

#endif
