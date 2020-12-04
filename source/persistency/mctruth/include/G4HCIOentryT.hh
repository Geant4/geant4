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
// G4HCIOentryT
//
// Class Description:
//
// Template class of HitsCollection I/O Manager for late binding.

// Author: Youhei Morita, 12.09.2001
// --------------------------------------------------------------------
#ifndef G4HCIOENTRYT_HH
#define G4HCIOENTRYT_HH 1

#include <string>
#include "G4Types.hh"
#include "G4VPHitsCollectionIO.hh"

#include "G4VHCIOentry.hh"

template <class T>
class G4HCIOentryT : public G4VHCIOentry
{
  public:

    G4HCIOentryT<T>(const G4String& n) : G4VHCIOentry(n)
    {
      if(m_verbose > 2)
      {
        G4cout << "G4HCIOentryT: Registering HitsCollection IO manager"
               << " for \"" << n << "\"" << G4endl;
      }
    }
      // Constructor

    ~G4HCIOentryT() {}
      // Destructor

    void CreateHCIOmanager(const G4String& detName, const G4String& colName)
    {
      if(f_manager == nullptr)
      {
        f_manager = new T(detName, colName);
        if(m_verbose > 2)
        {
          G4cout << "G4HCIOentryT: Constructing HitsCollection IO manager"
                 << " for \"" << detName << "\" " << f_manager << G4endl;
        }
        G4HCIOcatalog::GetHCIOcatalog()->RegisterHCIOmanager(f_manager);
        if(m_verbose > 2)
        {
          G4HCIOcatalog::GetHCIOcatalog()->PrintHCIOmanager();
        }
      }
    }
      // Create a new hits collection I/O manager

    void DeleteHCIOmanager()
    {
      if(f_manager != nullptr)
        delete f_manager;
    }
      // Delete a hits collection I/O manager

  private:

    G4VPHitsCollectionIO* f_manager = nullptr;
};

#endif
