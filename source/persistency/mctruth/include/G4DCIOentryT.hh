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
// File: G4DCIOentryT.hh
//
// Class Description:
//
// Template class of DigitsCollection I/O Manager for late binding

// Author: Youhei Morita, 12.09.2001
// --------------------------------------------------------------------
#ifndef G4DCIOENTRYT_HH
#define G4DCIOENTRYT_HH 1

#include <string>
#include "G4Types.hh"
#include "G4VPDigitsCollectionIO.hh"

// Class inherited:
#include "G4VDCIOentry.hh"

template <class T>
class G4DCIOentryT : public G4VDCIOentry
{
  public:

    G4DCIOentryT<T>(std::string n) : G4VDCIOentry(n)
    {
      if(m_verbose > 2)
      {
        G4cout << "G4DCIOentryT: Registering DigitsCollection IO manager"
               << " for \"" << n << "\"" << G4endl;
      }
    }
      // Constructor

    ~G4DCIOentryT() {}
      // Destructor

    void CreateDCIOmanager(const G4String& detName, const G4String& colName)
    {
      if(f_manager == nullptr)
      {
        f_manager = new T(detName, colName);
        if(m_verbose > 2)
        {
          G4cout << "G4DCIOentryT: Constructing DigitsCollection IO manager"
                 << " for \"" << detName << "\" " << f_manager << G4endl;
        }
        G4DCIOcatalog::GetDCIOcatalog()->RegisterDCIOmanager(f_manager);
        if(m_verbose > 2)
        {
          G4DCIOcatalog::GetDCIOcatalog()->PrintDCIOmanager();
        }
      }
    }
      // Create a new digits collection I/O manager

    void DeleteDCIOmanager()
    {
      if(f_manager != nullptr)  { delete f_manager; }
    }
      // Delete a digits collection I/O manager

  private:

    G4VPDigitsCollectionIO* f_manager = nullptr;

};

#endif
