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
// $Id$

/// \file hbook/include/ExG4HbookBaseHnManager.hh
/// \brief Definition of the ExG4HbookBaseHnManager class

// Author: Ivana Hrivnacova, 03/11/2014  (ivana@ipno.in2p3.fr)

#ifndef ExG4HbookBaseHnManager_h
#define ExG4HbookBaseHnManager_h 1

#include "globals.hh"

#include "tools/hbook/base_histo"
#include "tools/hbook/axis"

/// Manager class for tools::hbook::base_histo functions.
///
/// It implements functions common to all hbook histograms and profiles.

class ExG4HbookBaseHnManager //: public G4VH1Manager
{
  public:
    ExG4HbookBaseHnManager(const G4String& hnType);
    virtual ~ExG4HbookBaseHnManager();
    
    friend class ExG4HbookH1Manager;
    friend class ExG4HbookH2Manager;
    friend class ExG4HbookP1Manager;

  protected:
    typedef tools::hbook::base_histo G4HbookBaseHisto;

    // Access to datA parameters
    G4int    GetNbins(const tools::hbook::axis& axis) const;
    G4double GetMin(const tools::hbook::axis& axis) const;
    G4double GetMax(const tools::hbook::axis& axis) const;
    G4double GetWidth(const tools::hbook::axis& axis) const;

    // Attributes for plotting
    //

    // Setters
    //G4bool SetTitle(G4HbookBaseHisto& baseHisto, const G4String& title);
    G4bool SetAxisTitle(G4HbookBaseHisto& baseHisto, G4int dimension, const G4String& title);

    // Accessors
    //G4String GetTitle(const G4HbookBaseHisto& baseHisto) const;
    G4String GetAxisTitle(const G4HbookBaseHisto& baseHisto, G4int dimension) const;

    // Static data members
    static const G4int kX;
    static const G4int kY;
    static const G4int kZ;

  private:
    // Data members
    G4String  fHnType;
};

#endif

