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

// Manager class for tools::histo::base_histo functions.
// It implements functions common to all histograms and profiles.
//
// Author: Ivana Hrivnacova, 24/07/2014  (ivana@ipno.in2p3.fr)

#ifndef G4BaseToolsManager_h
#define G4BaseToolsManager_h 1

#include "globals.hh"

#include "tools/histo/base_histo"

class G4BaseToolsManager //: public G4VH1Manager
{
  public:
    G4BaseToolsManager(const G4String& hnType);
    virtual ~G4BaseToolsManager();
    
    friend class G4H1ToolsManager;
    friend class G4H2ToolsManager;
    friend class G4H3ToolsManager;
    friend class G4P1ToolsManager;
    friend class G4P2ToolsManager;

  protected:
    typedef tools::histo::base_histo<double, unsigned int, unsigned int, double, double> 
      G4ToolsBaseHisto;

    // Access to datA parameters
    G4int    GetNbins(const G4ToolsBaseHisto& baseHisto, G4int dimension) const;
    G4double GetMin(const G4ToolsBaseHisto& baseHisto, G4int dimension) const;
    G4double GetMax(const G4ToolsBaseHisto& baseHisto, G4int dimension) const;
    G4double GetWidth(const G4ToolsBaseHisto& baseHisto, G4int dimension) const;

    // Attributes for plotting
    //

    // Setters
    G4bool SetTitle(G4ToolsBaseHisto& baseHisto, const G4String& title);
    G4bool SetAxisTitle(G4ToolsBaseHisto& baseHisto, G4int dimension, const G4String& title);

    // Accessors
    G4String GetTitle(const G4ToolsBaseHisto& baseHisto) const;
    G4String GetAxisTitle(const G4ToolsBaseHisto& baseHisto, G4int dimension) const;

    // Static data members
    static const G4int kX;
    static const G4int kY;
    static const G4int kZ;

  private:
    // Data members
    G4String  fHnType;
};

#endif

