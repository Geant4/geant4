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

// The data class which defined the configuration parameters for plotting 
// which can be modified via UI commands.
//
// Author: Ivana Hrivnacova, 21/10/2015  (ivana@ipno.in2p3.fr)

#ifndef G4PlotParameters_h
#define G4PlotParameters_h 1

#include "globals.hh"
#include "G4PlotMessenger.hh"

#include <memory>

class G4PlotParameters 
{
  public:
    G4PlotParameters();

    // setters
    void SetLayout(G4int columns, G4int rows); 
    void SetDimensions(G4int width, G4int height);
    void SetStyle(const G4String& style);

    // getters
    // limits
    G4int GetMaxColumns(); 
    G4int GetMaxRows();
    G4String GetAvailableStyles(); 
    // data
    G4int GetColumns() const;
    G4int GetRows() const;
    G4int GetWidth() const;
    G4int GetHeight() const;
    G4float  GetScale() const;
    G4String GetStyle() const;

  private: 
    // data members
    std::unique_ptr<G4PlotMessenger> fMessenger;
    // defaults
    G4int fDefaultColumns; 
    G4int fDefaultRows; 
    G4int fDefaultWidth;
    G4int fDefaultHeight;
    G4String fDefaultStyle;
    G4float fDefaultScale;
    // limits
    G4int fMaxColumns; 
    G4int fMaxRows; 
    G4String fAvailableStyles; 
    // data
    G4int fColumns; 
    G4int fRows; 
    G4int fWidth;
    G4int fHeight;
    G4float fScale;
    G4String fStyle;
};  

// inline functions

inline G4int G4PlotParameters::GetMaxColumns()
{ return fMaxColumns; }

inline G4int G4PlotParameters::GetMaxRows()
{ return fMaxRows; }

inline G4String G4PlotParameters::GetAvailableStyles()
{ return fAvailableStyles; }

inline G4int G4PlotParameters::GetColumns() const
{ return fColumns; }

inline G4int G4PlotParameters::GetRows() const
{ return fRows; }

inline G4int G4PlotParameters::GetWidth() const
{ return fWidth; }

inline G4int G4PlotParameters::GetHeight() const
{ return fHeight; }

inline G4float G4PlotParameters::GetScale() const
{ return fScale; }

inline G4String G4PlotParameters::GetStyle() const
{ return fStyle; }

#endif
