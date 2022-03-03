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

// The data class which defined the configuration parameters for plotting
// which can be modified via UI commands.
//
// Author: Ivana Hrivnacova, 21/10/2015  (ivana@ipno.in2p3.fr)

#ifndef G4PlotParameters_h
#define G4PlotParameters_h 1

#include "globals.hh"
#include "G4PlotMessenger.hh"

#include <memory>
#include <string_view>

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
    // Static data members
    static constexpr std::string_view fkClass { "G4PlotParameters" };
    static constexpr G4int fkHRDefaultWidth { 2000 };
    static constexpr G4int fkLRDefaultWidth { 700 };
    static constexpr G4int fkDefaultColumns { 1 };
    static constexpr G4int fkDefaultRows    { 2 };
    static constexpr G4float fkDefaultScale { 0.9f };
#if defined(TOOLS_USE_FREETYPE)
   //Have vertical A4 :
    static constexpr G4int fkDefaultWidth   { fkHRDefaultWidth };
    static constexpr G4int fkDefaultHeight  { static_cast<G4int>(29.7f/21.0f*fkHRDefaultWidth) };
    // limits
    static constexpr G4int fkMaxColumns { 3 };
    static constexpr G4int fkMaxRows    { 5 };
#else
   //Have vertical A4 :
    static constexpr G4int fkDefaultWidth   { fkLRDefaultWidth };
    static constexpr G4int fkDefaultHeight  { static_cast<G4int>(29.7f/21.0f*fkLRDefaultWidth) };
    // limits
    static constexpr G4int fkMaxColumns { 2 };
    static constexpr G4int fkMaxRows    { 3 };
#endif

    // Data members
    std::unique_ptr<G4PlotMessenger> fMessenger;
    // defaults
    G4String fDefaultStyle;
    G4String fAvailableStyles;
    // data
    G4int fColumns  { fkDefaultColumns };
    G4int fRows     { fkDefaultRows };
    G4int fWidth    { fkDefaultWidth };
    G4int fHeight   { fkDefaultHeight };
    G4float fScale  { fkDefaultScale };
    G4String fStyle;
};

// inline functions

inline G4int G4PlotParameters::GetMaxColumns()
{ return fkMaxColumns; }

inline G4int G4PlotParameters::GetMaxRows()
{ return fkMaxRows; }

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
