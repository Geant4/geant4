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

// Author: Ivana Hrivnacova, 21/10/2015  (ivana@ipno.in2p3.fr)

#include "G4PlotParameters.hh"
#include "G4PlotMessenger.hh"
#include "G4AnalysisUtilities.hh"
#include "globals.hh"

using namespace G4Analysis;
using std::to_string;

//_____________________________________________________________________________
G4PlotParameters::G4PlotParameters()
 : fMessenger(nullptr),
#if defined(TOOLS_USE_FREETYPE)
   fDefaultStyle("ROOT_default"),
   fAvailableStyles("ROOT_default hippodrow inlib_default"),
#else
   fDefaultStyle("inlib_default"),
   fAvailableStyles("inlib_default"),
#endif
   fStyle(fDefaultStyle)
{
   fMessenger = std::make_unique<G4PlotMessenger>(this);
}

//
// private methods
//

//_____________________________________________________________________________
void G4PlotParameters::SetLayout(G4int columns, G4int rows)
{
  if ( columns > rows ||
       columns < 1 || columns > fkMaxColumns ||
       rows < 1 || rows > fkMaxRows ) {
    Warn("Layout: " + to_string(columns) + " x " + to_string(rows) +
         " was ignored.\n"
         "Supported layouts (columns <= rows): \n" +
         " columns = 1 .. " + to_string(fkMaxColumns) + "\n" +
         " rows    = 1 .. " + to_string(fkMaxRows),
         fkClass, "SetLayout");
    return;
  }
  fColumns = columns;
  fRows = rows;
}

//_____________________________________________________________________________
void G4PlotParameters::SetDimensions(G4int width, G4int height)
{
  fWidth = width;
  fHeight = height;
}

//_____________________________________________________________________________
void G4PlotParameters::SetStyle(const G4String& style)
{
// Set style and update scale according to the style selected

  if ( fAvailableStyles.find(style) == std::string::npos ) {
    Warn("Style: " + style + " was ignored.\n" +
         "Supported styles: " + fAvailableStyles,
         fkClass, "SetStyle");
    return;
  }

  fStyle = style;

  if ( fStyle == "ROOT_default" ) {
    fScale = fkDefaultScale;
  } else {
    fScale = 1.f;
  }
}
