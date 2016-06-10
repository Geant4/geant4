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

// Author: Ivana Hrivnacova, 21/10/2015  (ivana@ipno.in2p3.fr)

#include "G4PlotParameters.hh"
#include "G4PlotMessenger.hh"
#include "G4AnalysisUtilities.hh"
#include "globals.hh"

//_____________________________________________________________________________
G4PlotParameters::G4PlotParameters()
 : fMessenger(nullptr),
   fDefaultColumns(1),
   fDefaultRows   (2), 
#if defined(TOOLS_USE_FREETYPE)
   //Have vertical A4 :
   fDefaultWidth(2000), 
   fDefaultHeight((G4int)(29.7f/21.0f*fDefaultWidth)),
   fDefaultStyle("ROOT_default"),
   fDefaultScale(0.9f),
   fMaxColumns(3),
   fMaxRows   (5),
   fAvailableStyles("ROOT_default hippodrow inlib_default"), 
#else
   fDefaultWidth (700),
   fDefaultHeight((G4int)(29.7f/21.0f*fDefaultWidth)),
   fDefaultStyle("inlib_default"),
   fDefaultScale(0.9f),
   fMaxColumns(2),
   fMaxRows   (3),
   fAvailableStyles("inlib_default"),
#endif
   fColumns(fDefaultColumns), 
   fRows(fDefaultRows), 
   fWidth(fDefaultWidth), 
   fHeight(fDefaultHeight),
   fScale(fDefaultScale),
   fStyle(fDefaultStyle)
{
   fMessenger = G4Analysis::make_unique<G4PlotMessenger>(this);
}

//
// private methods
//

//_____________________________________________________________________________
void G4PlotParameters::SetLayout(G4int columns, G4int rows)
{
  if ( columns > rows ||
       columns < 1 || columns > fMaxColumns || 
       rows < 1 || rows > fMaxRows ) {
    G4ExceptionDescription description;
    description 
      << "Layout: " << columns << " x " << rows << " was ignored." << G4endl
      << "Suported layouts: " << G4endl
      << "  columns <= rows" << G4endl
      << "  columns = 1 .. " << fMaxColumns << G4endl
      << "  rows    = 1 .. " << fMaxRows << G4endl;
    G4Exception("G4PlotParameters::SetLayout",
                "Analysis_W013", JustWarning, description);
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
    G4ExceptionDescription description;
    description 
      << "Style: " << style << " was ignored." << G4endl
      << "Supported styles: " << fAvailableStyles << G4endl;
    G4Exception("G4PlotParameters::SetLayout",
                "Analysis_W013", JustWarning, description);
    return;
  }

  fStyle = style;

  if ( fStyle == "ROOT_default" ) {
    fScale = fDefaultScale;
  } else {
    fScale = 1.f;		
  }
}
