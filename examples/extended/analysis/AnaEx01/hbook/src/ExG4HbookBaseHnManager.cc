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

/// \file hbook/src/ExG4HbookBaseHnManager.cc
/// \brief Implementation of the ExG4HbookBaseHnManager class

// Author: Ivana Hrivnacova, 03/11/2014  (ivana@ipno.in2p3.fr)

#include "ExG4HbookBaseHnManager.hh"

// static data

const G4int ExG4HbookBaseHnManager::kX = 0;
const G4int ExG4HbookBaseHnManager::kY = 1;
const G4int ExG4HbookBaseHnManager::kZ = 2;

//
// Constructors, destructor
//

//_____________________________________________________________________________
ExG4HbookBaseHnManager::ExG4HbookBaseHnManager(const G4String& hnType)
 : fHnType(hnType)
{
}

//_____________________________________________________________________________
ExG4HbookBaseHnManager::~ExG4HbookBaseHnManager()
{  
}

//_____________________________________________________________________________
G4int ExG4HbookBaseHnManager::GetNbins(const tools::hbook::axis& axis) const
{
  return axis.bins();
}  

//_____________________________________________________________________________
G4double ExG4HbookBaseHnManager::GetMin(const tools::hbook::axis& axis) const
{
// Returns min data value

  return axis.lower_edge();
}  

//_____________________________________________________________________________
G4double ExG4HbookBaseHnManager::GetMax(const tools::hbook::axis& axis) const
{
// Returns max data value

  return axis.upper_edge();
}  

//_____________________________________________________________________________
G4double ExG4HbookBaseHnManager::GetWidth(const tools::hbook::axis& axis) const
{
  G4int nbins = axis.bins();
  if ( ! nbins ) {
    G4String functionName = "ExG4HbookBaseHnManager::Get";
    functionName += fHnType;
    functionName += "Width";
    G4ExceptionDescription description;
    description << "    nbins = 0 (for " << fHnType << ").";
    G4Exception(functionName, "Analysis_W014", JustWarning, description);
    return 0;
  }              
  
  return ( axis.upper_edge() - axis.lower_edge() )/nbins;
}  

/*
//_____________________________________________________________________________
G4bool ExG4HbookBaseHnManager::SetTitle(G4HbookBaseHisto& baseHisto, 
                                    const G4String& title)
{
  return baseHisto.set_title(title);
}  
*/
//_____________________________________________________________________________
G4bool ExG4HbookBaseHnManager::SetAxisTitle(G4HbookBaseHisto& baseHisto, 
                                        G4int dimension, const G4String& title)
{
  if ( dimension == 0 ) {
    baseHisto.add_annotation(tools::hbook::key_axis_x_title(), title);
  }
  else if ( dimension == 1 ) {  
    baseHisto.add_annotation(tools::hbook::key_axis_y_title(), title);
  }
  else if ( dimension == 2 ) {  
    baseHisto.add_annotation(tools::hbook::key_axis_z_title(), title);
  }
  
  return true;
}  
/*
//_____________________________________________________________________________
G4String ExG4HbookBaseHnManager::GetTitle(const G4HbookBaseHisto& baseHisto) const
{
  return baseHisto.title();
}  
*/

//_____________________________________________________________________________
G4String ExG4HbookBaseHnManager::GetAxisTitle(const G4HbookBaseHisto& baseHisto, 
                                          G4int dimension) const 
{
  G4String title;
  G4bool result = false;
  if ( dimension == 0 ) {
    result = baseHisto.annotation(tools::hbook::key_axis_x_title(), title);
  }  
  else if ( dimension == 1 ) {  
    result = baseHisto.annotation(tools::hbook::key_axis_y_title(), title);
  }
  else if ( dimension == 2 ) {  
    result = baseHisto.annotation(tools::hbook::key_axis_z_title(), title);
  }

  if ( ! result ) {
    G4String axes("xyz");
    G4String axis = axes(dimension, 1);
    G4String functionName = "ExG4HbookBaseHnManager::Get";
    functionName += fHnType;
    functionName += axis;
    functionName += "Title";
    G4ExceptionDescription description;
    description << "    Failed to get " << axis << " axis " << fHnType << " title.";
    G4Exception(functionName, "Analysis_W014", JustWarning, description);
    return "";
  }
  
  return title;              
}  
