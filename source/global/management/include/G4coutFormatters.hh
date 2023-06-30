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
// G4coutFormatters
//
// Description:
//
// Utilities to handle transformations of cout/cerr streams

//      ---------------- G4coutFormatters ----------------
//
// Author: A.Dotti (SLAC), April 2017
// --------------------------------------------------------------------
#ifndef G4COUTFORMATTERS_HH
#define G4COUTFORMATTERS_HH 1

#include <algorithm>
#include <ctime>
#include <functional>
#include <iomanip>
#include <sstream>
#include <unordered_map>
#include <vector>

#include "G4String.hh"
#include "G4ios.hh"
#include "G4coutDestination.hh"

namespace G4coutFormatters
{
  // Static definitions of provided formatters
  namespace ID
  {
    static const G4String SYSLOG  = "syslog";
    static const G4String DEFAULT = "default";
  }  // namespace ID

  using SetupStyle_f = std::function<G4int(G4coutDestination*)>;
  // A function that set ups a style for the destination
  // Example for a style that set to all capital the messages
  // to G4cerr:
  //  setupStyle_f myStyle = [](G4coutDestination* dest)->G4int {
  //            dest->SetCerrTransformer(
  //                         [](G4String& msg){
  //                            msg->toUpper();
  //                            return true; }
  //                    );
  //            };

  using String_V = std::vector<G4String>;

  String_V Names();
  // Return list of formatter names

  G4int HandleStyle(G4coutDestination* dest, const G4String& style);
  // Setup style (by name) to destination

  void SetMasterStyle(const G4String&);
  G4String GetMasterStyle();
  // Set/get name of the style for the master thread

  void SetupStyleGlobally(const G4String& news);
  // This function should be called in user application main function
  // to setup the style just after setting up RunManager

  void RegisterNewStyle(const G4String& name, SetupStyle_f& formatter);
  // To be used by user to register by name a new formatter.
  // So it can be used via one of the previous functions
}  // namespace G4coutFormatters

#endif
