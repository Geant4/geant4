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
//
// $Id: G4coutFormatters.hh 103582 2017-04-18 17:24:45Z adotti $
//
// 
// --------------------------------------------------------------------
// GEANT 4 header file 
//
// Class Description:
//
// Utilities to handle transformations of cout/cerr streams

//      ---------------- G4coutFormatters ----------------
//
// Author: A.Dotti (SLAC), April 2017
// --------------------------------------------------------------------
#ifndef G4COUTFORMATTERS_HH
#define G4COUTFORMATTERS_HH

#include <algorithm>
#include <sstream>
#include <vector>
#include <ctime>
#include <iomanip>
#include <functional>
#include <unordered_map>

#include "G4String.hh"
#include "G4ios.hh"
#include "G4strstreambuf.hh"

namespace G4coutFormatters
{
  // Static definitions of provided formatters
  namespace ID
  {
    static const G4String SYSLOG = "syslog";
    static const G4String DEFAULT= "default";
  }

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
  using SetupStyle_f = std::function<G4int(G4coutDestination*)>;

  using String_V=std::vector<G4String>;

  // Return list of formatter names
  String_V Names();

  // Setup style (by name) to destination
  G4int HandleStyle( G4coutDestination* dest , const G4String& style );

  // Set name of the style for the master thread
  void SetMasterStyle(const G4String& );
  G4String GetMasterStyle();

  // This function should be called in user application main function
  // to setup the style just after setting up RunManager
  void SetupStyleGlobally(const G4String& news);

  // To be used by user to register by name a new formatter.
  // So it can be used via one of the previous functions
  void RegisterNewStyle( const G4String& name , SetupStyle_f& formatter);
}

#endif // G4COUTFORMATTERS_HH
