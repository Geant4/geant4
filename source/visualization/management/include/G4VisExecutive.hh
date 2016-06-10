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
// $Id: G4VisExecutive.hh 66373 2012-12-18 09:41:34Z gcosmo $
//
// 
// John Allison 2nd February 2005 (based on MyVisManager, 24th January 1998).
//
// Class description
//
// Concrete Visualization Manager that implements the virtual
// functions RegisterGraphicsSystems and RegisterModelFactories.  This
// is executed when you Initialise() or Initialize() the vis manager.
// It exploits C-pre-processor variables G4VIS_USE_DAWN, etc., which
// are set by the GNUmakefiles if environment variables of the same
// name are set.
//
// Include this file and write code to instantiate G4VisExecutive just
// once at beginning of operations.  Before you compile, set
// appropriate environment variables (usually using "./Configure").
// If you change your environment you must force recompilation (the
// make files will not detect the need to do this).
//
// Typically, your main program file will contain:
//
// #ifdef G4VIS_USE
// #include "G4VisExecutive.hh"
// #endif
// ...
// int main() {
//   ...
// #ifdef G4VIS_USE
//   // Instantiate and initialise Visualization Manager.
//   G4VisManager* visManager = new G4VisExecutive;    // See Note (a).
//   visManager -> SetVerboseLevel (verbosityString);  // See Note (b).
//   visManager -> RegisterGraphicsSystem (new myGS);  // See Note (c).
//   visManager -> Initialize ();                      // See Note (d).
// #endif
//   ...
// #ifdef G4VIS_USE
//   G4cout << "Deleting vis manager..." << G4endl;
//   delete visManager;
//   G4cout << "Vis manager deleted." << G4endl;
// #endif
//
// Notes:
// (a) After instantiation, all references to this object should be as
//     a G4VisManager.  The functions RegisterGraphicsSystems and
//     RegisterModelFactories defined in G4VisExecutive.icc are
//     virtual functions of G4VisManager.  They are invoked by
//     G4VisManager::Initialise.  If you need to initialise in a
//     separate file, see advice below.
// (b) The verbosityString ("quiet", "errors", "warnings",
//     "confirmations", etc. - "help /vis/verbose" to see options) can be
//     set here or with /vis/verbose.  Alternatively, you can instantiate
//     with a verbosity string. e.g:
//       G4VisManager* visManager = new G4VisExecutive("quiet");
// (c) You can register your own graphics system like this.
// (d) Your can intialise like this with C++ code or use /vis/initialize.
//
// If you need to perform the instantiation and the initialisation in
// separate files, e.g., to establish the verbosity before
// initialisation, then the code that initialises must have access, of
// course, to the G4VisExecutive object, but this should be as a
// G4VisManager object, i.e., #include "G4VisManager.hh".
// RegisterGraphicsSystems and RegisterModelFactories are (pure)
// virtual methods of G4VisManager called from G4VisManager::Initialize.
// First file:
// #include "G4VisExecutive.hh"
// ...
//   fpVisManager = new G4VisExecutive;
// where fpVisManager is a G4VisManager*.
// Second file:
// #include "G4VisManager.hh"
// ...
//   fpVisManager -> Initialize ();
// where there is some mechanism for getting access to the pointer
// fpVisManager.
//
// The implementation is included as an .icc file because - for those
// graphics systems that need external libraries - only those systems
// that have been selected by the flags may be instantiated without
// causing unresolved references (only the user knows which libraries
// are available on his/her computer).  It also ensures that libraries
// can be linked in the right order, without circular dependencies.
// (Note that some graphics systems, notable those that write files
// for off-line viewing, do not suffer these restrictions and are
// always registered.)
//
// See class description of G4VisManager for more details.

#ifndef G4VISEXECUTIVE_HH
#define G4VISEXECUTIVE_HH

#include "G4VisManager.hh"

class G4VisExecutive: public G4VisManager {

public: // With description

  G4VisExecutive (const G4String& verbosityString = "warnings");

private:

  void RegisterGraphicsSystems();
  void RegisterModelFactories();

};

#include "G4VisExecutive.icc"

#endif
