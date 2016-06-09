//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4VisExecutive.hh,v 1.1 2005/02/04 16:42:37 johna Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// 
// John Allison 2nd February 2005 (based on MyVisManager, 24th January 1998).
//
// Class description
//
// Concrete Visualization Manager that implements the virtual function
//   RegisterGraphicsSystems.  This is executed when you Initialise()
//   or Initialize() the vis manager.  It exploits C-pre-processor
//   variables G4VIS_USE_DAWN, etc., which are set by the GNUmakefiles
//   if environment variables of the same name are set.
//
// Include this file and write code to instantiate G4VisExecutive just
//   once as beginning of operations.  Before you compile, set
//   appropriate environment variables.  If you change your
//   environment you must force recompilation (the make files will not
//   detect the need to do this).  Typically, your main program file
//   will contain:
//
// #ifdef G4VIS_USE
// #include "G4VisExecutive.hh"
// #endif
// ...
// int main() {
//   ...
// #ifdef G4VIS_USE
//   // Instantiate and initialise Visualization Manager.
//   G4VisManager* visManager = new G4VisExecutive;
//   //  visManager -> SetVerboseLevel (verbosityString);
//   //  visManager -> RegisterGraphicsSystem (new G4XXX);
//   visManager -> Initialize ();
// #endif
//   ...
// #ifdef G4VIS_USE
//   G4cout << "Deleting vis manager..." << G4endl;
//   delete visManager;
//   G4cout << "Vis manager deleted." << G4endl;
// #endif
//
// The implementation is included as an .icc file because - for those
//   graphics systems that need external libraries - only those
//   systems that have been selected by the flags may be instantiated
//   without causing unresolved references (only the user knows which
//   libraries are available on his/her computer).  It also ensures
//   that libraries can be linked in the right order, without circular
//   dependencies.  (Note that some graphics systems, notable those
//   that write files for off-line viewing, do not suffer these
//   restrictions and are always registered.)  Additional graphics
//   systems, XXX say, can be individually registered before
//   invocation of Initialise() with RegisterGraphicsSystem(new XXX).
//
// Alternatively, you can implement an empty function here and just
//   register the systems you want in your main(), e.g.:
//   G4VisManager* visManager = new G4VisExecutive;
//   visManager -> RegisterGraphicsSystem (new MyGraphicsSystem);
//
// See class description of G4VisManager for more details.

#ifndef G4VISEXECUTIVE_HH
#define G4VISEXECUTIVE_HH

#include "G4VisManager.hh"

class G4VisExecutive: public G4VisManager {

public: // With description

  G4VisExecutive ();

private:

  void RegisterGraphicsSystems ();

};

#include "G4VisExecutive.icc"

#endif

