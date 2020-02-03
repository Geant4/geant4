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
//
// 
// John Allison  5th April 2001
// A base class for a dummy viewer to dump geometry hierarchy.

#ifndef G4HepRepFileVIEWER_HH
#define G4HepRepFileVIEWER_HH

// #define G4HEPREPFILEDEBUG  // Comment this out to suppress debug code.

#include "G4VViewer.hh"

#include "G4HepRepFileXMLWriter.hh"
#include "G4HepRepFileSceneHandler.hh"

class G4HepRepFileViewer: public G4VViewer {
public:
  G4HepRepFileViewer(G4VSceneHandler&,const G4String& name);
  virtual ~G4HepRepFileViewer();
  void SetView();
  void ClearView();
  void DrawView();
  void ShowView();

private:
  G4HepRepFileXMLWriter *hepRepXMLWriter;
};

#endif
