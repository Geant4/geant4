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
// $Id: G4HepRepFileViewer.hh,v 1.2 2001-11-19 15:06:52 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  5th April 2001
// A base class for a dummy viewer to dump geometry hierarchy.

#ifndef G4HepRepFileVIEWER_HH
#define G4HepRepFileVIEWER_HH

// #define G4HEPREPFILEDEBUG  // Comment this out to suppress debug code.

#include "G4VViewer.hh"

#include "HepRepXMLWriter.hh"

class G4HepRepFileViewer: public G4VViewer {
public:
  G4HepRepFileViewer(G4VSceneHandler&,const G4String& name);
  virtual ~G4HepRepFileViewer();
  void SetView();
  void ClearView();
  void DrawView();
  void ShowView();

private:
  HepRepXMLWriter *hepRepXMLWriter;
};

#endif
