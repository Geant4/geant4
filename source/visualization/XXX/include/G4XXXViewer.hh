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
// $Id: G4XXXViewer.hh,v 1.7 2006-03-28 17:16:41 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  5th April 2001
// A template for a simplest possible graphics driver.
//?? Lines or sections marked like this require specialisation for your driver.

#ifndef G4XXXVIEWER_HH
#define G4XXXVIEWER_HH

//#define G4XXXDEBUG  // Comment this out to suppress debug code.

#include "G4VViewer.hh"

class G4XXXViewer: public G4VViewer {
public:
  G4XXXViewer(G4VSceneHandler&,const G4String& name);
  virtual ~G4XXXViewer();
  void SetView();
  void ClearView();
  void DrawView();
  void ShowView();
};

#endif
