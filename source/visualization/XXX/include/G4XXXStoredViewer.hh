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
// $Id: G4XXXStoredViewer.hh,v 1.1 2006-03-28 17:16:41 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  7th March 2006
// A template for a graphics driver with a store/database.
//?? Lines beginning like this require specialisation for your driver.

#ifndef G4XXXStoredVIEWER_HH
#define G4XXXStoredVIEWER_HH

#include "G4VViewer.hh"

class G4XXXStoredViewer: public G4VViewer {
public:
  G4XXXStoredViewer(G4VSceneHandler&,const G4String& name);
  virtual ~G4XXXStoredViewer();
  void SetView();
  void ClearView();
  void DrawView();
  void ShowView();
protected:
  void KernelVisitDecision ();
  G4bool CompareForKernelVisit(G4ViewParameters&);
  void DrawFromStore();
  G4ViewParameters fLastVP;  // Memory for making kernel visit decisions.
};

#endif
