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
// $Id: G4XXXSGViewer.hh,v 1.1 2006-03-28 17:16:41 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  10th March 2006
// A template for a sophisticated graphics driver with a scene graph.
//?? Lines beginning like this require specialisation for your driver.

#ifdef G4VIS_BUILD_XXXSG_DRIVER

#ifndef G4XXXSGVIEWER_HH
#define G4XXXSGVIEWER_HH

#include "G4VViewer.hh"

class G4XXXSGViewer: public G4VViewer {
public:
  G4XXXSGViewer(G4VSceneHandler&,const G4String& name);
  virtual ~G4XXXSGViewer();
  void SetView();
  void ClearView();
  void DrawView();
  void ShowView();
protected:
  void KernelVisitDecision ();
  G4bool CompareForKernelVisit(G4ViewParameters&);
  void DrawFromStore(const G4String& source);
  G4ViewParameters fLastVP;  // Memory for making kernel visit decisions.
};

#endif

#endif
