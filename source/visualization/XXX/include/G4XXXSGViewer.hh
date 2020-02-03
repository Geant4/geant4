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
// John Allison  10th March 2006
// A template for a sophisticated graphics driver with a scene graph.
//?? Lines beginning like this require specialisation for your driver.

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
