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
// $Id: G4VUserVisAction.hh 66373 2012-12-18 09:41:34Z gcosmo $
//
// 

// Class Description:
//
// G4VUserVisAction is added to the scene by the command
// /vis/scene/add/userAction, which instantiates a G4CallbackModel to
// add to the RunDurationModels.  Thus the pure virtual Draw() method
// is invoked whenever the viewer needs to "visit the kernel", e.g.,
// to remake its graphical database, if any, or simply to refresh the
// screen.  operator() is defined to satisfy the template
// G4CallbackModel<G4VUserVisAction>.
//
// See the User Guide for Application Developers, Section 8.8.7 and 8.8.8.


#ifndef G4VUSERVISACTION_HH
#define G4VUSERVISACTION_HH

#include "G4Transform3D.hh"

class G4VGraphicsScene;

class G4VUserVisAction
{
public: // With description
  G4VUserVisAction() {}
  virtual ~G4VUserVisAction() {}
  virtual void Draw() = 0;
  void operator()(G4VGraphicsScene&, const G4Transform3D&) {
    Draw();
  }
};

#endif
