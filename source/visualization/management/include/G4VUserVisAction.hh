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
// $Id: G4VUserVisAction.hh,v 1.1 2005-02-19 22:07:21 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

// Class Description:
//
// G4VUserVisAction is added to the scene by
// /vis/scene/add/userAction, which instantiates a G4UserActionModel
// to add to the RunDurationModels.  Thus the pure virtual Draw()
// method is invoked whenever the viewer need to "visit the kernel",
// i.e., to remake its graphical database, if any, or simply to
// refresh the screen.
// operator() allows G4CallbackModel<G4VUserVisAction>.


#ifndef G4VUSERVISACTION_HH
#define G4VUSERVISACTION_HH

class G4VUserVisAction {

public: // With description

  virtual void Draw() = 0;
  void operator()() {
    Draw();
  }

};

#endif
