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
// $Id: G4VUserVisAction.hh,v 1.2 2005-02-23 11:30:10 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

class G4VUserVisAction {

public: // With description

  virtual void Draw() = 0;
  void operator()() {
    Draw();
  }

};

#endif
