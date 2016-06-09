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
// $Id: G4CallbackModel.hh,v 1.3 2005/03/03 16:22:02 allison Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// 
// John Allison  31st December 1997.
//
// Class Description:
//
// G4CallbackModel calls a user-defined function containing, for
// example, calls to the Draw methods of G4VVisManager.  This allows
// the user to add his own graphical objects to a scene.  The idea is
// that the callback function is invoked by the vis system to refresh
// the screen or remake a graphical database whenever required.
//
// The user instantiates a function object containing, for example,
// calls to the Draw methods of G4VVisManager.  A base class
// G4VUserVisAction is provided in the visualisation category.  A
// G4CallbackModel is made by instantiating this template with a
// pointer to the function object.  The G4VisManager does this for the
// user if he/she registers it (SetUserAction) and issues the command
// /vis/scene/add/userAction.  See the User Guide for Application
// Developers, Section 8.8.7 and 8.8.8.

#ifndef G4CALLBACKMODEL_HH
#define G4CALLBACKMODEL_HH

#include "G4VModel.hh"

template <class F> class G4CallbackModel: public G4VModel {

 public:
  G4CallbackModel(F* function):
    fFunction(function) {}
  ~G4CallbackModel() {}
  void DescribeYourselfTo(G4VGraphicsScene&) {
    (*fFunction)(fTransform);
  }

protected:

  F* fFunction;

private:

  // Private copy constructor and assigment operator - copying and
  // assignment not allowed.  Keeps CodeWizard happy.
  G4CallbackModel (const G4CallbackModel&);
  G4CallbackModel& operator = (const G4CallbackModel&);
};

#endif
