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
// $Id: G4CallbackModel.hh,v 1.1 2005-02-19 22:07:21 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  31st December 1997.
//
// Class Description:
//
// G4CallbackModel calls a user-defined function containing, for example,
// calls to the Draw methods of G4VVisManager.  This allows the user to
// add his own graphical objects to a scene and they will be invoked
// appropriately to refresh the screen or remake a graphical database.

#ifndef G4CALLBACKMODEL_HH
#define G4CALLBACKMODEL_HH

// Doxygen generates separate HTML and LaTeX files for copying to the
// Geant4 User Guide documentation.
/**
   \page G4CallbackModelApplication G4CallbackModel Application Guide
   Here will be the documentation for the User Guide for
   Application Developers.
*/

#include "G4VModel.hh"

// Doxygen documentation...
/// A model that calls a user-defined function.
/** The user instantiates a function object containing, for example,
    calls to the Draw methods of G4VVisManager.  Typically, this would
    be invoked to remake a graphics database and/or refresh the
    viewer.  This allows the user to add his own graphical objects to
    a scene.  See \ref G4CallbackModelApplication for a more complete
    description and examples of its use. */

template <class F> class G4CallbackModel: public G4VModel {

 public:
  G4CallbackModel(F* function):
    fFunction(function) {}
  ~G4CallbackModel() {}
  void DescribeYourselfTo(G4VGraphicsScene&) {
    (*fFunction)();
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
