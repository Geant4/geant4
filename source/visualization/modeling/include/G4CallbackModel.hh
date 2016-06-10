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
// $Id: G4CallbackModel.hh 66373 2012-12-18 09:41:34Z gcosmo $
//
// 
// John Allison  31st December 1997.
//
// Class Description:
//
// G4CallbackModel calls a user-defined function containing, for
// example, calls to the Draw methods of G4VVisManager or the
// *Primitive* methods of G4VGraphicsScene.  This allows the user to
// add his own graphical objects to a scene.  The idea is that the
// callback function is invoked by the vis system to refresh the
// screen or remake a graphical database whenever required.  The
// function class must contain:
//
//   void operator()(G4VGraphicsScene&, const G4Transform3D&)
//
// A base class G4VUserVisAction is provided in the visualisation
// category.  A G4CallbackModel is made by instantiating this template
// with a pointer to the function object.  The G4VisManager does this
// for the user if he/she registers it (SetUserAction) and issues the
// command /vis/scene/add/userAction.  See the User Guide for
// Application Developers, Section 8.8.7 and 8.8.8.

#ifndef G4CALLBACKMODEL_HH
#define G4CALLBACKMODEL_HH

#include "G4VModel.hh"

template <class F> class G4CallbackModel: public G4VModel {

 public:
  G4CallbackModel(F* function):
    fFunction(function) {}
  ~G4CallbackModel() {}
  void DescribeYourselfTo(G4VGraphicsScene& sceneHandler) {
    (*fFunction)(sceneHandler, fTransform);
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
