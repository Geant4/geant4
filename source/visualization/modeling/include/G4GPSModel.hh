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
// $Id: G4GPSModel.hh 99642 2016-09-29 19:49:30Z allison $
//
// 
// John Allison  26th April 2017.
//
// Class Description:
//
// Model for a representation of the General Paricle Source.

#ifndef G4GPSMODEL_HH
#define G4GPSMODEL_HH

#include "G4VModel.hh"

#include "G4String.hh"
#include "G4Colour.hh"

class G4VGraphicsScene;

class G4GPSModel: public G4VModel {

public: // With description

  G4GPSModel (const G4Colour& colour = G4Colour(1.,0.,0.,0.3));
  // Create with colour (default red and transparent).

  virtual ~G4GPSModel ();

  void DescribeYourselfTo (G4VGraphicsScene&);
  // The main task of a model is to describe itself to the graphics scene
  // handler (a object which inherits G4VSceneHandler, which inherits
  // G4VGraphicsScene).

  G4String GetCurrentDescription () const;
  // A description which depends on the current state of the model.

  G4String GetCurrentTag () const;
  // A tag which depends on the current state of the model.

protected:

  G4Colour fColour;

private:
  // Private copy constructor and assigment operator - copying and
  // assignment not allowed.  Keeps CodeWizard happy.
  G4GPSModel (const G4GPSModel&);
  G4GPSModel& operator = (const G4GPSModel&);
};

#endif
