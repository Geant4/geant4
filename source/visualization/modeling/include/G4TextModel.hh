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
// $Id: G4TextModel.hh,v 1.6 2006/06/29 21:31:39 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// 
// John Allison  3rd April 2001
//
// Class Description:
//
// Model which knows how to draw text.
//
// For access to base class information, e.g., modeling parameters,
// use GetModelingParameters() inherited from G4VModel.  See Class
// Description of the base class G4VModel.

#ifndef G4TEXTMODEL_HH
#define G4TEXTMODEL_HH

#include "G4VModel.hh"
#include "G4Text.hh"

class G4TextModel: public G4VModel {

public: // With description

  G4TextModel (const G4Text&);
   
  virtual ~G4TextModel ();

  virtual void DescribeYourselfTo (G4VGraphicsScene&);
  // The main task of a model is to describe itself to the graphics scene.

  virtual G4String GetCurrentDescription () const;
  // A description which depends on the current state of the model.

  virtual G4String GetCurrentTag () const;
  // A tag which depends on the current state of the model.

private:

  // Private copy contructor and assignment to forbid use...
  G4TextModel (const G4TextModel&);
  G4TextModel& operator = (const G4TextModel&);

  G4Text fG4Text;

};

#include "G4TextModel.icc"

#endif
