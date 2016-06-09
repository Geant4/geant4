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
// $Id: G4PSHitsModel.hh,v 1.4 2009-11-04 12:44:39 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Created:  Mar. 31, 2009  Akinori Kimura  
//
// Model which knows how to draw GEANT4 scorers.
//

#ifndef G4SCORINGHITSMODEL_HH
#define G4SCORINGHITSMODEL_HH

#include "G4VModel.hh"

template <typename T> class G4THitsMap;

class G4PSHitsModel: public G4VModel {

public: // With description

  G4PSHitsModel (const G4String& requestedMapName = "all");
   
  virtual ~G4PSHitsModel ();

  virtual void DescribeYourselfTo (G4VGraphicsScene&);
  // The main task of a model is to describe itself to the graphics scene.

  const G4THitsMap<G4double>* GetCurrentHits() const
  {return fpCurrentHits;}

private:

  G4String fRequestedMapName;
  const G4THitsMap<G4double>* fpCurrentHits;
};

#endif
