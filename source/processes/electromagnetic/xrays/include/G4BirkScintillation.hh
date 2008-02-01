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
// $Id: G4BirkScintillation.hh,v 1.1 2008-02-01 14:56:06 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
////////////////////////////////////////////////////////////////////////
// Scintillation Light Class Definition 
////////////////////////////////////////////////////////////////////////
//
// File:        G4BirkScintillation.hh  
// Description:	Discrete Process - Generation of Birk's law Scintillation Photons
// Version:     1.0
// Created:     01.02.2008
// Author:      V. Grichine based on G4Scintillation of Peter Gumplinger
// Updated:     
// mail:        Vladimir.Grichine@cern.ch
//
////////////////////////////////////////////////////////////////////////

#ifndef G4BirkScintillation_h
#define G4BirkScintillation_h 1

/////////////
// Includes
/////////////

#include "globals.hh"
#include "templates.hh"
#include "Randomize.hh"
#include "G4Poisson.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleMomentum.hh"
#include "G4Step.hh"
#include "G4Scintillation.hh"
#include "G4OpticalPhoton.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh" 
#include "G4PhysicsTable.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4PhysicsOrderedFreeVector.hh"

// Class Description:
// RestDiscrete Process - Generation of Scintillation Photons accordind to Birk's law.
// Class inherits publicly from G4Scintillation.
// Class Description - End:

/////////////////////
// Class Definition
/////////////////////

class G4BirkScintillation : public G4Scintillation
{

public: 

	G4BirkScintillation(const G4String& processName = "BirkScintillation",
                                 G4ProcessType type = fElectromagnetic);


	~G4BirkScintillation();	
        
  // Methods

	G4VParticleChange* PostStepDoIt(const G4Track& aTrack, 
			                const G4Step&  aStep);
        G4VParticleChange* AtRestDoIt (const G4Track& aTrack,
                                       const G4Step& aStep);

        // These are the methods implementing the scintillation process.



};


#endif /* G4BirkScintillation_h */
