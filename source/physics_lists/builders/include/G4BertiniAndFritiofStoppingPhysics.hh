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
// $Id: G4BertiniAndFritiofStoppingPhysics.hh,v 1.6 2010-06-04 09:59:47 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:  G4BertiniAndFritiofStoppingPhysics
//
// Author:     Alberto Ribon
//
// Date:       27 July 2012
//
// Modified:  
//
// Class Description:
//
// This class provides the nuclear capture at rest of negatively charged
// particles, using: Bertini for pi-, K-, and Sigma-; 
//                   Fritiof/Precompound for anti-proton and anti-Sigma+;
//                   another (CHIPS-independent) model for mu-.
// No capture at rest is used for the other hadrons: Xi-, Omega-
// (this means that these hadrons can only decay).
//
//----------------------------------------------------------------------------

#ifndef G4BertiniAndFritiofStoppingPhysics_h
#define G4BertiniAndFritiofStoppingPhysics_h 1

#include "globals.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4MuonMinusCaptureAtRest.hh"


class G4HadronicAbsorptionBertini;
class G4HadronicAbsorptionFritiof;


class G4BertiniAndFritiofStoppingPhysics : public G4VPhysicsConstructor {

public: 

  G4BertiniAndFritiofStoppingPhysics( G4int ver = 1 );

  G4BertiniAndFritiofStoppingPhysics( const G4String& name,
		                      G4int ver = 1,
		                      G4bool UseMuonMinusCapture=true );

  virtual ~G4BertiniAndFritiofStoppingPhysics();

public: 

  // This method will be invoked in the Construct() method. 
  // each particle type will be instantiated
  virtual void ConstructParticle();
 
  // This method will be invoked in the Construct() method.
  // each physics process will be instantiated and
  // registered to the process manager of each particle type 
  virtual void ConstructProcess();

private:

  G4MuonMinusCaptureAtRest* muProcess;
  G4HadronicAbsorptionBertini* hBertiniProcess;
  G4HadronicAbsorptionFritiof* hFritiofProcess;
  
  G4int  verbose;
  G4bool wasActivated;
  G4bool useMuonMinusCaptureAtRest;
};


#endif
