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
// $Id: G4StoppingPhysics.hh 71041 2013-06-10 09:27:06Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:  G4StoppingPhysics
//
// Author:     Alberto Ribon
//
// Date:       27 July 2012
//
// Modified:  
// 20120921  M. Kelsey -- Move MuonMinusCapture.hh here; replace G4MMCAtRest
//		with new G4MuonMinusCapture.
// 16-Oct-2012 A. Ribon: renamed G4BertiniAndFritiofStoppingPhysics as 
//                       G4StoppingPhysics.
//
// Class Description:
//
// This class provides the nuclear capture at rest of negatively charged
// particles, using: Bertini for pi-, K-, Sigma-, Xi-, and Omega-. 
//                   Fritiof/Precompound for anti-proton and anti-Sigma+;
//                   another model for mu-.
//
//----------------------------------------------------------------------------

#ifndef G4StoppingPhysics_h
#define G4StoppingPhysics_h 1

#include "globals.hh"
#include "G4VPhysicsConstructor.hh"


class G4HadronicAbsorptionBertini;
class G4HadronicAbsorptionFritiof;
class G4MuonMinusCapture;


class G4StoppingPhysics : public G4VPhysicsConstructor {

public: 

  G4StoppingPhysics( G4int ver = 1 );

  G4StoppingPhysics( const G4String& name,
		                      G4int ver = 1,
		                      G4bool UseMuonMinusCapture=true );

  virtual ~G4StoppingPhysics();

public: 

  // This method will be invoked in the Construct() method. 
  // each particle type will be instantiated
  virtual void ConstructParticle();
 
  // This method will be invoked in the Construct() method.
  // each physics process will be instantiated and
  // registered to the process manager of each particle type 
  virtual void ConstructProcess();

private:

  //G4MuonMinusCapture* muProcess;
  //G4HadronicAbsorptionBertini* hBertiniProcess;
  //G4HadronicAbsorptionFritiof* hFritiofProcess;
  
  G4int  verbose;
  static G4ThreadLocal G4bool wasActivated;
  G4bool useMuonMinusCapture;
};


#endif
