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
//---------------------------------------------------------------------------
//
// ClassName:  G4StoppingPhysicsFritiofWithBinaryCascade
//
// Author:     Alberto Ribon
//
// Date:       July 2019
//
// Modified:  
//
// Class Description:
//
// This class provides the nuclear capture at rest of negatively charged
// particles, using: Bertini for pi-, K-, Sigma-, Xi-, and Omega-; 
//                   Fritiof/Binary for anti-proton and anti-neutron;
//                   Fritiof/Precompound for anti-lambda, anti-sigma0, 
//                   anti-sigma+, anti-xi0 and anti-nuclei;
//                   another model for mu-.
//
//----------------------------------------------------------------------------

#ifndef G4StoppingPhysicsFritiofWithBinaryCascade_h
#define G4StoppingPhysicsFritiofWithBinaryCascade_h 1

#include "globals.hh"
#include "G4VPhysicsConstructor.hh"


class G4StoppingPhysicsFritiofWithBinaryCascade : public G4VPhysicsConstructor {
  public: 
    G4StoppingPhysicsFritiofWithBinaryCascade( G4int ver = 1 );
    G4StoppingPhysicsFritiofWithBinaryCascade( const G4String& name, G4int ver = 1,
                                               G4bool UseMuonMinusCapture = true );
    virtual ~G4StoppingPhysicsFritiofWithBinaryCascade();

    // This method will be invoked in the Construct() method. 
    // each particle type will be instantiated
    virtual void ConstructParticle();
 
    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type 
    virtual void ConstructProcess();

  private:
    G4int  verbose;
    static G4ThreadLocal G4bool wasActivated;
    G4bool useMuonMinusCapture;
};


#endif
