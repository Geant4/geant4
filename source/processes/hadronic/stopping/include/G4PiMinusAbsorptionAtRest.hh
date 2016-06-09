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
//      File name:     G4PiMinusAbsorptionAtRest.hh 
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
// 
//      Creation date: 9 May 1998
//
// -------------------------------------------------------------------

#ifndef G4PIMINUSABSORPTIONATREST_HH
#define G4PIMINUSABSORPTIONATREST_HH

// Class Description:
//
// Alternative Process for absorption of pi- at rest. 
// To be used in your physics list in case you need this physics.


#include "globals.hh"
#include "G4VRestProcess.hh"
#include "G4ParticleTypes.hh"
#include "G4PiMinusStopMaterial.hh"
#include "G4PiMinusStopAbsorption.hh"
#include "G4StopDeexcitation.hh"
#include "G4StopDeexcitationAlgorithm.hh"
#include "G4HadronicProcessType.hh"

class G4DynamicParticle;

class G4PiMinusAbsorptionAtRest : public G4VRestProcess
{  

private:

  // Hide assignment operator as private 
  G4PiMinusAbsorptionAtRest& operator=(const G4PiMinusAbsorptionAtRest &right);

  // Copy constructor
  G4PiMinusAbsorptionAtRest(const G4PiMinusAbsorptionAtRest& );

public:

  // Constructor
  G4PiMinusAbsorptionAtRest(const G4String& processName ="PiMinusAbsorptionAtRest", 
                       G4ProcessType   aType = fHadronic );

  // Destructor
  ~G4PiMinusAbsorptionAtRest();

  G4bool IsApplicable(const G4ParticleDefinition& particle) 
    { return ( particle == *(G4PionMinus::PionMinus()) ); }

  void PreparePhysicsTable(const G4ParticleDefinition&);

  void BuildPhysicsTable(const G4ParticleDefinition&);

  G4VParticleChange* AtRestDoIt(const G4Track& aTrack, const G4Step& aStep); 

  void SetDeexcitationAlgorithm(G4int index);

protected:                         

  // zero mean lifetime
  G4double GetMeanLifeTime(const G4Track& aTrack,
			   G4ForceCondition* ) 
  {
     G4double result = 0;
     if(aTrack.GetMaterial()->GetNumberOfElements() == 1)
        if(aTrack.GetMaterial()->GetZ()<1.5) result = DBL_MAX;
     return result;
   }
private:

  //  G4PiMinusStopAbsorption* _stopAbsorption;
  //  G4StopDeexcitation* _stopDeexcitation;
  G4int _indexDeexcitation;

  G4PiMinusStopMaterial* LoadAlgorithm(int Z);
  G4StopDeexcitationAlgorithm* LoadNucleusAlgorithm();
};
 

#endif

