// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PiMinusAbsorptionAtRest.hh,v 1.1 1999-01-07 16:13:39 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      File name:     G4PiMinusAbsorptionAtRest.hh 
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
// 
//      Creation date: 9 May 1998
//
//      Modifications: 
// -------------------------------------------------------------------

#ifndef G4PIMINUSABSORPTIONATREST_HH
#define G4PIMINUSABSORPTIONATREST_HH

#include "globals.hh"
#include "G4VRestProcess.hh"
#include "G4ParticleTypes.hh"
#include "G4PiMinusStopMaterial.hh"
#include "G4PiMinusStopAbsorption.hh"
#include "G4StopDeexcitation.hh"
#include "G4StopDeexcitationAlgorithm.hh"

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
  G4PiMinusAbsorptionAtRest(const G4String& processName ="PiMinusAbsorptionAtRest");

  // Destructor
  ~G4PiMinusAbsorptionAtRest();

  G4bool IsApplicable(const G4ParticleDefinition& particle) 
    { return ( particle == *(G4PionMinus::PionMinus()) ); }

  G4VParticleChange* AtRestDoIt(const G4Track& aTrack, const G4Step& aStep); 

  void SetDeexcitationAlgorithm(G4int index);

protected:                         

  // zero mean lifetime
  G4double GetMeanLifeTime(const G4Track& aTrack,
			   G4ForceCondition* condition) 
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

