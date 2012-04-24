//
#ifndef BEAMTESTPRIMARYGENERATORACTION_HH
#define BEAMTESTPRIMARYGENERATORACTION_HH

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include <cmath>
//#include "BeamTestParameters.hh"
//#include "BeamTestConversion.hh"

class G4ParticleGun;
class G4Event;

class BeamTestPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction 
{

	public:

		// Constructor
		BeamTestPrimaryGeneratorAction(/*Parameters* parameter*/);    

		// Destructor
		virtual ~BeamTestPrimaryGeneratorAction();

		// Method
		void GeneratePrimaries(G4Event*);

	private:

		// Data member
		G4ParticleGun* particleGun;	 
		
		//void Position(Parameters* parameter);

};

#endif
