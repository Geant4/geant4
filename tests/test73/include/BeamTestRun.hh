//
#ifndef BEAMTESTRUN_HH
#define BEAMTESTRUN_HH

#include "G4Event.hh"
#include "G4Run.hh"
#include "G4THitsMap.hh"
#include <map>

class G4Event;

class BeamTestRun : public G4Run 
{

	public:

		// Constructor 
		BeamTestRun(const G4String& detectorName);

		// Destructor
		virtual ~BeamTestRun();

	public:


		// Dump all data
		void DumpData() const;



};

#endif
