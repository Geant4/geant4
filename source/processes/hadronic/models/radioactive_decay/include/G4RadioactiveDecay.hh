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
#ifndef G4RadioactiveDecay_h
#define G4RadioactiveDecay_h 1
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4RadioactiveDecay.hh
//
// Version:             0.b.4
// Date:                14/04/00
// Author:              F Lei & P R Truscott
// Organisation:        DERA UK
// Customer:            ESA/ESTEC, NOORDWIJK
// Contract:            12115/96/JG/NL Work Order No. 3
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 29 February 2000, P R Truscott, DERA UK
// 0.b.3 release.
//
// 13 April 2000, F Lei, DERA UK
// 0.b.4 release. No change to this file     
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
////////////////////////////////////////////////////////////////////////////////
//
#include "G4ios.hh"
#include "globals.hh"
#include "G4VRestDiscreteProcess.hh"
#include "G4ParticleChangeForRadDecay.hh"
#include "G4RadioactiveDecaymessenger.hh"  

#include "G4NucleusLimits.hh"
#include "G4RadioactiveDecayRate.hh"
#include "G4RadioactiveDecayRateVector.hh"
#include "G4RIsotopeTable.hh"
#include "G4RadioactivityTable.hh"

#include <vector>


//class G4UserlimitsForRD;

////////////////////////////////////////////////////////////////////////////////
//
class G4RadioactiveDecaymessenger;

typedef std::vector<G4RadioactiveDecayRateVector> G4RadioactiveDecayRateTable;
typedef std::vector<G4RadioactiveDecayRate> G4RadioactiveDecayRates;

class G4RadioactiveDecay : public G4VRestDiscreteProcess 
{
	// class description

	// Implementation of nuclear radioactive decay process. 
	// It simulates the decays of radioactive nuclei, which is submitted to RDM in the form
	// of G4Ions. Half-liffe and decay schemes are hold in the Radioactivity database
	// All decay products are submitted back to the particle tracking process through 
	// the G4ParticleChangeForRadDecay object.

	// class description - end 

public: // with description

	G4RadioactiveDecay (const G4String& processName="RadioactiveDecay");
	~G4RadioactiveDecay();
	// constructor and destructor
	//

	G4bool IsApplicable(const G4ParticleDefinition &);
	// Returns true if the specified isotope is
	//  1) defined as "nucleus" and
	//  2) it is within theNucleusLimit
	//
	G4bool IsLoaded (const G4ParticleDefinition &);
	// Returns true if the decay table of the specified nuclei is ready.
	//
	void SelectAVolume(const G4String aVolume);
	// Select a logical volume in which RDM applies.
	//
	void DeselectAVolume(const G4String aVolume);
	// remove a logical volume from the RDM applied list
	//
	void SelectAllVolumes();
	// Select all logical volumes for the application of RDM.
	//
	void DeselectAllVolumes();
	// Remove all logical volumes from RDM applications.
	//
	void SetDecayBias (G4String filename);
	//   Sets the decay biasing scheme using the data in "filename"
	//
	void SetHLThreshold (G4double hl) {halflifethreshold = hl;}
	// Set the half-life threshold for isomer production
	//
	void SetICM (G4bool icm) {applyICM = icm;}
	// Enable/disable ICM 
	//
	void SetARM (G4bool arm) {applyARM = arm;}
	// Enable/disable ARM
	// 
	void SetSourceTimeProfile (G4String filename) ;
	//  Sets source exposure function using histograms in "filename"
	//
	G4bool IsRateTableReady(const G4ParticleDefinition &);
	// Returns true if the coefficient and decay time table for all the
	// descendants of the specified isotope are ready.
	//
	// used in VR decay mode only 
	//
	void AddDecayRateTable(const G4ParticleDefinition &);
	// Calculates the coefficient and decay time table for all the descendants
	// of the specified isotope.  Adds the calculated table to the private data
	// member "theDecayRateTableVector".
	//
	//
	// used in VR decay mode only 
	//
	void GetDecayRateTable(const G4ParticleDefinition &);
	// Used to retrieve the coefficient and decay time table for all the
	// descendants of the specified isotope from "theDecayRateTableVector"
	// and place it in "theDecayRateTable".
	//
	//
	// used in VR decay mode only 
	//
	void SetDecayRate(G4int,G4int,G4double, G4int, std::vector<G4double>,
		std::vector<G4double>);
	// Sets "theDecayRate" with data supplied in the arguements.
	//
	//  //
	// used in VR decay mode only 
	//

	std::vector<G4RadioactivityTable*> GetTheRadioactivityTables() {return theRadioactivityTables;}
	// return the vector of G4Radioactivity map
	// should be used in VR mode only

	G4DecayTable *LoadDecayTable (G4ParticleDefinition & theParentNucleus);
	// Load the decay data of isotope theParentNucleus.
	//
	inline void  SetVerboseLevel(G4int value) {verboseLevel = value;}
	// Sets the VerboseLevel which controls duggering display.
	//
	inline G4int GetVerboseLevel() const {return verboseLevel;}
	// Returns the VerboseLevel which controls level of debugging output.
	//
	inline void SetNucleusLimits(G4NucleusLimits theNucleusLimits1)
	{theNucleusLimits = theNucleusLimits1 ;}
	//  Sets theNucleusLimits which specifies the range of isotopes
	//  the G4RadioactiveDecay applies.
	//
	inline G4NucleusLimits GetNucleusLimits() const
	{return theNucleusLimits;}
	//  Returns theNucleusLimits which specifies the range of isotopes
	//  the G4RadioactiveDecay applies.
	//
	inline void SetAnalogueMonteCarlo (G4bool r ) { 
	  AnalogueMC  = r; 
	  if (!AnalogueMC) halflifethreshold = 1e-6*s;
	}
	//   Controls whether G4RadioactiveDecay runs in analogue mode or
	//   variance reduction mode.
	inline void SetFBeta (G4bool r ) { FBeta  = r; }
	//   Controls whether G4RadioactiveDecay uses fast beta simulation mode
	//
	inline G4bool IsAnalogueMonteCarlo () {return AnalogueMC;}
	//   Returns true if the simulation is an analogue Monte Carlo, and false if
	//   any of the biassing schemes have been selected.
	//
	inline void SetBRBias (G4bool r ) { BRBias  = r;
	SetAnalogueMonteCarlo(0);}
	//   Sets whether branching ration bias scheme applies.
	//
	inline void SetSplitNuclei (G4int r) { NSplit = r; SetAnalogueMonteCarlo(0); }
	//  Sets the N number for the Nuclei spliting bias scheme
	//
	inline G4int GetSplitNuclei () {return NSplit;}
	//  Returns the N number used for the Nuclei spliting bias scheme
	//
public:

	void BuildPhysicsTable(const G4ParticleDefinition &);

protected:

	G4VParticleChange* DecayIt( const G4Track& theTrack,
		const G4Step&  theStep);

	G4DecayProducts* DoDecay(G4ParticleDefinition& theParticleDef);

	G4double GetMeanFreePath(const G4Track& theTrack,
		G4double previousStepSize,
		G4ForceCondition* condition);

	G4double GetMeanLifeTime(const G4Track& theTrack,
		G4ForceCondition* condition);

	G4double GetTaoTime(G4double,G4double);

	G4double GetDecayTime();

	G4int GetDecayTimeBin(const G4double aDecayTime);

  private:

    G4RadioactiveDecay(const G4RadioactiveDecay &right);
    G4RadioactiveDecay & operator=(const G4RadioactiveDecay &right);

  private:

    G4RadioactiveDecaymessenger* theRadioactiveDecaymessenger;
    G4PhysicsTable* aPhysicsTable;
    G4RIsotopeTable* theIsotopeTable;

    G4NucleusLimits theNucleusLimits;

    const G4double HighestBinValue;
    const G4double LowestBinValue;

    const G4int TotBin;

    G4bool AnalogueMC;
    G4bool BRBias;
    G4bool FBeta;
    G4int NSplit;

    G4double halflifethreshold;
    G4bool applyICM;
    G4bool applyARM;

    G4int NSourceBin;
    G4double SBin[100];
    G4double SProfile[100];
    G4int NDecayBin;
    G4double DBin[100];
    G4double DProfile[100];

	std::vector<G4String>         LoadedNuclei;
	std::vector<G4String>         ValidVolumes;

	G4RadioactiveDecayRate        theDecayRate;
	G4RadioactiveDecayRates       theDecayRateVector;
        G4RadioactiveDecayRateVector  theDecayRateTable;
        G4RadioactiveDecayRateTable   theDecayRateTableVector;

	// for the radioactivity tables
        std::vector<G4RadioactivityTable*>	theRadioactivityTables;
        G4int		              decayWindows[99];

	//
	static const G4double		  levelTolerance;

	// Remainder of life time at rest
	G4double                      fRemainderLifeTime;

	G4int                         verboseLevel;


	// ParticleChange for decay process
	G4ParticleChangeForRadDecay   fParticleChangeForRadDecay;

	inline G4double AtRestGetPhysicalInteractionLength
		(const G4Track& track, G4ForceCondition* condition)
	{fRemainderLifeTime = G4VRestDiscreteProcess::
	AtRestGetPhysicalInteractionLength(track, condition );
	return fRemainderLifeTime;}

	inline G4VParticleChange* AtRestDoIt
		(const G4Track& theTrack, const G4Step& theStep)
	{return DecayIt(theTrack, theStep);}

	inline G4VParticleChange* PostStepDoIt
		(const G4Track& theTrack, const G4Step& theStep)
	{return DecayIt(theTrack, theStep);}

};
#endif







