#include "globals.hh"
#include "G4ios.hh"
#include "G4qmdStringFragmentation.hh"

#include "G4KineticTrack.hh"
#include "G4KineticTrackVector.hh"

#include "G4ExcitedString.hh"
#include "G4ExcitedStringVector.hh"

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"

#include "G4ShortLivedConstructor.hh"
#include "G4ShortLivedTable.hh"

#include "G4VDecayChannel.hh"
#include "G4DecayTable.hh"


int main()
{
	G4ShortLivedConstructor ShortLived;
	ShortLived.ConstructParticle();

	G4ExcitedStringVector * qmdInitialState = new G4ExcitedStringVector();
	G4KineticTrackVector * ResultingHadrons = new G4KineticTrackVector();

//
// Generation of initial state
//

//--- pair 1:

	G4Parton* Quark = new G4Parton(1);
	G4Parton* AntiQuark = new G4Parton(-2);

	G4ThreeVector quarkposition(0,0,1);
	G4int quarkcolour = 1;

	Quark->SetPosition(quarkposition);
	Quark->SetColour(quarkcolour);
	AntiQuark->SetPosition(-quarkposition);
	AntiQuark->SetColour(-quarkcolour);

	Quark->GetDefinition()->DumpTable();
	AntiQuark->GetDefinition()->DumpTable();

	G4ExcitedString* TestString = new G4ExcitedString(Quark,AntiQuark);
	qmdInitialState->insert(TestString);

//--- pair 2:

	G4Parton* Quark2 = new G4Parton(3);
	G4Parton* AntiQuark2 = new G4Parton(-3);

	quarkposition = G4ThreeVector(0,1,0);
	quarkcolour = 2;

	Quark2->SetPosition(quarkposition);
	Quark2->SetColour(quarkcolour);
	AntiQuark2->SetPosition(-quarkposition);
	AntiQuark2->SetColour(-quarkcolour);

	Quark2->GetDefinition()->DumpTable();
	AntiQuark2->GetDefinition()->DumpTable();

	G4ExcitedString* TestString2 = new G4ExcitedString(Quark2,AntiQuark2);
	qmdInitialState->insert(TestString2);

//
// writing out initial state
//

	for (G4int StringCounter=0; StringCounter < qmdInitialState->length(); StringCounter++) {
 
		G4ExcitedString* ThisExcitedString = qmdInitialState->at(StringCounter);
		const G4PartonVector* ThePartons = ThisExcitedString->GetPartonList();

		G4cerr << "Extract partons from G4ExcitedString " << StringCounter << endl;

	  for (G4int QuarkCounter=0; QuarkCounter < ThePartons->length(); QuarkCounter++) {

			G4int flag = StringCounter*1000 + QuarkCounter;

 		  G4Parton* ThisParton = ThePartons->at(QuarkCounter);

      const G4ThreeVector quark_ThreeVector = ThisParton->GetPosition();
      const G4LorentzVector quark_4Momentum = ThisParton->Get4Momentum(); 
			G4int quark_PDGCode = ThisParton->GetPDGcode();
      G4int quark_Colour = ThisParton->GetColour();
      G4double quark_SpinZ = ThisParton->GetSpinZ();
      G4double quark_IsoSpinZ = ThisParton->GetIsoSpinZ();

			G4cerr << " Parton " << flag << " with properties " << G4endl;
      G4cerr << "    PDG-CODE:  " <<  quark_PDGCode << G4endl; 
      G4cerr << "  x-3-Vector:  " <<  quark_ThreeVector << G4endl; 
      G4cerr << "  p-4-Vector:  " <<  quark_4Momentum << G4endl; 
      G4cerr << "       color:  " <<  quark_Colour << G4endl; 
      G4cerr << "      Spin_Z:  " <<  quark_SpinZ << G4endl ;
      G4cerr << "   Isospin_Z:  " <<  quark_IsoSpinZ << G4endl;

		}

	}

//
// Start qMD
//

G4qmdStringFragmentation * theRun = new G4qmdStringFragmentation();
theRun->SetInputFile("Test");
//
// collect produced hadrons
//

// ResultingHadrons = theRun->FragmentStrings(qmdInitialState);

//

}



















