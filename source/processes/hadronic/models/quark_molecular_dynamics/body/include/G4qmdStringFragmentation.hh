// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// G4qmdStringFragmentation.hh  2001/03/29  Stefan Scherer
//
#ifndef G4qmdStringFragmentation_h
#define G4qmdStringFragmentation_h 1

#include "g4std/fstream"
#include "G4ios.hh"
#include <algo.h>
#include "newvector.hh"
#include "Random.hh"
#include "newBinning.hh"
#include "Geometry.hh"
#include "Arguments.hh"
#include "MathTools.hh"
#include "Error.hh"
#include "Propagation.hh"
#include "ParticleBase.hh"
#include "reactionChannels.hh"
#include "ParticleKinematics.hh"
#include "iso.hh"
#include "array.hh"
#include "Memory.hh"
#include "Collision.hh"
#include "genericRead.hh"
#include "StreamBuffer.hh"
#include "RunningVariable.hh"
#include "InputVariable.hh"
#include "InputReader.hh"
#include "Volume.hh"
#include "Potentials.hh"
#include "Quarkbox.hh"
#include "output.hh"

#include "globals.hh"
#include "G4KineticTrack.hh"
#include "G4KineticTrackVector.hh"

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"

#include "G4ShortLivedConstructor.hh"
#include "G4ShortLivedTable.hh"
#include "G4BosonConstructor.hh"

#include "G4VDecayChannel.hh"
#include "G4DecayTable.hh"

#include "G4VStringFragmentation.hh"

class Knot<ParticleType>;
class Colour;
class InverseFunction;
class Particle;
class CollisionTab;

class G4qmdStringFragmentation : public G4VStringFragmentation
{
  public:
      G4qmdStringFragmentation();
       ~G4qmdStringFragmentation();
      G4KineticTrackVector* FragmentStrings(const G4ExcitedStringVector* theStrings);
      G4KineticTrackVector* FragmentStringsFromFile(const G4String & anInputFile);

  private:
      G4qmdStringFragmentation(const G4qmdStringFragmentation &right);
      const G4qmdStringFragmentation & operator=(const G4qmdStringFragmentation &right);
      int operator==(const G4qmdStringFragmentation &right) const;
      int operator!=(const G4qmdStringFragmentation &right) const;

      void skipline(istream& in);
      double readEvent(istream& in);

  public:
      void SetInputFile(G4String & anInputFile);
      void SetColorStringDecay(G4String & aColorStringDecay);
      void SetColorCluster(G4String & aColorCluster);
      void SetDirectHadronFromString(G4String & aDirectHadronFromString);
      void SetFinalHadronDecay(G4String & aFinalHadronDecay);
      void SetForcedHadronDecay(G4String & aForcedHadronDecay);
      void SetColorPotential(G4String & aColorPotential);
      void SetParameterHadronizationCriterium(double aParameterHadronizationCriterium);
      void SetParameterKappa(double aParameterKappa);
      void SetFinalTime(double aFinalTime);
      void SetOutputTimestep(double aOutputTimestep);

      void SetupFromFile(const G4String & anInputFile="");
      void SetupFromG4ExcitedStringVector(const G4ExcitedStringVector * theInitalStrings);

			G4String GetInputFile();
      G4KineticTrackVector* TheHadrons() ;

  private:
      G4String  theInputFile;
      G4String  theColorStringDecay;
      G4String  theColorCluster;
      G4String  theDirectHadronFromString;
      G4String  theFinalHadronDecay;
      G4String  theForcedHadronDecay;
      G4String  theColorPotential;
      double  theParameterHadronizationCriterium;
      double  theParameterKappa;
      double  theFinalTime;
      double  theOutputTimestep;
      double  theInternalTimestep;

      Colour  * theQuarkSystem;

};



inline void G4qmdStringFragmentation::SetInputFile(G4String & anInputFile)
{
   theInputFile = anInputFile;
}

inline G4String G4qmdStringFragmentation::GetInputFile()
{
   return theInputFile;
}


inline void G4qmdStringFragmentation::SetColorStringDecay(G4String & aColorStringDecay)
{
   theColorStringDecay = aColorStringDecay;
   Colour::allowDecay = (theColorStringDecay == "yes");
}

inline void G4qmdStringFragmentation::SetColorCluster(G4String & aColorCluster)
{
   theColorCluster = aColorCluster;
   Colour::allowClustering = (theColorCluster == "yes");
}

inline void G4qmdStringFragmentation::SetDirectHadronFromString(G4String & aDirectHadronFromString)
{
   theDirectHadronFromString = aDirectHadronFromString;
}

inline void G4qmdStringFragmentation::SetFinalHadronDecay(G4String & aFinalHadronDecay)
{
   theFinalHadronDecay = aFinalHadronDecay;
}

inline void G4qmdStringFragmentation::SetForcedHadronDecay(G4String & aForcedHadronDecay)
{
   theForcedHadronDecay = aForcedHadronDecay;
}

inline void G4qmdStringFragmentation::SetColorPotential(G4String & aColorPotential)
{
   theColorPotential = aColorPotential;
}

inline void G4qmdStringFragmentation::SetParameterHadronizationCriterium(double aParameterHadronizationCriterium)
{
   theParameterHadronizationCriterium = aParameterHadronizationCriterium;
}

inline void G4qmdStringFragmentation::SetParameterKappa(double aParameterKappa)
{
   theParameterKappa = aParameterKappa;
   Colour::kappa = theParameterKappa;
}

inline void G4qmdStringFragmentation::SetFinalTime(double aFinalTime)
{
   theFinalTime = aFinalTime;
}

inline void G4qmdStringFragmentation::SetOutputTimestep(double aOutputTimestep)
{
   theOutputTimestep = aOutputTimestep;
}


#endif


