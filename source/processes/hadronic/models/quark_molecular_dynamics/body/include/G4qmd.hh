// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// G4qmd.hh  2000/08/03  Stefan Scherer
//
#ifndef G4qmd_h
#define G4qmd_h 1

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
#include "G4ExcitedStringDecay.hh"

class Knot<ParticleType>;
class Colour;
class InverseFunction;
class Particle;
class CollisionTab;

class G4qmd 
{
  public:
      G4qmd();
      G4qmd(const G4String & anInputFile="");
      G4qmd(const G4ExcitedStringVector * theStrings);
      ~G4qmd();

  private:
      G4qmd(const G4qmd &right);
      const G4qmd & operator=(const G4qmd &right);
      int operator==(const G4qmd &right) const;
      int operator!=(const G4qmd &right) const;

      void skipline(istream& in);
      double readEvent(istream& in);

  public:
      void SetInputFile(G4String & anInputFile);
      G4String GetInputFile();
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

      void SetupFromFile();
      void SetupFromG4ExcitedString();
 
      G4KineticTrackVector * TheHadrons() ;

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

      Colour  theQuarkSystem;

};



inline void G4qmd::SetInputFile(G4String & anInputFile)
{
   theInputFile = anInputFile;
}

inline G4String G4qmd::GetInputFile()
{
   return theInputFile;
}


inline void G4qmd::SetColorStringDecay(G4String & aColorStringDecay)
{
   theColorStringDecay = aColorStringDecay;
   Colour::allowDecay = (theColorStringDecay == "yes");
}

inline void G4qmd::SetColorCluster(G4String & aColorCluster)
{
   theColorCluster = aColorCluster;
   Colour::allowClustering = (theColorCluster == "yes");
}

inline void G4qmd::SetDirectHadronFromString(G4String & aDirectHadronFromString)
{
   theDirectHadronFromString = aDirectHadronFromString;
}

inline void G4qmd::SetFinalHadronDecay(G4String & aFinalHadronDecay)
{
   theFinalHadronDecay = aFinalHadronDecay;
}

inline void G4qmd::SetForcedHadronDecay(G4String & aForcedHadronDecay)
{
   theForcedHadronDecay = aForcedHadronDecay;
}

inline void G4qmd::SetColorPotential(G4String & aColorPotential)
{
   theColorPotential = aColorPotential;
}

inline void G4qmd::SetParameterHadronizationCriterium(double aParameterHadronizationCriterium)
{
   theParameterHadronizationCriterium = aParameterHadronizationCriterium;
}

inline void G4qmd::SetParameterKappa(double aParameterKappa)
{
   theParameterKappa = aParameterKappa;
   Colour::kappa = theParameterKappa;
}

inline void G4qmd::SetFinalTime(double aFinalTime)
{
   theFinalTime = aFinalTime;
}

inline void G4qmd::SetOutputTimestep(double aOutputTimestep)
{
   theOutputTimestep = aOutputTimestep;
}



//inline G4KineticTrackVector * TheHadrons() 
//{
//return 0;
//}

#endif


