// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// G4qmdDummy.hh  2000/08/03  Stefan Scherer
//
#ifndef G4qmdDummy_h
#define G4qmdDummy_h 1

#include "g4std/fstream"
#include "G4ios.hh"
#include "globals.hh"
#include "Propagation.hh"


class Colour;

class G4qmdDummy 
{
  public:
      G4qmdDummy();
      G4qmdDummy(const G4String & anInputFile="");
      ~G4qmdDummy();

  private:
      G4qmdDummy(const G4qmdDummy &right);
      const G4qmdDummy & operator=(const G4qmdDummy &right);
      int operator==(const G4qmdDummy &right) const;
      int operator!=(const G4qmdDummy &right) const;

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
      void justRun();

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

      Colour  theBox;

};



inline void G4qmdDummy::SetInputFile(G4String & anInputFile)
{
   theInputFile = anInputFile;
}

inline G4String G4qmdDummy::GetInputFile()
{
   return theInputFile;
}


inline void G4qmdDummy::SetColorStringDecay(G4String & aColorStringDecay)
{
   theColorStringDecay = aColorStringDecay;
}

inline void G4qmdDummy::SetColorCluster(G4String & aColorCluster)
{
   theColorCluster = aColorCluster;
}

inline void G4qmdDummy::SetDirectHadronFromString(G4String & aDirectHadronFromString)
{
   theDirectHadronFromString = aDirectHadronFromString;
}

inline void G4qmdDummy::SetFinalHadronDecay(G4String & aFinalHadronDecay)
{
   theFinalHadronDecay = aFinalHadronDecay;
}

inline void G4qmdDummy::SetForcedHadronDecay(G4String & aForcedHadronDecay)
{
   theForcedHadronDecay = aForcedHadronDecay;
}

inline void G4qmdDummy::SetColorPotential(G4String & aColorPotential)
{
   theColorPotential = aColorPotential;
}

inline void G4qmdDummy::SetParameterHadronizationCriterium(double aParameterHadronizationCriterium)
{
   theParameterHadronizationCriterium = aParameterHadronizationCriterium;
}

inline void G4qmdDummy::SetParameterKappa(double aParameterKappa)
{
   theParameterKappa = aParameterKappa;
}

inline void G4qmdDummy::SetFinalTime(double aFinalTime)
{
   theFinalTime = aFinalTime;
}

inline void G4qmdDummy::SetOutputTimestep(double anOutputTimestep)
{
   theOutputTimestep = anOutputTimestep;
}



#endif


