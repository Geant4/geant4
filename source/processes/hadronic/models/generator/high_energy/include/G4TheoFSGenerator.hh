// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4TheoFSGenerator.hh,v 1.2 1999/04/12 15:45:28 hpw Exp $
// GEANT4 tag $Name: geant4-02-00 $
//
#ifndef G4TheoFSGenerator_h
#define G4TheoFSGenerator_h 1

#include "G4VIntraNuclearTransportModel.hh"
#include "G4HadronicInteraction.hh"
#include "G4VHighEnergyGenerator.hh"

class G4TheoFSGenerator : public G4HadronicInteraction 

{
  public:
      G4TheoFSGenerator();
      ~G4TheoFSGenerator();

  private:
      G4TheoFSGenerator(const G4TheoFSGenerator &right);
      const G4TheoFSGenerator & operator=(const G4TheoFSGenerator &right);
      int operator==(const G4TheoFSGenerator &right) const;
      int operator!=(const G4TheoFSGenerator &right) const;

  public:
      G4VParticleChange * ApplyYourself(const G4Track & thePrimary, G4Nucleus & theNucleus);
      void SetTransport(G4VIntraNuclearTransportModel *const  value);
      void SetHighEnergyGenerator(G4VHighEnergyGenerator *const  value);

  private:
      const G4VIntraNuclearTransportModel * GetTransport() const;
      const G4VHighEnergyGenerator * GetHighEnergyGenerator() const;
      const G4ParticleChange * GetParticleChange() const;

  private: 
      G4VIntraNuclearTransportModel *theTransport;
      G4VHighEnergyGenerator *theHighEnergyGenerator;

      G4ParticleChange *theParticleChange;
};

inline const G4VIntraNuclearTransportModel * G4TheoFSGenerator::GetTransport() const
{
  return theTransport;
}

inline void G4TheoFSGenerator::SetTransport(G4VIntraNuclearTransportModel *const  value)
{
  theTransport = value;
}

inline const G4VHighEnergyGenerator * G4TheoFSGenerator::GetHighEnergyGenerator() const
{
  return theHighEnergyGenerator;
}

inline void G4TheoFSGenerator::SetHighEnergyGenerator(G4VHighEnergyGenerator *const  value)
{
  theHighEnergyGenerator= value;
}

inline const G4ParticleChange * G4TheoFSGenerator::GetParticleChange() const
{
  return theParticleChange;
}


#endif


