//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4TheoFSGenerator.hh,v 1.2 2005/06/04 13:50:03 jwellisc Exp $
// GEANT4 tag $Name: geant4-08-00 $
//

// Class Description
// Final state production model for theoretical models of hadron inelastic
// scattering in geant4;
// To be used in your physics list in case you need this physics.
// In this case you want to register an object of this class with 
// the corresponding process.
// Note: This class is part of an implementation framework. You need to
// register corresponding high energy generators and transport codes to 
// fill it with life; decay of strong resonances is done directly,
// in case there is no residual nucleus. 
// Class Description - End

#ifndef G4TheoFSGenerator_h
#define G4TheoFSGenerator_h 1

#include "G4VIntraNuclearTransportModel.hh"
#include "G4HadronicInteraction.hh"
#include "G4VHighEnergyGenerator.hh"
#include "G4DecayStrongResonances.hh"
#include "G4HadFinalState.hh"

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
      G4HadFinalState * ApplyYourself(const G4HadProjectile & thePrimary, G4Nucleus & theNucleus);
      void SetTransport(G4VIntraNuclearTransportModel *const  value);
      void SetHighEnergyGenerator(G4VHighEnergyGenerator *const  value);

  private:
      const G4VIntraNuclearTransportModel * GetTransport() const;
      const G4VHighEnergyGenerator * GetHighEnergyGenerator() const;
      const G4HadFinalState * GetFinalState() const;

  private: 
      G4VIntraNuclearTransportModel *theTransport;
      G4VHighEnergyGenerator *theHighEnergyGenerator;
      G4DecayStrongResonances theDecay;
      G4HadFinalState *theParticleChange;
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

inline const G4HadFinalState * G4TheoFSGenerator::GetFinalState() const
{
  return theParticleChange;
}


#endif


