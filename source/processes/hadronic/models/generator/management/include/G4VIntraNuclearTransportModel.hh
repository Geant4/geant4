// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VIntraNuclearTransportModel.hh,v 1.2 1998/11/27 08:26:52 pavliouk Exp $
// GEANT4 tag $Name: geant4-00 $
//
// $Id: G4IntraNuclearTransportMode.hh,v 1.0 1998/06/30
// -----------------------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, A. Feliciello, 30th June 1998
//      A.Pavliouk 26.11.98
//          In Set...() methods a pointer is deleted now before new
//          value will be asigned.
// -----------------------------------------------------------------------------

#ifndef G4VIntraNuclearTransportModel_h
#define G4VIntraNuclearTransportModel_h 1

#include "G4V3DNucleus.hh"
#include "G4VPreCompoundModel.hh"
#include "G4HadronicInteraction.hh"
#include "G4DynamicParticleVector.hh"
#include "G4Track.hh"



class G4VIntraNuclearTransportModel : public G4HadronicInteraction

{
  public:

      G4VIntraNuclearTransportModel();

      G4VIntraNuclearTransportModel(const G4VIntraNuclearTransportModel& right);

      virtual ~G4VIntraNuclearTransportModel();

      const G4VIntraNuclearTransportModel& operator=(const G4VIntraNuclearTransportModel &right);

      int operator==(const G4VIntraNuclearTransportModel& right) const;

      int operator!=(const G4VIntraNuclearTransportModel& right) const;

      virtual G4VParticleChange* ApplyYourself(const G4Track& aTrack, 
                                              G4Nucleus& theNucleus) = 0;

      virtual G4DynamicParticleVector* Propagate(G4KineticTrackVector* theSecondaries,
                                                 G4V3DNucleus* theNucleus) = 0;

      void SetDeExcitation(G4VPreCompoundModel* const  value);
      void Set3DNucleus(G4V3DNucleus* const value);


  protected:

      G4V3DNucleus* Get3DNucleus() const;

      G4VPreCompoundModel* GetDeExcitation() const;



  protected:

      G4V3DNucleus* the3DNucleus;

      G4VPreCompoundModel* theDeExcitation;
};



// Class G4VIntraNuclearTransportModel 


inline G4V3DNucleus* G4VIntraNuclearTransportModel::Get3DNucleus() const
{
  return the3DNucleus;
}

inline void G4VIntraNuclearTransportModel::Set3DNucleus(G4V3DNucleus* const value)
{
   delete the3DNucleus;  the3DNucleus = value;
}



inline G4VPreCompoundModel* G4VIntraNuclearTransportModel::GetDeExcitation() const
{
  return theDeExcitation;
}

inline void G4VIntraNuclearTransportModel::SetDeExcitation(G4VPreCompoundModel* const  value)
{
   delete theDeExcitation; theDeExcitation = value;
}

#endif


