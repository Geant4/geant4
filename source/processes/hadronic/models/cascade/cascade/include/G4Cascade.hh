//## begin module%3826C0E45798.cm preserve=no
//	  %X% %Q% %Z% %W%
//## end module%3826C0E45798.cm

//## begin module%3826C0E45798.cp preserve=no
//## end module%3826C0E45798.cp

//## Module: G4Cascade%3826C0E45798; Subprogram specification
//## Subsystem: cascade::include%3819C79D6220
//## Source file: /tmp/miheikki/geant4/source/processes/hadronic/models/cascade/cascade/include/G4Cascade.hh

#ifndef G4Cascade_h
#define G4Cascade_h 1

//## begin module%3826C0E45798.additionalIncludes preserve=no
//## end module%3826C0E45798.additionalIncludes

//## begin module%3826C0E45798.includes preserve=yes
//## end module%3826C0E45798.includes

#include "G4VIntraNuclearTransportModel.hh"
//## begin module%3826C0E45798.declarations preserve=no
//## end module%3826C0E45798.declarations

//## begin module%3826C0E45798.additionalDeclarations preserve=yes
//## end module%3826C0E45798.additionalDeclarations


//## begin G4Cascade%3819C64204D0.preface preserve=yes
//## end G4Cascade%3819C64204D0.preface

//## Class: G4Cascade%3819C64204D0; implementation
//## Category: GEANT4::\ncascade::cascade%3819C5CE6378
//## Subsystem: cascade::include%3819C79D6220
//## Persistence: Transient
//## Cardinality/Multiplicity: n

class G4Cascade : public G4VIntraNuclearTransportModel  //## Inherits: <unnamed>%38281FBBA9A0
{
  //## begin G4Cascade%3819C64204D0.initialDeclarations preserve=yes
  //## end G4Cascade%3819C64204D0.initialDeclarations

  public:
    //## Constructors (generated)
      G4Cascade();

    //## Destructor (generated)
      ~G4Cascade();


    //## Other Operations (specified)
      //## Operation: ApplyYourself%942063060
      G4VParticleChange* ApplyYourself (const G4Track& aTrack, G4Nucleus& theNucleus);

      //## Operation: Propagate%942064478
      G4ReactionProductVector* Propagate (G4KineticTrackVector* theSecondaries, G4V3DNucleus* theNucleus);

    // Additional Public Declarations
      //## begin G4Cascade%3819C64204D0.public preserve=yes
      //## end G4Cascade%3819C64204D0.public

    // Additional Protected Declarations
      //## begin G4Cascade%3819C64204D0.protected preserve=yes
      //## end G4Cascade%3819C64204D0.protected

  private:

    //## Other Operations (specified)
      //## Operation: operator==%942064479
      G4int operator == (G4Cascade& right);

      //## Operation: operator!=%942064480
      G4int operator != (G4Cascade& right);

    // Additional Private Declarations
      //## begin G4Cascade%3819C64204D0.private preserve=yes
      //## end G4Cascade%3819C64204D0.private

    // Additional Implementation Declarations
      //## begin G4Cascade%3819C64204D0.implementation preserve=yes
      //## end G4Cascade%3819C64204D0.implementation

};

//## begin G4Cascade%3819C64204D0.postscript preserve=yes
//## end G4Cascade%3819C64204D0.postscript

// Class G4Cascade 

//## begin module%3826C0E45798.epilog preserve=yes
//## end module%3826C0E45798.epilog


#endif
