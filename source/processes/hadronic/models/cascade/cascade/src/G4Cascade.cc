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
//## begin module%38280A1302D0.cm preserve=no
//	  %X% %Q% %Z% %W%
//## end module%38280A1302D0.cm

//## begin module%38280A1302D0.cp preserve=no
//## end module%38280A1302D0.cp

//## Module: G4Cascade%38280A1302D0; Subprogram body
//## Subsystem: cascade::src%3819C79F9ED0
//## Source file: /tmp/miheikki/geant4/source/processes/hadronic/models/cascade/cascade/src/G4Cascade.cc

//## begin module%38280A1302D0.additionalIncludes preserve=no
//## end module%38280A1302D0.additionalIncludes

//## begin module%38280A1302D0.includes preserve=yes
#include "G4Cascade.hh"
//## end module%38280A1302D0.includes

//## begin module%38280A1302D0.declarations preserve=no
//## end module%38280A1302D0.declarations

//## begin module%38280A1302D0.additionalDeclarations preserve=yes
//## end module%38280A1302D0.additionalDeclarations


// Class G4Cascade 

G4Cascade::G4Cascade()
  //## begin G4Cascade::G4Cascade%.hasinit preserve=no
  //## end G4Cascade::G4Cascade%.hasinit
  //## begin G4Cascade::G4Cascade%.initialization preserve=yes
  //## end G4Cascade::G4Cascade%.initialization
{
  //## begin G4Cascade::G4Cascade%.body preserve=yes
  //## end G4Cascade::G4Cascade%.body
}


G4Cascade::~G4Cascade()
{
  //## begin G4Cascade::~G4Cascade%.body preserve=yes
  //## end G4Cascade::~G4Cascade%.body
}



//## Other Operations (implementation)
G4VParticleChange* G4Cascade::ApplyYourself (const G4Track& aTrack, G4Nucleus& theNucleus)
{
  //## begin G4Cascade::ApplyYourself%942063060.body preserve=yes
 return new G4ParticleChange;
  //## end G4Cascade::ApplyYourself%942063060.body
}

G4ReactionProductVector* G4Cascade::Propagate (G4KineticTrackVector* theSecondaries, G4V3DNucleus* theNucleus)
{
  //## begin G4Cascade::Propagate%942064478.body preserve=yes
 return new G4ReactionProductVector;
  //## end G4Cascade::Propagate%942064478.body
}

G4int G4Cascade::operator == (G4Cascade& right)
{
  //## begin G4Cascade::operator==%942064479.body preserve=yes
return (this == &right);
  //## end G4Cascade::operator==%942064479.body
}

G4int G4Cascade::operator != (G4Cascade& right)
{
  //## begin G4Cascade::operator!=%942064480.body preserve=yes
return (this != &right);
  //## end G4Cascade::operator!=%942064480.body
}

// Additional Declarations
  //## begin G4Cascade%3819C64204D0.declarations preserve=yes
  //## end G4Cascade%3819C64204D0.declarations

//## begin module%38280A1302D0.epilog preserve=yes
//## end module%38280A1302D0.epilog
