// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPInelastic.hh,v 1.1 1999-01-07 16:13:01 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Process: High Precision low E neutron tracking
 // original by H.P. Wellisch, TRIUMF, 14-Feb-97
 // Builds and has the Cross-section data for one material.
 
#ifndef G4NeutronHPInelastic_h
#define G4NeutronHPInelastic_h 1

#include "globals.hh"
#include "G4NeutronHPChannel.hh"
#include "G4HadronicInteraction.hh"
#include "G4NeutronHPChannelList.hh"

#include "G4NeutronHP2AInelasticFS.hh"
#include "G4NeutronHP2N2AInelasticFS.hh"
#include "G4NeutronHP2NAInelasticFS.hh"
#include "G4NeutronHP2NDInelasticFS.hh"
#include "G4NeutronHP2NInelasticFS.hh"
#include "G4NeutronHP2NPInelasticFS.hh"
#include "G4NeutronHP2PInelasticFS.hh"
#include "G4NeutronHP3AInelasticFS.hh"
#include "G4NeutronHP3NAInelasticFS.hh"
#include "G4NeutronHP3NInelasticFS.hh"
#include "G4NeutronHP3NPInelasticFS.hh"
#include "G4NeutronHP4NInelasticFS.hh"
#include "G4NeutronHPAInelasticFS.hh"
#include "G4NeutronHPD2AInelasticFS.hh"
#include "G4NeutronHPDAInelasticFS.hh"
#include "G4NeutronHPDInelasticFS.hh"
#include "G4NeutronHPHe3InelasticFS.hh"
#include "G4NeutronHPN2AInelasticFS.hh"
#include "G4NeutronHPN2PInelasticFS.hh"
#include "G4NeutronHPN3AInelasticFS.hh"
#include "G4NeutronHPNAInelasticFS.hh"
#include "G4NeutronHPND2AInelasticFS.hh"
#include "G4NeutronHPNDInelasticFS.hh"
#include "G4NeutronHPNHe3InelasticFS.hh"
#include "G4NeutronHPNInelasticFS.hh"
#include "G4NeutronHPNPAInelasticFS.hh"
#include "G4NeutronHPNPInelasticFS.hh"
#include "G4NeutronHPNT2AInelasticFS.hh"
#include "G4NeutronHPNTInelasticFS.hh"
#include "G4NeutronHPNXInelasticFS.hh"
#include "G4NeutronHPPAInelasticFS.hh"
#include "G4NeutronHPPDInelasticFS.hh"
#include "G4NeutronHPPInelasticFS.hh"
#include "G4NeutronHPPTInelasticFS.hh"
#include "G4NeutronHPT2AInelasticFS.hh"
#include "G4NeutronHPTInelasticFS.hh"

class G4NeutronHPInelastic : public G4HadronicInteraction
{
  public: 
  
  G4NeutronHPInelastic();

  ~G4NeutronHPInelastic();
  
  G4VParticleChange * ApplyYourself(const G4Track & aTrack, G4Nucleus & aTargetNucleus);

  private:
  
  G4double * xSec;
  G4NeutronHPChannelList * theInelastic; // one List per element
  G4String dirName;
  G4int numEle;
  
  private:
  
   G4NeutronHP2AInelasticFS the2AFS;
   G4NeutronHP2N2AInelasticFS the2N2AFS;
   G4NeutronHP2NAInelasticFS the2NAFS;
   G4NeutronHP2NDInelasticFS the2NDFS;
   G4NeutronHP2NInelasticFS the2NFS;
   G4NeutronHP2NPInelasticFS the2NPFS;
   G4NeutronHP2PInelasticFS the2PFS;
   G4NeutronHP3AInelasticFS the3AFS;
   G4NeutronHP3NAInelasticFS the3NAFS;
   G4NeutronHP3NInelasticFS the3NFS;
   G4NeutronHP3NPInelasticFS the3NPFS;
   G4NeutronHP4NInelasticFS the4NFS;
   G4NeutronHPAInelasticFS theAFS;
   G4NeutronHPD2AInelasticFS theD2AFS;
   G4NeutronHPDAInelasticFS theDAFS;
   G4NeutronHPDInelasticFS theDFS;
   G4NeutronHPHe3InelasticFS theHe3FS;
   G4NeutronHPN2AInelasticFS theN2AFS;
   G4NeutronHPN2PInelasticFS theN2PFS;
   G4NeutronHPN3AInelasticFS theN3AFS;
   G4NeutronHPNAInelasticFS theNAFS;
   G4NeutronHPND2AInelasticFS theND2AFS;
   G4NeutronHPNDInelasticFS theNDFS;
   G4NeutronHPNHe3InelasticFS theNHe3FS;
   G4NeutronHPNInelasticFS theNFS;
   G4NeutronHPNPAInelasticFS theNPAFS;
   G4NeutronHPNPInelasticFS theNPFS;
   G4NeutronHPNT2AInelasticFS theNT2AFS;
   G4NeutronHPNTInelasticFS theNTFS;
   G4NeutronHPNXInelasticFS theNXFS;
   G4NeutronHPPAInelasticFS thePAFS;
   G4NeutronHPPDInelasticFS thePDFS;
   G4NeutronHPPInelasticFS thePFS;
   G4NeutronHPPTInelasticFS thePTFS;
   G4NeutronHPT2AInelasticFS theT2AFS;
   G4NeutronHPTInelasticFS theTFS;
};

#endif
