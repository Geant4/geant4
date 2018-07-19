//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
 // Hadronic Process: High Precision low E neutron tracking
 // original by H.P. Wellisch, TRIUMF, 14-Feb-97
 // Builds and has the Cross-section data for one material.
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#ifndef G4ParticleHPInelastic_h
#define G4ParticleHPInelastic_h 1

// Class Description
// Final state production model for a high precision (based on evaluated data
// libraries) description of neutron inelastic scattering below 20 MeV; 
// 36 exclusive final states are consideded.
// To be used in your physics list in case you need this physics.
// In this case you want to register an object of this class with 
// the corresponding process.
// Class Description - End

#include "globals.hh"
#include "G4ParticleHPChannel.hh"
#include "G4HadronicInteraction.hh"
#include "G4ParticleHPChannelList.hh"

/*
#include "G4ParticleHP2AInelasticFS.hh"
#include "G4ParticleHP2N2AInelasticFS.hh"
#include "G4ParticleHP2NAInelasticFS.hh"
#include "G4ParticleHP2NDInelasticFS.hh"
#include "G4ParticleHP2NInelasticFS.hh"
#include "G4ParticleHP2NPInelasticFS.hh"
#include "G4ParticleHP2PInelasticFS.hh"
#include "G4ParticleHP3AInelasticFS.hh"
#include "G4ParticleHP3NAInelasticFS.hh"
#include "G4ParticleHP3NInelasticFS.hh"
#include "G4ParticleHP3NPInelasticFS.hh"
#include "G4ParticleHP4NInelasticFS.hh"
#include "G4ParticleHPAInelasticFS.hh"
#include "G4ParticleHPD2AInelasticFS.hh"
#include "G4ParticleHPDAInelasticFS.hh"
#include "G4ParticleHPDInelasticFS.hh"
#include "G4ParticleHPHe3InelasticFS.hh"
#include "G4ParticleHPN2AInelasticFS.hh"
#include "G4ParticleHPN2PInelasticFS.hh"
#include "G4ParticleHPN3AInelasticFS.hh"
#include "G4ParticleHPNAInelasticFS.hh"
#include "G4ParticleHPND2AInelasticFS.hh"
#include "G4ParticleHPNDInelasticFS.hh"
#include "G4ParticleHPNHe3InelasticFS.hh"
#include "G4ParticleHPNInelasticFS.hh"
#include "G4ParticleHPNPAInelasticFS.hh"
#include "G4ParticleHPNPInelasticFS.hh"
#include "G4ParticleHPNT2AInelasticFS.hh"
#include "G4ParticleHPNTInelasticFS.hh"
#include "G4ParticleHPNXInelasticFS.hh"
#include "G4ParticleHPPAInelasticFS.hh"
#include "G4ParticleHPPDInelasticFS.hh"
#include "G4ParticleHPPInelasticFS.hh"
#include "G4ParticleHPPTInelasticFS.hh"
#include "G4ParticleHPT2AInelasticFS.hh"
#include "G4ParticleHPTInelasticFS.hh"
*/
#include "G4ParticleDefinition.hh"

class G4ParticleHPInelastic : public G4HadronicInteraction
{
  public: 

  G4ParticleHPInelastic(G4ParticleDefinition* projectile = G4Neutron::Neutron(), const char* name = "NeutronHPInelastic" );

  ~G4ParticleHPInelastic();
  
  G4HadFinalState * ApplyYourself(const G4HadProjectile & aTrack, G4Nucleus & aTargetNucleus);
  virtual const std::pair<G4double, G4double> GetFatalEnergyCheckLevels() const;

  public:
      G4int GetVerboseLevel() const;
      void SetVerboseLevel( G4int );
      void BuildPhysicsTable(const G4ParticleDefinition&);
      virtual void ModelDescription(std::ostream& outFile) const;

protected:
  
  //G4ParticleHPChannelList * theInelastic; // one List per element
  std::vector<G4ParticleHPChannelList*>* theInelastic; // one List per element
  G4String dataDirVariable;
  G4String dirName;
  G4int numEle;

  private:
 /* 
   G4ParticleHP2AInelasticFS the2AFS;
   G4ParticleHP2N2AInelasticFS the2N2AFS;
   G4ParticleHP2NAInelasticFS the2NAFS;
   G4ParticleHP2NDInelasticFS the2NDFS;
   G4ParticleHP2NInelasticFS the2NFS;
   G4ParticleHP2NPInelasticFS the2NPFS;
   G4ParticleHP2PInelasticFS the2PFS;
   G4ParticleHP3AInelasticFS the3AFS;
   G4ParticleHP3NAInelasticFS the3NAFS;
   G4ParticleHP3NInelasticFS the3NFS;
   G4ParticleHP3NPInelasticFS the3NPFS;
   G4ParticleHP4NInelasticFS the4NFS;
   G4ParticleHPAInelasticFS theAFS;
   G4ParticleHPD2AInelasticFS theD2AFS;
   G4ParticleHPDAInelasticFS theDAFS;
   G4ParticleHPDInelasticFS theDFS;
   G4ParticleHPHe3InelasticFS theHe3FS;
   G4ParticleHPN2AInelasticFS theN2AFS;
   G4ParticleHPN2PInelasticFS theN2PFS;
   G4ParticleHPN3AInelasticFS theN3AFS;
   G4ParticleHPNAInelasticFS theNAFS;
   G4ParticleHPND2AInelasticFS theND2AFS;
   G4ParticleHPNDInelasticFS theNDFS;
   G4ParticleHPNHe3InelasticFS theNHe3FS;
   G4ParticleHPNInelasticFS theNFS;
   G4ParticleHPNPAInelasticFS theNPAFS;
   G4ParticleHPNPInelasticFS theNPFS;
   G4ParticleHPNT2AInelasticFS theNT2AFS;
   G4ParticleHPNTInelasticFS theNTFS;
   G4ParticleHPNXInelasticFS theNXFS;
   G4ParticleHPPAInelasticFS thePAFS;
   G4ParticleHPPDInelasticFS thePDFS;
   G4ParticleHPPInelasticFS thePFS;
   G4ParticleHPPTInelasticFS thePTFS;
   G4ParticleHPT2AInelasticFS theT2AFS;
   G4ParticleHPTInelasticFS theTFS;
*/

   G4ParticleDefinition* theProjectile;

};

#endif
