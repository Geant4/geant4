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

// ====================================================================
//
//   G4HepMCInterface.cc
//   $Id: G4HepMCInterface.cc,v 1.1 2002-04-29 20:44:51 asaim Exp $
//
// ====================================================================
#include "G4HepMCInterface.hh"

#include "G4RunManager.hh"
#include "G4LorentzVector.hh"
#include "G4Event.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4TransportationManager.hh"

////////////////////////////////////
G4HepMCInterface::G4HepMCInterface()
  : hepmcEvent(0)
////////////////////////////////////
{
}

/////////////////////////////////////
G4HepMCInterface::~G4HepMCInterface()
/////////////////////////////////////
{
  delete hepmcEvent;
}

/////////////////////////////////////////////////////////
G4bool G4HepMCInterface::CheckVertexInsideWorld
                         (const G4ThreeVector& pos) const
/////////////////////////////////////////////////////////
{
  G4Navigator* navigator= G4TransportationManager::GetTransportationManager()
                                                 -> GetNavigatorForTracking();

  G4VPhysicalVolume* world= navigator-> GetWorldVolume();
  G4VSolid* solid= world-> GetLogicalVolume()-> GetSolid();
  EInside qinside= solid-> Inside(pos);

  if( qinside != kInside) return false;
  else return true;
}

////////////////////////////////////////////////////////////////
void G4HepMCInterface::HepMC2G4(const HepMC::GenEvent* hepmcevt, 
				G4Event* g4event)
////////////////////////////////////////////////////////////////
{
  for(HepMC::GenEvent::vertex_const_iterator vitr= hepmcevt->vertices_begin();
      vitr != hepmcevt->vertices_end(); ++vitr ) { // loop for vertex ...

    // real vertex?
    G4bool qvtx=false;
    for (HepMC::GenVertex::particle_iterator 
	   pitr= (*vitr)->particles_begin(HepMC::children);
	 pitr != (*vitr)->particles_end(HepMC::children); ++pitr) {

      if (!(*pitr)->end_vertex() && (*pitr)->status()==1) {
	qvtx=true;
	break;
      }
    }
    if (!qvtx) continue;

    // check world boundary
    G4LorentzVector xvtx= (*vitr)-> position();
    if (! CheckVertexInsideWorld(xvtx.vect()*mm)) continue;

    // create G4PrimaryVertex and associated G4PrimaryParticles
    G4PrimaryVertex* g4vtx= 
      new G4PrimaryVertex(xvtx.x()*mm, xvtx.y()*mm, xvtx.z()*mm, 
			  xvtx.t()*mm/c_light);

    for (HepMC::GenVertex::particle_iterator 
	   vpitr= (*vitr)->particles_begin(HepMC::children);
	 vpitr != (*vitr)->particles_end(HepMC::children); ++vpitr) {

      if( (*vpitr)->status() != 1 ) continue;

      G4int pdgcode= (*vpitr)-> pdg_id();
      G4LorentzVector p= (*vpitr)-> momentum();
      G4PrimaryParticle* g4prim= 
	new G4PrimaryParticle(pdgcode, p.x()*GeV, p.y()*GeV, p.z()*GeV);

      g4vtx-> SetPrimary(g4prim);
    }
    g4event-> AddPrimaryVertex(g4vtx);
  } 
} 


///////////////////////////////////////////////////////
HepMC::GenEvent* G4HepMCInterface::GenerateHepMCEvent()
///////////////////////////////////////////////////////
{
  HepMC::GenEvent* aevent= new HepMC::GenEvent();
  return aevent;
}

//////////////////////////////////////////////////////////////
void G4HepMCInterface::GeneratePrimaryVertex(G4Event* anEvent)
//////////////////////////////////////////////////////////////
{
  // delete previous event object
  delete hepmcEvent;

  // generate next event
  hepmcEvent= GenerateHepMCEvent();
  if(! hepmcEvent) {
    G4cout << "HepMCInterface: no generated particles. run terminated..." 
           << G4endl;
    G4RunManager::GetRunManager()-> AbortRun();
    return;
  }
  HepMC2G4(hepmcEvent, anEvent);
}

