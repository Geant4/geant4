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
// $Id$
//
//      ---------------- G4QCHIPSWorld ----------------
//             by Mikhail Kossov, Sept 1999.
//      class for the CHIPS World definition in CHIPS Model
// -------------------------------------------------------------------
// Short description: The CHIPS World is a world of elementary particles
// and nuclear fragments. This class is a singletone, but without fixed
// limits. E.g. the nuclear fragments as possible G4Candidates can be
// initialized in the CHIPS World only up to Be8 od C12 or other bigger
// nuclear fragment. If one need the heavy fragment production then the
// the CHIPS World must be initialized up to these fragments (see the
// CHIPS Manual), but the price in performans will be big, because in
// each act of the fragmentation competition these numerous candidates
// take place in the competition and the hadronization probability is
// calculated each time for each of them, so the Be8 limit (Be8->alpha+
// alpha decays very fast and contribute to the alpha-spectrum) is the
// most optimal.
// -------------------------------------------------------------------

//#define debug
//#define pdebug

#include "G4QCHIPSWorld.hh"

// Initialization of the CHIPSWorld Pointer
// G4QCHIPSWorld* G4QCHIPSWorld::aWorld =0;
// G4QParticleVector G4QCHIPSWorld::qWorld;

// Constructor
G4QCHIPSWorld::G4QCHIPSWorld()
{
}

G4QCHIPSWorld::~G4QCHIPSWorld()      // The CHIPS World is destructed only in the EndOfJob
{
  //G4int nP=GetQWorld().size();
  //G4cout<<"G4QCHIPSWorld::Destructor: Before nP="<<nP<<","<<GetQWorld()[0]<<G4endl; //TMP
  //if(nP) std::for_each(GetQWorld().begin(),GetQWorld().end(),DeleteQParticle());
  //G4cout<<"G4QCHIPSWorld::Destructor: After"<<G4endl; // TMP
  //GetQWorld().clear();
}

// Standard output for CHIPS World
std::ostream& operator<<(std::ostream& lhs, G4QCHIPSWorld& rhs)
{
  // @@ Later make a list of activated particles and clusters
  lhs << "[ Currently a#of particles in the CHIPS World = " << rhs.GetQPEntries() << "]";
  return lhs;
}

// Returns Pointer to the CHIPS World
G4QCHIPSWorld* G4QCHIPSWorld::Get()
{
  static G4QCHIPSWorld theWorld;                // *** Static body of the CHIPS World *** 
// Returns Pointer to the CHIPS World
//  if(!aWorld) aWorld=&theWorld;                 // Init the static pointer to CHIPS World
//  return aWorld;
  return &theWorld;
}

G4QParticleVector & G4QCHIPSWorld::GetQWorld()
{
  static G4QParticleVector theWorld;
  return theWorld;
}

// Return pointer to Particles of the CHIPS World
G4QParticleVector* G4QCHIPSWorld::GetParticles(G4int nOfParts)
{
  //static const G4int mnofParts = 486;           // max number of particles (up to A=80)
  //static const G4int mnofParts = 494;           // max number of particles (up to A=80)
  static const G4int mnofParts = 512; // max#of particles,A<57,G4QParticle::InitDecayVector
  static const G4bool cf = true;                // verbose=true G4QPDG construction flag
#ifdef debug
  G4cout<<"G4QCHIPSWorld::GetParticles: n="<<nOfParts<<" particles"<<G4endl;
#endif
  //if(GetQWorld().size())G4cout<<"G4QCHIPSWorld::GetPts:***Pt#0="<<GetQWorld()[0]<<G4endl;
  if(nOfParts>0)
  {
#ifdef debug
    G4cout<<"G4QCHIPSWorld::GetParticles: Creating CHIPS World of nP="<<nOfParts<<G4endl;
#endif
    G4int curNP=GetQWorld().size();
    //G4cout<<"G4QCHIPSWorld::GetParticles: Creating CHIPS World of curNP="<<curNP<<G4endl;
    if(nOfParts>curNP)                         // Initialization for increasing CHIPS World
    {
      if (nOfParts>mnofParts)
      {
        G4cerr<<"***G4QCHIPSWorld::GetPartics:nOfParts="<<nOfParts<<">"<<mnofParts<<G4endl;
        nOfParts=mnofParts;
      }
      if (nOfParts<10) nOfParts=10;            // Minimal number of particles for Vacuum(?)
#ifdef debug
      G4cout<<"G4QCHIPSWorld::GetParticles: n="<<nOfParts<<",c="<<curNP<<G4endl;
#endif
      for (G4int i=curNP; i<nOfParts; i++) 
      {
#ifdef debug
        G4cout<<"G4QCHIPSWorld::GetParticles: Create particle QCode="<<i<<G4endl;
#endif
        G4QParticle* curPart = new G4QParticle(cf,i); // Created with QCode=i
#ifdef debug
        G4cout<<"G4QCHIPSWorld::GetParticles: Particle QCode="<<i<<" is created."<<G4endl;
#endif
        //curPart->InitQParticle(i);             //
        //if(!i) G4cout<<"G4QCHIPSWorld::GetParticles:Pt#0="<<curPart<<G4endl;
        GetQWorld().push_back(curPart);           // Filled (forever but only once)
#ifdef debug
        G4cout<<"G4QCHIPSWorld::GetParticles: Pt#"<<i<<"("<<nOfParts<<") is done"<<G4endl;
#endif
      }
    }
    //else init--;//Recover theReInitializationCounter, if nothingWasAdded to theCHIPSWorld
  }
#ifdef debug
  G4cout<<"G4QCHIPSWorld::GetParticles: TotalPt#"<<GetQWorld().size()<<G4endl;
#endif
  return &GetQWorld();
}
