// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4QCHIPSWorld.cc,v 1.1 2000-09-04 07:46:39 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4QCHIPSWorld ----------------
//             by Mikhail Kossov, Sept 1999.
//      class for the CHIPS World definition in CHIPS Model
// -------------------------------------------------------------------

//#define debug
//#define pdebug

#include "G4QCHIPSWorld.hh"

G4QCHIPSWorld::G4QCHIPSWorld(G4int nOfParts)
{
  qWorld = InitCHIPSWorld(nOfParts); // Initialization of the CHIPS World with N particles
}

G4QCHIPSWorld::~G4QCHIPSWorld() {};

// Standard output for QPDGCode
ostream& operator<<(ostream& lhs, G4QCHIPSWorld& rhs)
//       ============================================
{
  // @@ Later make a list of activated particles and clusters
  //lhs << "[ PDG=" << rhs.GetPDGCode() << ", Q=" << rhs.GetQCode() << "]";
  return lhs;
}

//Initialize the CHIPS World of Quasmons
G4QParticleVector* G4QCHIPSWorld::InitCHIPSWorld(G4int nOfParts)
//                 =============================================
{
  static G4int mnofParts = 486;                 // max number of particles (up to A=80)
  static G4QParticleVector theWorld;            // *** A body of the CHIPS World *** 
  static int init = 0;                          // Initialization counter (private)
#ifdef debug
  G4cout<<"G4QCHIPSWorld::InitCHIPSWorld: n="<<nOfParts<<" particles, init="<<init<<G4endl;
#endif
  if(nOfParts>0)
  {
#ifdef debug
    G4cout<<"G4QCHIPSWorld::InitCHIPSWorld: Creating CHIPS World of nP="<<nOfParts<<G4endl;
#endif
    G4int curNP=theWorld.entries();
    if(curNP<0) curNP=0;
    if(!init++||nOfParts>curNP)                 // Initialize for increasing CHIPS World
    {
      if (nOfParts>mnofParts)
      {
        nOfParts=mnofParts;
        G4cerr<<"G4QCHIPSWorld::InitCHIPSWorld: nOfParts="<<nOfParts<<" >"<<mnofParts<<G4endl;
      }
      if (nOfParts<10) nOfParts=10;             // Minimal number of particles for Vacuum
      for (G4int i=curNP; i<nOfParts; i++) 
      {
        G4QParticle* curPart = new G4QParticle; // Created
        curPart->InitQParticle(i);              //   ||
        theWorld.insert(curPart);               // Filled
#ifdef debug
        G4cout<<"G4QCHIPSWorld::InitCHIPSWorld: Particle#"<<i<<"(of "<<nOfParts<<") done"<<endl;
#endif
      }
    }
    else if (nOfParts<0)
	{
      theWorld.clearAndDestroy();
      init=0;
	}
    else init--;
  }
  return &theWorld;
}





