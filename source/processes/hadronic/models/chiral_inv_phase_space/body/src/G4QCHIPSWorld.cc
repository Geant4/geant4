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
// $Id: G4QCHIPSWorld.cc,v 1.19 2002-12-12 19:14:33 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
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

G4QCHIPSWorld::G4QCHIPSWorld(const G4QCHIPSWorld& right)
{
  qWorld = right.qWorld;
}

G4QCHIPSWorld::G4QCHIPSWorld(G4QCHIPSWorld* right)
{
  qWorld = right->qWorld;
}

G4QCHIPSWorld::~G4QCHIPSWorld() {}

const G4QCHIPSWorld& G4QCHIPSWorld::operator=(const G4QCHIPSWorld &right)
{//   ===================================================================
  qWorld = right.qWorld;

  return *this;
}

// Standard output for CHIPS World
G4std::ostream& operator<<(G4std::ostream& lhs, G4QCHIPSWorld& rhs)
//       ============================================
{
  // @@ Later make a list of activated particles and clusters
  lhs << "[ a#of particles in the CHIPS World =" << rhs.GetQPEntries() << "]";
  return lhs;
}

//Initialize the CHIPS World of Quasmons
G4QParticleVector* G4QCHIPSWorld::InitCHIPSWorld(G4int nOfParts)
//                 =============================================
{
  static G4int mnofParts = 486;                 // max number of particles (up to A=80)
  static G4QParticleVector theWorld;            // *** A body of the CHIPS World *** 
  static G4int init = 0;                        // Initialization counter (private)
#ifdef debug
  G4cout<<"G4QCHIPSWorld::InitCHIPSWorld: n="<<nOfParts<<" particles, init="<<init<<G4endl;
#endif
  if(nOfParts>0)
  {
#ifdef debug
    G4cout<<"G4QCHIPSWorld::InitCHIPSWorld: Creating CHIPS World of nP="<<nOfParts<<G4endl;
#endif
    G4int curNP=theWorld.size();
    if(curNP<0) curNP=0;
    if(!init++||nOfParts>curNP)                 // Initialization for increasing CHIPS World
    {
      if (nOfParts>mnofParts)
      {
        G4cerr<<"***G4QCHIPSWorld::InitCHIPSWorld: nOfParts="<<nOfParts<<" >"<<mnofParts<<G4endl;
        nOfParts=mnofParts;
      }
      if (nOfParts<10) nOfParts=10;             // Minimal number of particles for Vacuum
#ifdef debug
      G4cerr<<"G4QCHIPSWorld::InitCHIPSWorld:Should be once, i="<<init<<",n="<<nOfParts<<",c="<<curNP<<G4endl;
#endif
      for (G4int i=curNP; i<nOfParts; i++) 
      {
        G4QParticle* curPart = new G4QParticle; // Created
        curPart->InitQParticle(i);              //   ||
        theWorld.push_back(curPart);            // Filled (forever but only once)
#ifdef debug
        G4cout<<"G4QCHIPSWorld::InitCHIPSWorld: Particle#"<<i<<"(of "<<nOfParts<<") done"<<endl;
#endif
      }
    }
    else init--;                                // Recover the init pointer if nothing was done
  }
  else if (nOfParts<0)                          // Cleaning up the CHIPS Word (a possibility)
  {
    G4std::for_each(theWorld.begin(), theWorld.end(), DeleteQParticle());
    theWorld.clear();
    init=0;
  }
  return &theWorld;
}





