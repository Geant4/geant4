// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4QDecayChan.cc,v 1.1 2000-08-17 13:55:49 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4QDecayChan ----------------
//             by Mikhail Kossov, Sept 1999.
//      class for Decay Channels of Hadrons in CHIPS Model
// -------------------------------------------------------------------
 
//#define debug
//#define pdebug

#include "G4QDecayChanVector.hh"

//G4QDecayChan::G4QDecayChan(){};

G4QDecayChan::G4QDecayChan(G4double pLev, G4int PDG1, G4int PDG2, G4int PDG3):
  aDecayChanLimit(pLev)
{
  G4QPDGCode* firstPDG = new G4QPDGCode(PDG1);
  theMinMass =firstPDG->GetMass();
  aVecOfSecHadrons.insert(firstPDG);
  G4QPDGCode* secondPDG = new G4QPDGCode(PDG2);
  theMinMass+=secondPDG->GetMass();
  aVecOfSecHadrons.insert(secondPDG);
  if(PDG3) 
  {
    G4QPDGCode* thirdPDG = new G4QPDGCode(PDG3);
    theMinMass+=thirdPDG->GetMass();
    aVecOfSecHadrons.insert(thirdPDG);
  }
#ifdef debug
  cout<<"G4QDecayChan is defined with pL="<<pLev<<",1="<<PDG1<<",2="<<PDG2<<",3="<<PDG3
      <<",m1="<<firstPDG->GetMass()<<",m2="<<secondPDG->GetMass()<<",minM="<<theMinMass<<endl;
#endif
}

G4QDecayChan::~G4QDecayChan() {aVecOfSecHadrons.clearAndDestroy();}

// Standard output for QDecayChan
ostream& operator<<(ostream& lhs, G4QDecayChan& rhs)
//       =========================================
{
  lhs << "[L=" << rhs.GetDecayChanLimit(); 
  G4QPDGCodeVector VSH = rhs.GetVecOfSecHadrons();
  G4int n = VSH.entries();
  lhs << ", N=" << n << ": ";
  for (int i=0; i<n; i++)
  {
    if(!i) lhs << ":";
    else   lhs << ",";
    lhs << VSH[i]->GetPDGCode();
  }
  lhs << "]";
  return lhs;
}

//G4int G4QDecayChan::() // Prototype
//    ===========================================
//{
//}





