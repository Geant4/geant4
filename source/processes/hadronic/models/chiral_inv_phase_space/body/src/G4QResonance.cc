// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4QResonance.cc,v 1.1 1999-11-17 11:04:17 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4QResonance ----------------
//             by Mikhail Kossov, Sept 1999.
//      class for Mesonic P-wave Resonances used by CHIPS Model
// ------------------------------------------------------------

//#define debug

#include "G4QResonance.hh"

G4QResonance::G4QResonance(G4int PDGcode, G4double maxM) 
{
#ifdef debug
  cout<<"G4QResonance is called with PDG="<<PDGcode<<", maxM="<<maxM<<endl;
#endif
}

G4double G4QResonance::CalculateMass(G4double maxM, G4int PDG)
{
  G4ParticleDefinition* curDefinition = GetDefinition();
  G4double meanM = curDefinition->GetPDGMass();
  G4double width = curDefinition->GetPDGWidth()/2.;
  if(width<=0.)
  {
	cerr<<"***G4QResonance::CalculateMass width="<<width<<" <= 0, PDG="<<PDG<<endl;
	G4Exception("G4QResonance::CalculateMass(): width of the Resonance <= 0");
  }
  G4int absPDG = abs(PDG);
  G4double minM=0.;
  if      (absPDG==223) minM=135.;   // (omega)=>PI0
  else if (absPDG==113) minM=279.12; // (rho0) =>2*mPI+
  else if (absPDG==333) minM=414.12; // (phi)  =>3PI
  else if (absPDG==213) minM=274.55; // (rho+) =>mPI0+mPI+
  else if (absPDG==313) minM=632.66; // (K0*)  =>K0+PI0
  else if (absPDG==323) minM=628.66; // (K+*)  =>PI0+K+
  else
  {
	cerr<<"***G4QResonance::CalculateMass unknown Resonance PDG="<<PDG<<endl;
	G4Exception("G4QResonance::CalculateMass(): unknown Resonance");
  }
  //Now calculate the Breit-Wigner distribution with two cuts
  G4double v1=atan((minM-meanM)/width);
  G4double v2=atan((maxM-meanM)/width);
  G4double dv=v2-v1;
#ifdef debug
  cout << "QRes::CalcMass: vMin=" << v1 << ", vMax=" << v2 << ", dv=" << dv << endl;
#endif
  return meanM+width*tan(v1+dv*G4UniformRand());
}

