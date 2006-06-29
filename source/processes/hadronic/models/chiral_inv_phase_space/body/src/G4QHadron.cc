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
// $Id: G4QHadron.cc,v 1.42 2006-06-29 20:07:01 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//      ---------------- G4QHadron ----------------
//             by Mikhail Kossov, Sept 1999.
//      class for Quasmon initiated Hadrons generated by CHIPS Model
// ------------------------------------------------------------------
//
//#define debug
//#define pdebug
//#define sdebug

#include "G4QHadron.hh"
#include <cmath>
using namespace std;

G4QHadron::G4QHadron() :
  theQPDG(0),theMomentum(0.,0.,0.,0.),valQ(0,0,0,0,0,0),nFragm(0)
{}

G4QHadron::G4QHadron(G4LorentzVector p) :
  theQPDG(0),theMomentum(p),valQ(0,0,0,0,0,0),nFragm(0)
{}

// For Chipolino or Quasmon doesn't make any sense
G4QHadron::G4QHadron(G4int PDGCode, G4LorentzVector p) :
  theQPDG(PDGCode),theMomentum(p),nFragm(0)
{
#ifdef debug
  G4cout<<"G4QHadron must be created with PDG="<<PDGCode<<", 4M="<<p<<G4endl;
#endif
  if(GetQCode()>-1)
  {
    if(theMomentum.e()==0.) theMomentum.setE(theQPDG.GetMass());
    valQ=theQPDG.GetQuarkContent();
  }
  else if(PDGCode>80000000) DefineQC(PDGCode);
  else G4cerr<<"***G4QHadron:(P) PDG="<<PDGCode<<", use other constructor"<<G4endl;
#ifdef debug
  G4cout<<"G4QHadron is created with QCode="<<GetQCode()<<", QC="<<valQ<<G4endl;
#endif
}

// For Chipolino or Quasmon doesn't make any sense
G4QHadron::G4QHadron(G4QPDGCode QPDG, G4LorentzVector p) :
  theQPDG(QPDG),theMomentum(p),nFragm(0)
{
  if(theQPDG.GetQCode()>-1)
  {
    if(theMomentum.e()==0.) theMomentum.setE(theQPDG.GetMass());
    valQ=theQPDG.GetQuarkContent();
  }
  else
  {
    G4int cPDG=theQPDG.GetPDGCode();
    if(cPDG>80000000) DefineQC(cPDG);
	   else G4cerr<<"***G4QHadr:(QP) PDG="<<cPDG<<" use other constructor"<<G4endl;
  }
}

// Make sense Chipolino or Quasmon
G4QHadron::G4QHadron(G4QContent QC, G4LorentzVector p) :
  theQPDG(0),theMomentum(p),valQ(QC),nFragm(0)
{
  G4int curPDG=valQ.GetSPDGCode();
  if(curPDG==10&&valQ.GetBaryonNumber()>0) curPDG=valQ.GetZNSPDGCode();
  if(curPDG&&curPDG!=10) theQPDG.SetPDGCode(curPDG);
  else theQPDG.InitByQCont(QC);
}

G4QHadron::G4QHadron(G4int PDGCode, G4double aMass, G4QContent QC) :
  theQPDG(PDGCode),theMomentum(0.,0.,0., aMass),valQ(QC),nFragm(0)
{}

G4QHadron::G4QHadron(G4QPDGCode QPDG, G4double aMass, G4QContent QC) :
  theQPDG(QPDG),theMomentum(0.,0.,0., aMass),valQ(QC),nFragm(0)
{}

G4QHadron::G4QHadron(G4int PDGCode, G4LorentzVector p, G4QContent QC) :
  theQPDG(PDGCode),theMomentum(p),valQ(QC),nFragm(0)
{}

G4QHadron::G4QHadron(G4QPDGCode QPDG, G4LorentzVector p, G4QContent QC) :
  theQPDG(QPDG),theMomentum(p),valQ(QC),nFragm(0)
{}

G4QHadron::G4QHadron(G4QParticle* pPart, G4double maxM) :
  theQPDG(pPart->GetQPDG()),theMomentum(0.,0.,0.,0.),nFragm(0)
{
#ifdef debug
  G4cout<<"G4QHadron is created & randomized with maxM="<<maxM<<G4endl;
#endif
  G4int PDGCode = theQPDG.GetPDGCode();
  if(PDGCode<2)G4cerr<<"***G4QHadron:(M) PDGC="<<PDGCode<<" use other constructor"<<G4endl;
  valQ=theQPDG.GetQuarkContent();
  theMomentum.setE(RandomizeMass(pPart, maxM));
}

G4QHadron::G4QHadron(const G4QHadron& right)
{
  theMomentum         = right.theMomentum;
  theQPDG             = right.theQPDG;
  valQ                = right.valQ;
  nFragm              = right.nFragm;
}

G4QHadron::G4QHadron(const G4QHadron* right)
{
  theMomentum         = right->theMomentum;
  theQPDG             = right->theQPDG;
  valQ                = right->valQ;
  nFragm              = right->nFragm;
}

const G4QHadron& G4QHadron::operator=(const G4QHadron &right)
{
  theMomentum         = right.theMomentum;
  theQPDG             = right.theQPDG;
  valQ                = right.valQ;
  nFragm              = right.nFragm;

  return *this;
}

G4QHadron::~G4QHadron() {}

// Decay of Hadron In2Particles f&s, f is in respect to the direction of HadronMomentumDir
G4bool G4QHadron::RelDecayIn2(G4LorentzVector& f4Mom, G4LorentzVector& s4Mom,
       G4LorentzVector& dir, G4double maxCost, G4double minCost)
{//    ===================================================================
  G4double fM2 = f4Mom.m2();
  G4double fM  = sqrt(fM2);              // Mass of the 1st Hadron
  G4double sM2 = s4Mom.m2();
  G4double sM  = sqrt(sM2);              // Mass of the 2nd Hadron
  G4double iM2 = theMomentum.m2();
  G4double iM  = sqrt(iM2);              // Mass of the decaying hadron
  G4double vP  = theMomentum.rho();      // Momentum of the decaying hadron
  G4double dE  = theMomentum.e();        // Energy of the decaying hadron
  if(dE<vP)
  {
    G4cerr<<"***G4QHad::RelDecIn2: Tachionic 4-mom="<<theMomentum<<", E-p="<<dE-vP<<G4endl;
    G4double accuracy=.000001*vP;
    G4double emodif=abs(dE-vP);
    if(emodif<accuracy)
				{
      G4cerr<<"G4QHadron::RelDecIn2: *Boost* E-p shift is corrected to "<<emodif<<G4endl;
      theMomentum.setE(vP+emodif);
    }
  }
  G4ThreeVector ltb = theMomentum.boostVector();// Boost vector for backward Lorentz Trans.
  G4ThreeVector ltf = -ltb;              // Boost vector for forward Lorentz Trans.
  G4LorentzVector cdir = dir;            // A copy to make a transformation to CMS
#ifdef debug
  if(cdir.e()+.001<cdir.rho()) G4cerr<<"*G4QH::RDIn2:*Boost* cd4M="<<cdir<<",e-p="
                                     <<cdir.e()-cdir.rho()<<G4endl;
#endif
  cdir.boost(ltf);                       // Direction transpormed to CMS of the Momentum
  G4ThreeVector vdir = cdir.vect();      // 3-Vector of the direction-particle
#ifdef debug
  G4cout<<"G4QHad::RelDI2:dir="<<dir<<",ltf="<<ltf<<",cdir="<<cdir<<",vdir="<<vdir<<G4endl;
#endif
  G4ThreeVector vx(0.,0.,1.);            // Ort in the direction of the reference particle
  G4ThreeVector vy(0.,1.,0.);            // First ort orthogonal to the direction
  G4ThreeVector vz(1.,0.,0.);            // Second ort orthoganal to the direction
  if(vdir.mag2() > 0.)                   // the refference particle isn't at rest in CMS
  {
    vx = vdir.unit();                    // Ort in the direction of the reference particle
    G4ThreeVector vv= vx.orthogonal();   // Not normed orthogonal vector (!)
    vy = vv.unit();                      // First ort orthogonal to the direction
    vz = vx.cross(vy);                   // Second ort orthoganal to the direction
  }
#ifdef debug
  G4cout<<"G4QHad::RelDecIn2:iM="<<iM<<"=>fM="<<fM<<"+sM="<<sM<<",ob="<<vx<<vy<<vz<<G4endl;
#endif
  if(maxCost> 1.) maxCost= 1.;
  if(minCost<-1.) minCost=-1.;
  if(maxCost<-1.) maxCost=-1.;
  if(minCost> 1.) minCost= 1.;
  if (fabs(iM-fM-sM)<.001)
  {
    G4double fR=fM/iM;
    G4double sR=sM/iM;
    f4Mom=fR*theMomentum;
    s4Mom=sR*theMomentum;
    return true;
  }
  else if (iM+.001<fM+sM || iM==0.)
  {//@@ Later on make a quark content check for the decay
    G4cerr<<"***G4QH::RelDecIn2: fM="<<fM<<"+sM="<<sM<<">iM="<<iM<<",d="<<iM-fM-sM<<G4endl;
    return false;
  }
  G4double d2 = iM2-fM2-sM2;
  G4double p2 = (d2*d2/4.-fM2*sM2)/iM2;    // Decay momentum(^2) in CMS of Quasmon
  if(p2<0.)
  {
#ifdef debug
    G4cout<<"**G4QH:RDIn2:p2="<<p2<<"<0,d2^2="<<d2*d2/4.<<"<4*fM2*sM2="<<4*fM2*sM2<<G4endl;
#endif
    p2=0.;
  }
  G4double p  = sqrt(p2);
  G4double ct = maxCost;
  if(maxCost>minCost)
  {
    G4double dcost=maxCost-minCost;
    ct = minCost+dcost*G4UniformRand();
  }
  G4double phi= 360.*deg*G4UniformRand();  // @@ Change 360.*deg to M_TWOPI (?)
  G4double ps = p * sqrt(1.-ct*ct);
  G4ThreeVector pVect=(ps*sin(phi))*vz+(ps*cos(phi))*vy+p*ct*vx;
#ifdef debug
  G4cout<<"G4QH::RelDIn2:ct="<<ct<<",p="<<p<<",ps="<<ps<<",ph="<<phi<<",v="<<pVect<<G4endl;
#endif

  f4Mom.setVect(pVect);
  f4Mom.setE(sqrt(fM2+p2));
  s4Mom.setVect((-1)*pVect);
  s4Mom.setE(sqrt(sM2+p2));
  
#ifdef debug
  G4cout<<"G4QHadr::RelDecIn2:p2="<<p2<<",v="<<ltb<<",f4M="<<f4Mom<<",s4M="<<s4Mom<<G4endl;
#endif
  if(f4Mom.e()+.001<f4Mom.rho())G4cerr<<"*G4QH::RDIn2:*Boost* f4M="<<f4Mom<<",e-p="
                                      <<f4Mom.e()-f4Mom.rho()<<G4endl;
  f4Mom.boost(ltb);                        // Lor.Trans. of 1st hadron back to LS
  if(s4Mom.e()+.001<s4Mom.rho())G4cerr<<"*G4QH::RDIn2:*Boost* s4M="<<s4Mom<<",e-p="
                                      <<s4Mom.e()-s4Mom.rho()<<G4endl;
  s4Mom.boost(ltb);                        // Lor.Trans. of 2nd hadron back to LS
#ifdef debug
  G4cout<<"G4QHadron::RelDecayIn2: ROOT OUTPUT f4Mom="<<f4Mom<<", s4Mom="<<s4Mom<<G4endl;
#endif
  return true;
} // End of "RelDecayIn2"

// Decay of the Hadron in 2 particles (f + s)
G4bool G4QHadron::DecayIn2(G4LorentzVector& f4Mom, G4LorentzVector& s4Mom)
{//    ===================================================================
  G4double fM  = f4Mom.m();                // Mass of the 1st Hadron
  G4double fM2 = f4Mom.m2();
  G4double sM  = s4Mom.m();                // Mass of the 2nd Hadron
  G4double sM2 = s4Mom.m2();
  G4double iM  = theMomentum.m();          // Mass of the decaying hadron
  G4double iM2 = theMomentum.m2();
#ifdef debug
  G4cout<<"G4QHadron::DecIn2: iM="<<iM<<" => fM="<<fM<<" + sM="<<sM<<" = "<<fM+sM<<G4endl;
#endif
  //@@ Later on make a quark content check for the decay
  if (fabs(iM-fM-sM)<.001)
  {
    G4double fR=fM/iM;
    G4double sR=sM/iM;
    f4Mom=fR*theMomentum;
    s4Mom=sR*theMomentum;
    return true;
  }
  else if (iM+.001<fM+sM || iM==0.)
  {
#ifdef pdebug
    G4cerr<<"***G4QHadron::DecayIn2*** fM="<<fM<<" + sM="<<sM<<"="<<fM+sM<<" > iM="<<iM
          <<", d="<<iM-fM-sM<<G4endl;
#endif
    return false;
  }

  G4double d2 = iM2-fM2-sM2;
  G4double p2 = (d2*d2/4.-fM2*sM2)/iM2;    // Decay momentum(^2) in CMS of Quasmon
  if (p2<0.)
  {
#ifdef debug
    G4cerr<<"***G4QH::DI2:p2="<<p2<<"<0,d2^2="<<d2*d2/4.<<"<4*fM2*sM2="<<4*fM2*sM2<<G4endl;
#endif
    p2=0.;
  }
  G4double p  = sqrt(p2);
  G4double ct = 1.-2*G4UniformRand();
#ifdef debug
  G4cout<<"G4QHadron::DecayIn2: ct="<<ct<<", p="<<p<<G4endl;
#endif
  G4double phi= 360.*deg*G4UniformRand();  // @@ Change 360.*deg to M_TWOPI (?)
  G4double ps = p * sqrt(1.-ct*ct);
  G4ThreeVector pVect(ps*sin(phi),ps*cos(phi),p*ct);

  f4Mom.setVect(pVect);
  f4Mom.setE(sqrt(fM2+p2));
  s4Mom.setVect((-1)*pVect);
  s4Mom.setE(sqrt(sM2+p2));

  if(theMomentum.e()<theMomentum.rho())
  {
	   G4cerr<<"*G4QH::DecIn2:*Boost* 4M="<<theMomentum<<",e-p="
          <<theMomentum.e()-theMomentum.rho()<<G4endl;
	//throw G4QException("G4QHadron::DecayIn2: Decay of particle with zero mass");
  }
  G4ThreeVector ltb = theMomentum.boostVector(); // Boost vector for backward Lor.Trans.
#ifdef pdebug
  G4cout<<"G4QHadron::DecIn2:LorTrans v="<<ltb<<",f4Mom="<<f4Mom<<",s4Mom="<<s4Mom<<G4endl;
#endif
  if(f4Mom.e()+.001<f4Mom.rho())G4cerr<<"*G4QH::DecIn2:*Boost* f4M="<<f4Mom<<G4endl;
  f4Mom.boost(ltb);                        // Lor.Trans. of 1st hadron back to LS
  if(s4Mom.e()+.001<s4Mom.rho())G4cerr<<"*G4QH::DecIn2:*Boost* s4M="<<s4Mom<<G4endl; 
  s4Mom.boost(ltb);                        // Lor.Trans. of 2nd hadron back to LS
#ifdef pdebug
  G4cout<<"G4QHadron::DecayIn2: ROOT OUTPUT f4Mom="<<f4Mom<<", s4Mom="<<s4Mom<<G4endl;
#endif
  return true;
} // End of "DecayIn2"

// Correction for the Hadron + fr decay in case of the new corM mass of the Hadron
G4bool G4QHadron::CorMDecayIn2(G4double corM, G4LorentzVector& fr4Mom)
{//    ===============================================================
  G4double fM  = fr4Mom.m();                // Mass of the Fragment
  G4LorentzVector comp=theMomentum+fr4Mom;  // 4Mom of the decaying compound system
  G4double iM  = comp.m();                  // mass of the decaying compound system
#ifdef pdebug
  G4cout<<"G4QH::CMDIn2: iM="<<iM<<comp<<"=>fM="<<fM<<"+corM="<<corM<<"="<<fM+corM<<G4endl;
#endif
  G4double dE=iM-fM-corM;
  //@@ Later on make a quark content check for the decay
  if (fabs(dE)<.001)
  {
    G4double fR=fM/iM;
    G4double cR=corM/iM;
    fr4Mom=fR*comp;
    theMomentum=cR*comp;
    return true;
  }
  else if (dE<-.001 || iM==0.)
  {
    G4cerr<<"***G4QH::CorMDIn2***fM="<<fM<<" + cM="<<corM<<" > iM="<<iM<<",d="<<dE<<G4endl;
    return false;
  }
  G4double corM2= corM*corM;
  G4double fM2 = fM*fM;
  G4double iM2 = iM*iM;
  G4double d2 = iM2-fM2-corM2;
  G4double p2 = (d2*d2/4.-fM2*corM2)/iM2;    // Decay momentum(^2) in CMS of Quasmon
  if (p2<0.)
  {
#ifdef pdebug
    G4cerr<<"**G4QH::CMDI2:p2="<<p2<<"<0,d="<<d2*d2/4.<<"<4*fM2*hM2="<<4*fM2*corM2<<G4endl;
#endif
    p2=0.;
  }
  G4double p  = sqrt(p2);
  if(comp.e()<comp.rho())G4cerr<<"*G4QH::CorMDecayIn2:*Boost* comp4M="<<comp<<",e-p="
                               <<comp.e()-comp.rho()<<G4endl;
  G4ThreeVector ltb = comp.boostVector();      // Boost vector for backward Lor.Trans.
  G4ThreeVector ltf = -ltb;                    // Boost vector for forward Lorentz Trans.
  G4LorentzVector cm4Mom=fr4Mom;               // Copy of fragment 4Mom to transform to CMS
  if(cm4Mom.e()+.001<cm4Mom.rho())G4cerr<<"*G4QH::CorMDecIn2:*Boost* c4M="<<cm4Mom<<G4endl;
  cm4Mom.boost(ltf);                           // Now it is in CMS (Forward Lor.Trans.)
  G4double pfx= cm4Mom.px();
  G4double pfy= cm4Mom.py();
  G4double pfz= cm4Mom.pz();
  G4double pt2= pfx*pfx+pfy*pfy;
  G4double tx=0.;
  G4double ty=0.;
  if(pt2<=0.)
  {
    G4double phi= 360.*deg*G4UniformRand();  // @@ Change 360.*deg to M_TWOPI (?)
    tx=sin(phi);
    ty=cos(phi);
  }
  else
  {
    G4double pt=sqrt(pt2);
    tx=pfx/pt;
    ty=pfy/pt;
  }
  G4double pc2=pt2+pfz*pfz;
  G4double ct=0.;
  if(pc2<=0.)
  {
    G4double rnd= G4UniformRand();
    ct=1.-rnd-rnd;
  }
  else
  {
    G4double pc=sqrt(pc2);
    ct=pfz/pc;
  }
#ifdef debug
  G4cout<<"G4QHadron::CorMDecayIn2: ct="<<ct<<", p="<<p<<G4endl;
#endif
  G4double ps = p * sqrt(1.-ct*ct);
  G4ThreeVector pVect(ps*tx,ps*ty,p*ct);
  fr4Mom.setVect(pVect);
  fr4Mom.setE(sqrt(fM2+p2));
  theMomentum.setVect((-1)*pVect);
  theMomentum.setE(sqrt(corM2+p2));
#ifdef pdebug
  G4LorentzVector dif2=comp-fr4Mom-theMomentum;
  G4cout<<"G4QH::CorMDIn2:c="<<comp<<"-f="<<fr4Mom<<"-4M="<<theMomentum<<"="<<dif2<<G4endl;
#endif
  if(fr4Mom.e()+.001<fr4Mom.rho())G4cerr<<"*G4QH::CorMDecIn2:*Boost*fr4M="<<fr4Mom<<G4endl;
  fr4Mom.boost(ltb);                        // Lor.Trans. of the Fragment back to LS
  if(theMomentum.e()+.001<theMomentum.rho())G4cerr<<"*G4QH::CMDI2:4="<<theMomentum<<G4endl;
  theMomentum.boost(ltb);                  // Lor.Trans. of the Hadron back to LS
#ifdef pdebug
  G4LorentzVector dif3=comp-fr4Mom-theMomentum;
  G4cout<<"G4QH::CorMDecIn2:OUTPUT:f4M="<<fr4Mom<<",h4M="<<theMomentum<<"d="<<dif3<<G4endl;
#endif
  return true;
} // End of "CorMDecayIn2"


// Fragment fr4Mom louse energy corE and transfer it to This Hadron 
G4bool G4QHadron::CorEDecayIn2(G4double corE, G4LorentzVector& fr4Mom)
{//    ===============================================================
  G4double fE  = fr4Mom.m();                // Energy of the Fragment
#ifdef pdebug
  G4cout<<"G4QH::CorEDecIn2:fE="<<fE<<fr4Mom<<">corE="<<corE<<",h4M="<<theMomentum<<G4endl;
#endif
  if (fE+.001<=corE)
  {
#ifdef pdebug
    G4cerr<<"***G4QHadron::CorEDecIn2*** fE="<<fE<<"<corE="<<corE<<", d="<<corE-fE<<G4endl;
#endif
    return false;
  }
  G4double fM2=fr4Mom.m2();                 // Squared Mass of the Fragment
  G4double iPx=fr4Mom.px();                 // Initial Px of the Fragment
  G4double iPy=fr4Mom.py();                 // Initial Py of the Fragment
  G4double iPz=fr4Mom.pz();                 // Initial Pz of the Fragment
  G4double fP2=iPx*iPx+iPy*iPy+iPz*iPz;     // Initial Squared 3-momentum of the Fragment
  G4double finE = fE - corE;                // Final energy of the fragment
  G4double rP = sqrt((finE*finE-fM2)/fP2);  // Reduction factor for the momentum
  G4double fPx=iPx*rP;
  G4double fPy=iPy*rP;
  G4double fPz=iPz*rP;
  fr4Mom= G4LorentzVector(fPx,fPy,fPz,finE);
  G4double Px=theMomentum.px()+iPx-fPx;
  G4double Py=theMomentum.py()+iPy-fPy;
  G4double Pz=theMomentum.pz()+iPz-fPz;
  G4double mE=theMomentum.e();
  ///////////G4double mM2=theMomentum.m2();
  theMomentum= G4LorentzVector(Px,Py,Pz,mE+corE);
#ifdef pdebug
  G4double difF=fr4Mom.m2()-fM2;
  G4cout<<"G4QH::CorEDecIn2: dF="<<difF<<",out:"<<theMomentum<<fr4Mom<<G4endl;
#endif
  return true;
} // End of "CorEDecayIn2"

// Decay of the hadron in 3 particles i=>r+s+t
G4bool G4QHadron::DecayIn3
                   (G4LorentzVector& f4Mom, G4LorentzVector& s4Mom, G4LorentzVector& t4Mom)
{//    ====================================================================================
#ifdef debug
  G4cout<<"G4QH::DIn3:"<<theMomentum<<"=>pf="<<f4Mom<<"+ps="<<s4Mom<<"+pt="<<t4Mom<<G4endl;
#endif
  G4double iM  = theMomentum.m();  // Mass of the decaying hadron
  G4double fM  = f4Mom.m();        // Mass of the 1st hadron
  G4double sM  = s4Mom.m();        // Mass of the 2nd hadron
  G4double tM  = t4Mom.m();        // Mass of the 3rd hadron
  G4double eps = 0.001;            // Accuracy of the split condition
  if (fabs(iM-fM-sM-tM)<=eps)
  {
    G4double fR=fM/iM;
    G4double sR=sM/iM;
    G4double tR=tM/iM;
    f4Mom=fR*theMomentum;
    s4Mom=sR*theMomentum;
    t4Mom=tR*theMomentum;
    return true;
  }
  if (iM+eps<fM+sM+tM)
  {
    G4cout<<"***G4QHadron::DecayIn3:fM="<<fM<<" + sM="<<sM<<" + tM="<<tM<<" > iM="<<iM
          <<",d="<<iM-fM-sM-tM<<G4endl;
    return false;
  }
  G4double fM2 = fM*fM;
  G4double sM2 = sM*sM;
  G4double tM2 = tM*tM;
  G4double iM2 = iM*iM;
  G4double m13sBase=(iM-sM)*(iM-sM)-(fM+tM)*(fM+tM);
  G4double m12sMin =(fM+sM)*(fM+sM);
  G4double m12sBase=(iM-tM)*(iM-tM)-m12sMin;
  G4double rR = 0.;
  G4double rnd= 1.;
#ifdef debug
  G4int    tr = 0;                 //@@ Comment if "cout" below is skiped @@
#endif
  G4double m12s = 0.;              // Fake definition before the Loop
  while (rnd > rR)
  {
    m12s = m12sMin + m12sBase*G4UniformRand();
    G4double e1=m12s+fM2-sM2;
    G4double e2=iM2-m12s-tM2;
    G4double four12=4*m12s;
    G4double m13sRange=0.;
    G4double dif=(e1*e1-four12*fM2)*(e2*e2-four12*tM2);
    if(dif<0.)
	   {
#ifdef debug
      if(dif<-.01) G4cerr<<"*G4QHadron::DecayIn3:iM="<<iM<<",tM="<<tM<<",sM="<<sM<<",fM="
                         <<fM<<",m12(s+f)="<<sqrt(m12s)<<", d="<<iM-fM-sM-tM<<G4endl;
#endif
    }
    else m13sRange=sqrt(dif)/m12s;
    rR = m13sRange/m13sBase;
    rnd= G4UniformRand();
#ifdef debug
    G4cout<<"G4QHadron::DecayIn3: try to decay #"<<++tr<<", rR="<<rR<<",rnd="<<rnd<<G4endl;
#endif
  }
  G4double m12 = sqrt(m12s);       // Mass of the H1+H2 system
  G4LorentzVector dh4Mom(0.,0.,0.,m12);
  
  if(!DecayIn2(t4Mom,dh4Mom))
  {
    G4cerr<<"***G4QHadron::DecayIn3: Exception1"<<G4endl;
	   //throw G4QException("G4QHadron::DecayIn3(): DecayIn2 did not succeed");
    return false;
  }
#ifdef debug
  G4cout<<"G4QHadron::DecayIn3: Now the last decay of m12="<<dh4Mom.m()<<G4endl;
#endif
  if(!G4QHadron(dh4Mom).DecayIn2(f4Mom,s4Mom))
  {
    G4cerr<<"***G4QHadron::DecayIn3: Error in DecayIn2 -> Exception2"<<G4endl;
	   //throw G4QException("G4QHadron::DecayIn3(): DecayIn2 did not succeed");
    return false;
  }
  return true;
}

// Randomize particle mass taking into account the width
G4double G4QHadron::RandomizeMass(G4QParticle* pPart, G4double maxM)
//       ===========================================================
{
  G4double meanM = theQPDG.GetMass();
  G4double width = theQPDG.GetWidth()/2.;
#ifdef debug
  G4cout<<"G4QHadron::RandomizeMass: meanM="<<meanM<<", halfWidth="<<width<<G4endl;
#endif
  if(maxM<meanM-3*width) 
  {
#ifdef debug
    G4cout<<"***G4QH::RandM:m=0 maxM="<<maxM<<"<meanM="<<meanM<<"-3*halfW="<<width<<G4endl;
#endif
    return 0.;
  }
  ///////////////G4double theMass  = 0.;
  if(width==0.)
  {
#ifdef debug
	   if(meanM>maxM) G4cerr<<"***G4QHadron::RandM:Stable m="<<meanM<<">maxM="<<maxM<<G4endl;
#endif
    return meanM;
    //return 0.;
  }
  else if(width<0.)
  {
	   G4cerr<<"***G4QHadron::RandM: width="<<width<<"<0,PDGC="<<theQPDG.GetPDGCode()<<G4endl;
	   throw G4QException("G4QHadron::RandomizeMass: with the width of the Hadron < 0.");
  }
  G4double minM = pPart->MinMassOfFragm();
  if(minM>maxM)
  {
#ifdef debug
	   G4cout<<"***G4QHadron::RandomizeMass:for PDG="<<theQPDG.GetPDGCode()<<" minM="<<minM
          <<" > maxM="<<maxM<<G4endl;
#endif
    return 0.;
  }
  //Now calculate the Breit-Wigner distribution with two cuts
  G4double v1=atan((minM-meanM)/width);
  G4double v2=atan((maxM-meanM)/width);
  G4double dv=v2-v1;
#ifdef debug
  G4cout<<"G4QHadr::RandM:Mi="<<minM<<",i="<<v1<<",Ma="<<maxM<<",a="<<v2<<","<<dv<<G4endl;
#endif
  return meanM+width*tan(v1+dv*G4UniformRand());
}



