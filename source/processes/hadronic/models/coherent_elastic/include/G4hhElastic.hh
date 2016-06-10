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
// $Id: G4hhElastic.hh,v 1.21 2010-11-09 09:04:29 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// G4 Model: hadron diffraction elastic scattering with 4-momentum balance 
//
// Class Description
// Final state production model for hadron-hadron elastic scattering 
// in the framework of quark-diquark model with springy Pomeron. 
// Projectiles are proton, neutron, pions, kaons. 
// Targets are proton (and neutron). 
// Class Description - End
//
// 02.05.14 V. Grichine - 1-st implementation
// 10.10.14 V. Grichine - change to combine with low mass diffraction

#ifndef G4hhElastic_h
#define G4hhElastic_h 1
 
#include "globals.hh"
#include <complex>
#include "G4Integrator.hh"

#include "G4HadronElastic.hh"
#include "G4HadProjectile.hh"
#include "G4Nucleus.hh"
#include "G4HadronNucleonXsc.hh"

#include "G4Exp.hh"
#include "G4Log.hh"


class G4ParticleDefinition;
class G4PhysicsTable;
class G4PhysicsLogVector;

class G4hhElastic  : public G4HadronElastic  
{
public:
  // PL constructor
  G4hhElastic();
  // test constructor
  G4hhElastic( G4ParticleDefinition* target, G4ParticleDefinition* projectile, 
                    G4double plab );
  // constructor used for low mass diffraction
  G4hhElastic( G4ParticleDefinition* target, G4ParticleDefinition* projectile);

  virtual ~G4hhElastic();

  virtual G4bool IsApplicable(const G4HadProjectile &/*aTrack*/, 
			      G4Nucleus & /*targetNucleus*/);


  void Initialise();

  void BuildTableT( G4ParticleDefinition* target, G4ParticleDefinition* projectile); // , G4double plab );
  void BuildTableTest( G4ParticleDefinition* target, G4ParticleDefinition* projectile, 
                    G4double plab );

  G4double SampleInvariantT( const G4ParticleDefinition* p,  G4double plab, G4int, G4int);
  G4double SampleBisectionalT( const G4ParticleDefinition* p, G4double plab);

  G4double SampleTest(G4double tMin ); // const G4ParticleDefinition* p, );

  G4double GetTransfer( G4int iMomentum, G4int iTransfer, G4double position );
 
private:


  G4ParticleDefinition* fTarget;
  G4ParticleDefinition* fProjectile;
  G4ParticleDefinition* theProton;
  G4ParticleDefinition* theNeutron;
  G4ParticleDefinition* thePionPlus;
  G4ParticleDefinition* thePionMinus;


  G4double lowEnergyRecoilLimit;  
  G4double lowEnergyLimitHE;  
  G4double lowEnergyLimitQ;  
  G4double lowestEnergyLimit;  
  G4double plabLowLimit;

  G4int fEnergyBin;
  G4int fBinT;

  G4PhysicsLogVector*           fEnergyVector;
  G4PhysicsTable*               fTableT;
  std::vector<G4PhysicsTable*>  fBankT;


  // Gauss model parameters


  G4double fMff2;
  G4double fMQ;
  G4double fMq;

  G4double fMassTarg;  // ~ A
  G4double fMassProj;  // ~ B
  G4double fMassSum2;
  G4double fMassDif2;


  G4double fRA; // hadron A
  G4double fRQ;
  G4double fRq;
  G4double fAlpha;
  G4double fBeta;

  G4double fRB; // hadron B
  G4double fRG;
  G4double fRg;
  G4double fGamma;
  G4double fDelta;

  G4double fAlphaP;
  G4double fLambdaFF;
  G4double fLambda;
  G4double fEta;
  G4double fImCof;
  G4double fCofF2;
  G4double fCofF3;
  G4double fRhoReIm;
  G4double fExpSlope;
  G4double fSo;

  G4double fSigmaTot;
  G4double fBq;
  G4double fBQ;
  G4double fBqQ;
  G4double fOptRatio;
  G4double fSpp;
  G4double fPcms;
  G4double fQcof;        // q prime when integrate

public:  // Gauss model methods

  void SetParameters();
  void SetSigmaTot(G4double stot){fSigmaTot = stot;};
  void SetSpp(G4double spp){fSpp = spp;};
  G4double GetSpp(){return fSpp;};
  void SetParametersCMS(G4double plab);

  G4double GetBq(){ return fBq;};
  G4double GetBQ(){ return fBQ;};
  G4double GetBqQ(){ return fBqQ;};
  void SetBq(G4double b){fBq = b;};
  void SetBQ(G4double b){fBQ = b;};
  void SetBqQ(G4double b){fBqQ = b;};
  G4double GetRhoReIm(){ return fRhoReIm;};

  void CalculateBQ(G4double b);
  void CalculateBqQ13(G4double b);
  void CalculateBqQ12(G4double b);
  void CalculateBqQ123(G4double b);
  void SetRA(G4double rn, G4double pq, G4double pQ);
  void SetRB(G4double rn, G4double pq, G4double pQ);

  void SetAlphaP(G4double a){fAlphaP = a;};
  void SetImCof(G4double a){fImCof = a;};
  G4double GetImCof(){return fImCof;};
  void SetLambda(G4double L){fLambda = L;};
  void SetEta(G4double E){fEta = E;};
  void SetCofF2(G4double f){fCofF2 = f;};
  void SetCofF3(G4double f){fCofF3 = f;};
  G4double GetCofF2(){return fCofF2;};
  G4double GetCofF3(){return fCofF3;};

  G4double GetRA(){ return fRA;};
  G4double GetRq(){ return fRq;};
  G4double GetRQ(){ return fRQ;};

  G4double GetRB(){ return fRB;};
  G4double GetRg(){ return fRg;};
  G4double GetRG(){ return fRG;};

  // FqQgG stuff

  G4complex Pomeron();

  G4complex Phi13();
  G4complex Phi14();
  G4complex Phi23();
  G4complex Phi24();

  G4complex GetF1qQgG(G4double qp);
  G4double GetdsdtF1qQgG(G4double s, G4double q);
  G4complex GetF2qQgG(G4double qp);
  G4double GetdsdtF12qQgG(G4double s, G4double q);
  G4complex GetF3qQgG(G4double qp);
  G4double GetdsdtF123qQgG(G4double q); // sampling ds/dt
  G4double GetdsdtF13qQG(G4double s, G4double q);


  // F123 stuff

  G4complex GetAqq();
  G4complex GetAQQ();
  G4complex GetAqQ();

  G4double GetCofS1();
  G4double GetCofS2();
  G4double GetCofS3();
  G4double GetOpticalRatio();

  G4complex GetF1(G4double qp);
  G4double GetdsdtF1(G4double s, G4double q);
  G4complex GetF2(G4double qp);
  G4double GetdsdtF12(G4double s, G4double q);
  G4complex GetF3(G4double qp);
  G4double GetdsdtF123(G4double q); // sampling ds/dt
  G4double GetExpRatioF123(G4double s, G4double q);

  // parameter arrays

private:

  G4int fInTkin;
  G4double fOldTkin;
  static const G4double theNuclNuclData[18][6];
  static const G4double thePiKaNuclData[8][6];
  G4HadronNucleonXsc* fHadrNuclXsc;
};

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////



inline G4bool G4hhElastic::IsApplicable(const G4HadProjectile & projectile, 
			      G4Nucleus & nucleus)
{
  if( ( projectile.GetDefinition() == G4Proton::Proton() ||
        projectile.GetDefinition() == G4Neutron::Neutron() ||
        projectile.GetDefinition() == G4PionPlus::PionPlus() ||
        projectile.GetDefinition() == G4PionMinus::PionMinus() ||
        projectile.GetDefinition() == G4KaonPlus::KaonPlus() ||
        projectile.GetDefinition() == G4KaonMinus::KaonMinus() ) &&

        nucleus.GetZ_asInt() < 2 ) return true;
  else                              return false;
}


inline void G4hhElastic::SetParameters()
{
  // masses

  fMq     = 0.36*CLHEP::GeV; // 0.441*GeV; // 0.36*GeV;
  fMQ     = 0.441*CLHEP::GeV;
  fMff2   = 0.26*CLHEP::GeV*CLHEP::GeV; // 0.25*GeV*GeV; // 0.5*GeV*GeV;

  fAlpha = 1./3.;
  fBeta  = 1. - fAlpha;

  fGamma = 1./2.; // 1./3.; // 
  fDelta = 1. - fGamma; // 1./2.;

  // radii and exp cof

  fRA       = 6.5/CLHEP::GeV; // 7.3/GeV;  //  3.25/GeV; // 7./GeV; // 2./GeV; // 1./GeV;
  fRq       = 0.173*fRA; // 2.4/GeV;
  fRQ       = 0.316*fRA; // 1./GeV; // 2./GeV; // 1./GeV;
  fRB       = 6.5/CLHEP::GeV; // 7.3/GeV;  //  3.25/GeV; // 7./GeV; // 2./GeV; // 1./GeV;
  fRg       = 0.173*fRA; // 2.4/GeV;
  fRG       = 0.173*fRA; // 1./GeV; // 2./GeV; // 1./GeV;

  fAlphaP   = 0.15/CLHEP::GeV/CLHEP::GeV; // 0.15/GeV/GeV;
  fLambda   = 0.25*fRA*fRA; // 0.25
  fEta      = 0.25*fRB*fRB; // 0.25
  fImCof    = 6.5;
  fCofF2 = 1.;
  fCofF3 = 1.;

  fBq = 0.02; // 0.21; // 1./3.;
  fBQ = 1. + fBq - 2*std::sqrt(fBq); // 1 - fBq; // 2./3.;
  fBqQ = std::sqrt(fBq*fBQ);

  fLambdaFF = 1.5/CLHEP::GeV/CLHEP::GeV; // 0.15/GeV/GeV;
  fSo     = 1.*CLHEP::GeV*CLHEP::GeV;
  fQcof   = 0.009*CLHEP::GeV;
  fExpSlope = 19.9/CLHEP::GeV/CLHEP::GeV;
}


////////////////////////////////////////////////////////////////////////
//
// Set target and projectile masses and calculate mass sum and difference squared for Pcms

inline void G4hhElastic::SetParametersCMS(G4double plab)
{
  G4int i;
  G4double trMass = 900.*CLHEP::MeV, Tkin;
  G4double sl, sh, ds, rAl, rAh, drA, rBl, rBh, drB, bql, bqh, dbq, bQl, bQh, dbQ, cIl, cIh, dcI;

  Tkin = std::sqrt(fMassProj*fMassProj + plab*plab) - fMassProj;

  G4DynamicParticle*  theDynamicParticle = new G4DynamicParticle(fProjectile,
                                              G4ParticleMomentum(0.,0.,1.),
                                              Tkin);
  fSigmaTot = fHadrNuclXsc->GetHadronNucleonXscNS( theDynamicParticle, fTarget );

  delete theDynamicParticle;

  fSpp  = fMassProj*fMassProj + fMassTarg*fMassTarg + 2.*fMassTarg*std::sqrt(plab*plab + fMassProj*fMassProj);
  fPcms = std::sqrt( (fSpp - fMassSum2)*(fSpp - fMassDif2)/4./fSpp);

  G4double sCMS = std::sqrt(fSpp);

  if( fMassProj > trMass ) // p,n,pb on p
  {
    this->SetCofF2(1.);
    this->SetCofF3(1.);
    fGamma = 1./3.; // 1./3.; // 
    fDelta = 1. - fGamma; // 1./2.;

    if( sCMS <= theNuclNuclData[0][0]*CLHEP::GeV ) // low edge, as s=2.76754
    {
      this->SetRA(theNuclNuclData[0][1]/CLHEP::GeV,0.173,0.316);
      this->SetRB(theNuclNuclData[0][2]/CLHEP::GeV,0.173,0.316); 

      this->SetBq(theNuclNuclData[0][3]);
      this->SetBQ(theNuclNuclData[0][4]);
      this->SetImCof(theNuclNuclData[0][5]);

      this->SetLambda(0.25*this->GetRA()*this->GetRA());
      this->SetEta(0.25*this->GetRB()*this->GetRB());
    }
    else if( sCMS >= theNuclNuclData[17][0]*CLHEP::GeV ) // high edge, as s=7000 ???
    {
      this->SetRA(theNuclNuclData[17][1]/CLHEP::GeV,0.173,0.316);
      this->SetRB(theNuclNuclData[17][2]/CLHEP::GeV,0.173,0.316); 

      this->SetBq(theNuclNuclData[17][3]);
      this->SetBQ(theNuclNuclData[17][4]);
      this->SetImCof(theNuclNuclData[17][5]);

      this->SetLambda(0.25*this->GetRA()*this->GetRA());
      this->SetEta(0.25*this->GetRB()*this->GetRB());
    }
    else // in approximation between array points
    {
      for( i = 0; i < 18; i++ ) if( sCMS <= theNuclNuclData[i][0]*CLHEP::GeV ) break;
      if( i == 0 ) i++;
      if( i == 18 ) i--;

      sl = theNuclNuclData[i-1][0]*CLHEP::GeV;
      sh = theNuclNuclData[i][0]*CLHEP::GeV;
      ds = (sCMS - sl)/(sh - sl);

      rAl = theNuclNuclData[i-1][1]/CLHEP::GeV;
      rAh = theNuclNuclData[i][1]/CLHEP::GeV;
      drA = rAh - rAl;

      rBl = theNuclNuclData[i-1][2]/CLHEP::GeV;
      rBh = theNuclNuclData[i][2]/CLHEP::GeV;
      drB = rBh - rBl;

      bql = theNuclNuclData[i-1][3];
      bqh = theNuclNuclData[i][3];
      dbq = bqh - bql;

      bQl = theNuclNuclData[i-1][4];
      bQh = theNuclNuclData[i][4];
      dbQ = bQh - bQl;

      cIl = theNuclNuclData[i-1][5];
      cIh = theNuclNuclData[i][5];
      dcI = cIh - cIl;

      this->SetRA(rAl+drA*ds,0.173,0.316);
      this->SetRB(rBl+drB*ds,0.173,0.316); 

      this->SetBq(bql+dbq*ds);
      this->SetBQ(bQl+dbQ*ds);
      this->SetImCof(cIl+dcI*ds);

      this->SetLambda(0.25*this->GetRA()*this->GetRA());
      this->SetEta(0.25*this->GetRB()*this->GetRB());
    }
  }
  else // pi, K
  {
    this->SetCofF2(1.);
    this->SetCofF3(-1.);
    fGamma = 1./2.; // 1./3.; // 
    fDelta = 1. - fGamma; // 1./2.;
  
    if( sCMS <= thePiKaNuclData[0][0]*CLHEP::GeV ) // low edge, as s=2.76754
    {
      this->SetRA(thePiKaNuclData[0][1]/CLHEP::GeV,0.173,0.316);
      this->SetRB(thePiKaNuclData[0][2]/CLHEP::GeV,0.173,0.173); 

      this->SetBq(thePiKaNuclData[0][3]);
      this->SetBQ(thePiKaNuclData[0][4]);
      this->SetImCof(thePiKaNuclData[0][5]);

      this->SetLambda(0.25*this->GetRA()*this->GetRA());
      this->SetEta(this->GetRB()*this->GetRB()/6.);
    }
    else if( sCMS >= thePiKaNuclData[7][0]*CLHEP::GeV ) // high edge, as s=7000 ???
    {
      this->SetRA(thePiKaNuclData[7][1]/CLHEP::GeV,0.173,0.316);
      this->SetRB(thePiKaNuclData[7][2]/CLHEP::GeV,0.173,0.173); 

      this->SetBq(thePiKaNuclData[7][3]);
      this->SetBQ(thePiKaNuclData[7][4]);
      this->SetImCof(thePiKaNuclData[7][5]);

      this->SetLambda(0.25*this->GetRA()*this->GetRA());
      this->SetEta(this->GetRB()*this->GetRB()/6.);
    }
    else // in approximation between array points
    {
      for( i = 0; i < 8; i++ ) if( sCMS <= thePiKaNuclData[i][0]*CLHEP::GeV ) break;
      if( i == 0 ) i++;
      if( i == 8 ) i--;

      sl = thePiKaNuclData[i-1][0]*CLHEP::GeV;
      sh = thePiKaNuclData[i][0]*CLHEP::GeV;
      ds = (sCMS - sl)/(sh - sl);

      rAl = thePiKaNuclData[i-1][1]/CLHEP::GeV;
      rAh = thePiKaNuclData[i][1]/CLHEP::GeV;
      drA = rAh - rAl;

      rBl = thePiKaNuclData[i-1][2]/CLHEP::GeV;
      rBh = thePiKaNuclData[i][2]/CLHEP::GeV;
      drB = rBh - rBl;

      bql = thePiKaNuclData[i-1][3];
      bqh = thePiKaNuclData[i][3];
      dbq = bqh - bql;

      bQl = thePiKaNuclData[i-1][4];
      bQh = thePiKaNuclData[i][4];
      dbQ = bQh - bQl;

      cIl = thePiKaNuclData[i-1][5];
      cIh = thePiKaNuclData[i][5];
      dcI = cIh - cIl;

      this->SetRA(rAl+drA*ds,0.173,0.316);
      this->SetRB(rBl+drB*ds,0.173,0.173); 

      this->SetBq(bql+dbq*ds);
      this->SetBQ(bQl+dbQ*ds);
      this->SetImCof(cIl+dcI*ds);

      this->SetLambda(0.25*this->GetRA()*this->GetRA());
      this->SetEta(this->GetRB()*this->GetRB()/6.);
    }
   }
  return;
}

/////////////////////////////////////////////////////
//
// RA for qQ

inline void G4hhElastic::SetRA(G4double rA, G4double pq, G4double pQ)
{
  fRA = rA;
  fRq = fRA*pq;
  fRQ = fRA*pQ;
}

/////////////////////////////////////////////////////
//
// RB for gG

inline void G4hhElastic::SetRB(G4double rB, G4double pg, G4double pG)
{
  fRB = rB;
  fRg = fRB*pg;
  fRG = fRB*pG;
}

////////////////////////////////////////////////////
//
// Returns Pomeron parametrization with Im part modified, *= fImCof

inline G4complex G4hhElastic::Pomeron()
{
  G4double re, im;

  re = fAlphaP*G4Log(fSpp/fSo);
  im = -0.5*fAlphaP*fImCof*CLHEP::pi;
  return G4complex(re,im);
} 

//////////////////////////////////////////////////

inline G4complex G4hhElastic::Phi13()
{
  G4double re = (fRq*fRq + fRg*fRg)/16.;
  G4complex result(re,0.);
  result += Pomeron();
  return result;
}

//////////////////////////////////////////////////

inline G4complex G4hhElastic::Phi14()
{
  G4double re = (fRq*fRq + fRG*fRG)/16.;
  G4complex result(re,0.);
  result += Pomeron();
  return result;
}

//////////////////////////////////////////////////

inline G4complex G4hhElastic::Phi23()
{
  G4double re = (fRQ*fRQ + fRg*fRg)/16.;
  G4complex result(re,0.);
  result += Pomeron();
  return result;
}

//////////////////////////////////////////////////

inline G4complex G4hhElastic::Phi24()
{
  G4double re = (fRQ*fRQ + fRG*fRG)/16.;
  G4complex result(re,0.);
  result += Pomeron();
  return result;
}

/////////////////////////////////////////////////////
//
// F1, case qQ-gG

inline G4complex G4hhElastic::GetF1qQgG(G4double t)
{
  G4double p = std::sqrt((fSpp - fMassSum2)*(fSpp - fMassDif2)/4./fSpp);
  G4double k = p/CLHEP::hbarc;

  G4complex exp13 = fBq*std::exp(-(Phi13() + fBeta*fBeta*fLambda + fDelta*fDelta*fEta)*t);
  G4complex exp14 = fBq*std::exp(-(Phi14() + fBeta*fBeta*fLambda + fGamma*fGamma*fEta)*t);
  G4complex exp23 = fBQ*std::exp(-(Phi23() + fAlpha*fAlpha*fLambda + fDelta*fDelta*fEta)*t);
  G4complex exp24 = fBQ*std::exp(-(Phi24() + fAlpha*fAlpha*fLambda + fGamma*fGamma*fEta)*t);

  G4complex res  = exp13 + exp14 + exp23 + exp24;

  res *=  0.25*k*fSigmaTot/CLHEP::pi;
  res *= G4complex(0.,1.);

  return res;
}

/////////////////////////////////////////////////////
//
// 

inline G4double G4hhElastic::GetdsdtF13qQG(G4double spp, G4double t)
{
  fSpp = spp;
  G4double p = std::sqrt((fSpp - fMassSum2)*(fSpp - fMassDif2)/4./fSpp);
  G4double k = p/CLHEP::hbarc;

  G4complex exp14 = fBqQ*std::exp(-(Phi14() + fBeta*fBeta*fLambda + fGamma*fGamma*fEta)*t);
  G4complex exp24 = fBQ*std::exp(-(Phi24() + fAlpha*fAlpha*fLambda + fGamma*fGamma*fEta)*t);

  G4complex F1  = exp14 + exp24;

  F1 *=  0.25*k*fSigmaTot/CLHEP::pi;
  F1 *= G4complex(0.,1.);

  // 1424

  G4complex z1424  = -(Phi24() + fAlpha*fLambda)*(Phi24() + fAlpha*fLambda);
            z1424 /= Phi14() + Phi24() + fLambda;
            z1424 += Phi24() + fAlpha*fAlpha*fLambda + fGamma*fGamma*fEta;

 
  G4complex exp1424  = std::exp(-z1424*t);
            exp1424 /= Phi14() + Phi24() + fLambda;

  G4complex F3  = fBqQ*fBQ*exp1424;


  F3 *=  0.25*k/CLHEP::pi;
  F3 *= G4complex(0.,1.);
  F3 *= fSigmaTot*fSigmaTot/(8.*CLHEP::pi*CLHEP::hbarc*CLHEP::hbarc);

  G4complex F13  = F1 - F3;

  G4double dsdt = CLHEP::pi/p/p;
  dsdt *= real(F13)*real(F13) + imag(F13)*imag(F13);  

  return dsdt;
}

//////////////////////////////////////////////////////////////
//
// dsigma/dt(s,t) F1qQgG

inline  G4double G4hhElastic::GetdsdtF1qQgG(G4double spp, G4double t)
{
  fSpp = spp;
  G4double p = std::sqrt((fSpp - fMassSum2)*(fSpp - fMassDif2)/4./fSpp);

  G4complex F1 = GetF1qQgG(t);

  G4double dsdt = CLHEP::pi/p/p;
  dsdt *= real(F1)*real(F1) + imag(F1)*imag(F1);  
  return dsdt;
}

/////////////////////////////////////////////////////
//
// 

inline G4complex G4hhElastic::GetF2qQgG(G4double t)
{
  G4double p = std::sqrt((fSpp - fMassSum2)*(fSpp - fMassDif2)/4./fSpp);
  G4double k = p/CLHEP::hbarc;

  G4complex z1324  = -(Phi24() + fAlpha*fLambda + fGamma*fEta)*(Phi24() + fAlpha*fLambda + fGamma*fEta);
            z1324 /= Phi13() + Phi24() + fLambda + fEta;
            z1324 += Phi24() + fAlpha*fAlpha*fLambda + fGamma*fGamma*fEta;

  G4complex exp1324  = std::exp(-z1324*t);
            exp1324 /= Phi13() + Phi24() + fLambda + fEta;

  G4complex z1423 = -(Phi23() + fAlpha*fLambda + fDelta*fEta)*(Phi24() + fAlpha*fLambda + fDelta*fEta);;
            z1423 /= Phi14() + Phi23() + fLambda + fEta;
            z1423 += Phi23() + fAlpha*fAlpha*fLambda + fDelta*fDelta*fEta;

  G4complex exp1423 = std::exp(-z1423*t);
            exp1423 /= Phi14() + Phi23() + fLambda + fEta; 

  G4complex res  = exp1324 + exp1423;


  res *=  0.25*k/CLHEP::pi;
  res *= G4complex(0.,1.);
  res *= fBq*fBQ*fSigmaTot*fSigmaTot/(8.*CLHEP::pi*CLHEP::hbarc*CLHEP::hbarc); // or 4. ???

  return res;
}

//////////////////////////////////////////////////////////////
//
// dsigma/dt(s,t) F12

inline  G4double G4hhElastic::GetdsdtF12qQgG( G4double spp, G4double t)
{
  fSpp = spp;
  G4double p = std::sqrt((fSpp - fMassSum2)*(fSpp - fMassDif2)/4./fSpp);

  G4complex F12 = GetF1qQgG(t) - GetF2qQgG(t);

  G4double dsdt = CLHEP::pi/p/p;
  dsdt *= real(F12)*real(F12) + imag(F12)*imag(F12);  
  return dsdt;
}

/////////////////////////////////////////////////////
//
// 

inline G4complex G4hhElastic::GetF3qQgG(G4double t)
{
  G4double p = std::sqrt( (fSpp - fMassSum2)*(fSpp - fMassDif2)/4./fSpp);
  G4double k = p/CLHEP::hbarc;

            // 1314

  G4complex z1314  = -(Phi14() + fGamma*fEta)*(Phi14() + fGamma*fEta);
            z1314 /= Phi13() + Phi14() + fEta;
            z1314 += Phi14() + fBeta*fBeta*fLambda + fGamma*fGamma*fEta;

  G4complex exp1314  = std::exp(-z1314*t);
            exp1314 /= Phi13() + Phi14() + fEta;

	    // 2324

  G4complex z2324 = -(Phi24() + fGamma*fEta)*(Phi24() + fGamma*fEta);;
            z2324 /= Phi24() + Phi23() + fEta;
            z2324 += Phi24() + fAlpha*fAlpha*fLambda + fGamma*fGamma*fEta;

  G4complex exp2324 = std::exp(-z2324*t);
            exp2324 /= Phi24() + Phi23() + fEta; 

	    // 1323

  G4complex z1323  = -(Phi23() + fAlpha*fLambda)*(Phi23() + fAlpha*fLambda);
            z1323 /= Phi13() + Phi23() + fLambda;
            z1323 += Phi23() + fAlpha*fAlpha*fLambda + fDelta*fDelta*fEta;

  G4complex exp1323  = std::exp(-z1323*t);
            exp1323 /= Phi13() + Phi23() + fLambda;

	    // 1424

  G4complex z1424  = -(Phi24() + fAlpha*fLambda)*(Phi24() + fAlpha*fLambda);
            z1424 /= Phi14() + Phi24() + fLambda;
            z1424 += Phi24() + fAlpha*fAlpha*fLambda + fGamma*fGamma*fEta;

  G4complex exp1424  = std::exp(-z1424*t);
            exp1424 /= Phi14() + Phi24() + fLambda;

	    G4complex res  = fBq*fBq*exp1314 + fBQ*fBQ*exp2324 + fBq*fBQ*exp1323 + fBq*fBQ*exp1424;

  res *=  0.25*k/CLHEP::pi;
  res *= G4complex(0.,1.);
  res *= fSigmaTot*fSigmaTot/(8.*CLHEP::pi*CLHEP::hbarc*CLHEP::hbarc);

  return res;
}

//////////////////////////////////////////////////////////////
//
// dsigma/dt(s,t) F123 sampling ds/dt

inline  G4double G4hhElastic::GetdsdtF123qQgG(G4double t)
{
  G4double p = std::sqrt( (fSpp - fMassSum2)*(fSpp - fMassDif2)/4./fSpp );

  G4complex F123  = GetF1qQgG(t); // - fCofF2*GetF2qQgG(t) - fCofF3*GetF3qQgG(t);
            F123 -= fCofF2*GetF2qQgG(t);
            F123 -= fCofF3*GetF3qQgG(t);

  G4double dsdt = CLHEP::pi/p/p;
  dsdt *= real(F123)*real(F123) + imag(F123)*imag(F123);  
  return dsdt;
}

/////////////////////////////////////////////////////
//
// Set fBqQ at a given fBQ=b2 according to the optical theorem,qQ-G

inline void G4hhElastic::CalculateBqQ13(G4double b2)
{
  fBQ = b2;

  G4complex z1424 = G4complex(1./8./CLHEP::pi,0.);
            z1424 /= Phi14() + Phi24() + fAlpha;
	    G4double c1424 = real(z1424)/(CLHEP::hbarc*CLHEP::hbarc);

  fBqQ  = 1. - fBQ;
  fBQ /= 1. - fSigmaTot*fBQ*c1424; 

  G4cout<<"fSigmaTot*fBQ*c1424 = "<<fSigmaTot*fBQ*c1424<<G4endl;

  G4double ratio  = fBqQ + fBQ - fSigmaTot*fBqQ*fBQ*c1424; 
  G4cout<<"ratio = "<<ratio<<G4endl;

  return ;
}

/////////////////////////////////////////////////////
//
// Set fBQ at a given fBq=b according to the optical theorem, F1-F2

inline void G4hhElastic::CalculateBqQ12(G4double b1)
{
  fBq = b1;

  G4complex z1324 = G4complex(1./8./CLHEP::pi,0.);
            z1324 /= Phi13() + Phi24() + fLambda + fEta;
	    G4double c1324 = real(z1324)/(CLHEP::hbarc*CLHEP::hbarc);

  G4complex z1423 = G4complex(1./8./CLHEP::pi,0.);
            z1423 /= Phi14() + Phi23() + fLambda + fEta;
	    G4double c1423 = real(z1423)/(CLHEP::hbarc*CLHEP::hbarc);

  fBQ  = 1. - 2.*fBq;
  fBQ /= 2. - fSigmaTot*fBq*(c1324+1423); 

  G4double ratio  = 2.*(fBq + fBQ) - fSigmaTot*fBq*fBQ*(c1324 + c1423); 
  G4cout<<"ratio = "<<ratio<<G4endl;

  return ;
}

/////////////////////////////////////////////////////
//
// Set fBQ at a given fBq=b according to the optical theorem, F1-F2-F3, 
// simplified meson-barion case g=G=q

inline void G4hhElastic::CalculateBqQ123(G4double b1)
{
  fBq = b1;

  G4complex z1324 = fCofF2*G4complex(1./8./CLHEP::pi,0.);
            z1324 /= Phi13() + Phi24() + fLambda + fEta;
	    G4double c1324 = real(z1324)/(CLHEP::hbarc*CLHEP::hbarc);

  G4complex z1423 = fCofF2*G4complex(1./8./CLHEP::pi,0.);
            z1423 /= Phi14() + Phi23() + fLambda + fEta;
	    G4double c1423 = real(z1423)/(CLHEP::hbarc*CLHEP::hbarc);

  G4complex z1314 = fCofF3*G4complex(1./8./CLHEP::pi,0.);
            z1314 /= Phi13() + Phi14() + fEta;
	    G4double c1314 = real(z1314)/(CLHEP::hbarc*CLHEP::hbarc);

  G4complex z2324 = fCofF3*G4complex(1./8./CLHEP::pi,0.);
            z2324 /= Phi23() + Phi24() + fEta;
	    G4double c2324 = real(z2324)/(CLHEP::hbarc*CLHEP::hbarc);

  G4complex z1323 = fCofF3*G4complex(1./8./CLHEP::pi,0.);
            z1323 /= Phi13() + Phi23() + fLambda;
	    G4double c1323 = real(z1323)/(CLHEP::hbarc*CLHEP::hbarc);

  G4complex z1424 = fCofF3*G4complex(1./8./CLHEP::pi,0.);
            z1424 /= Phi14() + Phi24() + fLambda;
	    G4double c1424 = real(z1424)/(CLHEP::hbarc*CLHEP::hbarc);

  G4double A = fSigmaTot*c2324;
  G4double B = fSigmaTot*fBq*(c1324 + c1423 + c1323 + c1424) - 2.;
  G4double C = 1. + fSigmaTot*fBq*fBq*c1314 - 2*fBq;
  G4cout<<"A = "<<A<<"; B = "<<B<<"; C = "<<C<<G4endl;
  G4cout<<"determinant = "<<B*B-4.*A*C<<G4endl;

  G4double x1 = ( -B - std::sqrt(B*B-4.*A*C) )/2./A;
  G4double x2 = ( -B + std::sqrt(B*B-4.*A*C) )/2./A;
  G4cout<<"x1 = "<<x1<<"; x2 = "<<x2<<G4endl;

  if( B*B-4.*A*C < 1.e-6 ) fBQ = std::abs(-B/2./A);
  else if ( B < 0.)        fBQ = std::abs( ( -B - std::sqrt(B*B-4.*A*C) )/2./A);
  else                     fBQ  = std::abs( ( -B + std::sqrt(B*B-4.*A*C) )/2./A);
 
  fOptRatio  = 2*(fBq+fBQ) - fSigmaTot*fBq*fBQ*(c1324 + c1423 + c1323 + c1424);
  fOptRatio -= fSigmaTot*fBq*fBq*c1314 + fSigmaTot*c2324*fBQ*fBQ;
  G4cout<<"BqQ123, fOptRatio = "<<fOptRatio<<G4endl;

  return ;
}



/////////////////////  F123 stuff hh-elastic, qQ-qQ ///////////////////
//////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

inline G4complex G4hhElastic::GetAqq()
{
  G4double re, im;

  re = fRq*fRq/8. + fAlphaP*G4Log(fSpp/fSo) + 8.*fLambda/9.;
  im = -0.5*fAlphaP*fImCof*CLHEP::pi;
  return G4complex(re,im);
}

/////////////////////////////////////////////////////
//
//

inline G4complex G4hhElastic::GetAQQ()
{
  G4double re, im;

  re = fRQ*fRQ/8. + fAlphaP*G4Log(fSpp/fSo) + 2.*fLambda/9.;
  im = -0.5*fAlphaP*fImCof*CLHEP::pi;
  return G4complex(re,im);
}

/////////////////////////////////////////////////////
//
//

inline G4complex G4hhElastic::GetAqQ()
{
  G4complex z = 0.5*( GetAqq() + GetAQQ() );
  return z;
}

/////////////////////////////////////////////////////
//
//

inline G4double G4hhElastic::GetCofS1()
{
  G4complex z = 1./( GetAqQ() + 4.*fLambda/9. );
  G4double result = real(z);
  result /= 4.*CLHEP::pi*CLHEP::hbarc*CLHEP::hbarc;
  result *= fSigmaTot*fCofF2;
  return result;
}

/////////////////////////////////////////////////////
//
//

inline G4double G4hhElastic::GetCofS2()
{
  G4complex z = 1./( GetAqq() + GetAqQ() - 4.*fLambda/9. );
  G4double result = real(z);
  result /= 4.*CLHEP::pi*CLHEP::hbarc*CLHEP::hbarc;
  result *= fSigmaTot*fCofF3;
  return result;
}

/////////////////////////////////////////////////////
//
//

inline G4double G4hhElastic::GetCofS3()
{
  G4complex z = 1./( GetAQQ() + GetAqQ() + 2.*fLambda/9. );
  G4double result = real(z);
  result /= 4.*CLHEP::pi*CLHEP::hbarc*CLHEP::hbarc;
  result *= fSigmaTot*fCofF3;
  return result;
}

/////////////////////////////////////////////////////
//
//

inline G4double G4hhElastic::GetOpticalRatio()
{
  return fOptRatio;
  // G4double sqrtBqBQ = std::sqrt(fBq*fBQ);
  // G4double result = fBq + fBQ + 2.*sqrtBqBQ - 1.;
  // result /= sqrtBqBQ*( GetCofS1()*sqrtBqBQ + GetCofS2()*fBq + GetCofS3()*fBQ );
  // return result;
}



/////////////////////////////////////////////////////
//
// Set fBQ at a given fBq=b according to the optical theorem

inline void G4hhElastic::CalculateBQ(G4double b1)
{
  fBq = b1;
  G4double s1 = GetCofS1();
  G4double s2 = GetCofS2();
  G4double s3 = GetCofS3();
  G4double sqrtBq = std::sqrt(fBq);

  // cofs of the fBQ 3rd equation

  G4double a = s3*sqrtBq;
  G4double b = s1*fBq - 1.;
  G4double c = (s2*fBq - 2.)*sqrtBq;
  G4double d = 1. - fBq;

  // cofs of the incomplete 3rd equation

  G4double p  = c/a;
           p -= b*b/a/a/3.;
  G4double q  = d/a;
           q -= b*c/a/a/3.;
           q += 2*b*b*b/a/a/a/27.;

  // cofs for the incomplete colutions

  G4double D  = p*p*p/3./3./3.;
           D += q*q/2./2.;
  G4complex A1 = G4complex(- q/2., std::sqrt(-D) );
  G4complex A  = std::pow(A1,1./3.);

  G4complex B1 = G4complex(- q/2., -std::sqrt(-D) );
  G4complex B  = std::pow(B1,1./3.);

  // roots of the incomplete 3rd equation

  G4complex y1 = A + B;
  G4complex y2 = -0.5*(A + B) + 0.5*std::sqrt(3.)*(A - B)*G4complex(0.,1.);
  G4complex y3 = -0.5*(A + B) - 0.5*std::sqrt(3.)*(A - B)*G4complex(0.,1.);
 
  G4complex x1 = y1 - b/a/3.;
  G4complex x2 = y2 - b/a/3.;
  G4complex x3 = y3 - b/a/3.;

  G4cout<<"re_x1 = "<<real(x1)<<"; re_x2 = "<<real(x2)<<"; re_x3 = "<<real(x3)<<G4endl;
  G4cout<<"im_x1 = "<<imag(x1)<<"; im_x2 = "<<imag(x2)<<"; im_x3 = "<<imag(x3)<<G4endl;

  G4double r1 = real(x1)*real(x1);
  G4double r2 = real(x2)*real(x2);
  G4double r3 = real(x3)*real(x3);

  if( r1 <= 1. && r1 >= 0. )      fBQ = r1;
  else if( r2 <= 1. && r2 >= 0. ) fBQ = r2;
  else if( r3 <= 1. && r3 >= 0. )      fBQ = r3;
  else fBQ = 1.;
  // fBQ = real(x3)*real(x3);
  G4double sqrtBqBQ = std::sqrt(fBq*fBQ);
  fOptRatio  = fBq + fBQ + 2.*sqrtBqBQ - 1.;
  fOptRatio /= sqrtBqBQ*( GetCofS1()*sqrtBqBQ + GetCofS2()*fBq + GetCofS3()*fBQ );
  G4cout<<"F123, fOptRatio = "<<fOptRatio<<G4endl;
  
  return ;
}

/////////////////////////////////////////////////////
//
//

inline G4complex G4hhElastic::GetF1(G4double t)
{
  G4double p = std::sqrt(0.25*fSpp - CLHEP::proton_mass_c2*CLHEP::proton_mass_c2);
  G4double  k    = p/CLHEP::hbarc;
  G4complex exp1 = fBq*std::exp(-GetAqq()*t);
  G4complex exp2 = fBQ*std::exp(-GetAQQ()*t);
  G4complex exp3 = 2.*std::sqrt(fBq*fBQ)*std::exp(-GetAqQ()*t);

  G4complex res  = exp1 + exp2 + exp3;
  res *=  0.25*k*fSigmaTot/CLHEP::pi;
  res *= G4complex(0.,1.);
  return res;
}

//////////////////////////////////////////////////////////////
//
// dsigma/dt(s,t) F1

inline  G4double G4hhElastic::GetdsdtF1(G4double spp, G4double t)
{
  fSpp = spp;
  G4double p = std::sqrt(0.25*spp - CLHEP::proton_mass_c2*CLHEP::proton_mass_c2);
  G4complex F1 = GetF1(t);

  G4double dsdt = CLHEP::pi/p/p;
  dsdt *= real(F1)*real(F1) + imag(F1)*imag(F1);  
  return dsdt;
}

/////////////////////////////////////////////////////
//
//

inline G4complex G4hhElastic::GetF2(G4double t)
{
  G4double p = std::sqrt(0.25*fSpp - CLHEP::proton_mass_c2*CLHEP::proton_mass_c2);
  G4double  k   = p/CLHEP::hbarc;
  G4complex z1  = GetAqq()*GetAQQ() - 16.*fLambda*fLambda/81.;
            z1 /= 2.*(GetAqQ() + 4.*fLambda/9.);
  G4complex exp1 = std::exp(-z1*t);

  G4complex z2 = 0.5*( GetAqQ() - 4.*fLambda/9.);

  G4complex exp2 = std::exp(-z2*t);

  G4complex res  = exp1 + exp2;

  G4complex z3 = GetAqQ() + 4.*fLambda/9.;

  res *=  0.25*k/CLHEP::pi;
  res *= G4complex(0.,1.);
  res /= z3;
  res *= fBq*fBQ*fSigmaTot*fSigmaTot/(8.*CLHEP::pi*CLHEP::hbarc*CLHEP::hbarc);

  return res;
}

//////////////////////////////////////////////////////////////
//
// dsigma/dt(s,t) F12

inline  G4double G4hhElastic::GetdsdtF12(G4double spp, G4double t)
{
  fSpp = spp;
  G4double p = std::sqrt(0.25*spp - CLHEP::proton_mass_c2*CLHEP::proton_mass_c2);
  G4complex F1 = GetF1(t) - GetF2(t);

  G4double dsdt = CLHEP::pi/p/p;
  dsdt *= real(F1)*real(F1) + imag(F1)*imag(F1);  
  return dsdt;
}

/////////////////////////////////////////////////////
//
//

inline G4complex G4hhElastic::GetF3(G4double t)
{
  G4double p = std::sqrt(0.25*fSpp - CLHEP::proton_mass_c2*CLHEP::proton_mass_c2);
  G4double  k   = p/CLHEP::hbarc;
  G4complex z1  = GetAqq()*GetAqQ() - 4.*fLambda*fLambda/81.;
            z1 /= GetAqq() + GetAqQ() - 4.*fLambda/9.;

  G4complex exp1 = std::exp(-z1*t)*fBq/(GetAqq() + GetAqQ() - 4.*fLambda/9.);

  G4complex z2 = GetAqQ()*GetAQQ() - 1.*fLambda*fLambda/81.;
            z2 /= GetAQQ() + GetAqQ() + 2.*fLambda/9.;

  G4complex exp2 = std::exp(-z2*t)*fBQ/(GetAQQ() + GetAqQ() + 2.*fLambda/9.);

  G4complex res  = exp1 + exp2;


  res *=  0.25*k/CLHEP::pi;
  res *= G4complex(0.,1.);
  res *= std::sqrt(fBq*fBQ)*fSigmaTot*fSigmaTot/(4.*CLHEP::pi*CLHEP::hbarc*CLHEP::hbarc);

  return res;
}

//////////////////////////////////////////////////////////////
//
// dsigma/dt(s,t) F123, sampling ds/dt

inline  G4double G4hhElastic::GetdsdtF123(G4double t)
{
  G4double  p   = std::sqrt(0.25*fSpp - CLHEP::proton_mass_c2*CLHEP::proton_mass_c2);
  G4complex F1  = GetF1(t);
            F1 -= fCofF2*GetF2(t);
            F1 -= fCofF3*GetF3(t);
  G4double dsdt = CLHEP::pi/p/p;
  dsdt         *= real(F1)*real(F1) + imag(F1)*imag(F1);  
  return dsdt;
}

//////////////////////////////////////////////////////////////
//
// dsigma/dt(s,t) F123

inline  G4double G4hhElastic::GetExpRatioF123(G4double spp, G4double t)
{
  fSpp = spp;
  G4double p = std::sqrt(0.25*spp - CLHEP::proton_mass_c2*CLHEP::proton_mass_c2);

  // qQ-ds/dt

  G4complex F1 = GetF1(t) - fCofF2*GetF2(t) - fCofF3*GetF3(t);

  G4double dsdt = CLHEP::pi/p/p;
  dsdt *= real(F1)*real(F1) + imag(F1)*imag(F1); 

  // exponent ds/dt

  G4complex F10 = GetF1(0.) - fCofF2*GetF2(0.) - fCofF3*GetF3(0.);

  fRhoReIm = real(F10)/imag(F10);

  G4double dsdt0 = CLHEP::pi/p/p;
  dsdt0 *= real(F10)*real(F10) + imag(F10)*imag(F10); 

  dsdt0 *= G4Exp(-fExpSlope*t);

  G4double ratio = dsdt/dsdt0;
 
  return ratio;
}


//
//
////////////////////////////////////////////////////////////////////////

#endif






