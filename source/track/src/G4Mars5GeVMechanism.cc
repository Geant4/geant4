// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Mars5GeVMechanism.cc,v 1.4 1999-12-15 14:53:55 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//	For information related to this code contact:
//	CERN, CN Division, ASD group
// 
// ------------------------------------------------------------
//   First Implemention    17 Nov. 1998  M.Asai, H.Kurahige
// 
// ------------------------------------------------------------
//  This is a Event Biasing mechanism based on MARS code
//   This model is applicable to 
//   proton/neutron/pi+-/K+-/gamma/anti_proton
//   with energy < 5.0GeV
//
//  Original code is MARS13 written by Nikolai Mokhov (FNAL)
//**************************************************************
//*   MARS13: 9. hA EVENT GENERATOR:
//*     Copyright Nikolai Mokhov (Fermilab)
//*
//*     LAST CHANGE: 14-NOV-1998
//**************************************************************
//*     Copyright Nikolai Mokhov (Fermilab)
//*
//*     MARS13(98)
//*
//*     INCLUSIVE HADRON(photon)-NUCLEUS VERTEX AT E < 5 GEV !!!
//*     THREE WEIGHTED HADRONS IN FINAL STATE:     !!!
//*     IP+A -> N/P(CASC)+ PI+/PI-(K+/K-) + PI0
//

#include "G4Mars5GeVMechanism.hh"
#include "globals.hh"
#include "G4ios.hh"
#include "Randomize.hh"
#include "G4VParticleChange.hh"
#include "G4Material.hh"
#include "G4Track.hh"
#include "G4Step.hh"
//-------------------------------------------------------

G4Mars5GeVMechanism::G4Mars5GeVMechanism(const G4String& name):
   G4VEvtBiasMechanism(name),
   EthForIncident(5.0*GeV)
{
   theParticleTable = G4ParticleTable::GetParticleTable();
   G4ParticleDefinition* pProton = theParticleTable->FindParticle("proton");
   if(pProton) ProtonMass = pProton->GetPDGMass();

   // set some constants 
   selec3.Eth = 1.0*MeV;
}  

G4Mars5GeVMechanism::G4Mars5GeVMechanism(const G4Mars5GeVMechanism& right):
   G4VEvtBiasMechanism(right),
   EthForIncident(right.EthForIncident)
{
  theParticleTable = G4ParticleTable::GetParticleTable();
  ProtonMass = right.ProtonMass;
  // set some constants 
  selec3.Eth = right.selec3.Eth;
}

G4Mars5GeVMechanism::~G4Mars5GeVMechanism()
{
}
 
void G4Mars5GeVMechanism::GetTargetNuclei(const G4Material* material)
{
  // get elements in the actual material,
  const G4ElementVector* theElementVector = material->GetElementVector();
  const G4double* theAtomicNumDensityVector = material->GetAtomicNumDensityVector();
  const G4int numberOfElements = material->GetNumberOfElements() ;
#ifdef G4VERBOSE
  if (GetVerboseLevel() > 2) {
    G4cout << " G4Mars5GeVMechanism::GetTargetNuclei" << G4endl;
 }
#endif 
  
  fANucl = 0.0;
  fZNucl = 0.0;
  G4double totNumAtoms = 0.0; 
  for (G4int iel=0; iel < numberOfElements; iel +=1) {
    totNumAtoms +=  theAtomicNumDensityVector[iel];
	fZNucl += theAtomicNumDensityVector[iel]*((*theElementVector)(iel)->GetZ());
	fANucl += theAtomicNumDensityVector[iel]*((*theElementVector)(iel)->GetN());
#ifdef G4VERBOSE
    if (GetVerboseLevel() > 2) {
      G4cout << iel << ": " << theAtomicNumDensityVector[iel];
      G4cout << "  Z=" << (*theElementVector)(iel)->GetZ() << "  A=" << (*theElementVector)(iel)->GetN();
      G4cout << G4endl; 
   }
#endif 
  }
  fANucl /= totNumAtoms;
  fZNucl /= totNumAtoms;
#ifdef G4VERBOSE
  if (GetVerboseLevel() > 2) {
    G4cout << "<Z>=" << fZNucl;
    G4cout << "<A>=" << fANucl;
    G4cout << G4endl; 
  }
#endif 
} 


void G4Mars5GeVMechanism::Treem5()
{
 
  G4double pMass = incidentParticle->GetDefinition()->GetPDGMass();
  G4double pE    = incidentParticle->GetKineticEnergy();
  G4int    pType = incidentMarsEncoding;

#ifdef G4VERBOSE
  if (GetVerboseLevel() > 2) {
    G4cout << " G4Mars5GeVMechanism::Treem5() ";
    G4cout << "       Incident Particle: " << incidentParticle->GetDefinition()->GetParticleName();
    G4cout << "  : energy = " <<   pE/GeV << "[GeV]" << G4endl; 
  }
#endif

  // CoulombBarrier
  if (CoulombBarrier(pType, pE)) return;

  G4int ib;
  if (pType==MarsAP) {
    ib = MarsP;
  } else if (pType==MarsGAM){
    if ( G4UniformRand() >0.5) {
      ib = MarsPIplus;
    } else {
      ib = MarsPIminus;
    }
  } else {
    ib = pType;
  }


  selec1.Einc = pE;
  if (pE < 0.5*MeV) pE =  0.5*MeV;
  selec3.Emax = pE;
  selec3.X  = 0.0;
  selec3.Pt = 0.0;
  selec3.P  = 0.0;

  // Nucleons at E < 5GeV
  CreateNucleon(ib, pType, pE);

  // Pion+- or Kaon+- at E < 5GeV
  CreatePion(ib, pType, pE);

  // Pi0 at E < 5GeV
  CreatePionZero(ib, pType, pE);
}

G4bool G4Mars5GeVMechanism::CoulombBarrier(G4int pType, G4double  pE){
  static const G4double EthCoulombBarrier = 20.0* MeV;
  static const G4double AvCoulomb = 1.11*MeV;
  static const G4double RCoulombTh = 1.0e-5;
  // CoulombBarrier
  if ( (  pType == MarsP) || (  pType ==MarsPIplus) || ( pType ==MarsKplus) ) { 
    if ( ( pE < EthCoulombBarrier ) && (fANucl >=1.5) ) {   
      G4double pMass =  GetParticleDefinition(pType)->GetPDGMass();
      G4double vCoulomb = AvCoulomb*pow(fZNucl/fANucl, 1./3.);
      G4double tc = pE*(fANucl*ProtonMass)/(pMass+(fANucl*ProtonMass));
      G4double rCoulomb = 1.0-vCoulomb/tc;
      if ( rCoulomb < RCoulombTh ) {
#ifdef G4VERBOSE
	if (GetVerboseLevel() > 2) {
	  G4cout << " Can not interact because of Coulomb Barrier " << G4endl; 
	}
#endif 
	return true;
      }
    }
  }
  return false; 
}

void G4Mars5GeVMechanism::CreateNucleon(G4int ib, G4int pType, G4double  pE)
{
#ifdef G4VERBOSE
  if (GetVerboseLevel() > 2) {
    G4cout << " G4Mars5GeVMechanism::CreateNucleon()" << G4endl;
  }
#endif
 if ( pType == MarsGAM) {
    selec1.Treac = MarsPIplus;
    selec1.Tprod = MarsN;
    selec1.V10   = 2.5;
  } else {
    if ( ib == MarsP ) {
      selec1.Treac = MarsP;
    } else if  ( ib == MarsN ) {
      selec1.Treac = MarsN;
    } else if  ( ib == MarsPIplus ) {
      selec1.Treac = MarsPIplus;
    } else if  ( ib == MarsPIminus ) {
      selec1.Treac = MarsPIminus;
    } else if  ( ib == MarsKplus ) {
      selec1.Treac = MarsPIplus;
    } else if  ( ib == MarsKminus ) {
      selec1.Treac = MarsPIminus;
    } else {
      selec1.Treac = MarsPIminus;
    }
    if (G4UniformRand()<0.5) {
      selec1.Tprod = MarsN;
    } else {
      selec1.Tprod = MarsP;
    }
   selec1.V10 = 2.0;
 }

  if ( SelBS(pType, fANucl, fZNucl) >0.0 ) AddSecondary();
  
}

void G4Mars5GeVMechanism::CreatePion(G4int ib, G4int pType, G4double  pE)
{
#ifdef G4VERBOSE
  if (GetVerboseLevel() > 2) {
    G4cout << " G4Mars5GeVMechanism::CreatePion()" << G4endl;
  }
#endif
 static const G4double PionProductionEth = 0.28*GeV; 
  static const G4double KaonProductionEth = 2.0*GeV; 

  if ( pE<PionProductionEth ) {
    if ((ib==MarsP)||(ib==MarsN)) return;
    pE +=  GetParticleDefinition(MarsPIminus)->GetPDGMass();
  }
  selec1.Einc = pE;		
	       
  if ( ib == MarsP ) {
    selec1.Treac = MarsP;
  } else if  ( ib == MarsN ) {
    selec1.Treac = MarsN;
  } else if  ( ib == MarsPIplus ) {
    selec1.Treac = MarsPIplus;
  } else if  ( ib == MarsPIminus ) {
    selec1.Treac = MarsPIminus;
  } else if  ( ib == MarsKplus ) {
    selec1.Treac = MarsPIplus;
  } else if  ( ib == MarsKminus ) {
    selec1.Treac = MarsPIminus;
  } else {
    selec1.Treac = MarsPIminus;
  }
  if (G4UniformRand()<0.5) {
    selec1.Tprod = MarsPIplus;
  } else {
    selec1.Tprod = MarsPIminus;
  }
  selec1.V10 = 2.1;
  if ( SelBS(pType, fANucl, fZNucl) >0.0 ){
    // change secondary into Kaon 
    if ( pE  > PionProductionEth ) {
      if ( Rkaon(ib,selec1.Tprod,pE) > G4UniformRand()) {
	if (selec1.Tprod==MarsPIminus) {
	  selec1.Tprod=MarsKminus;
	} else {
	  selec1.Tprod=MarsKplus;
	}
      }
    }
    AddSecondary();
  }
}

void G4Mars5GeVMechanism::CreatePionZero(G4int ib, G4int pType, G4double  pE)
{
#ifdef G4VERBOSE
  if (GetVerboseLevel() > 2) {
    G4cout << " G4Mars5GeVMechanism::CreatePionZero()" << G4endl;
  }
#endif
  static const G4double PionProductionEth = 0.28*GeV; 

  if ( pE<PionProductionEth ) {
    if ((ib==MarsP)||(ib==MarsN)) return;
  }
  			       
  if ( ib == MarsP ) {
    selec1.Treac = MarsP;
  } else if  ( ib == MarsN ) {
    selec1.Treac = MarsN;
  } else if  ( ib == MarsPIplus ) {
    selec1.Treac = MarsPIplus;
  } else if  ( ib == MarsPIminus ) {
    selec1.Treac = MarsPIminus;
  } else if  ( ib == MarsKplus ) {
    selec1.Treac = MarsPIplus;
  } else if  ( ib == MarsKminus ) {
    selec1.Treac = MarsPIminus;
  } else {
    selec1.Treac = MarsPIminus;
  }
  selec1.Tprod = MarsKplus;
  selec1.V10 = 1.0;
  if ( SelBS(pType, fANucl, fZNucl) >0.0 ) {
    selec1.Tprod = MarsPI0;
    AddSecondary();
  }
}

void G4Mars5GeVMechanism::AddSecondary()
{
#ifdef G4VERBOSE
  if (GetVerboseLevel() > 2) {
    G4cout << " G4Mars5GeVMechanism::AddSecondary()" << G4endl;
    G4cout << " Particle :" << selec1.Tprod;
    G4cout << ":" << GetParticleName(selec1.Tprod) <<G4endl;
    G4cout << " Energy   :" << selec1.EN <<G4endl;
    G4cout << "Weight    :" <<  selec1.V * incidentWeight << G4endl;
  }
#endif
  // determine direction cosine
  G4double g  = 1.0; 
  while (g>=1.0) {
    G4double g1 = G4UniformRand();
    G4double g2 = G4UniformRand();
    G4double gg = 2.0*g1 - 1.0;
    g = gg*gg + g2*g2;
    selec2.Ch = (gg*gg - g2*g2)/g;
    selec2.Sh = 2.0*gg*g2/g;
  }
  G4ThreeVector pin = incidentParticle->GetMomentumDirection();
  G4ThreeVector pout;
  Trans(&pin, &pout);

  if (numberOfSecondaries>=FastVectorSize) {
    G4Exception(" G4Mars5GeVMechanism::AddSecondary() too many secondaries");
  }

  // create seconday Dynamic Particle
  G4DynamicParticle* secondary = 
        new G4DynamicParticle(GetParticleDefinition(selec1.Tprod),
			      pout.unit(),
			      selec1.EN);
  // add secondary into list
  secondaries.SetElement(numberOfSecondaries, secondary);
  weightOfSecondaries[numberOfSecondaries] = selec1.V * incidentWeight;
  numberOfSecondaries +=1;
}

G4double G4Mars5GeVMechanism::SelBS(G4int pType, G4double aNucl, G4double zNucl)
{
  static const G4double Atau= 0.2;
  static const G4double Btau= 0.5*GeV;

  G4int nc = 0;
  G4int ip = selec1.Treac;  // reaction particle type  
  G4int jp = selec1.Tprod;  // procduction particle type
  G4int jj = pType;         // incident particle type 
  G4double e0 = selec1.Einc;
  G4double en;
  G4double v2 = 0.0;

#ifdef G4VERBOSE
  if (GetVerboseLevel() > 2) {
    G4cout << " G4Mars5GeVMechanism::SelBS" << G4endl;
    G4cout << "   pType = " << pType <<  "   e0    = " << e0 << G4endl;
    G4cout << "   aNucl = " << aNucl <<  "   zNucl = " << zNucl << G4endl;
    G4cout << "   Treac = " <<selec1.Treac;
    G4cout << "   Tprod = " <<selec1.Tprod << G4endl;
 }
#endif 

  while(1){

    G4double g1 =  G4UniformRand();
    G4double g2 =  G4UniformRand();
    
    // calculate energy
    G4double dw = 0.0;
    if (ip==jp) {
      G4double ea = e0 * 0.01;
      if (ea < selec3.Eth) {
	dw = selec3.Emax-selec3.Eth;
	en = selec3.Eth + g1*dw;
      } else {
	G4double cb = log(ea/selec3.Eth);
	G4double ca = cb + 99.0;
	if (g1<cb/ca) {
	  en = selec3.Eth*exp(g1*ca);
	  dw = en*ca;
	} else {
	  en = ea*(g1*ca + 1.0 - cb);
	  dw = ea*ca;
	}
      }
    } else {
      en = selec3.Eth*pow(selec3.Emax/selec3.Eth, g1);
      dw = en*log(selec3.Emax/selec3.Eth);
    }

    selec1.EN = en;

#ifdef G4VERBOSE
    if (GetVerboseLevel() > 2) {
      G4cout << "selec1.EN = " << en << G4endl;
    }
#endif 
    
    if (en<0.5*MeV) {
      selec1.V  = 0.0;
      return selec1.V;
    }
    
    // calculate direction cosine
    G4double tau = en/Atau/e0*(Btau+e0);
    G4double c5 = 1.0-exp(-pi*tau);
    G4double c4 = 1.0-g2*c5; 
    G4double t1 = -log(c4)/tau;
    G4double rcs = cos(t1);
    G4double rss = sqrt(1.0-rcs*rcs);
    G4double da  = 2.0*pi*rss*c5/(tau+c4);
    selec2.Cs = rcs;
    selec2.Ss = rss;
    
    // select particle type
    G4int ib = ip;
    if (ip == MarsP) {
      ib = MarsN;
    } else if (ip == MarsN) {
      ib = MarsP;
    }
    G4int jb = jp;
    if ( ( jj==MarsGAM ) && ((jp!=MarsP)||(jp!=MarsN)) ){
      jb = MarsKplus;
    } else if (jp == MarsP) {
      jb = MarsN;
    } else if (jp == MarsN) {
      jb = MarsP;
    }
    
    // calculate V
    nc +=1;
    v2 = dw*D2N2(jj, e0, en, t1, ib, jb, aNucl, zNucl)*da*(selec1.V10);
#ifdef G4VERBOSE
    if (GetVerboseLevel() > 2) {
      G4cout << " D2N2 = " <<  v2/(dw*da*(selec1.V10));
      G4cout << " v2 = " << v2 << G4endl;
    }
#endif 

    if (v2>0.0) break;

    if (nc >=3) {
      selec1.V  = 0.0;
#ifdef G4VERBOSE
      if (GetVerboseLevel() > 2) {
	G4cout << "exceed retry limit !!" << G4endl;
      }
#endif 
      return selec1.V;
    }
  }
  selec1.V  = v2;
  return v2;
}


G4double G4Mars5GeVMechanism::D2N2(G4int    pType,       G4double incidentE, 
	      G4double prodE,       G4double tin,
	      G4int    reacType,    G4int    proType,
              G4double ai,          G4double z)
{
  // Hadron inclusive yield at E0 < 5 GeV
  // All parametrizations are based on 
  // the energy unit of MeV
  //
  // Original code is written by Nikolai Mokhov (Fermilab)
  //C     Copyright Nikolai Mokhov (Fermilab)
  //C
  //C     MARS13(98)
  //C
  //C     HADRON INCLUSIVE YIELD AT E0 < 5 GEV
  //C-----
  //C     CREATED:     1979        BY B.SYCHEV
  //C     MODIFIED:    1979-1998   BY NVM
  //C     LAST CHANGE: 16-JUL-1998 BY NVM

  static const G4double o2pi = 1./twopi;
  static const G4double ospi = 1./sqrt(pi);

  static G4double abu = 0.0;
  static G4double alga = 0.0;
  static G4double a13 = 1.0;
  static G4double a23 = 1.0;
  static G4double a125 = 1.0;
  static G4double am25 = 0.0;
  static G4double sqa = 1.0;
  static G4double sqa1 = 0.0;
  static G4double bm = 2.0;
  static G4double sl;
  static G4double sa;

  // input of this method
  G4double e0 =  incidentE/MeV; // SHOULD BE GIVEN BY MEV !!!
  G4double e  =  prodE/MeV;     // SHOULD BE GIVEN BY MEV !!!
  G4double t  =  tin;
  G4int i = reacType;
  G4int j = proType;
  G4int jj = pType;

  // output of this method
  G4double d2n = 0.0;
  G4double dnde = 0.0;  // this value is not used anywhere else

  if(ai<1.0) return 0.0;

  G4double a = ai;
  if(abu!=a)
  {
    abu = a;
    if(a<=2.0)
    {
      alga = 0.0;
      a13 = 1.0;
      a23 = 1.0;
      a125 = 1.0;
      am25 = 0.0;
      sqa = 1.0;
      sqa1 = 0.0;
      bm = 2.0;
    }
    else
    {
      alga = log(a);
      a13 = pow(a,1./3.);
      a23 = a13*a13;
      a125 = pow(a,-1.25);
      am25 = pow(a-1.0,0.25);
      sqa = sqrt(a);
      sqa1 = sqrt(a-1.);
      bm = 1.0 + sqa;
    }
    sl = 0.72/pow(1.+alga,0.4);
    sa = 0.087*a23 + 4.15;
  }

  G4double bn;
  if(a<=2.0)
  {
    if(i*j<9 && t>halfpi) return 0.0;
    bn = 2.0;
    if(i>=3) bn = 1.0;
  }
  else
  { bn = bm*exp(-sa*pow(3.68/e0,sl)); }
    
  G4double emm = e0;
  G4double e1ge = 0.001*e0;
  G4double e2ge = e1ge*e1ge;
  G4double f21 = 0.04/(e2ge*e2ge);
  G4double f31 = 0.38*pow(e1ge,-0.65);
  G4double f22 = 0.25/e2ge;
  G4double f32 = am25*0.7/(e1ge+1.);
  G4double ei2 = 0.0;
  G4double x1 = f21 + f31;
  if(x1<60.) ei2 = 0.8*exp(-x1);
  G4double ei1 = ei2;
  if(a>2.0)
  {
    ei1 = 0.;
    x1 = f22 + f32;
    if(x1<60.) ei1 = exp(-x1);
  }

  G4double ew1 = ei2;
  G4double dnl = 0.0;
  G4double dnl1 = 0.0;
  G4double eli = 0.0;
  if(ew1>=1.e-19)
  {
    G4double dli = 35.0*ew1/(a+69.0);
    eli = 0.5*dli*e0;
    if(i==j)
    { dnl1 = dli*(2./3.)/e0; }
    else
    { dnl1 = dli*(1.-e/e0)/e0; }
    dnl = dli*(5./3.-e/e0)/e0;
  }

  G4double qel = 1.0 - ei1;
  if(a>2.0)
  {
    G4double e02 = pow(e0/350.,1.5);
    G4double ex2 = 1.0;
    if(e02<60.) ex2 = 1.0-exp(-e02);
    qel = 0.0;
    if(t<halfpi) qel = 1.17*ex2*exp(-0.08*sqa1)*(1.-ew1);
  }

  G4int in = i;       // save i
  G4double sw2 = (0.5+10.*e2ge/(2.+e1ge))*(4.+e0/470.);
  G4double sql = sw2/(2.+e0/940.)-1.0;
  G4double eql = 0.0;
  G4double dnq = 0.0;

  if(jj!=MarsGAM || j!=5)
  {
    if(qel>1.e-25)
    {
      eql = qel*e0*(sql+1.)/(sql+2.);
      G4double bp1x = -60./log(e/e0);
      if(sql<=bp1x)
      { 
        bp1x = pow(e/e0,sql);
        dnq = qel/e0*(sql+1.)*bp1x;
      }
    }
  }

  G4double bp1 = e0;
  if(e0<1.e9) bp1 = sqrt(e0*e0+1880.*e0);
  G4double pul = 1.e-3*bp1;
  bp1 = 3.*pow(pul,0.25) - 2.0;
  if(bp1<1.) bp1 = 1.;
  G4double bpi = 0.0;
  if(ei1>0.) bpi = bp1*exp(0.075*sqa1)*ei1;
  G4double ec = 0.0;
  if(a>2.0) ec = 10.5 - 0.02*a;
  G4double g = 0.1*alga + 0.2;
  G4double eog = pow(e0,g);
  G4double f1 = 1./3.*ec*a/(1.8*eog);
  x1 = 1.0;
  if(f1<60.) x1 = 1.0 - exp(-f1);
  G4double fm = 1.8*eog*x1;
  G4double ez = ec + fm;
  if(fm>=e0) ez = ec + e0;
  G4double d = 1.0;
  if(i>=3) d = 0.0;
  x1 = 1.0;
  if(a<=44.) x1 = exp(-exp(4.-a));
  G4double ez2 = 33.5*a125*x1*(ez-ec)*(1.-ez/(ec-e0));
  G4double epw = e0 - ez - (bn-d)*ec - 140.*bpi - ez2;
  G4double e2 = epw - eli - eql;

  G4double ak1 = 3.0;
  G4double ak20 = 5.e-4*(1.+a13)*e0;
  x1 = 1.0;
  if(ak20<60.) x1 = 1.0 - exp(-ak20);
  G4double ga = pow(e1ge,0.06)*ak1*x1;
  G4double egr = e0/(ga+1.);
  G4double d2 = 250.*(1.+2.5*e0*exp(-0.02*a)/(e0+1.e3))/sqa;
  G4double aea = e2/(e0*(1./(1.+ga)-d2*log(1.+egr/d2)/e0));
  aea *= 1./(1.+d2*(ga+1.)/(3.*e0*(d2/e0+1.75)));
  if(i<=2 && j>=3)
  {
    emm = e0 - 140.;
    if(j!=5 && a!=1. && (i+j)!=5) emm = e0 - 280;
  }
  G4double dn = 0.0;
  if(e<=emm) dn = aea*(e0/emm)*pow(1.-e/emm,ga)/(e+d2);
  if(i>=3) bpi += 1;
  dnde = dn + dnq + dnl;
  // In original code, check nupr. But in this code, nupr is aliways set to 0
  // if(nupr==1) return;

  G4double pna = bpi/bp1;
  G4double pns = bn+bpi;

  // Angular distribution
  G4double qe = 0.0;
  if((jj!=MarsGAM || j!=5)
   &&(t<halfpi)
   &&(i<3||i==j||j>4)
   &&(i>2||j<3)
   &&(a>2.0||i!=2||j!=1)
   &&(qel>=1.e-26))
  {
    G4double d1 = 25.*(1.+0.008*e0*t);
    G4double dp;
    if(i==j) 
    { dp = 0.8; }
    else
    { dp = 0.2; }
    if(a<=2.&&i==1&&j<=2) dp = 0.5;
    if(a<=2.&&i==2&&j==2) dp = 1.0;
    G4double eq = e0*sqr(cos(t))/(1.+e0*sqr(sin(t))/1880.) - 25.0;
    G4double exq = sqr((e-eq)/d1) + 0.5*sw2*t*t;
    if(exq<60.) qe = qel*sw2*exp(-exq)*dp*ospi/d1;
  }

  G4int iold = i;
  if(i==3) i = 2;
  if(i==4) i = 1;
  G4bool condA = a<2. && i==2;
  G4bool condB = a<2.;
  G4double az;
  G4double pn;
  G4double pr;
  if(!condB)
  {
    if(i==2)
    { az = (z+1.)/(a-z); }
    else
    { az = (a+1.-z)/z; }
    x1 = 0.5*e1ge;
    pn = az;
    if(x1<60.) pn *= 1. + exp(-x1);
    if(i==j)
    { pr = bn*pn/(1.+pn); }
    else
    { pr = bn/(1.+pn); }
  }
  if(condA || !condB)
  {
    if(i==1) az = z/(a-z);
    if(i==2) az = (a-z)/z;
    G4double bp = 1.0;
    G4double e0g = e1ge*e2ge;
    if(e0g<60.) bp -= 0.5*exp(-e0g);
    G4double ap = az * bp;
    bp = 6.*(1.+ap);
    if((i==1&&j==3)||(i==2&&j==4)) pr = pna*(bp1/3.-(2.+ap)/bp);
    if((i==2&&j==3)||(i==1&&j==4)) pr = pna*(bp1/3.+(3.-ap)/bp);
    if(j==5) pr = pna*(bp1/3.+(2.*ap-1.)/bp);
  }
  if(condB)
  {
    switch(i)
    {
    case 1:
      switch(j)
      {
      case 1:
      case 2:
        pr = bn/2.; break;
      case 3:
      case 4:
        pr = pna*(bp1/3.-1./6.); break;
      case 5:
        pr = pna*(bp1+1.)/3.; break;
      }
      break;
    case 2:
      switch(j)
      {
      case 1: 
        pr = 0.33*ew1/bn; break;
      case 2:
        pr = (1.-0.33*ew1)/bn; break;
      }
      break;
    }
  }

  G4double ek3 = 0.01*sqrt(e1ge)*(1.+alga/4.);
  G4double tay = 200.*e0/(e0+560.);
  G4double w = e/tay;
  G4double ek4 = 1.21*e0*w/(sqrt(1.+alga)*(e0+2000.));
  if(j>=3) ek4 = 0.3*w*(e0-1000.)/(e0+1000.);
  G4double wpic = w*pi;
  G4double w2 = w*w;
  G4double ex8 = 1.0;
  if(wpic<60.) ex8 /= 1.0 + exp(-wpic);
  G4double ek = (1.+w2)*(1.+5.2*ek4/(2.+w2))*ex8;
  G4double wtw = 2.*(sqrt(1.+ek3*e*t*1.e-3)-1.)/(tay*ek3*1.e-3)+ek4*t*t;
  G4double sm = 0.0;
  if(dn>=1.e-26 && wtw<60.) sm = pr*dn*ek*exp(-wtw)/pns;

  G4double dl = 0.0;
  i = iold;
  if((dnl1>=1.e-20)
   &&(i<3||i==j||j>4)
   &&(i>2||j<3))
  {
    tay = 200.*e0/(e0+2600.);
    w = e/tay;
    ek4 = 1.21*e0*w/(sqrt(1.+alga)*(e0+2000.));
    i = in;
    wpic = w*pi;
    w2 = w*w;
    ex8 = 1.0;
    if(wpic<60.) ex8 /= 1.0 + exp(-wpic);
    wtw = w*t + ek4*t*t;
    if(wtw<60.)
    {
      G4double ft = (1.+w2)*(1.+5.2*ek4/(2.+w2))*exp(-wtw)*ex8;
      if(ft>=1.-16)
      {
        G4double dp;
        if(i==j)
        { dp = 2./3.; }
        else
        { dp = 1./3.; }
        if(a<=2.&&i==1&&j==2) dp = 0.5;
        dl = dp*dnl1*ft;
      }
    }
  }

  if(jj==MarsGAM && j==5)
  {
    sm *= 0.6;
    dl *= 2.0;
  }
  d2n = o2pi*(qe+sm+dl);
  // d2n value is calculated in unit of [1/MeV] 
  
  d2n *= (1./MeV);
  return d2n;
}


G4double G4Mars5GeVMechanism::Rkaon(G4int ib, G4int jp, G4double eRaw)
{
  // Energy dependent K/pi ratio
  // Parametrizations are valid for energy range of
  // incident particle as 2.0 GeV to 100 GeV
  // All parametrizations in this method are based on
  // the energy unit of GeV.
  //
  // Original code is written by Nikolai Mokhov (Fermilab)
  //C     Copyright Nikolai Mokhov (Fermilab)
  //C
  //C     MARS13(98)
  //C     ENERGY DEPENDENT K/PI RATIO
  //C     FOR GIVEN TREEM AND SELMO PARAMETERS
  //C-----
  //C     CREATED:     1996        BY N.MOKHOV (NVM)
  //C     LAST CHANGE: 12-FEB-1996 BY NVM

  static const G4double rkp = 0.071;
  static const G4double rkm = 0.083;
  static const G4double al2 = 0.69314718;
  static const G4double al100 = 4.6051702;
  static const G4double al21 = 3.0445224;
  static const G4double al51 = 3.9318256;

  G4double eGeV = eRaw / GeV;
  G4double rK = 0.;
  if(eGeV < 2.1) return rK;
  G4double ale = log(eGeV);

  // No.1
  rK = rkp;
  if(jp == MarsPIminus) rK = rkm;
  if(ib == MarsPIplus || ib == MarsPIminus) rK *= 1.3;
  else if(ib == MarsKplus || ib == MarsKminus) rK *= 2.0;
  G4double rK1 = rK;
  if(eGeV<100.)
  { 
    G4double rmi = 0.03;
    if(ib >= MarsPIplus) rmi = 0.08;
    rK1 = rmi + (rK-rmi)*(ale-al2)/(al100-al2);
  }

  // No.2
  if(eGeV<=5.2 || eGeV>=51.0) { 
    rK = 1.3*rK1; 
  } else if(eGeV<7.2) { 
    rK = rK1*(1.3+0.15*(eGeV-5.2)); 
  } else if(eGeV<21.) {
    rK = 1.6*rK1; 
  } else { 
    rK = rK1*(1.3+0.3*(al51-ale)/(al51-al21)); 
  }

  return rK;
}

void G4Mars5GeVMechanism::Trans(G4ThreeVector* d1,G4ThreeVector* d2)
{
 
#ifdef G4VERBOSE
  if (GetVerboseLevel() > 2) {
    G4cout << " G4Mars5GeVMechanism::Trans() " << G4endl;
  }
#endif
 // Direction cosine transformation
  // using (cs,ss,ch,sh)
 
 // inputs
  G4double cs = selec2.Cs;
  G4double ss = selec2.Ss;
  G4double ch = selec2.Ch; 
  G4double sh = selec2.Sh;

  G4double sss, ttt, uuu;
  G4double dx1 = d1->x();
  G4double dy1 = d1->y();
  G4double dz1 = d1->z();
  G4double sz = dx1*dx1 + dy1*dy1;
  if(sz > 1.e-50)
  {
    sz = sqrt(sz);
    sss = ss*(ch*dz1*dx1-sh*dy1)/sz + cs*dx1;
    ttt = ss*(ch*dz1*dy1+sh*dx1)/sz + cs*dy1;
    uuu = - ss*ch*sz + cs*dz1;
  }
  else
  {
    sss = ss*ch + dx1;
    ttt = ss*sh + dy1;
    uuu = cs*dz1;
  }
  G4double den = sqrt(sss*sss+uuu*uuu+ttt*ttt);
  d2->setX(sss/den);
  d2->setY(ttt/den);
  d2->setZ(uuu/den);
  return;
}


G4VParticleChange* G4Mars5GeVMechanism::ApplyMath( 
                          G4VParticleChange* pVPChange, 
						  const G4Step& aStep )
{
  G4VParticleChange* pChange = (pVPChange);
#ifdef G4VERBOSE
  if (GetVerboseLevel() > 2) {
    G4cout << " G4Mars5GeVMechanism::ApplyMath" << G4endl;
  }
#endif
  G4Track* incidentTrack = aStep.GetTrack();
  incidentParticle = incidentTrack->GetDynamicParticle();
  incidentWeight = pChange->GetParentWeight();

  // check incident energy below this model is active 
  if (incidentParticle->GetKineticEnergy() > EthForIncident) return pChange;  
  if (incidentParticle->GetKineticEnergy() < 1.0*MeV) return pChange;  
  // check the incident particle type
  incidentMarsEncoding = GetMarsEncoding(incidentParticle->GetDefinition());
  if ( !IsApplicable(incidentMarsEncoding) ) return pChange; 
#ifdef G4VERBOSE
  if (GetVerboseLevel() > 2) {
    G4cout << " OK the particle is applicable" << G4endl;
  }
#endif

  // examine whether the current process is "hadronic interaction"
  //   by checking secondaries 
  G4int idx;
  G4bool flag = false;
  for (idx=0; (!flag) && (idx<pChange->GetNumberOfSecondaries()); idx+=1){
    G4String type = pChange->GetSecondary(idx)->GetDefinition()->GetParticleType();
    flag |= (type=="baryon")||(type=="nucleus")||(type=="meson");
  }
  if (!flag) return pChange;

  // Atomic and charge number 
  GetTargetNuclei( incidentTrack->GetMaterial() );
  
  // initialize secondary information 
  numberOfSecondaries = 0;
  secondaries.Initialize(FastVectorSize);

  // clean up ParticleChange
  for (idx=0; idx<pChange->GetNumberOfSecondaries(); idx+=1){
    delete pChange->GetSecondary(idx);
  } 
  pChange->Clear();
 
  // invoke MARS
  Treem5();

  // 
  pChange->SetNumberOfSecondaries(numberOfSecondaries);
  G4Track* track;
  for (idx=0; idx<numberOfSecondaries; idx+=1){
    track = new G4Track(
			       secondaries[idx],
			       incidentTrack->GetGlobalTime(),
			       incidentTrack->GetPosition()
                             );
    track->SetWeight( weightOfSecondaries[idx] );
	pChange->AddSecondary(track);
  }
  pChange->SetStatusChange(fStopAndKill);
  pChange->SetLocalEnergyDeposit (0.);
  return pChange;
}

