// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Quasmon.cc,v 1.2 1999-12-15 14:52:11 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4Quasmon ----------------
//             by Mikhail Kossov, July 1999.
//  class for an excited hadronic state used by the CHIPS Model
// ------------------------------------------------------------

//#define debug
#define pdebug

#include "G4Quasmon.hh"

G4Quasmon::G4Quasmon(G4int Z, G4int N, G4int S, G4LorentzVector FourMomentum) :
  qZ(Z),qN(N),qS(S),q4Mom(FourMomentum)
{
  if (Z>0)
  {
    valQ.IncU(Z+Z);
    valQ.IncD(Z);
  }
  else if (Z<0)
  {
    valQ.IncAU(Z+Z);
    valQ.IncAD(Z);
  }

  if (N>0)
  {
    int t=valQ.GetAU();
    if(t>=N) valQ.DecAU(N);
    else if(t<=0) valQ.IncU(N);
    else
    {
      valQ.SetAU();
      valQ.IncU(N-t);
    }

    int dN=N+N;
    t=valQ.GetAD();
    if(t>=dN) valQ.DecAD(dN);
    else if(t<=0) valQ.IncD(dN);
    else
    {
      valQ.SetAD();
      valQ.IncD(dN-t);
    }
  }
  else if (N<0)
  {
    int t=valQ.GetU();
    if(t>=N) valQ.DecU(N);
    else if(t<=0) valQ.IncAU(N);
    else
    {
      valQ.SetU();
      valQ.IncAU(N-t);
    }

    int dN=N+N;
    t=valQ.GetD();
    if(t>=dN) valQ.DecD(dN);
    else if(t<=0) valQ.IncAD(dN);
    else
    {
      valQ.SetD();
      valQ.IncAD(dN-t);
    }
  }

  if (S>0)
  {
    valQ.SetS(S);
 
    int t=valQ.GetAU();
    if(t>=S) valQ.DecAU(S);
    else if(t<=0) valQ.IncU(S);
    else
    {
      valQ.SetAU();
      valQ.IncU(S-t);
    }

    t=valQ.GetAD();
    if(t>=S) valQ.DecAD(S);
    else if(t<=0) valQ.IncD(S);
    else
    {
      valQ.SetAD();
      valQ.IncD(S-t);
    }
  }
  else if (S<0)
  {
    valQ.SetAS(S);
 
    int t=valQ.GetU();
    if(t>=S) valQ.DecU(S);
    else if(t<=0) valQ.IncAU(S);
    else
    {
      valQ.SetU();
      valQ.IncAU(S-t);
    }

    t=valQ.GetD();
    if(t>=S) valQ.DecD(S);
    else if(t<=0) valQ.IncAD(S);
    else
    {
      valQ.SetD();
      valQ.IncAD(S-t);
    }
  };
  //@@ Temporary debugging test ==========
  if (valQ.CheckNegative())
  {
    G4cerr << "G4Quasmon:NegativeQContent(" << valQ.GetU() << "," << valQ.GetD() << ","
         << valQ.GetS() << "," << valQ.GetAU() << "," << valQ.GetAD() << "," 
         << valQ.GetAS() << ")" << " for Z=" << Z << ",N=" << N << ",S=" << S << G4endl;
    //G4cerr<<"G4Quasmon:NegativeQContent(" << valQ <<"  for Z="<< Z <<",N="<< N <<",S="<< S << G4endl;
    abort();
  };
  InitCandidateVector();
}

G4Quasmon::G4Quasmon(G4int projPDG, G4int targPDG, G4LorentzVector proj4Mom,
                     G4LorentzVector targ4Mom)
{
  q4Mom = proj4Mom + targ4Mom;
  G4ParticleDefinition * projDefinition =
              G4ParticleTable::GetParticleTable()->FindParticle(projPDG);
  G4ParticleDefinition * targDefinition =
              G4ParticleTable::GetParticleTable()->FindParticle(targPDG);
  valQ.SetU (projDefinition->GetQuarkContent(2)     + targDefinition->GetQuarkContent(2));
  valQ.SetD (projDefinition->GetQuarkContent(1)     + targDefinition->GetQuarkContent(1));
  valQ.SetS (projDefinition->GetQuarkContent(3)     + targDefinition->GetQuarkContent(3));
  valQ.SetAU(projDefinition->GetAntiQuarkContent(2) + targDefinition->GetAntiQuarkContent(2));
  valQ.SetAD(projDefinition->GetAntiQuarkContent(1) + targDefinition->GetAntiQuarkContent(1));
  valQ.SetAS(projDefinition->GetAntiQuarkContent(3) + targDefinition->GetAntiQuarkContent(3));
#ifdef debug
  cout << "QInit: u=" << valQ.GetU() << ", d=" << valQ.GetD() << ", s=" << valQ.GetS()
       << ", au="<<valQ.GetAU() << ", ad="<<valQ.GetAD() << ", as="<<valQ.GetAS() << G4endl;
#endif
  InitCandidateVector();
}

G4Quasmon::G4Quasmon(const G4Quasmon &right) {}

G4Quasmon::~G4Quasmon()
{
  theQCandidates.clearAndDestroy();
  theQResonances.clearAndDestroy();
  //theQHadrons.clearAndDestroy();
}

// Additional Declarations

void G4Quasmon::InitCandidateVector()
{
  //@@ Temporary to find out the nature of "PDG" codes
  //for (G4int k=0; k<G4ParticleTable::GetParticleTable()->entries(); k++)
  //{
  //  G4String pName = G4ParticleTable::GetParticleTable()->GetParticleName(k);
  //  G4ParticleDefinition* pdf = G4ParticleTable::GetParticleTable()->GetParticle(k);
  //  G4int codePDG = pdf->GetPDGEncoding();
  //  cout << k << ". Name=" << pName << ", code=" << codePDG << G4endl;
  //}
  // Originaly the only particles which were in the Table are:
  //0. Name=proton, code=2212
  //1. Name=deuteron, code=0 (???)
  //2. Name=neutron, code=2112
  //3. Name=lambda, code=3122
  //4. Name=e-, code=11
  //5. Name=alpha, code=0 (???)
  //@@ It is necessary to define all mesons in the main!! (temporary)
  // Initialization if S-Meson Candidates starting a number of the Candidates
  const G4int nOfSMesons = 9;   // a#of S-wave mesons - candidates to output hadrons
  // @@Selfdescribing but not economic______________________________________________________
  //G4int cEta=221,cPiZ=111,cPiP=211,cPiM=-211,cKaZ=311,cKaP=321,cAKZ=-311,cKaM=-321,cEtP=331;
  //G4int sMezonPDG[nOfSMesons]={cEta,cPiZ,cPiP,cPiM,cKaZ,cKaP,cAKZ,cKaM,cEtP};
  //@@Alternative possibility:   Eta,Pi0,Pi+,APi-,Ka0,Ka+,AKa0,AKa-,Eta*
  G4int sMesonPDG[nOfSMesons] = {221,111,211,-211,311,321,-311,-321,331};
  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  for (int i=0; i<nOfSMesons; i++) 
  {
    //cout << "Meson # " << i << "before initialization" << G4endl;
    theQCandidates.insert(new G4QCandidate(sMesonPDG[i]));
#ifdef debug
    cout << "Meson # " << i << " with code = " << sMesonPDG[i]
         << ", QC=" << theQCandidates[i]->GetQC() << " is initialized" << G4endl;
#endif
  }

  const G4int nOfPResons = 9;  // a#of P-wave resonances for the final decay
  //@@ If the set is changed => change the initialization value of nOfPMesonResonances
  // @@Selfdescribing but not economic________________________
  //G4int cOmeg=223, cRhoZ=113, cRhoP=213, cRhoM=-213,
  //      cKZSt=313, cKPSt=323, cAKZS=-313,cKMSt=-323, cPhi=333;
  //G4int pMesonPDG[nOfPResons]={cOmeg,cRhoZ,cRhoP,cRhoM,cKZSt,cKPSt,cAKZS,cKMSt,cPhi};
  //@@Alternative possibility: Omega,Rh0,Rh+,Rho-,K0*,K+*,AK0*,AK-*,Phi
  G4int pMesonPDG[nOfPResons] = {223,113,213,-213,313,323,-313,-323,333};
  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  for (int j=0; j<nOfPResons; j++) 
  {
    //cout << "Resonance # " << j << "before initialization" << G4endl;
    theQResonances.insert(new G4QHadron(pMesonPDG[j]));
#ifdef debug
    cout << "Resonance # " << j << " with code = " << pMesonPDG[j] << ", QC="
         << theQResonances[j]->GetQC() << " is initialized" << G4endl;
#endif
  }
}

G4QHadronVector G4Quasmon::HadronizeQuasmon()
{
  static const G4double mPi0 = 134.9764;
  static const G4double mPi0SqDiv2 = mPi0*mPi0/2.;
  int j=0;
  G4int sPDG=0;                               // Fake definition of a Selected Candidate
  G4double MaxSqThresh = 0.;
  G4double selMass = 0.;
  while (q4Mom.m2()>MaxSqThresh)
  {
    G4double quaM2= q4Mom.m2();
#ifdef debug
    cout<<"QFragm: Yes. RQ4mom2="<<q4Mom.m2()<<" > Threshold2="<<MaxSqThresh<< G4endl;
#endif
    G4double quasM= q4Mom.m();                // Current mass of Quasmon 
    G4double kMin = mPi0SqDiv2/quasM;
#ifdef debug
    cout << "QFragm: kMax=" << quasM/2. << ", kMin=" << kMin << G4endl;
#endif
    CalculateContentOfQPartons(quasM);
#ifdef debug
    cout << "QFragm: Now call GetQPartonMomentum QC of Quasmon=" << valQ << G4endl;
#endif
    G4double kMom = GetQPartonMomentum(kMin); // Call it only after CalculateContentOfQPartons
#ifdef debug
    cout << "QFragm: =======>" << ++j << ": kMom=" << kMom << G4endl;
#endif
    CalculateHadronizationProbabilities(kMom);
    G4int nCandid = theQCandidates.entries();
#ifdef debug
    cout << "QFragm: NofCandidates=" << nCandid << G4endl;
#endif
    G4double totP = theQCandidates[nCandid-1]->GetIntegProbability() * G4UniformRand();
    int i=0;
    while(theQCandidates[i]->GetIntegProbability() < totP) i++;
#ifdef debug
    cout << "QFragm: i=" << i << " for p=" << totP << G4endl;
#endif
    //@@(?) now "i" points to the selected candidate
    if (i>=nCandid) G4cerr << "Candidate # " << i << " >= Total # " << nCandid << G4endl;
    G4double selM2 = theQCandidates[i]->GetMass2();
    sPDG  = theQCandidates[i]->GetPDGcode();
	//@@ Temporary
    selMass = theQCandidates[i]->GetMass();
    if (selMass+mPi0>quasM) G4cerr << "***QFrag:Too big candidate Mh=" << selMass
                                 << ", PDG=" << sPDG << ", Mq=" << quasM << G4endl;
    //
    //---Now the residual Quasmon mass cann't be smaller than a mass of corresponding S-meson 
    G4QContent curQ = valQ;                //Temporary copy of valQ to estimate the MinM2
    curQ   -= theQCandidates[i]->GetQC();  //Subtract outHadron from Quarc Content of Quasmon
#ifdef debug
    cout << "QFragm: -Q=" << theQCandidates[i]->GetQC()<< " = curQ=" << curQ  << G4endl;
#endif
    G4double minSqT = GetMinSqThresh(curQ);//Calculates mS+nP+gammaP Squared Threshold
    G4double dk = kMom+kMom;
    G4double kt = (quasM-dk)*(quasM-selM2/dk);
    G4double np = nOfQ - 3;                //Power for an integrated residual quasmon mass
#ifdef debug
    cout << "QFragm: kt=" << kt << ", minSqT=" << minSqT << ", np=" << np << G4endl;
#endif
    G4double cM = pow(minSqT/kt,np);       //Bottom cut for a possible Quasmon residual mass
    G4double uR = G4UniformRand();
    G4double rn = pow(cM +(1.-cM)*uR, 1./np);
#ifdef debug
    cout << "QFragm: cM=" << cM << ", uR=" << uR << ", rn=" << rn << G4endl;
#endif
    G4double m2 = kt*rn;                   //Squared Mass of residual quasmon
    //---Now check that the residual quasmon can decay in something (M_Q > mSMes1+mSMes2)
    // Residual S limit (absolute low limit for P-resonance) for S+S Decay
    //G4double minM2=GetMinSqThresh(curQ);
    //cout << "QFragm: If m2="<<m2<<"> minThresh2="<<minM2<<" to avoid S+S decay?"<< G4endl;
    //if (m2<=minM2)                         //In this case decay in ThisHadron+minSqTHadron
    // Residual S+S limit (absolute low limit for P-resonance) for S+S Decay
    G4double midM2=GetMidSqThresh(curQ);
#ifdef debug
    cout << "QFragm: Is m2="<<m2<< "> midThresh2="<<midM2<<" to avoid S+S decay?"<< G4endl;
#endif
    if (m2<=midM2)      //======>>>>>      In this case decay in ThisHadron+minSqTHadron
	{
      G4int flag=0;                              // Formal flag to get rPDG
      G4int rPDG=GetNofSqMassSum(valQ, flag);    // PDG of residual hadron (@@ temporary)
#ifdef debug
      cout << "QFragm:>>>NO final Q->hr+hs: rPDG=" << rPDG << ", sPDG=" << sPDG << G4endl;
#endif
      if(!Quasmon2HDecay(rPDG, sPDG)) G4cerr<<"***QFragm: Exception1"<<G4endl;
      return theQHadrons;                        // This is the last decay of the quasmon...
	}
    MaxSqThresh = GetMaxHSqThresh(valQ);         // Calculates mS+nP+GamP Squared Threshold
#ifdef debug
    cout<<"QFragm: YES. Is m2="<<m2<<" > maxTh2="<<MaxSqThresh<<" to stay in LOOP?"<< G4endl;
#endif
    if (m2<=MaxSqThresh) break;                  // In this case decay in Hadron+Resonance
    G4LorentzVector resQ4Mom(0.,0.,0.,sqrt(m2)); // 4-mom of residual Quasmon in CMS
    if(!Quasmon2HDecay(sPDG, resQ4Mom)) G4cerr << "***QFragm: Exception2" << G4endl;
    q4Mom = resQ4Mom;                            // Update the residual quasmon Lor.Vect.
    valQ    = curQ;                              // Update the Quark Content of Quasmon
#ifdef debug
    cout<<"QFragm: residual Quasmon's Lorents Vector="<<q4Mom << ", QC="<<valQ<< G4endl;
#endif
  }
  //---Now the {S-meson + P-meson} decay should take place (with Bright-Wigner probabilities)
#ifdef debug
  cout<<"QFragm:>>>NO S+P decay starts resM="<<q4Mom.m()<<",LorV="<<q4Mom<<",QC="<<valQ<< G4endl;
#endif
  G4int rPDG=GetNofSqMassSum(valQ, sPDG);  // PDG of residual hadron
#ifdef debug
  cout << "QFragm: P-resonance rPDG=" << rPDG << ", the hadron sPDG=" << sPDG << G4endl;
#endif
  G4double maxM = q4Mom.m()-selMass;
  //@@G4QResonance* curRes = new G4QResonance(rPDG, maxM); // The Resonance
  G4QHadron* curRes = new G4QHadron(rPDG, maxM); // The Resonance (as aHadron@@)
  G4LorentzVector tmp4Mom = curRes->Get4Momentum();
  if(!Quasmon2HDecay(sPDG, tmp4Mom)) G4cerr << "***QFragm: Exception3" << G4endl;
  //---Decay of residual P-resonance to the main channels listed in the PDG Decay Table
  q4Mom = tmp4Mom;                               // Update the residual quasmon = P-resonance
  //@@Later all the following part shold be in the G4QResonance class
  G4double rM = q4Mom.m();                       // Mass of decaying resonance
  G4int thePDG = curRes->GetPDGcode();           // Get PDG code of P-res. to switch
#ifdef debug
  cout<<"QFragm: decay Resonance with LV="<<q4Mom<<", m="<<rM<<", PDG="<<thePDG<< G4endl;
#endif
  G4double rnd = G4UniformRand();          // Random value to select brunching
  if     (thePDG== 223) // omega decays ===============================================
  {
    if(rnd<.891&&rM>414.12)
      {if(!Quasmon3HDecay(211,-211,111))G4cerr<<"***QFragm:1 Rez. decay"<<G4endl; ;}//3pi
    else if (rnd<.914&&rM>279.12)
      {if(!Quasmon2HDecay(211,-211))G4cerr<<"***QFragm:2 Rez. decay"<<G4endl; ;}//2pi
	else if(!Quasmon2HDecay(111,  22))G4cerr<<"***QFragm:3 Rez. decay"<<G4endl; ;//pi0,gam

  } 
  else if(thePDG== 113) // rho0 decay ============================================
         {if(!Quasmon2HDecay(211,-211))G4cerr<<"***QFragm:4 Rez. decay"<<G4endl; ;}//rho0=>pi+,pi-
  else if(thePDG== 213) // rho+ decay ============================================
         {if(!Quasmon2HDecay(211, 111))G4cerr<<"***QFragm:5 Rez. decay"<<G4endl; ;}//rho+=>pi+,pi0
  else if(thePDG==-213) // rho- decay ============================================
         {if(!Quasmon2HDecay(111,-211))G4cerr<<"***QFragm:6 Rez. decay"<<G4endl; ;}//rho-=>pi0,pi-
  else if(thePDG== 313) // K0* decays  ===========================================
  {
    if(rnd<.5&&rM>633.25)
      {if(!Quasmon2HDecay(321,-211))G4cerr<<"***QFragm:7 Rez. decay"<<G4endl; ;}//K+,pi-
    else if(!Quasmon2HDecay(311, 111))G4cerr<<"***QFragm:8 Rez. decay"<<G4endl; ;//K0,pi0
  } 
  else if (thePDG== 323) // K+* decays  ======================================
  {
    if(rnd<.5&&rM>637.25)
      {if(!Quasmon2HDecay(311,211))G4cerr<<"***QFragm:9 Rez. decay"<<G4endl; ;}//K0,pi+
    else if(!Quasmon2HDecay(321,111))G4cerr<<"***QFragm:A Rez. decay"<<G4endl; ;//K+,pi0
  }
  else if (thePDG==-313)// AK0* decays =========================================
  {
    if(rnd<.5&&rM>633.25)
      {if(!Quasmon2HDecay(-321,211))G4cerr<<"***QFragm:B Rez. decay"<<G4endl; ;}// K-,pi+
    else if(!Quasmon2HDecay(-311,111))G4cerr<<"***QFragm:C Rez. decay"<<G4endl; ;//AK0,pi0
  }
  else if (thePDG==-323)// K-* decays  ==========================================
  {
    if(rnd<.5&&rM>637.25)
      {if(!Quasmon2HDecay(-311,-211))G4cerr<<"***QFragm:D Rez. decay"<<G4endl; ;}//AK0,pi-
    else if(!Quasmon2HDecay(-321, 111))G4cerr<<"***QFragm:E Rez. decay"<<G4endl; ;// K-,pi0
  } 
  else if (thePDG==333) // phi decays  ===============================================
  {
	if(rnd<.343&&rM>995.36)
           {if(!Quasmon2HDecay(311,-311))G4cerr<<"***QFragm:F Rez. decay"<<G4endl; ;}//K0,AK0
	else if(rnd<.846&&rM>987.36)
           {if(!Quasmon2HDecay(321,-321))G4cerr<<"***QFragm:G Rez. decay"<<G4endl; ;}//K+,K-
    else if(rnd<.012&&rM>547.5)
           {if(!Quasmon2HDecay(221,  22))G4cerr<<"***QFragm:H Rez. decay"<<G4endl; ;}//eta,gam
    else if(!Quasmon3HDecay(211,-211,111))G4cerr<<"***QFragm:I Rez. decay"<<G4endl; ;//pi+,pi-,pi0
  }
  else G4cerr << "***QFragm: Unknown P-resonance PDGcode=" << thePDG << G4endl;
  return theQHadrons;
}

// Decay of current quasmon q4mom -> 2 particles: rPDG-hadron & sPDG-hadron
G4bool G4Quasmon::Quasmon2HDecay(const G4int& rPDG, const G4int& sPDG)
{
  G4QHadron* curHadr = new G4QHadron(rPDG);  // creation of the 1st Hadron
  G4LorentzVector tmp4Mom = curHadr->Get4Momentum();
  if(!Quasmon2HDecay(sPDG, tmp4Mom))
  {
    G4cerr<<"***Q22H(r/s): Exception"<<G4endl;
    return false;
  }
  curHadr->Set4Momentum(tmp4Mom);
  theQHadrons.insert(curHadr);               // Fill 1st Hadron to Output Hadrons Vector
  return true;
}

// Decay of current quasmon q4mom -> 2 particles: r4Mom & sPDG-hadron
G4bool G4Quasmon::Quasmon2HDecay(const G4int& sPDG, G4LorentzVector& r4Mom)
{
#ifdef debug
  cout << "Q22H is called with sPDJ="<<sPDG<<", r4Mom="<<r4Mom<< G4endl;
#endif
  G4QHadron* curHadr = new G4QHadron(sPDG);    // Creation of the 2nd Hadron
  G4double sM    = curHadr->GetMass();
  G4double sM2   = curHadr->GetMass2();
  G4double rM    = r4Mom.m();                 // Mass of (res.Quasmon/P-res./1st Hadron)
  G4double rM2   = r4Mom.m2();
  G4double quasM = q4Mom.m();                  // Mass of the current quasmon
  G4double quaM2 = q4Mom.m2();
#ifdef debug
  cout << "Q22H: sM="<<sM<<"(sPDG="<<sPDG<<"), rM="<<rM<<", qM="<<quasM<< G4endl;
#endif
  //@@ Later on make a quark content check for the decay
  if (quasM<rM+sM)
  {
#ifdef debug
    cout << "***Q22H*** sM="<<sM<<" + rM="<<rM<<" > M="<<quasM<< G4endl;
#endif
    delete curHadr; // The secondary hadron should be cleaned up if failed
    return false;
  }
  G4double d2 = quaM2-rM2-sM2;
  G4double p2 = (d2*d2/4.-rM2*sM2)/quaM2;      // Decay momentum(^2) in CMS of Quasmon
  if (p2<0.) G4cerr << "*?*p2="<<p2<< "d2*d2="<<d2*d2/4.<< "4m2selM2="<<rM2*sM2 << G4endl;
  G4double p  = sqrt(p2);
  G4double ct = 1.-2*G4UniformRand();
#ifdef debug
  cout << "Q22H: ct=" << ct << G4endl;
#endif
  G4double phi= 360.*deg*G4UniformRand();      // @@ Change 360.*deg to M_TWOPI (?)
  G4double ps = p * sqrt(1.-ct*ct);
  G4ThreeVector pVect(ps*sin(phi),ps*cos(phi),p*ct);
  G4LorentzVector had4Mom(pVect,sqrt(sM2+p2));
  r4Mom.setVect((-1)*pVect);
  r4Mom.setE(sqrt(rM2+p2));
  
  G4ThreeVector ltb = q4Mom.boostVector();// Boost vector for backward Lor.Trans.
#ifdef debug
  cout << "Q22H: Now making LorTrans with "<<ltb << ", 4Mom="<<q4Mom << G4endl;
#endif
  had4Mom.boost(ltb);                          // Lor.Trans. of hadron back to LS
  r4Mom.boost(ltb);                            // Lor.Trans. of residual Quazmon back to LS
#ifdef debug
  cout<<"Q22H: the hadron LorVector="<<had4Mom<<", the residual LorVector="<<r4Mom<<G4endl;
#endif
  curHadr->Set4Momentum(had4Mom);              // Define the LS 4-mom of the outgoing hadron
  theQHadrons.insert(curHadr);                 // Fill 2nd Hadron to Output Hadrons Vector
  return true;
}

// Decay of current quasmon q4mom (can be a P-res.) to 3 particles: rPDG, sPDG, & tPDG
G4bool G4Quasmon::Quasmon3HDecay(const G4int& rPDG, const G4int& sPDG, const G4int& tPDG)
{
#ifdef debug
  cout << "Q23H is called: rPDG="<<rPDG << ", sPDG="<<sPDG << ", tPDG="<<tPDG << G4endl;
#endif
  G4QHadron* curHadr1 = new G4QHadron(rPDG);   // 1-st hadron
  G4double rM  = curHadr1->GetMass();
  G4QHadron* curHadr2 = new G4QHadron(sPDG);   // 2-st hadron
  G4double sM  = curHadr2->GetMass();
  G4QHadron* curHadr3 = new G4QHadron(tPDG);   // 3-st hadron
  G4double tM  = curHadr3->GetMass();
  G4double qM  = q4Mom.m();                    // Mass of the current (decaying) quasmon
  //@@ Later on make a quark content check for the decay
  if (qM < rM + sM + tM)
  {
    G4cerr << "***Q23H: rM="<<rM<<" + sM="<<sM<<" + tM="<<tM<<" > qM="<<qM<< G4endl;
    delete curHadr1;                           // Should be cleaned up if failed
    delete curHadr2;                           // Should be cleaned up if failed
    delete curHadr3;                           // Should be cleaned up if failed
    return false;
  }
  G4double rM2 = rM*rM;
  G4double sM2 = sM*sM;
  G4double tM2 = tM*tM;
  G4double qM2 = qM*qM;
  G4double m13sBase=(qM-sM)*(qM-sM)-(rM+tM)*(rM+tM);
  G4double m12sMin =(rM+sM)*(rM+sM);
  G4double m12sBase=(qM-tM)*(qM-tM)-m12sMin;
  G4double rR = 0.;
  G4double rnd= 1.;
  G4int    tr = 0;                  // !Comment if "cout" belw is skiped !
  G4double m12s = 0.;               // Fake definition before the Loop
  while (rnd > rR)
  {
    m12s = m12sMin + m12sBase*G4UniformRand();
    G4double e1=m12s-rM2-sM2;
    G4double e2=qM2-m12s-tM2;
    G4double four12=4*m12s;
    G4double m13sRange=sqrt((e1*e1-four12*rM2)*(e2*e2-four12*tM2))/m12s;
    rR = m13sRange/m13sBase;
    rnd= G4UniformRand();
#ifdef debug
    cout << "Q23H: try to decay #"<<++tr << ", rR="<<rR << ", rnd="<<rnd << G4endl;
#endif
  }
  delete curHadr3;                  // Hadron3(tPDG) is created later in Q22H
  G4ThreeVector v0(0.,0.,0.);       // Zero Vector...
  G4double m12 = sqrt(m12s);        // Mass of H1+H2 system
  G4LorentzVector hadQ4Mom(v0,m12); // 4-mom of H1+H2 system in CMS
  if(!Quasmon2HDecay(tPDG, hadQ4Mom)) G4cerr << "***Q23H: Exception1" << G4endl;
  q4Mom = hadQ4Mom;                 // Update the residual quasmon = H1+H2 system
#ifdef debug
  cout << "Q23H: Now the last decay for mSS=" << q4Mom.m() << G4endl;
#endif
  delete curHadr2;                  // Hadron2(sPDG) is created later in Q22H
  hadQ4Mom.setVect(v0);             // Redefine hadQ4Mom for Hadron1(rPDG)
  hadQ4Mom.setE(rM);
  if(!Quasmon2HDecay(sPDG, hadQ4Mom)) G4cerr << "***Q23H: Exception2" << G4endl;
  curHadr1->Set4Momentum(hadQ4Mom); // Define the LS 4-mom of the outgoing hadron
  theQHadrons.insert(curHadr1);     // Fill the Output Hadrons Vector
  return true;
}

// mS
G4double G4Quasmon::GetMinSqThresh(const G4QContent& qCon) //@@ Temporary not used
{
  static const G4double sqm[5] = {0., 19479.771, 243716.98, 247677.42, 18218.629};

  G4int ind = (qCon.GetCharge()+1)*3 + qCon.GetStrangeness() + 1;
  if (ind<1 || ind>7) G4cerr << "G4Quasmon::GetMinSqThresh: ind=" << ind << G4endl;
  if (ind>5) ind = 8 - ind;

  return sqm[ind];
}

// mS+mS
G4double G4Quasmon::GetMidSqThresh(const G4QContent& qCon)
{
  // Squared MassSum (m1+m2)^2: [0] pi0+pi0,[1]pi0+pi,[2]pi0+K, pi0+K0,  pi0+eta, [5]pi+pi,
  static const G4double mm[15] = {72874.514,75375.698,395205.1,400244.,465501.09,77919.084,
  401001.7,406077.3,471790.33,974867.92,982772.84,1083633.1,990709.68,1091966.5,1198149.2};
  //  pi+K,[7]pi+K0,[8]pi+eta,  [9] K+K, [10]K+K0,[11]K+eta,[12]K0+K0,   K0+eta,  eta+eta

  G4int flag=10;
  return mm[GetNofSqMassSum(qCon,flag)];
}

// mP-res - GammaP-res
G4double G4Quasmon::GetMaxLSqThresh(const G4QContent& qCon) //@@ Not used(?!) Get rid of L/H?
{
  // Squared Sum (m1+m2-Gam)^2:[0]pi0+rho,[1]pi0+rho,[2]rho+K,[3]pi0+K0*,[4]pi0+phi,[5]pi+rho,
  static const G4double mm[15]={568480.41, 568480.41,1238050.1,961530.08, 1322406.6,575428.37,
   961242.89,970559.83,1332992.5,1780989.,1791667.9,1926988.2, 1804379.7, 1940170.4, 2440728.2};
  //[6]pi+K*,[7]pi+K0*,[8]pi+phi,[9]K+K*,[10]K*+K0,[11]K*+eta,[12]K0+K0*,[13]K0*+eta,[14]eta+phi

  G4int flag=10;
  return mm[GetNofSqMassSum(qCon,flag)];
}

// mP-res + GammaP-res
G4double G4Quasmon::GetMaxHSqThresh(const G4QContent& qCon)
{
  // Squared Sum (m1+m2+Gam)^2:[0]pi0+rho,[1]pi+rho,[2]rho+K,[3]rho+K0,[4]rho+eta,[5]pi+rho,
  static const G4double mm[15]={1115086.2,1124808.6,2001311.,2012630.2, 2155904.9,1124808.6,
    2001311.,2012630.2,2155904.9,2062489.5,2074397.8,2302867., 2085921.6, 2315007.9,2468490.3};
  //[6]rho+K,[7]rho+K0,[8]rho+eta,[9]K+K*,[10]K+K0*,[11]K+phi,[12]K0+K0*,[13]K0+phi,[14]eta+phi

  G4int flag=10;
  return mm[GetNofSqMassSum(qCon,flag)];
}

// Gives a number (0-14) corresponding to the decay pair (H+H or H+R)
G4int G4Quasmon::GetNofSqMassSum(const G4QContent& qCon, G4int& sPDG)
{
  G4QContent cQ = qCon;   // To modify if it's not enough or too many valence quarks
  if(sPDG==0 || sPDG==10) // Fake magic numbers for S+S & thresholds calculations
  {
    G4int     tot = cQ.GetTot();
    if (tot<4)      cQ.IncQAQ();
    else if (tot>4)      // Too many valent quarks to estimate the MidThreshold
    {
      G4int dec = (tot-4)/2;
#ifdef debug
      cout << "NofSMS: tot=" << tot << ", dec=" << dec << G4endl;
#endif
      G4int r   = cQ.DecQAQ(dec);
      if (r<0) G4cerr<<"***NofSMS:4 fail to decrement r="<<r<<", sPDG="<<sPDG<< G4endl;
    }
  }
  // Now we have 1 diquark and 1 anti-diquark
  //              0 (+4/3)     1 (-2/3)     2 (-2/3)     3 (+1/3)     4 (+1/3)      5 (-2/3)
  G4int dq[6] = {cQ.GetUU(),  cQ.GetDD(),  cQ.GetSS(),  cQ.GetUD(),  cQ.GetUS(),   cQ.GetDS()};
  G4int ad[6] = {cQ.GetAUAU(),cQ.GetADAD(),cQ.GetASAS(),cQ.GetAUAD(),cQ.GetAUAS(), cQ.GetADAS()};
  G4int dn=0;
  G4int an=0;
  for (int i=1; i<6; i++)
  {
    if (dq[i]) dn=i;
    if (ad[i]) an=i;
  }
  if (an<dn) // swap to make an>=dn
  {
    G4int r = an;
    an      = dn;
    dn      =  r;
  }

  if (sPDG==10) // fake sPDG just to get the threshold indices
  {
    if     (dn==an && an==2) return 14; // eta+eta // eta+phi-G // eta+phi+G
    else if (dn==an && an<4) return  0; // pi0+pi0 // pi0+rho-G // pi0+rho+G
    else if (dn==0 && an==1) return  5; // pi +pi  // pi +rho-G // pi +rho+G
    else if (dn<2  && an==3) return  1; // pi0+pi  // pi +rho-G // pi +rho+G
    else if (dn==0 && an==4) return  2; // pi0+K   // rho+K  -G // rho+K  +G
    else if (dn==3 && an==5) return  2; // pi0+K   // rho+K  -G // rho+K  +G
    else if (dn==1 && an==5) return  3; // pi0+K0  // pi0+K0*-G // rho+K0 +G
    else if (dn==3 && an==4) return  3; // pi0+K0  // pi0+K0*-G // rho+K0 +G
    else if (dn==0 && an==5) return  6; // pi +K   // pi +K* -G // rho+K  +G
    else if (dn==1 && an==4) return  7; // pi +K0  // pi +K0*-G // rho+K0 +G
    else if (dn==an)         return  4; // pi0+eta // pi0+phi-G // rho+eta+G
    else if (dn==4 && an==5) return  8; // pi +eta // pi +phi-G // rho+eta+G
    else if (dn==0 && an==2) return  9; // K  +K   // K  +K* -G // K  +K* +G
    else if (dn==2 && an==3) return 10; // K  +K0  // K0 +K* -G // K  +K0*+G
    else if (dn==1 && an==2) return 12; // K0 +K0  // K0 +K0*-G // K0 +K0*+G
    else if (dn==2 && an==4) return 11; // K  +eta // eta+K* -G // K  +phi+G
    else if (dn==2 && an==5) return 13; // K0 +eta // eta+K0*-G // K0 +phi+G
    else G4cerr << "Impossible combination dn=" << dn << ", an=" << G4endl;
  }
  else if (sPDG==0) // fake sPDG to get PDGs for S+S decay
  {
    //@@ Should be modified when compared with data...
    G4int charge  = cQ.GetCharge();
    G4int strange = cQ.GetStrangeness();
    if     (dn==an && an==2) // eta+eta
	{
      sPDG = 221;
      return 221;
    }
    else if (dn==an && an<4) // pi0+pi0
	{
      sPDG = 111;
      return 111;
    }
    else if (dn==0 && an==1) // pi +pi
	{
      if (charge>0)
	  {
        sPDG = 211;
        return 211;
	  }
      else
	  {
        sPDG = -211;
        return -211;
	  }
    }
    else if (dn<2  && an==3) // pi0+pi
	{
      if (charge>0)
	  {
        sPDG = 111;
        return 211;
	  }
      else
	  {
        sPDG = 111;
        return -211;
	  }
    }
    else if ((dn==0 && an==4) ||(dn==3 && an==5))  // pi0+K
	{
      if (charge>0 && strange<0)
	  {
        sPDG = 111;
        return 321;
	  }
      else if (charge<0 && strange>0)
	  {
        sPDG = 111;
        return -321;
	  }
      else
	  {
        G4cerr << cQ << " in pi0+K decay" << G4endl;
	  }
    }
    else if ((dn==1 && an==5) || (dn==3 && an==4)) // pi0+K0
	{
      if (strange>0)
	  {
        sPDG = 111;
        return 311;
	  }
      else
	  {
        sPDG = 111;
        return -311;
	  }
    }
    else if (dn==0 && an==5) // pi +K
	{
      if (charge==0 && strange<0)
	  {
        sPDG = -211;
        return 321;
	  }
      else if (charge==0 && strange>0)
	  {
        sPDG = 211;
        return -321;
	  }
      else if (charge==2 && strange<0)
	  {
        sPDG = 211;
        return 321;
	  }
      else if (charge==-2 && strange>0)
	  {
        sPDG = -211;
        return -321;
	  }
      else
	  {
        G4cerr << cQ << " in pi+K decay" << G4endl;
	  }
    }
    else if (dn==1 && an==4) // pi +K0
	{
      if (charge==1 && strange<0)
	  {
        sPDG = 211;
        return 311;
	  }
      else if (charge==1 && strange>0)
	  {
        sPDG = 211;
        return -311;
	  }
      else if (charge==-1 && strange<0)
	  {
        sPDG = -211;
        return 311;
	  }
      else if (charge==-1 && strange>0)
	  {
        sPDG = -211;
        return -311;
	  }
      else
	  {
        G4cerr << cQ << " in pi+K decay" << G4endl;
	  }
    }
    else if (dn==an)         // pi0+eta
	{
      sPDG = 111;
      return 221;
    }
    else if (dn==4 && an==5) // pi +eta
	{
      if (charge>0)
	  {
        sPDG = 221;
        return 211;
	  }
      else
	  {
        sPDG = 221;
        return -211;
	  }
    }
    else if (dn==0 && an==2) // K  +K
	{
      if (charge==2 && strange==-2)
	  {
        sPDG = 321;
        return 321;
	  }
      else if (charge==-2 && strange==2)
	  {
        sPDG = -321;
        return -321;
	  }
      //else if (charge==0 && strange==0) // It can be in eta+pi0 channel which is minimal
	  //{
      //  sPDG = 321;
      //  return -321;
	  //}
      else
	  {
        G4cerr << cQ << " in pi0+K decay" << G4endl;
	  }
    }
    else if (dn==2 && an==3) // K  +K0
	{
      if (charge>0)
	  {
        sPDG = -321;
        return 311;
	  }
      else
	  {
        sPDG = 321;
        return -311;
	  }
    }
    else if (dn==1 && an==2) // K0 +K0
	{
      sPDG = 311;
      return -311;
    }
    else if (dn==2 && an==4) // K  +eta
	{
      if (charge>0 && strange<0)
	  {
        sPDG = 221;
        return 321;
	  }
      else if (charge<0 && strange>0)
	  {
        sPDG = 221;
        return -321;
	  }
      else
	  {
        G4cerr << cQ << " in eta+K decay" << G4endl;
	  }
    }
    else if (dn==2 && an==5) // K0 +eta
	{
      if (strange>0)
	  {
        sPDG = 221;
        return 311;
	  }
      else
	  {
        sPDG = 221;
        return -311;
	  }
    }
    else G4cerr << "Impossible combination dn=" << dn << ", an=" << G4endl;
  }
  else  // For real PDG of outgoing hadron we should find the corresponding P-resonance
  {
    G4ParticleDefinition* def=G4ParticleTable::GetParticleTable()->FindParticle(sPDG);
    G4int nU = cQ.GetU()  - def->GetQuarkContent(2);
    G4int nD = cQ.GetD()  - def->GetQuarkContent(1);
    G4int nS = cQ.GetS()  - def->GetQuarkContent(3);
    G4int nAU= cQ.GetAU() - def->GetAntiQuarkContent(2);
    G4int nAD= cQ.GetAD() - def->GetAntiQuarkContent(1);
    G4int nAS= cQ.GetAS() - def->GetAntiQuarkContent(3);
    cQ.SetU(nU);
    cQ.SetD(nD);
    cQ.SetS(nS);
    cQ.SetAU(nAU);
    cQ.SetAD(nAD);
    cQ.SetAS(nAS);
    G4int     tot = cQ.GetTot();
    if (tot<2)      cQ.IncQAQ();
    else if (tot>2)      // Too many valent quarks to estimate the MidThreshold
    {
      G4int dc = (tot-2)/2;
#ifdef debug
      cout << "NofSMS: SP case: tot=" << tot << ", dc=" << dc << G4endl;
#endif
      G4int r  = cQ.DecQAQ(dc);
      if (r<0) G4cerr<<"***NofSMS:2 fail to decrement r="<<r<<", sPDG="<<sPDG<<G4endl;
      nU = cQ.GetU();
      nD = cQ.GetD();
      nS = cQ.GetS();
      nAU= cQ.GetAU();
      nAD= cQ.GetAD();
      nAS= cQ.GetAS();
    }
	if(nU==1&&nAU==1||nD==1&&nAD==1) // rho0/omega case
	{
      if (G4UniformRand()>0.5) return 223; // omega
      else                     return 113; // rho0
	}
    else if(nU==1&&nAD==1)     return 213; // rho+
    else if(nD==1&&nAU==1)     return -213;// rho-
    else if(nD==1&&nAS==1)     return 313; // K0*
    else if(nU==1&&nAS==1)     return 323; // K+*
    else if(nS==1&&nAD==1)     return -313;// AK0*
    else if(nS==1&&nAU==1)     return -323;// K-*
    else if(nS==1&&nAS==1)     return 333; // phi
	else G4cerr<<"***NofSMS:U"<<nU<<",D"<<nD<<",S"<<nS<<",AU"<<nAU<<",AD"<<nAD<<",AS"<<nAS<<G4endl;
  }
  return 0;
}

// Calculate a momentum of quark-parton greater then minimum value kMin
G4double G4Quasmon::GetQPartonMomentum(G4double kMin)
{
  //gives k>kMin QParton Momentum for the current Quasmon
#ifdef debug
  cout << "G4Quasmon::GetQPartonMomentum is called with kMin = " << kMin << G4endl;
#endif
  G4double qMass = q4Mom.m();
  if (kMin+kMin>qMass || qMass<=0. || nOfQ<2) return 0.;

  G4double kMax=qMass/2.;
#ifdef debug
  cout << "GetQPartonMomentum: kMax = " << kMax << ", nQ=" << nOfQ << G4endl;
#endif
  if (nOfQ==2) return kMax;
  G4int n=nOfQ-2;
  G4double fn=n;
  G4double vRndm = G4UniformRand();
  if (kMin>0.) // If there is a minimum cut for the QpMomentum 
  {
    G4double xR=kMin/kMax;
    vRndm = vRndm*pow((1.-xR),n)*(1.+n*xR);
  }
  if (n==1) return kMax*sqrt(1.-vRndm); // Direct solution
  else                                  // Needs iterations
  {
    G4double x = 1. - pow(vRndm*(1+n*vRndm)/(fn+1.),1/fn);
    int it=0;
    //while (it<7)
    G4double d = 1.;
    G4double df= 1./static_cast<double>(nOfQ);
    G4double f = df*(static_cast<int>(nOfQ*nOfQ*n*x/5.)+(nOfQ/2));
#ifdef debug
    cout << "First x=" << x << ", f=" << f << G4endl;
#endif
    while ( abs(d/vRndm) > 0.001 )
	{
      G4double r = 1. - x;
      G4double p = r;
      if (n>2) p = pow(r,n-1);
      G4double nx = n*x;
      G4double c = p*r*(1.+nx);
      d = c - vRndm;
#ifdef debug
	  cout << "Iteration#" << it << ": (c=" << c << " - R=" << vRndm 
           << ") =" << d/vRndm << ", x=" << x << ", f=" << f << G4endl;
#endif
      x = x + f*d/(r*nx*(fn+1.));
      if (f>1.0001)    f-=df;
      else if (f<.999) f+=df;
      it++;
	};
    return kMax*x;
  }
  return 0.;
} 

void G4Quasmon::CalculateContentOfQPartons(G4double qMass)
{
  // @@ Temporary here. To have 3 quarks in Nucleon Temperature should be < M_N/4 (234 MeV)
  static const G4double Temperature = 200.; // "Vacuum Temperature" parameter of the CHIPS model
  // M^2=4*n*(n-1)*T^2 => n = M/2T + 1/2 + T/4M + O (T^3/16M^3)
#ifdef debug
  cout << "#ofQP is called: qMass = " << qMass << ", T=" << Temperature<< G4endl;
#endif
  G4double qMOver2T = qMass/(Temperature+Temperature);
  G4double est = qMOver2T+1.+0.125/qMOver2T;
  G4int valc = valQ.GetTot();
  // @@ two edges integer values can be randomized later on ...
  nOfQ = est;          // nOfQ canbe an odd number !
  if (nOfQ > valc)
  {
    if ((nOfQ-valc)%2) // Randomize the case of an odd N partons
    {
      if (G4UniformRand()>0.5) nOfQ++;
      else                     nOfQ--;
    }
    G4int nSeaPairs = (nOfQ-valc)/2;
#ifdef debug
    cout<<"#ofQP: #Initial="<<valc<<", #Final="<<nOfQ<<", #Sea="<<nSeaPairs<<G4endl;
#endif
	if (nSeaPairs) valQ.IncQAQ(nSeaPairs);
  }
  else nOfQ=valc;
} 

void G4Quasmon::CalculateHadronizationProbabilities(G4double kVal)
{
  G4double mQ = q4Mom.m();
  G4double accumulatedProbability = 0.;
#ifdef debug
  cout << "QHadProb: Quasmon mass="<<mQ << ", Quark Content (QC) ="<<valQ << G4endl;
#endif
  for (G4int index=0; index<theQCandidates.entries(); index++)
  {
    G4double mH2 = theQCandidates[index]->GetMass2();
    G4double probability = 0.;
    G4double possibility = mH2/kVal/mQ;
    if (possibility < 2.)
	{
      if (nOfQ>3) probability = pow ( 1.- 0.5*possibility , nOfQ-3 );
      else        probability = 1.;
	}
#ifdef debug
    cout << "QHadProb: #"<<index <<", m2="<<mH2 << ", n="<<nOfQ << ", p="<<probability << G4endl;
#endif
    G4QContent candQC = theQCandidates[index]->GetQC();
    probability *= valQ.NOfCombinations(candQC);
    theQCandidates[index]->SetRelProbability(probability);
    accumulatedProbability += probability;
    theQCandidates[index]->SetIntegProbability(accumulatedProbability);
#ifdef debug
    cout << "QHadProb: Final p="<<probability << ", P_sum="<<accumulatedProbability << G4endl;
#endif
  }
}



