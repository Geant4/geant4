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
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      ---------------- G4QString ----------------
//      by Mikhail Kossov Oct, 2006
//      class for excited string used by Parton String Models
//   For comparison mirror member functions are taken from G4 classes:
//   G4FragmentingString
//   G4ExcitedStringDecay
// ---------------------------------------------------------------------
// Short description: If partons from the G4QPartonPair are close in
// rapidity, they create Quasmons, but if they are far in the rapidity
// space, they can not interact directly. Say the bottom parton (quark)
// has rapidity 0, and the top parton (antiquark) has rapidity 8, then
// the top quark splits in two by radiating gluon, and each part has
// rapidity 4, then the gluon splits in quark-antiquark pair (rapidity
// 2 each), and then the quark gadiates anothe gluon and reachs rapidity
// 1. Now it can interact with the bottom antiquark, creating a Quasmon
// or a hadron. The intermediate partons is the string ladder.
// ---------------------------------------------------------------------

//#define debug
//#define pdebug
//#define edebug

#include <algorithm>

#include "G4QString.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

// Static parameters definition
G4double G4QString::MassCut=350.*MeV;     // minimum mass cut for the string
G4double G4QString::SigmaQT=0.5*GeV;      // quarkTransverseMomentum distribution parameter
G4double G4QString::DiquarkSuppress=0.1;  // is Diquark suppression parameter  
G4double G4QString::DiquarkBreakProb=0.1; // is Diquark breaking probability 
G4double G4QString::SmoothParam=0.9;      // QGS model parameter
G4double G4QString::StrangeSuppress=0.435;// Strangeness suppression (u:d:s=1:1:0.3 ?M.K.)
G4double G4QString::widthOfPtSquare=-0.72*GeV*GeV; // pt -width2 forStringExcitation

G4QString::G4QString() : theDirection(0), thePosition(G4ThreeVector(0.,0.,0.)),
                         theStableParton(0), theDecayParton(0){}

G4QString::G4QString(G4QParton* Color, G4QParton* AntiColor, G4int Direction)
  : SideOfDecay(0)
{
#ifdef debug
  G4cout<<"G4QString::PPD-Constructor: Direction="<<Direction<<G4endl;
#endif
  ExciteString(Color, AntiColor, Direction);
#ifdef debug
  G4cout<<"G4QString::PPD-Constructor: ------>> String is excited"<<G4endl;
#endif
}

G4QString::G4QString(G4QPartonPair* CAC): SideOfDecay(0)
{
#ifdef debug
  G4cout<<"G4QString::PartonPair-Constructor: Is CALLED"<<G4endl;
#endif
  ExciteString(CAC->GetParton1(), CAC->GetParton2(), CAC->GetDirection());
#ifdef debug
  G4cout<<"G4QString::PartonPair-Constructor: ------>> String is excited"<<G4endl;
#endif
}

G4QString::G4QString(G4QParton* QCol,G4QParton* Gluon,G4QParton* QAntiCol,G4int Direction)
 : theDirection(Direction), thePosition(QCol->GetPosition()), SideOfDecay(0)
{
  thePartons.push_back(QCol);
  G4LorentzVector sum=QCol->Get4Momentum();
  thePartons.push_back(Gluon);
  sum+=Gluon->Get4Momentum();
  thePartons.push_back(QAntiCol);
  sum+=QAntiCol->Get4Momentum();
  Pplus =sum.e() + sum.pz();
  Pminus=sum.e() - sum.pz();
  decaying=None;
}

G4QString::G4QString(const G4QString &right) : theDirection(right.GetDirection()),
					       thePosition(right.GetPosition()), SideOfDecay(0)
{
  //LeftParton=right.LeftParton;
  //RightParton=right.RightParton;
  Ptleft=right.Ptleft;
  Ptright=right.Ptright;
  Pplus=right.Pplus;
  Pminus=right.Pminus;
  decaying=right.decaying;
}

G4QString::~G4QString()
{if(thePartons.size()) std::for_each(thePartons.begin(),thePartons.end(),DeleteQParton());}

G4LorentzVector G4QString::Get4Momentum() const
{
  G4LorentzVector momentum(0.,0.,0.,0.);
  for(unsigned i=0; i<thePartons.size(); i++) momentum += thePartons[i]->Get4Momentum();
  return momentum;
}

void G4QString::LorentzRotate(const G4LorentzRotation & rotation)
{
  for(unsigned i=0; i<thePartons.size(); i++)
     thePartons[i]->Set4Momentum(rotation*thePartons[i]->Get4Momentum());
}

//void G4QString::InsertParton(G4QParton* aParton, const G4QParton* addafter)
//{
//  G4QPartonVector::iterator insert_index;   // Begin by default (? M.K.)
//  if(addafter != NULL) 
//  {
//    insert_index=std::find(thePartons.begin(), thePartons.end(), addafter);
//    if (insert_index == thePartons.end())  // No object addafter in thePartons
//    {
//      G4cerr<<"***G4QString::InsertParton: Addressed Parton is not found"<<G4endl;
//      G4Exception("G4QString::InsertParton:","72",FatalException,"NoAddressParton");
//    }
//  }
//  thePartons.insert(insert_index+1, aParton);
//} 

void G4QString::Boost(G4ThreeVector& Velocity)
{
  for(unsigned cParton=0; cParton<thePartons.size(); cParton++)
  {
    G4LorentzVector Mom = thePartons[cParton]->Get4Momentum();
    Mom.boost(Velocity);
    thePartons[cParton]->Set4Momentum(Mom);
  }
}

//G4QParton* G4QString::GetColorParton(void) const
//{
//  G4QParton* start = *(thePartons.begin());
//  G4QParton* end = *(thePartons.end()-1);
//  G4int Encoding = start->GetPDGCode();
//  if (Encoding<-1000 || (Encoding  < 1000 && Encoding > 0)) return start;
//  return end; 
//}

//G4QParton* G4QString::GetAntiColorParton(void) const
//{
//  G4QParton* start = *(thePartons.begin());
//  G4QParton* end = *(thePartons.end()-1);
//  G4int Encoding = start->GetPDGCode();
//  if (Encoding < -1000 || (Encoding  < 1000 && Encoding > 0)) return end; 
//  return start; 
//}

// Fill parameters
void G4QString::SetParameters(G4double mCut, G4double sigQT, G4double DQSup, G4double DQBU,
                              G4double smPar, G4double SSup, G4double SigPt)
{
  MassCut         = mCut;           // minimum mass cut for the string
  SigmaQT         = sigQT;          // quark transverse momentum distribution parameter 
  DiquarkSuppress = DQSup;          // is Diquark suppression parameter  
  DiquarkBreakProb= DQBU;           // is Diquark breaking probability 
  SmoothParam     = smPar;          // QGS model parameter
  StrangeSuppress = SSup;           // Strangeness suppression parameter
  widthOfPtSquare = -2*SigPt*SigPt; // width^2 of pt for string excitation
}

// Pt distribution @@ one can use 1/(1+A*Pt^2)^B
G4ThreeVector G4QString::GaussianPt(G4double widthSquare, G4double maxPtSquare) const
{
  G4double pt2; do{pt2=widthSquare*std::log(G4UniformRand());} while (pt2>maxPtSquare);
  pt2=std::sqrt(pt2);
  G4double phi=G4UniformRand()*twopi;
  return G4ThreeVector(pt2*std::cos(phi),pt2*std::sin(phi),0.);    
} // End of GaussianPt

// Diffractively excite the string
//void G4QString::DiffString(G4QHadron* hadron, G4bool isProjectile)
//{
//  hadron->SplitUp();
//  G4QParton* start = hadron->GetNextParton();
//  if( start==NULL) 
//  {
//    G4cerr<<"***G4QString::DiffString: No start parton found, nothing is done"<<G4endl;
//    return;
//  }
//  G4QParton* end = hadron->GetNextParton();
//  if( end==NULL) 
//  {
//    G4cerr<<"***G4QString::DiffString: No end parton found, nothing is done"<<G4endl;
//    return;
//  }
//  if(isProjectile) ExciteString(end, start, 1); //  1 = Projectile
//  else             ExciteString(start, end,-1); // -1 = Target
//  SetPosition(hadron->GetPosition());
//  // momenta of string ends
//  G4double ptSquared= hadron->Get4Momentum().perp2();
//  G4double hmins=hadron->Get4Momentum().minus();
//  G4double hplus=hadron->Get4Momentum().plus();
//  G4double transMassSquared=hplus*hmins;
//  G4double maxMomentum = std::sqrt(transMassSquared) - std::sqrt(ptSquared);
//  G4double maxAvailMomentumSquared = maxMomentum * maxMomentum;
//  G4ThreeVector pt=GaussianPt(widthOfPtSquare,maxAvailMomentumSquared);
//  G4LorentzVector Pstart(G4LorentzVector(pt,0.));
//  G4LorentzVector Pend(hadron->Get4Momentum().px(), hadron->Get4Momentum().py(), 0.);
//  Pend-=Pstart;
//  G4double tm1=hmins+(Pend.perp2()-Pstart.perp2())/hplus;
//  G4double tm2=std::sqrt( std::max(0., tm1*tm1-4*Pend.perp2()*hmins/hplus ) );
//  G4int Sign= isProjectile ? TARGET : PROJECTILE;
//  G4double endMinus  = 0.5 * (tm1 + Sign*tm2);
//  G4double startMinus= hmins - endMinus;
//  G4double startPlus = Pstart.perp2() / startMinus;
//  G4double endPlus   = hplus - startPlus;
//  Pstart.setPz(0.5*(startPlus - startMinus));
//  Pstart.setE (0.5*(startPlus + startMinus));
//  Pend.setPz  (0.5*(endPlus - endMinus));
//  Pend.setE   (0.5*(endPlus + endMinus));
//  start->Set4Momentum(Pstart);
//  end->Set4Momentum(Pend);
//#ifdef debug
//  G4cout<<"G4QString::DiffString: StartQ="<<start->GetPDGCode()<<start->Get4Momentum()<<"("
//        <<start->Get4Momentum().mag()<<"), EndQ="<<end->GetPDGCode()<<end ->Get4Momentum()
//        <<"("<<end->Get4Momentum().mag()<<"), sumOfEnds="<<Pstart+Pend<<", H4M(original)="
//        <<hadron->Get4Momentum()<<G4endl;
//#endif
//} // End of DiffString (The string is excited diffractively)

// Excite the string by two partons
void G4QString::ExciteString(G4QParton* Color, G4QParton* AntiColor, G4int Direction)
{
#ifdef debug
  G4cout<<"G4QString::ExciteString: **Called**, direction="<<Direction<<G4endl;
#endif
  if(thePartons.size()) std::for_each(thePartons.begin(),thePartons.end(),DeleteQParton());
  thePartons.clear();
  theDirection = Direction;
  thePosition  = Color->GetPosition();
#ifdef debug
  G4cout<<"G4QString::ExciteString: ColourPosition = "<<thePosition<<", beg="<<Color->GetPDGCode()
          <<",end="<<AntiColor->GetPDGCode()<<G4endl;
#endif
  thePartons.push_back(Color);
  G4LorentzVector sum=Color->Get4Momentum();
  thePartons.push_back(AntiColor); // @@ Complain of Valgrind
  sum+=AntiColor->Get4Momentum();
  Pplus =sum.e() + sum.pz();
  Pminus=sum.e() - sum.pz();
  decaying=None;
#ifdef debug
  G4cout<<"G4QString::ExciteString: ***Done***, beg="<<(*thePartons.begin())->GetPDGCode()
          <<",end="<<(*(thePartons.end()-1))->GetPDGCode()<<G4endl;
#endif
} // End of ExciteString

// LUND Longitudinal fragmentation
G4double G4QString::GetLundLightConeZ(G4double zmin, G4double zmax, G4int , // @@ ? M.K.
                                      G4QHadron* pHadron, G4double Px, G4double Py)
{
  static const G4double  alund = 0.7/GeV/GeV; 
  // If blund get restored, you MUST adapt the calculation of zOfMaxyf.
  //static const G4double  blund = 1;
  G4double z, yf;
  G4double Mt2 = Px*Px + Py*Py + pHadron->GetMass2();
  G4double zOfMaxyf=alund*Mt2/(alund*Mt2+1.);
  G4double maxYf=(1.-zOfMaxyf)/zOfMaxyf * std::exp(-alund*Mt2/zOfMaxyf);
  do
  {
     z = zmin + G4UniformRand()*(zmax-zmin);
     // yf = std::pow(1. - z, blund)/z*std::exp(-alund*Mt2/z);
     yf = (1-z)/z * std::exp(-alund*Mt2/z);
  } while (G4UniformRand()*maxYf>yf); 
  return z;
} // End of GetLundLightConeZ


// QGSM Longitudinal fragmentation
G4double G4QString::GetQGSMLightConeZ(G4double zmin, G4double zmax, G4int PartonEncoding,
                                      G4QHadron* , G4double, G4double) // @@ M.K.
{
  static const G4double arho=0.5;
  static const G4double aphi=0.;
  static const G4double an=-0.5;
  static const G4double ala=-0.75;
  static const G4double aksi=-1.;
  static const G4double alft=0.5;
  G4double z;    
  G4double theA(0), d2, yf;
  G4int absCode = std::abs(PartonEncoding);
  // @@ Crazy algorithm is used for simple power low...
  if (absCode < 10)                            // Quarks (@@ 9 can be a gluon)
  { 
    if(absCode == 1 || absCode == 2) theA = arho;
    else if(absCode == 3)            theA = aphi;
    else
    {
     G4cerr<<"***G4QString::GetQGSMLightConeZ: CHIPS is SU(3), quakCode="<<absCode<<G4endl;
     G4Exception("G4QString::GetQGSMLightConeZ:","72",FatalException,"WrongQuark");
    }
    do  
    {
      z  = zmin + G4UniformRand()*(zmax - zmin);
      d2 = alft - theA;
      yf = std::pow( 1.-z, d2);
    } 
    while (G4UniformRand()>yf);
  }
  else                                     // Di-quarks (@@ Crazy codes are not checked)
  {       
    if (absCode==3101||absCode==3103||absCode==3201||absCode==3203)   d2=alft-ala-ala+arho;
    else if(absCode==1103||absCode==2101||absCode==2203||absCode==2103) d2=alft-an-an+arho;
    else                                                            d2=alft-aksi-aksi+arho;
    do  
    {
      z = zmin + G4UniformRand()*(zmax - zmin);
      yf = std::pow(1.-z, d2);
    } 
    while (G4UniformRand()>yf); 
  }
  return z;
} // End of GetQGSMLightConeZ

// Diffractively excite the string (QL=true - QGS Light Cone, =false - Lund Light Cone)
G4QHadronVector* G4QString::FragmentString(G4bool QL)
{
  // Can no longer modify Parameters for Fragmentation.
#ifdef edebug
  G4LorentzVector string4M=Get4Momentum();         // Just for Energy-Momentum ConservCheck
#endif 
#ifdef debug
  G4cout<<"G4QString::FragmentString:-->Called,QL="<<QL<<", M="<<Get4Momentum().m()<<", L="
        <<GetLeftParton()->Get4Momentum()<<",R="<<GetRightParton()->Get4Momentum()<<G4endl;
#endif 
  //  check if string has enough mass to fragment. If not, convert to one or two hadrons
  G4QHadronVector* LeftVector = LightFragmentationTest();
  if(LeftVector)
  {
#ifdef edebug
    G4LorentzVector sL=string4M;
    for(unsigned L = 0; L < LeftVector->size(); L++)
    {
      G4QHadron* LH = (*LeftVector)[L];
      G4LorentzVector L4M=LH->Get4Momentum();
      sL-=L4M;
      G4cout<<"-EMC-G4QStr::FragStr:L#"<<L<<",PDG="<<LH->GetPDGCode()<<",4M="<<L4M<<G4endl;
    }
    G4cout<<"-EMC-G4QString::FragmentString:---LightFragmentation---> Res4M="<<sL<<G4endl;
#endif 
    return LeftVector; //@@ Just decay in 2 or 1 (?) hadron, if below theCut
  }
#ifdef debug
  G4cout<<"G4QString::FragmentString:OUTPUT is not yet defined, define Left/Right"<<G4endl;
#endif  
  LeftVector = new G4QHadronVector; // Valgrind complain to LeftVector
  G4QHadronVector* RightVector = new G4QHadronVector;
  // Remember 4-momenta of the string ends (@@ only for the two-parton string, no gluons)
  G4LorentzVector left4M=GetLeftParton()->Get4Momentum(); // For recovery when failed
  G4LorentzVector right4M=GetRightParton()->Get4Momentum();
#ifdef debug
  G4cout<<"G4QString::FragmString: ***Remember*** L4M="<<left4M<<", R4M="<<right4M<<G4endl;
#endif
  G4int leftPDG=GetLeftParton()->GetPDGCode();
  G4int rightPDG=GetRightParton()->GetPDGCode();
  // Transform string to CMS
  G4LorentzVector momentum=Get4Momentum();
  G4LorentzRotation toCms(-(momentum.boostVector()));
  momentum= toCms * thePartons[0]->Get4Momentum();
  toCms.rotateZ(-1*momentum.phi());
  toCms.rotateY(-1*momentum.theta());
  for(unsigned index=0; index<thePartons.size(); index++)
  {
    momentum = toCms * thePartons[index]->Get4Momentum(); // @@ reuse of the momentum
    thePartons[index]->Set4Momentum(momentum);
  }
  // Copy the string for independent attempts
  G4QParton* LeftParton = new G4QParton(GetLeftParton());
  G4QParton* RightParton= new G4QParton(GetRightParton());
  G4QString* theStringInCMS = new G4QString(LeftParton,RightParton,GetDirection());
#ifdef debug
  G4cout<<"G4QString::FragmentString: Copy with nP="<<theStringInCMS->thePartons.size()
        <<", beg="<<(*(theStringInCMS->thePartons.begin()))->GetPDGCode()
        <<", end="<<(*(theStringInCMS->thePartons.end()-1))->GetPDGCode()<<G4endl;
#endif 
  G4bool success=false;
  G4bool inner_sucess=true;
  G4int attempt=0;
  G4int StringLoopInterrupt=27;  // String fragmentation LOOP limit 
#ifdef debug
  G4cout<<"G4QString::FragmentString: BeforeWhileLOOP, max = "<<StringLoopInterrupt
        <<", nP="<<thePartons.size()<<", beg="<<(*thePartons.begin())->GetPDGCode()
        <<",end="<<(*(thePartons.end()-1))->GetPDGCode()<<G4endl;
#endif
#ifdef edebug
  G4LorentzVector cmS4M=theStringInCMS->Get4Momentum();
  G4cout<<"-EMC-G4QString::FragmString: c4M="<<cmS4M<<",Max="<<StringLoopInterrupt<<G4endl;
#endif
  while (!success && attempt++ < StringLoopInterrupt) // Try fragment String till success
  {
    // Recover the CMS String
    LeftParton = new G4QParton(theStringInCMS->GetLeftParton());
    RightParton= new G4QParton(theStringInCMS->GetRightParton());
    ExciteString(LeftParton, RightParton, theStringInCMS->GetDirection());
#ifdef edebug
    G4LorentzVector cm4M=cmS4M;    // Copy the full momentum for the reduction and check
    G4cout<<"-EMC-.G4QString::FragmentString: CHEK "<<cm4M<<" ?= "<<Get4Momentum()<<G4endl;
#endif 
#ifdef debug
    G4cout<<"G4QString::FragmentString:===>LOOP, attempt = "<<attempt<<", nP="
          <<thePartons.size()<<", beg="<<(*thePartons.begin())->GetPDGCode()
          <<",end="<<(*(thePartons.end()-1))->GetPDGCode()<<G4endl;
#endif 
    // Now clean up all hadrons in the Left and Right vectors for the new attempt
    if(LeftVector->size())
    {
      std::for_each(LeftVector->begin(), LeftVector->end(), DeleteQHadron());
      LeftVector->clear();
    }
    //delete LeftVector; // @@ Valgrind ?
    if(RightVector->size())
    {
      std::for_each(RightVector->begin(), RightVector->end(), DeleteQHadron());
      RightVector->clear();
    }
    //delete RightVector; // @@ Valgrind ?
    inner_sucess=true;                                           // set false on failure
    while (!StopFragmentation())                                 // LOOP with break
    {  // Split current string into hadron + new string state
#ifdef debug
      G4cout<<"G4QString::FragmentString:-->Begin LOOP/LOOP, decaying="<<decaying<<G4endl;
#endif 
      G4QHadron* Hadron=Splitup(QL); // MAIN: where the hadron is split from the string
#ifdef edebug
      cm4M-=Hadron->Get4Momentum();
      G4cout<<"-EMC-G4QString::FragmentString:LOOP/LOOP,d4M="<<cm4M-Get4Momentum()<<G4endl;
#endif
      G4bool canBeFrag=IsFragmentable();
#ifdef debug
      G4cout<<"G4QString::FragmentString: LOOP/LOOP, canBeFrag="<<canBeFrag<<", decay="
            <<decaying<<", H="<<Hadron<<", newStringMass="<<Get4Momentum().m()<<G4endl;
#endif 
      if(Hadron && canBeFrag)
      {
#ifdef debug
        G4cout<<">>G4QString::FragmentString: LOOP/LOOP-NO FRAGM-, dec="<<decaying<<G4endl;
#endif 
        if(GetDecayDirection()>0) LeftVector->push_back(Hadron);
        else RightVector->push_back(Hadron);
      }
      else
      {
        // @@ Try to convert to the resonance and decay, abandon if M<mGS+mPI0
        // abandon ... start from the beginning
#ifdef debug
        G4cout<<"G4QString::FragmentString: LOOP/LOOP, Start from scratch"<<G4endl;
#endif 
        if (Hadron) delete Hadron;
        inner_sucess=false;
        break;
      }
#ifdef debug
      G4cout<<"G4QString::FragmentString: LOOP/LOOP End, nP="
            <<thePartons.size()<<", beg="<<(*thePartons.begin())->GetPDGCode()
            <<",end="<<(*(thePartons.end()-1))->GetPDGCode()<<G4endl;
#endif 
    } 
#ifdef edebug
    G4LorentzVector fLR=cmS4M-Get4Momentum();
    for(unsigned L = 0; L < LeftVector->size(); L++)
    {
      G4QHadron* LH = (*LeftVector)[L];
      G4LorentzVector L4M=LH->Get4Momentum();
      fLR-=L4M;
      G4cout<<"-EMC-G4QStr::FrStr:L#"<<L<<",PDG="<<LH->GetPDGCode()<<",4M="<<L4M<<G4endl;
    }
    for(unsigned R = 0; R < RightVector->size(); R++)
    {
      G4QHadron* RH = (*RightVector)[R];
      G4LorentzVector R4M=RH->Get4Momentum();
      fLR-=R4M;
      G4cout<<"-EMC-G4QStr::FrStr:R#"<<R<<",PDG="<<RH->GetPDGCode()<<",4M="<<R4M<<G4endl;
    }
    G4cout<<"-EMC-G4QString::FragmentString:L/R_BeforLast->r4M/M2="<<fLR<<fLR.m2()<<G4endl;
#endif 
    // Split current string into 2 final Hadrons
#ifdef debug
    G4cout<<"G4QString::FragmentString: inner_success = "<<inner_sucess<<G4endl;
#endif
    if(inner_sucess)
    {
      success=true;                                      // Default prototype
      //... perform last cluster decay
      G4LorentzVector tot4M = Get4Momentum();
      G4double totM    = tot4M.m(); 
#ifdef debug
      G4cout<<"G4QString::FragmString: string4M="<<tot4M<<totM<<G4endl;
#endif 
      G4QHadron* LeftHadron;
      G4QHadron* RightHadron;
      G4QParton* RQuark = 0;
      SetLeftPartonStable();              // to query quark contents
      if(DecayIsQuark() && StableIsQuark()) // There're quarks on clusterEnds
      {
#ifdef debug
        G4cout<<"G4QString::FragmentString: LOOP Quark Algorithm"<<G4endl;
#endif 
        LeftHadron= QuarkSplitup(GetLeftParton(), RQuark);
      }
      else
      {
#ifdef debug
        G4cout<<"G4QString::FragmentString: LOOP Di-Quark Algorithm"<<G4endl;
#endif 
        //... there is a Diquark on cluster ends
        G4int IsParticle;
        if(StableIsQuark()) IsParticle=(GetLeftParton()->GetPDGCode()>0)?-1:1;
        else                IsParticle=(GetLeftParton()->GetPDGCode()>0)?1:-1;
        G4QPartonPair QuarkPair = CreatePartonPair(IsParticle,false); // no diquarks
        RQuark = QuarkPair.GetParton2();
        G4QParton* LQuark = QuarkPair.GetParton1();
        LeftHadron = CreateHadron(LQuark, GetLeftParton()); // Create Left Hadron
        delete LQuark;                                      // Delete the temporaryParton
      }
      RightHadron = CreateHadron(GetRightParton(), RQuark); // Create Right Hadron
      delete RQuark;                                        // Delete the temporaryParton
      //... repeat procedure, if mass of cluster is too low to produce hadrons
      G4double LhM=LeftHadron->GetMass();
      G4double RhM=RightHadron->GetMass();
#ifdef debug
      G4cout<<"G4QStr::FrSt:L="<<LeftHadron->GetPDGCode()<<",R="<<RightHadron->GetPDGCode()
            <<",ML="<<LhM<<",MR="<<RhM<<",SumM="<<LhM+RhM<<",tM="<<totM<<G4endl;
#endif
      if(totM < LhM + RhM) success=false;
      //... compute hadron momenta and energies   
      if(success)
      {
        G4ThreeVector    Pos=GetPosition();
        G4LorentzVector  Lh4M(0.,0.,0.,LhM);
        G4LorentzVector  Rh4M(0.,0.,0.,RhM);
        if(G4QHadron(tot4M).DecayIn2(Lh4M,Rh4M))
        {
          LeftVector->push_back(new G4QHadron(LeftHadron, 0, Pos, Lh4M));
          delete LeftHadron;
          RightVector->push_back(new G4QHadron(RightHadron, 0, Pos, Rh4M));
          delete RightHadron;
        }
#ifdef debug
        G4cout<<"->>G4QStr::FragString:HFilled (L) PDG="<<LeftHadron->GetPDGCode()<<", 4M="
              <<Lh4M<<", (R) PDG="<<RightHadron->GetPDGCode()<<", 4M="<<Rh4M<<G4endl;
#endif
#ifdef edebug
	  G4cout<<"-EMC-G4QString::FragmentString: Residual4M="<<tot4M-Lh4M-Rh4M<<G4endl;
#endif
      }
      else
      {
        if(LeftHadron)  delete LeftHadron;
        if(RightHadron) delete RightHadron;
      }
    } // End of inner success
  } // End of while
  delete theStringInCMS;
#ifdef debug
  G4cout<<"G4QString::FragmentString: LOOP/LOOP, success="<<success<<G4endl;
#endif 
  if (!success)
  {
    if(RightVector->size())
    {
      std::for_each(RightVector->begin(), RightVector->end(), DeleteQHadron());
      RightVector->clear();
    }
    delete RightVector;
    if(LeftVector->size())
    {
      std::for_each(LeftVector->begin(), LeftVector->end(), DeleteQHadron());
      LeftVector->clear();
    }
    delete LeftVector;
#ifdef debug
    G4cout<<"G4QString::FragmString:StringNotFragm,L4M="<<left4M<<",R4M="<<right4M<<G4endl;
#endif
    // Recover the Left/Right partons 4-moms of the String in ZLS
    GetLeftParton()->SetPDGCode(leftPDG);
    GetRightParton()->SetPDGCode(rightPDG);
    GetLeftParton()->Set4Momentum(left4M);
    GetRightParton()->Set4Momentum(right4M);
    return 0;                         // The string can not be fragmented
  }
  // @@@@@ Print collected Left and Right Hadrons (decay resonances!)
#ifdef edebug
  G4LorentzVector sLR=cmS4M;
  for(unsigned L = 0; L < LeftVector->size(); L++)
  {
    G4QHadron* LH = (*LeftVector)[L];
    G4LorentzVector L4M=LH->Get4Momentum();
    sLR-=L4M;
    G4cout<<"-EMC-G4QStr::FragmStri:L#"<<L<<",PDG="<<LH->GetPDGCode()<<",4M="<<L4M<<G4endl;
  }
  for(unsigned R = 0; R < RightVector->size(); R++)
  {
    G4QHadron* RH = (*RightVector)[R];
    G4LorentzVector R4M=RH->Get4Momentum();
    sLR-=R4M;
    G4cout<<"-EMC-G4QStr::FragmStri:R#"<<R<<",PDG="<<RH->GetPDGCode()<<",4M="<<R4M<<G4endl;
  }
  G4cout<<"-EMC-G4QString::FragmentString:---L/R_BeforeMerge---> Res4M="<<sLR<<G4endl;
#endif 
  // Join Left and Right Vectors into LeftVector in correct order.
  while(!RightVector->empty())
  {
     LeftVector->push_back(RightVector->back());
     RightVector->erase(RightVector->end()-1);
  }
  delete RightVector;
  // @@ A trick, the real bug should be found !!
  G4QHadronVector::iterator ilv;                                           // @@
  for(ilv = LeftVector->begin(); ilv < LeftVector->end(); ilv++)           // @@
  {
    G4ThreeVector CV=(*ilv)->Get4Momentum().vect();                        // @@
    if(CV.x()==0. && CV.y()==0. && CV.z()==0.) LeftVector->erase(ilv);     // @@
  }
  // Calculate time and position of hadrons with @@ very rough formation time
  G4double StringMass=Get4Momentum().mag();
  static const G4double dkappa = 2.0 * GeV/fermi; // @@ 2*kappa kappa=1 GeV/fermi (?)
  for(unsigned c1 = 0; c1 < LeftVector->size(); c1++)
  {
    G4double SumPz = 0; 
    G4double SumE  = 0;
    for(unsigned c2 = 0; c2 < c1; c2++)
    {
      G4LorentzVector hc2M=(*LeftVector)[c2]->Get4Momentum();
      SumPz += hc2M.pz();
      SumE  += hc2M.e();   
    }
    G4QHadron* hc1=(*LeftVector)[c1];
    G4LorentzVector hc1M=hc1->Get4Momentum();
    G4double HadronE = hc1M.e();
    G4double HadronPz= hc1M.pz();
    hc1->SetFormationTime((StringMass-SumPz-SumPz+HadronE-HadronPz)/dkappa);
    hc1->SetPosition(G4ThreeVector(0,0,(StringMass-SumE-SumE-HadronE+HadronPz)/dkappa));
  } 
  G4LorentzRotation toObserverFrame(toCms.inverse());
#ifdef debug
  G4cout<<"G4QString::FragmentString: beforeLoop LVsize = "<<LeftVector->size()<<G4endl;
#endif 
  for(unsigned C1 = 0; C1 < LeftVector->size(); C1++)
  {
    G4QHadron* Hadron = (*LeftVector)[C1];
    G4LorentzVector Momentum = Hadron->Get4Momentum();
    Momentum = toObserverFrame*Momentum;
    Hadron->Set4Momentum(Momentum);
    G4LorentzVector Coordinate(Hadron->GetPosition(), Hadron->GetFormationTime());
    Momentum = toObserverFrame*Coordinate;
    Hadron->SetFormationTime(Momentum.e());
    Hadron->SetPosition(GetPosition()+Momentum.vect());
  }
#ifdef edebug
  G4LorentzVector sLA=string4M;
  for(unsigned L = 0; L < LeftVector->size(); L++)
  {
    G4QHadron* LH = (*LeftVector)[L];
    G4LorentzVector L4M=LH->Get4Momentum();
    sLA-=L4M;
    G4cout<<"-EMC-G4QStr::FragmStri:L#"<<L<<",PDG="<<LH->GetPDGCode()<<",4M="<<L4M<<G4endl;
  }
  G4cout<<"-EMC-G4QString::FragmentString:---LSAfterMerge&Conv---> Res4M="<<sLA<<G4endl;
#endif 
#ifdef debug
  G4cout<<"G4QString::FragmentString: *** Done *** "<<G4endl;
#endif 
  return LeftVector; // Should be deleted by user (@@ Valgrind complain ?)
} // End of FragmentString

// Simple decay of the string if the excitation mass is too small for HE fragmentation
// !! If the mass is below the single hadron threshold, make warning (improve) and convert
// the string to the single S-hadron breaking energy conservation (temporary) and improve,
// taking the threshold into account on the level of the String creation (merge strings) !!
G4QHadronVector* G4QString::LightFragmentationTest()
{
  // Check string decay threshold
  G4LorentzVector tot4M=Get4Momentum();
#ifdef debug
  G4cout<<"G4QString::LightFragmentationTest: ***Called***, string4M="<<tot4M<<G4endl;
#endif 
  G4QHadronVector* result=0;  // return 0 when string exceeds the mass cut or below mh1+mh2
 
  G4QHadronPair hadrons((G4QHadron*)0, (G4QHadron*)0); // pair of hadrons for output of FrM
  G4double fragMass = FragmentationMass(0,&hadrons);   // Minimum mass to decay the string
#ifdef debug
  G4cout<<"G4QString::LightFragTest: before check nP="<<thePartons.size()<<", MS2="
        <<Mass2()<<", MCut="<<MassCut<<", beg="<<(*thePartons.begin())->GetPDGCode()
        <<",end="<<(*(thePartons.end()-1))->GetPDGCode()<<", fM="<<fragMass<<G4endl;
#endif  
  if(Mass2() > sqr(fragMass+MassCut))// Big enough to fragment in a lader (avoid the decay)
  {
    if(hadrons.first) delete hadrons.first;
    if(hadrons.second) delete hadrons.second;
#ifdef debug
    G4cout<<"G4QString::LightFragTest:NO,M2="<<Mass2()<<">"<<sqr(fragMass+MassCut)<<G4endl;
#endif  
    return result;                          // =0. Depends on the parameter of the Mass Cut
  }
  G4double totM= tot4M.m();
  G4QHadron* h1=hadrons.first;
  G4QHadron* h2=hadrons.second;
  if(h1 && h2)
  {
    G4double h1M = h1->GetMass();
    G4double h2M = h2->GetMass();
#ifdef debug
    G4cout<<"G4QString::LightFragTest:tM="<<totM<<","<<h1M<<"+"<<h2M<<"+"<<h1M+h2M<<G4endl;
#endif
    if(h1M + h2M <= totM)                   // The string can decay in these two hadrons
    {  
      // Create two stable hadrons
      G4LorentzVector  h4M1(0.,0.,0.,h1M);
      G4LorentzVector  h4M2(0.,0.,0.,h2M);
      if(G4QHadron(tot4M).DecayIn2(h4M1,h4M2))
      {
        result = new G4QHadronVector;        
        result->push_back(new G4QHadron(hadrons.first, 0, GetPosition(), h4M1));
        result->push_back(new G4QHadron(hadrons.second,0, GetPosition(), h4M2));
      }
#ifdef edebug
      G4int L=result->size(); if(L) for(G4int i=0; i<L; i++) 
      {
        tot4M-=(*result)[i]->Get4Momentum();
        G4cout<<"-EMC-G4QString::LightFragTest: i="<<i<<", residual4M="<<tot4M<<G4endl;
      }
#endif
    }
#ifdef debug
    else G4cout<<"-Warning-G4QString::LightFragTest: TooBigHadronMasses to decay"<<G4endl;
#endif
  }
#ifdef debug
  else G4cout<<"-Warning-G4QString::LightFragTest: No Hadrons have been proposed"<<G4endl;
#endif
  delete hadrons.first;
  delete hadrons.second;
  return result;
} // End of LightFragmentationTest

// Calculate Fragmentation Mass (if pdefs # 0, returns two hadrons)
G4double G4QString::FragmentationMass(G4int HighSpin, G4QHadronPair* pdefs)
{
  G4double mass=0.;
#ifdef debug
  G4cout<<"G4QString::FragmMass: ***Called***, s4M="<<Get4Momentum()<<G4endl;
#endif 
  // Example how to use an interface to different member functions
  G4QHadron* Hadron1 = 0;
  G4QHadron* Hadron2 = 0;
#ifdef debug
  G4cout<<"G4QString::FragmentationMass: Create spin-0 or spin-1/2 hadron: nP="
        <<thePartons.size()<<", beg="<<(*thePartons.begin())->GetPDGCode()
        <<",end="<<(*(thePartons.end()-1))->GetPDGCode()<<G4endl;
#endif 
  G4int iflc = (G4UniformRand() < 0.5) ? 1 : 2; // Create additional Q-antiQ pair @@ No S
  G4int LPDG= GetLeftParton()->GetPDGCode();
  G4int LT  = GetLeftParton()->GetType();
  if ( (LPDG > 0 && LT == 1)  || (LPDG < 0 && LT == 2) ) iflc = -iflc; // anti-quark
  G4QParton* piflc = new G4QParton( iflc);
  G4QParton* miflc = new G4QParton(-iflc);
  if(HighSpin)
  {
    Hadron1 = CreateHighSpinHadron(GetLeftParton(),piflc);
    Hadron2 = CreateHighSpinHadron(GetRightParton(),miflc);
#ifdef debug
    G4cout<<"G4QString::FragmentationMass: High, PDG1="<<Hadron1->GetPDGCode()
          <<", PDG2="<<Hadron2->GetPDGCode()<<G4endl;
#endif 
  }
  else
  {
    Hadron1 = CreateLowSpinHadron(GetLeftParton(),piflc);
    Hadron2 = CreateLowSpinHadron(GetRightParton(),miflc);
#ifdef debug
    G4cout<<"G4QString::FragmentationMass: Low, PDG1="<<Hadron1->GetPDGCode()
          <<", PDG2="<<Hadron2->GetPDGCode()<<G4endl;
#endif 
  }
  mass    = Hadron1->GetMass() + Hadron2->GetMass();
  if(pdefs) // need to return hadrons as well as the mass estimate 
  {
    pdefs->first  = Hadron1;  // To be deleted by the calling program if not zero
    pdefs->second = Hadron2;  // To be deleted by the calling program if not zero
  }
  else      // Forget about the hadrons
  {
    if(Hadron1) delete Hadron1;
    if(Hadron2) delete Hadron2;
  }
  delete piflc;
  delete miflc;
#ifdef debug
  G4cout<<"G4QString::FragmentationMass: ***Done*** mass="<<mass<<G4endl;
#endif 
  return mass;
} // End of FragmentationMass

void G4QString::SetLeftPartonStable()
{
  theStableParton=GetLeftParton();
  theDecayParton=GetRightParton();
  decaying=Right;
}

void G4QString::SetRightPartonStable()
{
  theStableParton=GetRightParton();
  theDecayParton=GetLeftParton();
  decaying=Left;
}

G4int G4QString::GetDecayDirection() const
{
  if      (decaying == Left ) return +1;
  else if (decaying == Right) return -1;
  else
  {
    G4cerr<<"***G4QString::GetDecayDirection: wrong DecayDirection="<<decaying<<G4endl;
    G4Exception("G4QString::GetDecayDirection:","72",FatalException,"WrongDecayDirection");
  }
  return 0;
}

//G4ThreeVector G4QString::StablePt()
//{
//  if (decaying == Left ) return Ptright;
//  else if (decaying == Right ) return Ptleft;
//  else
//  {
//    G4cerr<<"***G4QString::StablePt: wrong DecayDirection="<<decaying<<G4endl;
//    G4Exception("G4QString::StablePt:","72",FatalException,"WrongDecayDirection");
//  }
//  return G4ThreeVector();
//}

G4ThreeVector G4QString::DecayPt()
{
  if      (decaying == Left  ) return Ptleft;
  else if (decaying == Right ) return Ptright;
  else
  {
    G4cerr<<"***G4QString::DecayPt: wrong DecayDirection="<<decaying<<G4endl;
    G4Exception("G4QString::DecayPt:","72",FatalException,"WrongDecayDirection");
  }
  return G4ThreeVector();
}

// Random choice of string end to use it for creating the hadron (decay)   
G4QHadron* G4QString::Splitup(G4bool QL)
{
  SideOfDecay = (G4UniformRand() < 0.5) ? 1: -1;
#ifdef debug
  G4cout<<"G4QString::Splitup:**Called**,s="<<SideOfDecay<<",s4M="<<Get4Momentum()<<G4endl;
#endif
  if(SideOfDecay<0) SetLeftPartonStable();  // Decay Right parton
  else              SetRightPartonStable(); // Decay Left parton
  G4QParton* newStringEnd;
  G4QHadron* Hadron;
  if(DecayIsQuark()) Hadron= QuarkSplitup(theDecayParton, newStringEnd);   // Split Quark
  else               Hadron= DiQuarkSplitup(theDecayParton, newStringEnd); // Split DiQuark
#ifdef debug
  G4cout<<"G4QString::Splitup: newStringEndPDG="<<newStringEnd->GetPDGCode()<<", nP="
          <<thePartons.size()<<", beg="<<(*thePartons.begin())->GetPDGCode()
          <<",end="<<(*(thePartons.end()-1))->GetPDGCode()<<G4endl;
#endif
  // create new String from the old one: keep Left and Right order, but replace decay
  G4LorentzVector* HadronMomentum=SplitEandP(Hadron, QL);//The decayed parton isn't changed
#ifdef debug
  G4cout<<"G4QString::Splitup: HM="<<HadronMomentum<<", nP="
          <<thePartons.size()<<", beg="<<(*thePartons.begin())->GetPDGCode()
          <<",end="<<(*(thePartons.end()-1))->GetPDGCode()<<G4endl;
#endif
  if(HadronMomentum) // The decay succeeded, now the new 4-mon can be set to NewStringEnd
  {    
#ifdef pdebug
    G4cout<<"---->>G4QString::Splitup: HFilled 4M="<<*HadronMomentum<<",PDG="
          <<Hadron->GetPDGCode()<<",s4M-h4M="<<Get4Momentum()-*HadronMomentum<<G4endl;
#endif 
    newStringEnd->Set4Momentum(theDecayParton->Get4Momentum()-*HadronMomentum);
    Hadron->Set4Momentum(*HadronMomentum);
    Hadron->SetPosition(GetPosition());
    if(decaying == Left)
    {
      G4QParton* theFirst = thePartons[0];            // Substitute for the First Parton
      delete theFirst;                                // The OldParton instance is deleted
      thePartons[0] = newStringEnd;                   // Delete equivalent for newStringEnd
#ifdef debug
      G4cout<<"G4QString::Splitup:  theFirstPDG="<<theFirst->GetPDGCode()<<G4endl;
#endif 
      Ptleft  -= HadronMomentum->vect();
      Ptleft.setZ(0.);                                // @@ (Z is anyway ignored) M.K. (?)
    }
    else if (decaying == Right)
    {
      G4QParton* theLast = thePartons[thePartons.size()-1]; // Substitute for theLastParton
      delete theLast;                                 // The OldParton instance is deleted
      thePartons[thePartons.size()-1] = newStringEnd; // Delete equivalent for newStringEnd
#ifdef debug
      G4cout<<"G4QString::Splitup:  theLastPDG="<<theLast->GetPDGCode()<<", nP="
            <<thePartons.size()<<", beg="<<thePartons[0]->GetPDGCode()<<",end="
            <<thePartons[thePartons.size()-1]->GetPDGCode()<<",P="<<theLast
            <<"="<<thePartons[thePartons.size()-1]<<G4endl;
#endif 
      Ptright -= HadronMomentum->vect();
      Ptright.setZ(0.);                               // @@ (Z is anyway ignored) M.K. (?)
    }
    else
    {
      G4cerr<<"***G4QString::Splitup: wrong oldDecay="<<decaying<<G4endl;
      G4Exception("G4QString::Splitup","72",FatalException,"WrongDecayDirection");
    }
    Pplus  -= HadronMomentum->e() + HadronMomentum->pz();// Reduce Pplus ofTheString (Left)
    Pminus -= HadronMomentum->e() - HadronMomentum->pz();// Reduce Pminus ofTheString(Rite)
#ifdef debug
    G4cout<<"G4QString::Splitup:  P+="<<Pplus<<",P-="<<Pminus<<", nP="
          <<thePartons.size()<<", beg="<<(*thePartons.begin())->GetPDGCode()
          <<",end="<<(*(thePartons.end()-1))->GetPDGCode()<<G4endl;
    G4cout<<">...>G4QString::Splitup: NewString4M="<<Get4Momentum()<<G4endl;
#endif 
    delete HadronMomentum;
  }      
#ifdef debug
  G4cout<<"G4QString::Splitup: ***Done*** H4M="<<Hadron->Get4Momentum()<<", nP="
          <<thePartons.size()<<", beg="<<(*thePartons.begin())->GetPDGCode()
          <<",end="<<(*(thePartons.end()-1))->GetPDGCode()<<G4endl;
#endif 
  return Hadron;
} // End of Splitup

// QL=true for QGSM and QL=false for Lund fragmentation
G4LorentzVector* G4QString::SplitEandP(G4QHadron* pHadron, G4bool QL)
{
  G4double HadronMass = pHadron->GetMass();
#ifdef debug
  G4cout<<"G4QString::SplitEandP: ***Called*** HMass="<<HadronMass<<G4endl;
#endif 
  // calculate and assign hadron transverse momentum component HadronPx andHadronPy
  G4ThreeVector HadronPt = SampleQuarkPt() + DecayPt(); // @@ SampleQuarkPt & DecayPt once
  HadronPt.setZ(0.);
  //...  sample z to define hadron longitudinal momentum and energy
  //... but first check the available phase space
  G4double HadronMass2T = HadronMass*HadronMass + HadronPt.mag2();
  if (HadronMass2T >= SmoothParam*Mass2() )  return 0;  // restart!
  //... then compute allowed z region  z_min <= z <= z_max 
  G4double zMin = HadronMass2T/Mass2();
  G4double zMax = 1.;
#ifdef debug
  G4cout<<"G4QString::SplitEandP: zMin="<<zMin<<", zMax"<<zMax<<G4endl;
#endif 
  if (zMin >= zMax) return 0;  // have to start all over!  
  G4double z=0;
  if(QL) z = GetQGSMLightConeZ(zMin, zMax, theDecayParton->GetPDGCode(), pHadron,
                               HadronPt.x(), HadronPt.y());      
  else   z = GetLundLightConeZ(zMin, zMax, theDecayParton->GetPDGCode(), pHadron,
                               HadronPt.x(), HadronPt.y());      
  //... now compute hadron longitudinal momentum and energy
  // longitudinal hadron momentum component HadronPz
  G4double zl= z;
  if      (decaying == Left  ) zl*=Pplus;
  else if (decaying == Right ) zl*=Pminus;
  else                                                // @@ Is that possible?
  {
    G4cerr<<"***G4QString::SplitEandP: wrong DecayDirection="<<decaying<<G4endl;
    G4Exception("G4QString::SplitEandP:","72",FatalException,"WrongDecayDirection");
  }
  G4double HadronE = (zl + HadronMass2T/zl)/2;
  HadronPt.setZ( GetDecayDirection() * (zl - HadronE) );
  G4LorentzVector* a4Momentum= new G4LorentzVector(HadronPt,HadronE);
  return a4Momentum;
}

G4ThreeVector G4QString::SampleQuarkPt()
{
   G4double Pt = SigmaQT * std::sqrt( -std::log(G4UniformRand()));
   G4double phi = twopi*G4UniformRand();
   return G4ThreeVector(Pt * std::cos(phi),Pt * std::sin(phi),0);
}

G4QHadron* G4QString::QuarkSplitup(G4QParton* decay, G4QParton* &created)// VGComplTo decay
{
  G4int IsParticle=(decay->GetPDGCode()>0) ? -1 : +1; // a quark needs antiquark or diquark
  G4QPartonPair QuarkPair = CreatePartonPair(IsParticle);
  created = QuarkPair.GetParton2();                   // New Parton after splitting
#ifdef debug
  G4cout<<"G4QString::QuarkSplitup: ***Called*** crP="<<created->GetPDGCode()<<G4endl;
#endif
  G4QParton* P1=QuarkPair.GetParton1();
  G4QHadron* result=CreateHadron(P1, decay);         // New Hadron after splitting
  delete P1;                                         // Clean up the temporary parton
  return result; 
} // End of QuarkSplitup

//
G4QHadron* G4QString::DiQuarkSplitup(G4QParton* decay, G4QParton* &created)
{
  //... can Diquark break or not? 
  if (G4UniformRand() < DiquarkBreakProb )
  {
    //... Diquark break
    G4int stableQuarkEncoding = decay->GetPDGCode()/1000;
    G4int decayQuarkEncoding = (decay->GetPDGCode()/100)%10;
    if (G4UniformRand() < 0.5)
    {
      G4int Swap = stableQuarkEncoding;
      stableQuarkEncoding = decayQuarkEncoding;
      decayQuarkEncoding = Swap;
    }
    G4int IsParticle=(decayQuarkEncoding>0) ? -1 : +1;// Diquark is equivalent to antiquark
    G4QPartonPair QuarkPair = CreatePartonPair(IsParticle,false);  // no diquarks wanted
    G4QParton* P2=QuarkPair.GetParton2();
    G4int QuarkEncoding=P2->GetPDGCode();
    delete P2;
    G4int i10  = std::max(std::abs(QuarkEncoding), std::abs(stableQuarkEncoding));
    G4int i20  = std::min(std::abs(QuarkEncoding), std::abs(stableQuarkEncoding));
    G4int spin = (i10 != i20 && G4UniformRand() <= 0.5) ? 1 : 3;
    G4int NewDecayEncoding = -1*IsParticle*(i10 * 1000 + i20 * 100 + spin);
    created = new G4QParton(NewDecayEncoding);
#ifdef debug
    G4cout<<"G4QString::DiQuarkSplitup: inside, crP="<<created->GetPDGCode()<<G4endl;
#endif 
    G4QParton* decayQuark= new G4QParton(decayQuarkEncoding);
    G4QParton* P1=QuarkPair.GetParton1();
    G4QHadron* newH=CreateHadron(P1, decayQuark);
    delete P1;
    delete decayQuark;
    return newH;
  }
  else
  {
    //... Diquark does not break
    G4int IsParticle=(decay->GetPDGCode()>0) ? +1 : -1; 
    G4QPartonPair QuarkPair = CreatePartonPair(IsParticle,false);  // no diquarks wanted
    created = QuarkPair.GetParton2();
#ifdef debug
    G4cout<<"G4QString::DiQuarkSplitup: diQ not break, crP="<<created->GetPDGCode()<<G4endl;
#endif 
    G4QParton* P1=QuarkPair.GetParton1();
    G4QHadron* newH=CreateHadron(P1, decay);
    delete P1;
    return newH;
 }
} // End of DiQuarkSplitup

G4QPartonPair G4QString::CreatePartonPair(G4int NeedParticle, G4bool AllowDiquarks)
{
#ifdef debug
		G4cout<<"G4QString::CreatePartonPair: ***Called***, P="<<NeedParticle<<", ALLOWdQ="
        <<AllowDiquarks<<G4endl;
#endif 
  //  NeedParticle = {+1 for Particle, -1 for AntiParticle}
  if(AllowDiquarks && G4UniformRand() < DiquarkSuppress)
  {
    // Create a Diquark - AntiDiquark pair , first in pair is anti to IsParticle
    G4int q1  = SampleQuarkFlavor();
    G4int q2  = SampleQuarkFlavor();
    G4int spin = (q1 != q2 && G4UniformRand() <= 0.5) ? 1 : 3; // @@ 0.5 M.K.?
    // Convention: quark with higher PDG number is first
    G4int PDGcode = (std::max(q1,q2) * 1000 + std::min(q1,q2) * 100 + spin) * NeedParticle;
#ifdef debug
		  G4cout<<"G4QString::CreatePartonPair: Created dQ-AdQ, PDG="<<PDGcode<<G4endl;
#endif 
    return G4QPartonPair(new G4QParton(-PDGcode), new G4QParton(PDGcode));
  }
  else
  {
    // Create a Quark - AntiQuark pair, first in pair is a Particle
    G4int PDGcode=SampleQuarkFlavor()*NeedParticle;
#ifdef debug
  		G4cout<<"G4QString::CreatePartonPair: Created Q-aQ, PDG="<<PDGcode<<G4endl;
#endif 
    return G4QPartonPair(new G4QParton(PDGcode), new G4QParton(-PDGcode));
  }
} // End of CreatePartonPair

// Creation of the Meson out of two partons (q, anti-q) 
G4QHadron* G4QString::CreateMeson(G4QParton* black, G4QParton* white, Spin theSpin)
{
  static G4double scalarMesonMixings[6]={0.5, 0.25, 0.5, 0.25, 1.0, 0.5}; 
  static G4double vectorMesonMixings[6]={0.5, 0.0, 0.5, 0.0, 1.0, 1.0};
  G4int id1= black->GetPDGCode();
  G4int id2= white->GetPDGCode();
#ifdef debug
  G4cout<<"G4QString::CreateMeson: bT="<<black->GetType()<<"("<<id1<<"), wT="
        <<white->GetType()<<"("<<id2<<")"<<G4endl;
#endif 
  if (std::abs(id1) < std::abs(id2))     // exchange black and white
  {
    G4int xchg = id1; 
    id1 = id2;  
    id2 = xchg;
  }
  if(std::abs(id1)>3) 
  {
    G4cerr<<"***G4QString::CreateMeson: q1="<<id1<<", q2="<<id2
          <<" while CHIPS is only SU(3)"<<G4endl;
    G4Exception("G4QString::CreateMeson:","72",FatalException,"HeavyQuarkFound");
  } 
  G4int PDGEncoding=0;
  if(!(id1+id2))                      // annihilation case (neutral)
  {     
    G4double rmix = G4UniformRand();
    G4int    imix = 2*std::abs(id1) - 1;
    if(theSpin == SpinZero)
       PDGEncoding = 110*(1 + G4int(rmix + scalarMesonMixings[imix - 1])
                            + G4int(rmix + scalarMesonMixings[imix]    ) ) +  theSpin;
    else
       PDGEncoding = 110*(1 + G4int(rmix + vectorMesonMixings[imix - 1])
                            + G4int(rmix + vectorMesonMixings[imix]    ) ) +  theSpin;
  }
  else
  {
    PDGEncoding = 100 * std::abs(id1) + 10 * std::abs(id2) +  theSpin;  
    G4bool IsUp = (std::abs(id1)&1) == 0; // quark 1 is up type quark (u or c?)
    G4bool IsAnti = id1 < 0;              // quark 1 is an antiquark?
    if( (IsUp && IsAnti) || (!IsUp && !IsAnti) )  PDGEncoding = - PDGEncoding;
    // Correction for the true neutral mesons
    if( PDGEncoding == -111 || PDGEncoding == -113 || PDGEncoding == -223 || 
        PDGEncoding == -221 || PDGEncoding == -331 || PDGEncoding == -333 )
                                                               PDGEncoding = - PDGEncoding;
  }
  G4QHadron* Meson= new G4QHadron(PDGEncoding);
#ifdef debug
  G4cout<<"G4QString::CreateBaryon: Meson is created with PDG="<<PDGEncoding<<G4endl;
#endif
  //delete black;                          // It is better to delete here and consider  
  //delete white;                          // the hadron creation as a delete equivalent
  return Meson;
}

// Creation of the Baryon out of two partons (q, di-q), (anti-q, anti-di-q) 
G4QHadron* G4QString::CreateBaryon(G4QParton* black, G4QParton* white, Spin theSpin)
{
  G4int id1= black->GetPDGCode();
  G4int id2= white->GetPDGCode();
#ifdef debug
  G4cout<<"G4QString::CreateBaryon: bT="<<black->GetType()<<"("<<id1<<"), wT="
        <<white->GetType()<<"("<<id2<<")"<<G4endl;
#endif 
  if(std::abs(id1) < std::abs(id2))
  {
    G4int xchg = id1; 
    id1 = id2;  
    id2 = xchg;
  }
  if(std::abs(id1)<1000 || std::abs(id2)> 3) 
  {
    G4cerr<<"***G4QString::CreateBaryon: q1="<<id1<<", q2="<<id2
          <<" can't create a Baryon"<<G4endl;
    G4Exception("G4QString::CreateBaryon:","72",FatalException,"WrongQdQSequence");
  }
  G4int ifl1= std::abs(id1)/1000;
  G4int ifl2 = (std::abs(id1) - ifl1 * 1000)/100;
  G4int diquarkSpin = std::abs(id1)%10; 
  G4int ifl3 = id2;
  if (id1 < 0) {ifl1 = - ifl1; ifl2 = - ifl2;}
  //... Construct baryon, distinguish Lambda and Sigma baryons.
  G4int kfla = std::abs(ifl1);
  G4int kflb = std::abs(ifl2);
  G4int kflc = std::abs(ifl3);
  G4int kfld = std::max(kfla,kflb);
        kfld = std::max(kfld,kflc);
  G4int kflf = std::min(kfla,kflb);
        kflf = std::min(kflf,kflc);
  G4int kfle = kfla + kflb + kflc - kfld - kflf;
  //... baryon with content uuu or ddd or sss has always spin = 3/2
  if(kfla==kflb && kflb==kflc) theSpin=SpinThreeHalf;   

  G4int kfll = 0;
  if(theSpin == SpinHalf && kfld > kfle && kfle > kflf)
  { 
    // Spin J=1/2 and all three quarks different
    // Two states exist: (uds -> lambda or sigma0)
    //   -  lambda: s(ud)0 s : 3122; ie. reverse the two lighter quarks
    //   -  sigma0: s(ud)1 s : 3212
    if(diquarkSpin == 1 )
    {
       if ( kfla == kfld) kfll = 1; // heaviest quark in diquark
       else kfll = G4int(0.25 + G4UniformRand());
    }   
    if(diquarkSpin==3 && kfla!=kfld) kfll = G4int(0.75+G4UniformRand());
  }
  G4int PDGEncoding=0;
  if (kfll == 1) PDGEncoding = 1000 * kfld + 100 * kflf + 10 * kfle + theSpin;
  else           PDGEncoding = 1000 * kfld + 100 * kfle + 10 * kflf + theSpin;
  if (id1 < 0) PDGEncoding = -PDGEncoding;
  G4QHadron* Baryon= new G4QHadron(PDGEncoding);
#ifdef debug
  G4cout<<"G4QString::CreateBaryon: Baryon is created with PDG="<<PDGEncoding<<G4endl;
#endif 
  //delete black;                          // It is better to delete here and consider  
  //delete white;                          // the hadron creation as a delete equivalent
  return Baryon;
} // End of CreateBaryon

G4QHadron* G4QString::CreateHadron(G4QParton* black, G4QParton* white)
{
  //static G4double mesonLowSpin = 0.25;      // probability to create scalar meson (2s+1) 
  //static G4double baryonLowSpin= 1./3.;     // probability to create 1/2 baryon  (2s+1)
  static G4double mesonLowSpin = 0.5;      // probability to create scalar meson (spFlip)
  static G4double baryonLowSpin= 0.5;      // probability to create 1/2 baryon (spinFlip)
  G4int bT=black->GetType();
  G4int wT=white->GetType();
#ifdef debug
  G4cout<<"G4QString::CreateHadron: bT="<<bT<<"("<<black->GetPDGCode()<<"), wT="<<wT<<"("
        <<white->GetPDGCode()<<")"<<G4endl;
#endif 
  if(bT==2 || wT==2)
  {
    // Baryon consists of quark and at least one di-quark
    Spin spin = (G4UniformRand() < baryonLowSpin) ? SpinHalf : SpinThreeHalf;
#ifdef debug
    G4cout<<"G4QString::CreateHadron: ----> Baryon is under creation"<<G4endl;
#endif 
    return CreateBaryon(black, white, spin);
  }
  else
  { 
    // Meson consists of quark and abnti-quark
    Spin spin = (G4UniformRand() < mesonLowSpin) ? SpinZero : SpinOne;
#ifdef debug
    G4cout<<"G4QString::CreateHadron: ----> Meson is under creation"<<G4endl;
#endif 
    return CreateMeson(black, white, spin);
  }
} // End of Create Hadron

// Creation of only High Spin (2,3/2) hadrons
G4QHadron* G4QString::CreateLowSpinHadron(G4QParton* black, G4QParton* white)
{
  G4int bT=black->GetType();
  G4int wT=white->GetType();
#ifdef debug
  G4cout<<"G4QString::CreateLowSpinHadron: ***Called***, bT="<<bT<<"("<<black->GetPDGCode()
        <<"), wT="<<wT<<"("<<white->GetPDGCode()<<")"<<G4endl;
#endif 
  if(bT == 1 && wT == 1)
  {
#ifdef debug
    G4cout<<"G4QString::CreateLowSpinHadron: ----> Meson is under creation"<<G4endl;
#endif 
    return CreateMeson(black, white, SpinZero);
  }
  else                  // returns a SpinThreeHalf Baryon if all quarks are the same
  {
#ifdef debug
    G4cout<<"G4QString::CreateLowSpinHadron: ----> Baryon is under creation"<<G4endl;
#endif 
    return CreateBaryon(black, white, SpinHalf);
  }
} // End of CreateLowSpinHadron

// Creation of only High Spin (2,3/2) hadrons
G4QHadron* G4QString::CreateHighSpinHadron(G4QParton* black, G4QParton* white)
{
  G4int bT=black->GetType();
  G4int wT=white->GetType();
#ifdef debug
  G4cout<<"G4QString::CreateHighSpinHadron:***Called***, bT="<<bT<<"("<<black->GetPDGCode()
        <<"), wT="<<wT<<"("<<white->GetPDGCode()<<")"<<G4endl;
#endif 
  if(bT == 1 && wT == 1)
  {
#ifdef debug
    G4cout<<"G4QString::CreateHighSpinHadron: ----> Meson is created"<<G4endl;
#endif 
    return CreateMeson(black,white, SpinOne);
  }
  else
  {
#ifdef debug
    G4cout<<"G4QString::CreateHighSpinHadron: ----> Baryon is created"<<G4endl;
#endif 
    return CreateBaryon(black,white,SpinThreeHalf);
  }
} // End of CreateHighSpinHadron
