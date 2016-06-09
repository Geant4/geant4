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
//      ---------------- G4QPDGCode ----------------
//             by Mikhail Kossov, Sept 1999.
//      class for Hadron definitions in CHIPS Model
// -------------------------------------------------------------------
// Short description: The PDG Code is made on the basis of the Quark
// Content (G4QuarkContent) of the hadronic state (including nuclear
// fragments). The PDG code of the ground state (e.g. pi, N, etc.) is
// calculated. It includes a complicated algortithm of the G.S. mass
// calculation for nuclear fragments (now it is synchronised with the
// G4 nuclear massess).
// -------------------------------------------------------------------

//#define debug
//#define pdebug
//#define qdebug
//#define idebug
//#define sdebug

#include <cmath>
#include <cstdlib>

#include "G4QPDGCodeVector.hh"

using namespace std;

G4QPDGCode::G4QPDGCode(G4int PDGCode): thePDGCode(PDGCode)
{
#ifdef sdebug
  G4cout<<"G4QPDGCode:Constructer is called with PDGCode="<<PDGCode<<G4endl;  
#endif
  if(PDGCode==130) PDGCode= 311; // Safety. Should not happen.
  if(PDGCode==310) PDGCode=-311; // Safety. Should not happen.
  if(PDGCode==90000000)
  {
    thePDGCode=22;
    theQCode=6;
  }
  else if(PDGCode) theQCode=MakeQCode(PDGCode);
  else        
  {
#ifdef sdebug
    G4cout<<"***G4QPDGCode: Constructed with PDGCode=0, QCode=-2"<<G4endl;  
#endif
    theQCode=-2;
  }
#ifdef debug
  if(PDGCode==3222)G4cout<<"G4QPDGCd:Con(PDG) PDG="<<PDGCode<<", QCode="<<theQCode<<G4endl;
#endif
}

G4QPDGCode::G4QPDGCode(G4bool f, G4int QCode): theQCode(QCode)
{
  if(f&&QCode<0)G4cerr<<"***G4QPDGCode::Constr. QCode="<<QCode<<G4endl;
  thePDGCode = MakePDGCode(QCode);
#ifdef debug
  G4cout<<"G4QPDGCode::Constr: PDGCode="<<thePDGCode<<G4endl;
#endif
}

G4QPDGCode::G4QPDGCode(G4QContent QCont)
{
  InitByQCont(QCont);
}

G4QPDGCode::G4QPDGCode(const G4QPDGCode& rhs)
{
  thePDGCode =rhs.thePDGCode;
  theQCode   =rhs.theQCode;
}

G4QPDGCode::G4QPDGCode(G4QPDGCode* rhs)
{
  thePDGCode =rhs->thePDGCode;
  theQCode   =rhs->theQCode;
}

const G4QPDGCode& G4QPDGCode::operator=(const G4QPDGCode& rhs)
{
  if(this != &rhs)                          // Beware of self assignment
  {
    thePDGCode =rhs.thePDGCode;
    theQCode   =rhs.theQCode;
  }
  return *this;
}

G4QPDGCode::~G4QPDGCode() {}

// Standard output for QPDGCode
ostream& operator<<(ostream& lhs, G4QPDGCode& rhs)
{
  lhs << "[ PDG=" << rhs.GetPDGCode() << ", Q=" << rhs.GetQCode() << "]";
  return lhs;
}

// Standard output for const QPDGCode
ostream& operator<<(ostream& lhs, const G4QPDGCode& rhs)
{
  lhs << "[ PDG=" << rhs.GetPDGCode() << ", Q=" << rhs.GetQCode() << "]";
  return lhs;
}

// Overloading of QPDGCode addition
G4int operator+(const G4QPDGCode& lhs, const G4QPDGCode& rhs)
{
  G4int  s_value  = lhs.GetPDGCode();
  return s_value += rhs.GetPDGCode();
}
G4int operator+(const G4QPDGCode& lhs, const G4int& rhs)
{
  G4int  s_value  = lhs.GetPDGCode();
  return s_value += rhs;
}
G4int operator+(const G4int& lhs, const G4QPDGCode& rhs)
{
  G4int  s_value  = lhs;
  return s_value += rhs.GetPDGCode();
}

// Overloading of QPDGCode subtraction
G4int operator-(const G4QPDGCode& lhs, const G4QPDGCode& rhs)
{
  G4int  s_value  = lhs.GetPDGCode();
  return s_value -= rhs.GetPDGCode();
}
G4int operator-(const G4QPDGCode& lhs, const G4int& rhs)
{
  G4int  s_value  = lhs.GetPDGCode();
  return s_value -= rhs;
}
G4int operator-(const G4int& lhs, const G4QPDGCode& rhs)
{
  G4int  s_value  = lhs;
  return s_value -= rhs.GetPDGCode();
}

// Overloading of QPDGCode multiplication
G4int operator*(const G4QPDGCode& lhs, const G4QPDGCode& rhs)
{
  G4int  s_value  = lhs.GetPDGCode();
  return s_value *= rhs.GetPDGCode();
}

G4int operator*(const G4QPDGCode& lhs, const G4int& rhs)
{
  G4int  s_value  = lhs.GetPDGCode();
  return s_value *= rhs;
}

G4int operator*(const G4int& lhs, const G4QPDGCode& rhs)
{
  G4int  s_value  = lhs;
  return s_value *= rhs.GetPDGCode();
}

// Overloading of QPDGCode division
G4int operator/(const G4QPDGCode& lhs, const G4QPDGCode& rhs)
{
  G4int  s_value  = lhs.GetPDGCode();
  return s_value /= rhs.GetPDGCode();
}

G4int operator/(const G4QPDGCode& lhs, const G4int& rhs)
{
  G4int  s_value  = lhs.GetPDGCode();
  return s_value /= rhs;
}

G4int operator/(const G4int& lhs, const G4QPDGCode& rhs)
{
  G4int  s_value  = lhs;
  return s_value /= rhs.GetPDGCode();
}

// Overloading of QPDGCode residual
G4int operator%(const G4QPDGCode& lhs, const G4int& rhs)
{
  G4int  s_value  = lhs.GetPDGCode();
  return s_value %= rhs;
}

// TRUE if it is not RealNeutral (111,221,331 etc), FALSE if it is.
G4bool G4QPDGCode::TestRealNeutral(const G4int& PDGCode)
{
  if(PDGCode>0 && PDGCode<999)    // RealNeutral are always positive && mesons
  {
    if(PDGCode==22) return false; // Photon
    G4int p=PDGCode/10;
    if(p/10==p%10) return false; // This is a RealNeutral
  }
  return true;
}

// Make a Q Code out of the PDG Code
G4int G4QPDGCode::MakePDGCode(const G4int& QCode)
{
  //static const G4int modi = 81;  // Q Codes for more than di-baryon nuclei
  //static const G4int modi = 89;  // Q Codes for more than di-baryon nuclei "IsoNuclei"
  //static const G4int modi = 122; // Q Codes: more than quarta-baryon nuclei "Lept/Hyper"
  static const G4int modi = 85;    // Reduced Q Codes: > quarta-baryon nuclei "Lept/Hyper"
  //static G4int qC[modi]={ 11,   12,   13,   14,   15,   16,   22,   23,   24,   25, // 10
  //                        37,  110,  220,  330,  111,  211,  221,  311,  321,  331, // 20
  //                      2112, 2212, 3122, 3112, 3212, 3222, 3312, 3322,  113,  213, // 30
  //                       223,  313,  323,  333, 1114, 2114, 2214, 2224, 3124, 3114, // 40
  //                      3214, 3224, 3314, 3324, 3334,  115,  215,  225,  315,  325, // 50
  //                       335, 2116, 2216, 3126, 3116, 3216, 3226, 3316, 3326,  117, // 60
  //                       217,  227,  317,  327,  337, 1118, 2118, 2218, 2228, 3128, // 70
  //                      3118, 3218, 3228, 3318, 3328, 3338,  119,  219,  229,  319, // 80
  //                       329,  339, 90002999 , 89999003 , 90003998 ,
  //                       89998004 , 90003999 , 89999004 , 90004998 , 89998005 ,     // 90
  //                       90000001 , 90001000 , 91000000 , 90999001 , 91000999 ,
  //                       91999000 , 91999999 , 92999000 , 90000002 , 90001001 ,     //100
  //                       90002000 , 91000001 , 91001000 , 92000000 , 90999002 ,
  //                       91001999 , 90001002 , 90002001 , 91000002 , 91001001 ,     //110
  //                       91002000 , 92000001 , 92001000 , 90999003 , 90001003 ,
  //                       90002002 , 90003001 , 91001002 , 91002001 , 92000002 ,     //120
  //                       92001001 , 92002000};                                      //122
  static G4int qC[modi] ={  11,   12,   13,   14,   15,   16,   22,   23,   24,   25, // 10
                            37,  110,  220,  330,  111,  211,  221,  311,  321,  331, // 20
                          2112, 2212, 3122, 3112, 3212, 3222, 3312, 3322,  113,  213, // 30
                           223,  313,  323,  333, 1114, 2114, 2214, 2224, 3124, 3114, // 40
				  3214, 3224, 3314, 3324, 3334,                               // 45
				    90002999 , 89999003 , 90003998 , 89998004 , 90003999 ,    // 50
				    89999004 , 90004998 , 89998005 , 90000001 , 90001000 ,    // 55
				    91000000 , 90999001 , 91000999 , 91999000 , 91999999 ,    // 60
				    92999000 , 90000002 , 90001001 , 90002000 , 91000001 ,    // 65
				    91001000 , 92000000 , 90999002 , 91001999 , 90001002 ,    // 70
				    90002001 , 91000002 , 91001001 , 91002000 , 92000001 ,    // 75
				    92001000 , 90999003 , 90001003 , 90002002 , 90003001 ,    // 80
				    91001002 , 91002001 , 92000002 , 92001001 , 92002000};    // 85
  static G4int aC[15] = {1,1000,999001,1000000,1000999,1999000,1999999,        // sum 1
                         2,1001,2000,1000001,1001000,1999001,2000000,2000999}; // sum 2
  if      (QCode<0)
  {
    G4cerr<<"***G4QPDGCode::MakePDGCode: negative Q Code ="<<QCode<<G4endl;
    return 0;
  }
  else if (QCode>=modi)
  {
    //G4int q=QCode-modi;              // Starting BarNum=3
    //G4int a=q/15+1;                  // BarNum/2
    //G4int b=q%15;
    G4int q=QCode-modi;              // Starting BarNum=5
    G4int a=q/15+2;                  // BarNum/2
    G4int b=q%15;
    return 90000000+a*1001+aC[b];
  }
  return qC[theQCode];
}

// Hadronic masses synhronized with the Geant4 hadronic masses
G4double G4QPDGCode:: QHaM(G4int nQ)
{
  static G4bool iniFlag=true;
  //static G4double mass[nQHM]={.511,0.,105.65837, 0., 1777., 0.,   0., 91188., 80403., 140.00
  //,120.000,    800.,     980.,    1370.,  134.98,  139.57, 547.51, 497.65, 493.68, 957.78
  //,939.5654,938.272, 1115.683,  1197.45, 1192.64, 1189.37,1321.31,1314.83,  775.5,  775.5
  //, 782.65,   896.0,   891.66,  1019.46,   1232.,   1232.,  1232.,  1232., 1519.5, 1387.2
  //, 1383.7,  1382.8,    1535.,   1531.8, 1672.45,  1318.3, 1318.3, 1275.4, 1432.4, 1425.6
  //,  1525.,   1680.,    1680.,    1820.,   1915.,   1915.,  1915.,  2025.,  2025.,  1691.
  //,  1691.,   1667.,    1776.,    1776.,   1854.,   1950.,  1950.,  1950.,  1950.,  2100.
  //,  2030.,   2030.,    2030.,    2127.,   2127.,   2252.,  2020.,  2020.,  2044.,  2045.
  //, 2045., 2297., 2170.272, 2171.565, 2464., 2464., 3108.544, 3111.13,3402.272,3403.565};
  // Reduced:
  static G4double mass[nQHM]={.511, 0., 105.65837, 0., 1777., 0.,   0., 91188., 80403., 140.00
    ,120.000,    800.,     980.,    1370.,  134.98,  139.57, 547.51, 497.65, 493.68, 957.78
    ,939.5654,938.272, 1115.683,  1197.45, 1192.64, 1189.37,1321.31,1314.83,  775.5,  775.5
    , 782.65,   896.0,   891.66,  1019.46,   1232.,   1232.,  1232.,  1232., 1519.5, 1387.2
    , 1383.7,  1382.8,    1535.,   1531.8, 1672.45,2170.272,2171.565, 2464., 2464.,3108.544
    ,3111.13,3402.272, 3403.565};
  if(iniFlag) // Initialization of the Geant4 hadronic masses
  {
    mass[ 0]=      G4Electron::Electron()->GetPDGMass();
    mass[ 1]=    G4NeutrinoE::NeutrinoE()->GetPDGMass();
    mass[ 2]=    G4MuonMinus::MuonMinus()->GetPDGMass();
    mass[ 3]=  G4NeutrinoMu::NeutrinoMu()->GetPDGMass();
    mass[ 4]=      G4TauMinus::TauMinus()->GetPDGMass();
    mass[ 5]=G4NeutrinoTau::NeutrinoTau()->GetPDGMass();
    mass[14]=      G4PionZero::PionZero()->GetPDGMass();
    mass[15]=    G4PionMinus::PionMinus()->GetPDGMass();
    mass[16]=                G4Eta::Eta()->GetPDGMass();
    mass[17]=      G4KaonZero::KaonZero()->GetPDGMass();
    mass[18]=    G4KaonMinus::KaonMinus()->GetPDGMass();
    mass[19]=      G4EtaPrime::EtaPrime()->GetPDGMass();
    mass[20]=        G4Neutron::Neutron()->GetPDGMass();
    mass[21]=          G4Proton::Proton()->GetPDGMass();
    mass[22]=          G4Lambda::Lambda()->GetPDGMass();
    mass[23]=  G4SigmaMinus::SigmaMinus()->GetPDGMass();
    mass[24]=    G4SigmaZero::SigmaZero()->GetPDGMass();
    mass[25]=    G4SigmaPlus::SigmaPlus()->GetPDGMass();
    mass[26]=        G4XiMinus::XiMinus()->GetPDGMass();
    mass[27]=          G4XiZero::XiZero()->GetPDGMass();
    mass[44]=  G4OmegaMinus::OmegaMinus()->GetPDGMass();
    iniFlag=false;
  }
  if(nQ<0 || nQ>=nQHM)
  {
    G4cout<<"***G4QPDGCode::QHaM: negative Q-code or Q="<<nQ<<" >= nQmax = "<<nQHM<<G4endl;
    return 0.;
  }
  return mass[nQ];
}

// Make a Q Code out of the PDG Code
G4int G4QPDGCode::MakeQCode(const G4int& PDGCode)
{
  static const G4int qr[10]={0,13,19,27,33,44,50,58,64,75};
  G4int PDGC=abs(PDGCode);        // Qcode is always not negative
  G4int s_value=0;
  G4int z=0;
  G4int n=0;
  if (PDGC>100000000)             // Not supported
  {
#ifdef debug
    G4cout<<"***G4QPDGCode::MakeQCode: Unknown in Q-System code: "<<PDGCode<<G4endl;
#endif
    return -2;
  }
  else if (PDGC>80000000 && PDGC<100000000) // Try to convert the NUCCoding to PDGCoding
  {
    //if(PDGC==90000000) return 6;            // @@ already done in the constructor
    ConvertPDGToZNS(PDGC, z, n, s_value);
    G4int b=n+z+s_value;                           // Baryon number
#ifdef debug
    G4cout<<"***G4QPDGCode::Z="<<z<<",N="<<n<<",S="<<s_value<<G4endl;
#endif
    if(b<0)                                        // ---> Baryons & Fragments
    {
      b=-b;
      n=-n;
      z=-z;
      s_value=-s_value;
      PDGC=90000000+s_value*1000000+z*1000+n;      // New PDGC for anti-baryons
    }
    else if(!b)                                    // --> Mesons
    {
      //G4bool anti=false;                           // For the PDG conversion
      if(z<0)                                      // --> Mesons conversion
      {
        n=-n;
        z=-z;
        s_value=-s_value;
        //anti=true;                                 // For the PDG conversion
      }
      if(!z)
      {
        if(s_value>0)
        {
          n=-n;
          s_value=-s_value;
          //anti=true;                               // For the PDG conversion
        }
        if     (s_value==-1) return 17;            // K0
        else if(s_value==-2) return -1;            // K0+K0 chipolino
        else           return -2;                  // Not supported by Q Code
      }
      else                                         // --> z>0
      {
        if(z==1)
        {
          if   (s_value==-1) return 18;            // K+
          else         return 15;                  // pi+
        }
        else if(z==2)  return -1;                  // Chipolino
        else           return -2;                  // Not supported by Q Code
      }
    } // End of meson case
    if(b>0)                                        // --> Baryoniums case
    {
      if(b==1)                                     // --> Baryons+Hyperons
      {
        if(PDGC>80000000)
        {
          if(!s_value)                               // --> Baryons
          {
            if     (!z)  return 53;                  // neutron
            else if(z==1)return 54;                  // proton
            else         return -2;                  // Not supported by Q Code
          }
          else if(s_value==1)                        // --> Hyperons
          {
            if(z==-1)    return 56;                  // Sigma-
            else if(!z)  return 55;                  // Lambda
            else if(z==1)return 57;                  // Sigma+
            else         return -2;                  // Not supported by Q Code
          }
          else if(s_value==2)                        // --> Xi Hyperons
          {
            if(z==-1)    return 58;                  // Xi-
            else if(!z)  return 59;                  // Xi0
            else         return -2;                  // Not supported by Q Code
          }
          else if(s_value==3)                        // --> Xi Hyperons
          {
            if(z==-1)    return 60;                  // Omega-
            else         return -2;                  // Not supported by Q Code
          }
        }
        else
        {
          if(!s_value)                               // --> Baryons
          {
            if(z==-1)    return 34;                  // Delta-
            else if(!z)  return 20;                  // neutron
            else if(z==1)return 21;                  // proton
            else if(z==2)return 37;                  // Delta++
            else if(z==3||z==-2)return -1;           // Delta+pi Chipolino
            else         return -2;                  // Not supported by Q Code
          }
          else if(s_value==1)                        // --> Hyperons
          {
            if(z==-1)    return 23;                  // Sigma-
            else if(!z)  return 22;                  // Lambda (@@ 24->Sigma0)
            else if(z==1)return 25;                  // Sigma+
            else if(z==2||z==-2) return -1;          // Sigma+pi Chipolino
            else         return -2;                  // Not supported by Q Code
          }
          else if(s_value==2)                        // --> Xi Hyperons
          {
            if(z==-1)    return 26;                  // Xi-
            else if(!z)  return 27;                  // Xi0
            else if(z==1||z==-2)return -1;           // Xi+pi Chipolino
            else         return -2;                  // Not supported by Q Code
          }
          else if(s_value==3)                        // --> Xi Hyperons
          {
            if(z==-1)    return 44;                  // Omega-
            else if(!z||z==-2)  return -1;           // Omega+pi Chipolino
            else         return -2;                  // Not supported by Q Code
          }
        }
      }
      else
      {
        if(b==2)
        {
          if     (PDGC==90002999) return 45;       // p DEL++
          else if(PDGC==89999003) return 46;       // n DEL-
          else if(PDGC==90003998) return 47;       // DEL++ DEL++
          else if(PDGC==89998004) return 48;       // DEL-  DEL-
          else if(PDGC==90999002) return 67;      // n Sigma-
          else if(PDGC==91001999) return 68;      // p Sigma+
        }
        if(b==3)
        {
          if     (PDGC==90003999) return 49;       // p p DEL++
          else if(PDGC==89999004) return 50;       // n n DEL-
          else if(PDGC==90004998) return 51;       // p DEL++ DEL++
          else if(PDGC==89998005) return 52;       // n DEL-  DEL-
          else if(PDGC==90999003) return 76;      // n n Sigma-
        }
      }
    }
  }
  if (PDGC<80000000)                // ----> Direct Baryons & Mesons
  {
    if     (PDGC<100)               // => Leptons and field bosons
    {
      if     (PDGC==10)  return -1; // Chipolino
      else if(PDGC==11)  return  0; // e-
      else if(PDGC==12)  return  1; // nu_e
      else if(PDGC==13)  return  2; // mu-
      else if(PDGC==14)  return  3; // nu_mu
      else if(PDGC==15)  return  4; // tau-
      else if(PDGC==16)  return  5; // nu_tau
      else if(PDGC==22)  return  6; // Photon
      else if(PDGC==23)  return  7; // Z0 boson
      else if(PDGC==24)  return  8; // W- boson
      else if(PDGC==25)  return  9; // H0 (neutral Higs boson)
      else if(PDGC==37)  return 10; // H- (charged Higs boson)
    }
    G4int r=PDGC%10;                // 2s+1
    G4int         Q= 0;
    if     (!r)
    {
      // Internal CHIPS codes for the wide f_0 states must be 9000221, 9010221, 10221
      if     (PDGC==110) return 11; // Low R-P: Sigma (pi,pi S-wave)
      else if(PDGC==220) return 12; // Midle Regeon-Pomeron
      else if(PDGC==330) return 13; // High Regeon-Pomeron
#ifdef debug
      G4cout<<"***G4QPDGCode::MakeQCode: (0) Unknown in Q-System code: "<<PDGCode<<G4endl;
#endif
      return -2;
    }
    else Q=qr[r];
    G4int p=PDGC/10;                // Quark Content
    if(r%2)                         // (2s+1 is odd) Mesons
    {
      if     (p==11) return Q+=1;
      else if(p==21) return Q+=2;
      else if(p==22) return Q+=3;
      else if(p==31) return Q+=4;
      else if(p==32) return Q+=5;
      else if(p==33) return Q+=6;
      else
      {
#ifdef debug
        G4cout<<"*Warning*G4QPDGCode::MakeQCode:(1)UnknownQCode for PDG="<<PDGCode<<G4endl;
#endif
        return -2;
      }
    }
    else                    // (2s+1 is even) Baryons
    {
      s_value=r/2;
      if(s_value%2)         // ((2s+1)/2 is odd) N Family
      {
        if     (p==211) return Q+=1;
        else if(p==221) return Q+=2;
        else if(p==312) return Q+=3;
        else if(p==311) return Q+=4;
        else if(p==321) return Q+=5;
        else if(p==322) return Q+=6;
        else if(p==331) return Q+=7;
        else if(p==332) return Q+=8;
        else
        {
#ifdef debug
          G4cout<<"*Warning*G4QPDGCode::MakeQCode:(2) UnknownQCode, PDG="<<PDGCode<<G4endl;
#endif
          return -2;
        }
      }
      else                  // ((2s+1)/2 is odd) Delta Family
      {
        if     (p==111) return Q+= 1;
        else if(p==211) return Q+= 2;
        else if(p==221) return Q+= 3;
        else if(p==222) return Q+= 4;
        else if(p==312) return Q+= 5;
        else if(p==311) return Q+= 6;
        else if(p==321) return Q+= 7;
        else if(p==322) return Q+= 8;
        else if(p==331) return Q+= 9;
        else if(p==332) return Q+=10;
        else if(p==333) return Q+=11;
        else
        {
#ifdef debug
          G4cout<<"**G4QPDGCode::MakeQCode:(3) Unknown in Q-System code:"<<PDGCode<<G4endl;
#endif
          return -2;
        }
      }
    }
  }
  else                        // Nuclear Fragments
  {
    G4int d=n+n+z+s_value;    // a#of d quarks
    G4int u=n+z+z+s_value;    // a#of u quarks
    G4int t=d+u+s_value;      // tot#of quarks
    if(t%3 || abs(t)<3)       // b=0 are in mesons
    {
#ifdef debug
      G4cout<<"***G4QPDGCode::MakeQCode: Unknown PDGCode="<<PDGCode<<", t="<<t<<G4endl;
#endif
      return -2;
    }
    else
    {
      G4int b=t/3;            // baryon number
      if(b==1)                // baryons
      {
        if     (s_value==0&&u==1&&d==2) return 53; // n
        else if(s_value==0&&u==2&&d==1) return 54; // p
        else if(s_value==1&&u==1&&d==1) return 55; // Lambda
        else if(s_value==1&&u==0&&d==2) return 56; // Sigma-
        else if(s_value==1&&u==2&&d==0) return 57; // Sigma+
        else if(s_value==2&&u==0&&d==1) return 58; // Xi-
        else if(s_value==2&&u==1&&d==0) return 59; // Xi0
        else if(s_value==3&&u==0&&d==0) return 60; // Omega-
        else
        {
#ifdef debug
          G4cout<<"**G4QPDGCode::MakeQCode:(5) Unknown in Q-System code:"<<PDGCode<<G4endl;
#endif
          return -2;
        }
      }
      else if(b==2)           // di-baryons
      {
        if     (s_value==0&&u==2&&d==4) return 61; // nn
        else if(s_value==0&&u==3&&d==3) return 62; // np
        else if(s_value==0&&u==4&&d==2) return 63; // pp
        else if(s_value==1&&u==2&&d==3) return 64; // nLambda
        else if(s_value==1&&u==3&&d==2) return 65; // pLambda
        else if(s_value==2&&u==2&&d==2) return 66; // LambdaLambda
        else
        {
#ifdef debug
          G4cout<<"**G4QPDGCode::MakeQCode:(6) Unknown in Q-System code:"<<PDGCode<<G4endl;
#endif
          return -2;
        }
      }
      else if(b==3)           // tri-baryons
      {
        if     (s_value==0&&u==4&&d==5) return 69; // pnn
        else if(s_value==0&&u==5&&d==4) return 70; // npp
        else if(s_value==1&&u==3&&d==5) return 71; // Lnn
        else if(s_value==1&&u==4&&d==4) return 72; // Lnp
        else if(s_value==1&&u==5&&d==3) return 73; // Lpp
        else if(s_value==2&&u==3&&d==4) return 74; // LLn
        else if(s_value==2&&u==4&&d==3) return 75; // LLp
        else if(s_value==1&&u==2&&d==6) return 76; // nnSigma-
        else
        {
#ifdef debug
          G4cout<<"**G4QPDGCode::MakeQCode:(7) Unknown in Q-System code:"<<PDGCode<<G4endl;
#endif
          return -2;
        }
      }
      G4int c=b/2;              // From here b>3: (4,5):c=2,g=0,1; (6,7):c=3,g=0,1; ...
      G4int g_value=b%2;
      G4int h=c*3;
      //G4int Q=57+c*15;
      //G4int Q=65+c*15;           // "IsoNuclei"
      //G4int Q=83+c*15;           // "Leptons/Hyperons"
      G4int Q=46+c*15;           // Reduced "Leptons/Hyperons"
      u-=h;
      d-=h;
      if(g_value)
      {
        if     (s_value==0&&u==1&&d==2) return Q+= 9;
        else if(s_value==0&&u==2&&d==1) return Q+=10;
        else if(s_value==1&&u==0&&d==2) return Q+=11;
        else if(s_value==1&&u==1&&d==1) return Q+=12;
        else if(s_value==1&&u==2&&d==0) return Q+=13;
        else if(s_value==2&&u==0&&d==1) return Q+=14;
        else if(s_value==2&&u==1&&d==0) return Q+=15;
        else
        {
#ifdef debug
          G4cout<<"**G4QPDGCode::MakeQCode:(8) Unknown in Q-System code:"<<PDGCode<<G4endl;
#endif
          return -2;
        }
      }
      else
      {
        if     (s_value==0&&u==-1&&d== 1) return Q+=1;
        else if(s_value==0&&u== 0&&d== 0) return Q+=2;
        else if(s_value==0&&u== 1&&d==-1) return Q+=3;
        else if(s_value==1&&u==-1&&d== 0) return Q+=4;
        else if(s_value==1&&u== 0&&d==-1) return Q+=5;
        else if(s_value==2&&u==-2&&d== 0) return Q+=6;
        else if(s_value==2&&u==-1&&d==-1) return Q+=7;
        else if(s_value==2&&u== 0&&d==-2) return Q+=8;
        else
        {
#ifdef debug
          G4cout<<"**G4QPDGCode::MakeQCode:(9) Unknown in Q-System code:"<<PDGCode<<G4endl;
#endif
          return -2;
        }
      }
    }
  }
  G4cout<<"*Warning*G4QPDGCode::MakeQCode:() Unknown Q Code for PDG = "<<PDGCode<<G4endl;
  return -2; //not reachable statement (fake, only for compiler)  
}

// Get the mean mass value for the PDG
G4double G4QPDGCode::GetMass()
{
  G4int ab=theQCode;
#ifdef pdebug
  G4bool pPrint = thePDGCode == 3222 || ab == 25;
  if(pPrint)
  G4cout<<"G4QPDGCode::GetMass: Mass for Q="<<ab<<",PDG="<<thePDGCode<<",N="<<nQHM<<G4endl;
#endif
  if ( (ab < 0 && thePDGCode < 80000000) || !thePDGCode) {
#ifdef pdebug
    if(thePDGCode!=10 && pPrint)
      G4cout<<"**G4QPDGCode::GetMass:m=100000.,QC="<<theQCode<<",PDG="<<thePDGCode<<G4endl;
#endif
    return 100000.;
  }
  else if(ab>-1 && ab<nQHM)
  {
    G4double mass = QHaM(ab);
#ifdef pdebug
    //if(pPrint)
    if(thePDGCode == 3222 || ab == 25)
    G4cout<<"G4QPDGCode::GetMa:m="<<mass<<",Q="<<theQCode<<",PDG="<<thePDGCode<<G4endl;
#endif
    return mass;                                // Get mass from the table
  }
  //if(szn==0) return m[15];                    // @@ mPi0   @@ MK ?
  if(thePDGCode==90000000)
  {
#ifdef pdebug
    if(pPrint)
    G4cout<<"G4QPDGCode::GetMass:***m=0, QC="<<theQCode<<",PDG="<<thePDGCode<<G4endl;
#endif
    return 0.;
  }
  G4int s_value=0;
  G4int z=0;
  G4int n=0;
  ConvertPDGToZNS(thePDGCode, z, n, s_value);
  G4double m_value=GetNuclMass(z,n,s_value);
#ifdef pdebug
  if(pPrint)
  G4cout<<"G4QPDG::GetM:PDG="<<thePDGCode<<"=>Z="<<z<<",N="<<n<<",S="<<s_value<<",M="<<m_value<<G4endl;
#endif
  return m_value;
}

// Get the width value for the PDG
G4double G4QPDGCode::GetWidth()
{
  //static const int nW = 72;
  //static const int nW = 80; // "Isobars"
  //static G4double width[nQHM] = {0.,0.,0.,0.,0.,0.,0.,2.495,2.118,10.
  //  ,  10., 800.,  75., 350.,   0.,   0., .00118,  0.,   0., .203
  //  ,   0.,   0.,   0.,   0.,   0.,   0.,   0.,    0., 160., 160.
  //  , 8.41, 50.5, 50.8, 4.43, 120., 120., 120.,  120., 15.6,  39.
  //  ,  36., 35.8,  10.,   9.,   0., 107., 107., 185.5, 109., 98.5
  //  ,  76., 130., 130.,  80., 120., 120., 120.,   20.,  20., 160.
  //  , 160., 168., 159., 159.,  87., 300., 300.,  300., 300., 200.
  //  , 180., 180., 180.,  99.,  99.,  55., 387.,  387., 208., 198.
  //  , 198., 149., 120., 120., 170., 170., 120.,  120., 170., 170.};
  // Reduced:
  static G4double width[nQHM] = {0.,0.,0.,0.,0.,0.,0.,2.495,2.118,10.
    ,  10., 800.,  75., 350.,   0.,   0., .00118,  0.,   0., .203
    ,   0.,   0.,   0.,   0.,   0.,   0.,   0.,    0., 160., 160.
    , 8.41, 50.5, 50.8, 4.43, 120., 120., 120.,  120., 15.6,  39.
    ,  36., 35.8,  10.,   9.,   0., 120., 120.,  170., 170., 120.
    , 120., 170., 170.};
  G4int ab=abs(theQCode);
  if(ab<nQHM) return width[ab];
  return 0.;             // @@ May be real width should be implemented for nuclei e.g. pp
}

// Get a nuclear mass for Z (a#of protons), N (a#of neutrons), & S (a#of lambdas) 
G4double G4QPDGCode::GetNuclMass(G4int z, G4int n, G4int s_value)
{
  static const G4double anb = .01; // Antibinding for Ksi-n,Sig-n,Sig+p,Sig-nn,
  static const G4double mNeut= QHaM(20);
  static const G4double mProt= QHaM(21);
  static const G4double mLamb= QHaM(22);
  static const G4double mPiC = QHaM(15);
  static const G4double mKZ  = QHaM(17);
  static const G4double mKM  = QHaM(18);
  static const G4double mSiM = QHaM(23);
  static const G4double mSiP = QHaM(25);
  static const G4double mKsZ = QHaM(27);
  static const G4double mKsM = QHaM(26);
  static const G4double mOmM = QHaM(44);
  static const G4double mKZa = mKZ +anb;
  static const G4double mKMa = mKM +anb;
  static const G4double mSigM= mSiM+anb;
  static const G4double mSigP= mSiP+anb;
  static const G4double mKsiZ= mKsZ+anb;
  static const G4double mKsiM= mKsM+anb;
  static const G4double mOmeg= mOmM+anb;
  static const G4double mDiPi= mPiC+mPiC+anb;
  static const G4double mDiKZ= mKZa+mKZ;
  static const G4double mDiKM= mKMa+mKM;
  static const G4double mDiPr= mProt+mProt;
  static const G4double mDiNt= mNeut+mNeut;
  static const G4double mSmPi= mSiM+mDiPi;
  static const G4double mSpPi= mSiP+mDiPi;
  static const G4double mOmN = mOmeg+mNeut;
  static const G4double mSpP = mSigP+mProt;
  static const G4double mSpPP= mSpP +mProt;
  static const G4double mSmN = mSigM+mNeut;
  static const G4double mSmNN= mSmN +mNeut;
  // -------------- DAM Arrays ----------------------
  static const G4int iNR=76;    // Neutron maximum range for each Z
  static const G4int nEl = 105; // Maximum Z of the associative memory is "nEl-1=104"
  static const G4int iNF[nEl]={0,0,0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, // 14
                         1  ,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, // 29
                         16 , 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, // 44
                         31 , 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 53, 54, 55, // 59
                         56 , 56, 57, 57, 58, 60, 61, 63, 64, 65, 66, 67, 68, 69, 70, // 74
                         71 , 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, // 89
                         86 , 87, 88, 89, 91, 94, 98,103,109,115,122,128,134,140,146};//104
#ifdef qdebug
  static G4int iNmin[nEl]={0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, // 14
                         1  ,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, // 29
                         16 , 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, // 44
                         31 , 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 53, 54, 55, // 59
                         56 , 56, 57, 57, 58, 60, 61, 63, 64, 65, 66, 67, 68, 69, 70, // 74
                         71 , 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, // 89
                         86 , 87, 88, 89, 91, 94, 98,103,109,115,122,128,134,140,146};//104
  static G4int iNmax=iNR;
  static G4int iNran[nEl]={19,20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, // 14
                         34 , 35, 36, 37, 38, 39, 40, 48, 48, 48, 48, 50, 50, 50, 52, // 29
                         53 , 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, // 44
                         68 , 69, 70, 70, 70, 71, 71, 71, 71, 71, 72, 72, 72, 72, 72, // 59
                         73 , 73, 73, 73, 74, 74, 74, 74, 74, 74, 74, 74, 74, 75, 76, // 74
                         76 , 76, 76, 76, 76, 75, 74, 73, 72, 71, 70, 70, 69, 69, 69, // 89
                         68 , 68, 68, 67, 63, 59, 55, 51, 47, 43, 39, 35, 31, 27, 23};//104
#endif
  static const G4int iNL[nEl]={19,20,21,22,23,24, 25, 26, 27, 28, 29, 30, 31, 32, 33, // 14
                         34 , 35, 36, 37, 38, 39, 40, 48, 48, 48, 48, 50, 50, 50, 52, // 29
                         53 , 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, // 44
                         68 , 69, 70, 70, 70, 71, 71, 71, 71, 71, 72, 72, 72, 72, 72, // 59
                         73 , 73, 73, 73, 74, 74, 74, 74, 74, 74, 74, 74, 74, 75, 76, // 74
                         76 , 76, 76, 76, 76, 75, 74, 73, 72, 71, 70, 70, 69, 69, 69, // 89
                         68 , 68, 68, 67, 63, 59, 55, 51, 47, 43, 39, 35, 31, 27, 23};//104
   // ********* S=-4 vectors *************
  static G4bool iNin6[nEl]={false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false};
  static G4double VZ6[nEl][iNR];
  //********* S=-3 vectors *************
  static G4bool iNin7[nEl]={false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false};
  static G4double VZ7[nEl][iNR];
  // ********* S=-2 vectors *************
  static G4bool iNin8[nEl]={false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false};
  static G4double VZ8[nEl][iNR];
  // ********* S=-1 vectors *************
  static G4bool iNin9[nEl]={false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false};
  static G4double VZ9[nEl][iNR];
  // ********* S=0 vectors *************
  static G4bool iNin0[nEl]={false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false};
 static G4double VZ0[nEl][iNR];
  // ********* S=1 vectors *************
  static G4bool iNin1[nEl]={false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false};
  static G4double VZ1[nEl][iNR];
  // ********* S=2 vectors *************
  static G4bool iNin2[nEl]={false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false};
  static G4double VZ2[nEl][iNR];
  // ********* S=3 vectors *************
  static G4bool iNin3[nEl]={false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false};
  static G4double VZ3[nEl][iNR];
  // ********* S=2 vectors *************
  static G4bool iNin4[nEl]={false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,false,false,false,false,false,false};
  static G4double VZ4[nEl][iNR];
  //
#ifdef qdebug
  static G4int Smin=-1; // Dynamic Associated memory is implemented for S>-2 nuclei
  static G4int Smax= 2; // Dynamic Associated memory is implemented for S< 3 nuclei
  static G4int NZmin= 0; // Dynamic Associated memory is implemented for Z>-1 nuclei
  static G4int NNmin= 0; // Dynamic Associated memory is implemented for N>-1 nuclei
  static G4int NZS1max= 0; // Dynamic Associated memory is implemented for S<3, Z=-1, N<2
  static G4int NNS1max= 0; // Dynamic Associated memory is implemented for S<3, Z=-1, N<2 
#endif
  // -------------------------------------------------------------------------------------
  G4double rm=0.;
  G4int nz=n+z;
  G4int zns=nz+s_value;
  if(nz+s_value<0)
  {
    z=-z;
    n=-n;
    s_value=-s_value;
    nz=-nz;
    zns=-zns;
  }
  if(z<0)
  {
    if(z==-1)
    {
      if(!s_value)
      {
        if(n==1)      return mPiC;              // pi-
        else          return mPiC+(n-1)*mNeut;  // Delta- + (N-1)*n
      }
      else if(s_value==1)                       // Strange negative hadron
      {
        if(!n)        return mKM;               // K-
        else if(n==1) return mSiM;              // Sigma-
        else if(n==2) return mSmN ;             // Sigma- + n DiBaryon
        else if(n==3) return mSmNN;             // Sigma- +2n TriBaryon
        else          return mSigM+mNeut*(n-1); // Sigma- + (N-1)*n
      }
      else if(s_value==2)                       // --> Double-strange negative hadrons
      {
        if(!n)        return mKsM;              // Ksi-
        else if(n==1) return mKsiM+mNeut;       // Ksi- + n
        else if(n==2) return mKsiM+mNeut+mNeut; // Ksi- + 2n
        else          return mKsiM+mNeut*n;     // Ksi- + Z*n
      }
      else if(s_value==-2)
      {
        if     (nz==2)         return mDiKZ+mPiC; // 2K0 + Pi-
        else                   return mDiKZ+mPiC+(nz-2)*mProt;
      }
      else if(s_value==3)                       // --> Triple-strange negative hadrons
      {
        if     (n==-1) return mOmM;       // Triple-strange Omega minus
        else if(!n   ) return mOmN;       // Triple-strange (Omega-) + n DiBaryon
        else if(n==-2) return mDiKZ+mKM;  // Triple-strange K- + 2*K0
        else           return mOmeg+mNeut*(n+2);
      }
      else if(s_value==4)
      {
        if(n==-2)      return mOmeg+mKM;  // Omega- + K-
        else if(n==-1) return mOmeg+mLamb;// Omega- + Lambda
        else           return mOmeg+mLamb+(n+1)*mNeut; // Omega- + Lambda + (n+1)*Neutrons
      }
      else if(!n)            return mOmeg+(s_value-2)*mLamb; // Multy-Lambda + Omega minus
      else
      {
#ifdef qdebug
        if(s_value>NZS1max)
        {
          NZS1max=s_value;
          G4cout<<"-------->>G4QPDGCode::GetNucMass: Z=-1, S="<<s_value<<">2 with N="<<n<<G4endl;
        }
#endif
          return CalculateNuclMass(z,n,s_value);
      }
    }
    else if(!s_value)
    {
      if(!nz)
      {
        if(n==2)             return mDiPi;
        else                 return mDiPi+(n-2)*mPiC;
      }
      else                   return mNeut*nz-z*mPiC+anb;
    }
    else if(!zns)           // !!! s=0 is treated above !!
    {
      if(s_value>0)          return anb+s_value*mKM+n*mPiC; // s*K- + n*Pi-
      else                   return anb-s_value*mKZ-z*mPiC; // (-s)*aK0 + (-z)*Pi-
    }
    else if(s_value==1)
    {
      if(!nz)
      {
        if(n==2)             return mSmPi;
        else                 return mSmPi+(n-2)*mPiC;
      }
      else                   return mSigM+nz*mNeut-(z+1)*mPiC;
    }
    else if(s_value==-1)     return mKZa-z*mPiC+(nz-1)*mNeut; // aK0+(nz-1)n+(-z)*Pi-
    else if(s_value==2)
    {
      if     (nz==-1)        return mKsiM+n*mPiC;
      else if(!nz)           return mKsiM+mNeut-(z+1)*mPiC;
      else                   return mKsiM+(nz+1)*mNeut-(z+1)*mPiC;
    }
    else if(s_value==-2)     return mDiKZ-z*mPiC+(nz-2)*mNeut;
    else if(s_value==3)
    {
      if     (nz==-2)        return mOmeg+(n+1)*mPiC; // Omega- + (n+1)* Pi-
      else if(nz==-1)        return mOmeg+mNeut+n*mPiC; // Omega- + n * Pi-
      else if(!nz)           return mOmeg+mDiNt+(n-1)*mPiC; // Omega- + 2N + (n-1)*Pi-
      else                   return mOmeg+(nz+2)*mProt-(z+1)*mPiC;
    }
    else if(s_value<-2)      return anb-s_value*mKZ-z*mPiC+(nz+s_value)*mNeut;
    else if(s_value==4)
    {
      if     (nz==-3)        return mOmeg+mKM+(n+1)*mPiC;        // Om- + K- + (n+1)*Pi-
      else if(nz==-2)        return mOmeg+mSigM+n*mPiC;          // Om- + Sig- + n*Pi-
      else if(nz==-1)        return mOmeg+mSigM+mNeut+(n-1)*mPiC;//Om- + Sig- +N +(n-1)*Pi-
      else if(!nz)           return mOmeg+mSigM+mDiNt+(n-2)*mPiC;//Om- +Sig- +2N +(n-2)*Pi-
      else                   return mOmeg+mSigM+(nz+2)*mDiNt-(z+2)*mPiC;//+(nz+2)N-(z+2)Pi-
    }
    // s=5: 2*K-, Ksi-; s=6: 3*K-, Omega-; s>6 adds K- and Sigma- instead of Protons
    else                     // !!All s<0 are done and s>4 can be easy extended if appear!!
    {
#ifdef qdebug
      if(z<NZmin)
      {
        NZmin=z;
        G4cout<<"---->>G4QPDGCode::GetNucMass: Z="<<z<<"<-1 with N="<<n<<", S="<<s_value<<G4endl;
      }
#endif
      return CalculateNuclMass(z,n,s_value);
    }
  }
  else if(n<0)
  {
    if(n==-1)
    {
      if(!s_value)
      {
        if(z==1)      return mPiC;              // pi+
        else          return mPiC+(z-1)*mProt;  // Delta++ + (Z-1)*p
      }
      else if(s_value==1)                       // --> Strange neutral hadrons
      {
        if(!z)        return mKZ;               // K0
        else if(z==1) return mSiP;              // Sigma+
        else if(z==2) return mSpP ;             // Sigma+ + p DiBaryon
        else if(z==3) return mSpPP;             // Sigma+ +2p TriBaryon
        else          return mSigP+mProt*(z-1); // Sigma+ + (Z-1)*p
      }
      else if(s_value==2)                       // --> Double-strange negative hadrons
      {
        if(!z)        return mKsZ;              // Ksi0
        else if(z==1) return mKsiZ+mProt;       // Ksi- + p
        else if(z==2) return mKsiZ+mProt+mProt; // Ksi- + 2p
        else          return mKsiZ+mProt*z;     // Ksi- + Z*p
      }
      else if(s_value==-2)
      {
        if     (nz==2)         return mDiKM+mPiC; // 2K+ + Pi+
        else                   return mDiKM+mPiC+(nz-2)*mProt;
      }
      else if(s_value==3)
      {
        if(z==1) return mOmeg+mDiPr;
        else     return mOmeg+(z+1)*mProt;
      }
      else if(s_value==4) return mOmeg+mLamb+(z+1)*mProt;
      else if(!z) return mKZa+(s_value-1)*mLamb;      // Multy-Lambda + K0
      else
      {
#ifdef qdebug
        if(s_value>NNS1max)
        {
          NNS1max=s_value;
          G4cout<<"-------->>G4QPDGCode::GetNucMass: N=-1, S="<<s_value<<">2 with Z="<<z<<G4endl;
        }
#endif
          return CalculateNuclMass(z,n,s_value);
      }
    }
    else if(!s_value)
    {
      if(!nz)
      {
        if(z==2)             return mDiPi;
        else                 return mDiPi+(z-2)*mPiC;
      }
      else                   return mProt*nz-n*mPiC+anb;
    }
    else if(!zns)           // !!! s=0 is treated above !!
    {
      if(s_value>0)          return anb+s_value*mKZ+z*mPiC; // s*K0 + n*Pi+
      else                   return anb-s_value*mKM-n*mPiC; // (-s)*aK+ + (-n)*Pi+
    }
    else if(s_value==1)
    {
      if(!nz)
      {
        if(z==2)             return mSpPi;
        else                 return mSpPi+(z-2)*mPiC;
      }
      else                   return mSigP+nz*mProt-(n+1)*mPiC;
    }
    else if(s_value==-1)     return mKMa-n*mPiC+(nz-1)*mProt; // K+ + (nz-1)*P + (-n)*Pi+
    else if(s_value==2)
    {
      if     (nz==-1)        return mKsiZ+z*mPiC;
      else if(!nz)           return mKsiZ+mProt-(n+1)*mPiC;
      else                   return mKsiZ+(nz+1)*mProt-(n+1)*mPiC;
    }
    else if(s_value==-2)     return mDiKM-n*mPiC+(nz-2)*mProt;
    else if(s_value==3)
    {
      if     (nz==-2)        return mOmeg+(z+1)*mPiC;       // Omega + (z+1)*Pi+
      else if(nz==-1)        return mOmeg+mProt+z*mPiC;     // Omega- + P +z*Pi+
      else if(!nz)           return mOmeg+mDiPr+(z-1)*mPiC; // Omega- + 2P + (z-1)* Pi+
      else                   return mOmeg+(nz+2)*mProt-(n+1)*mPiC;
    }
    else if(s_value<-2)      return anb-s_value*mKM-n*mPiC+(nz+s_value)*mProt;
    else if(s_value==4)
    {
      if     (nz==-3)        return mOmeg+mKZ+(z+1)*mPiC;        // Om- + K0 + (z+1)*Pi+
      else if(nz==-2)        return mOmeg+mSigP+z*mPiC;          // Om- + Sig+ + z*Pi+
      else if(nz==-1)        return mOmeg+mSigP+mProt+(z-1)*mPiC;// Om- +Sig+ +P +(z-1)*Pi+
      else if(!nz)           return mOmeg+mSigP+mDiPr+(z-2)*mPiC;//Om- +Sig+ +2P +(z-2)*Pi+
      else                   return mOmeg+mSigP+(nz+2)*mProt-(n+2)*mPiC;//+(nz+2)P-(n+2)Pi+
    }
    // s=5: 2*KZ, Ksi0; s=6: 3*KZ, Omega-; s>6 adds K0 and Sigma+ instead of Protons
    else                     // !!All s<0 are done and s>4 can be easy extended if appear!!
    {
#ifdef qdebug
      if(n<NNmin)
      {
        NNmin=n;
        G4cout<<"---->>G4QPDGCode::GetNucMass: N="<<n<<"<-1 with Z="<<z<<", S="<<s_value<<G4endl;
      }
#endif
      return CalculateNuclMass(z,n,s_value);
    }
  }
  //return CalculateNuclMass(z,n,s_value); // @@ This is just to compare the calculation time @@
  if(nz >= 256 || z >= nEl) return 256000.;
  if(!s_value)             // **************> S=0 nucleus
  {
    G4int Nmin=iNF[z];     // Minimun N(Z) for the Dynamic Associative Memory (DAM)
    if(!iNin0[z])          // =--=> This element is already initialized
    {
#ifdef idebug
      G4cout<<"**>G4QPDGCode::GetNucM:Z="<<z<<", S=0 is initialized. F="<<iNin0[z]<<G4endl;
#endif
      G4int iNfin=iNL[z];
      if(iNfin>iNR) iNfin=iNR;
      for (G4int in=0; in<iNfin; in++) VZ0[z][in] = CalculateNuclMass(z,in+Nmin,s_value);
      iNin0[z]=true;
    }
    G4int dNn=n-Nmin;
    if(dNn<0)              // --> The minimum N must be reduced 
    {
#ifdef qdebug
      if(n<iNmin[z])
      {
        G4cout<<"-->>G4QPDGCode::GetNucM:Z="<<z<<", S=0 with N="<<n<<"<"<<iNmin[z]<<G4endl;
        iNmin[z]=n;
      }
#endif
      return CalculateNuclMass(z,n,s_value);
    }
    else if(dNn<iNL[z]) return VZ0[z][dNn]; // Found in DAM
    else                   // --> The maximum N must be increased
    {
#ifdef qdebug
      if(dNn>iNmax)
      {
        G4cout<<"**>>G4QPDGCode::GetNucM:Z="<<z<<", S=0 with dN="<<dNn<<">"<<iNmax<<G4endl;
        iNmax=dNn;
      }
      if(dNn>iNran[z])
      {
        G4cout<<">G4QPDGCode::GetNucM:Z="<<z<<", S=0 with dN="<<dNn<<">"<<iNran[z]<<G4endl;
        iNran[z]=dNn;
      }
#endif
      return CalculateNuclMass(z,n,s_value);
    }
  }
  else if(s_value==1)      // ******************> S=1 nucleus
  {
    G4int Nmin=iNF[z];     // Minimun N(Z) for the Dynamic Associative Memory (DAM)
#ifdef pdebug
    G4bool pPrint = !z && !n;
    if(pPrint)
      G4cout<<"G4QPDGC::GetNucM:Nmin="<<Nmin<<",iNin1="<<iNin1[0]<<",iNL="<<iNL[0]<<G4endl;
#endif
    if(!iNin1[z])          // =--=> This element is already initialized
    {
#ifdef pdebug
      if(pPrint)
      G4cout<<"**>G4QPDGCode::GetNucM:Z="<<z<<", S=1 is initialized. F="<<iNin1[z]<<G4endl;
#endif
      G4int iNfin=iNL[z];
      if(iNfin>iNR) iNfin=iNR;
      for (G4int in=0; in<iNfin; in++) VZ1[z][in] = CalculateNuclMass(z,in+Nmin,s_value);
      iNin1[z]=true;
    }
    G4int dNn=n-Nmin;
    if(dNn<0)              // --> The minimum N must be reduced 
    {
#ifdef qdebug
      if(n<iNmin[z])
      {
        G4cout<<"-->>G4QPDGCode::GetNucM:Z="<<z<<", S=1 with N="<<n<<"<"<<iNmin[z]<<G4endl;
        iNmin[z]=n;
      }
#endif
      return CalculateNuclMass(z,n,s_value);
    }
    else if(dNn<iNL[z]) return VZ1[z][dNn]; // Found in DAM
    else                   // --> The maximum N must be increased
    {
#ifdef qdebug
      if(dNn>iNmax)
      {
        G4cout<<"**>>G4QPDGCode::GetNucM:Z="<<z<<", S=1 with dN="<<dNn<<">"<<iNmax<<G4endl;
        iNmax=dNn;
      }
      if(dNn>iNran[z])
      {
        G4cout<<">G4QPDGCode::GetNucM:Z="<<z<<", S=1 with dN="<<dNn<<">"<<iNran[z]<<G4endl;
        iNran[z]=dNn;
      }
#endif
      return CalculateNuclMass(z,n,s_value);
    }
  }
  else if(s_value==-1)     // ******************> S=-1 nucleus
  {
    G4int Nmin=iNF[z];     // Minimun N(Z) for the Dynamic Associative Memory (DAM)
    if(!iNin9[z])          // =--=> This element is already initialized
    {
#ifdef idebug
      G4cout<<"*>G4QPDGCode::GetNucM:Z="<<z<<", S=-1 is initialized. F="<<iNin9[z]<<G4endl;
#endif
      G4int iNfin=iNL[z];
      if(iNfin>iNR) iNfin=iNR;
      for (G4int in=0; in<iNfin; in++) VZ9[z][in] = CalculateNuclMass(z,in+Nmin,s_value);
      iNin9[z]=true;
    }
    G4int dNn=n-Nmin;
    if(dNn<0)              // --> The minimum N must be reduced 
    {
#ifdef qdebug
      if(n<iNmin[z])
      {
        G4cout<<"->>G4QPDGCode::GetNucM:Z="<<z<<" ,S=-1 with N="<<n<<"<"<<iNmin[z]<<G4endl;
        iNmin[z]=n;
      }
#endif
      return CalculateNuclMass(z,n,s_value);
    }
    else if(dNn<iNL[z]) return VZ9[z][dNn]; // Found in DAM
    else                   // --> The maximum N must be increased
    {
#ifdef qdebug
      if(dNn>iNmax)
      {
        G4cout<<"**>G4QPDGCode::GetNucM:Z="<<z<<", S=-1 with dN="<<dNn<<">"<<iNmax<<G4endl;
        iNmax=dNn;
      }
      if(dNn>iNran[z])
      {
        G4cout<<"G4QPDGCode::GetNucM:Z="<<z<<", S=-1 with dN="<<dNn<<">"<<iNran[z]<<G4endl;
        iNran[z]=dNn;
      }
#endif
      return CalculateNuclMass(z,n,s_value);
    }
  }
  else if(s_value==2)      // ******************> S=2 nucleus
  {
    G4int Nmin=iNF[z];     // Minimun N(Z) for the Dynamic Associative Memory (DAM)
    if(!iNin2[z])          // =--=> This element is already initialized
    {
#ifdef idebug
      G4cout<<"**>G4QPDGCode::GetNucM:Z="<<z<<", S=2 is initialized. F="<<iNin2[z]<<G4endl;
#endif
      G4int iNfin=iNL[z];
      if(iNfin>iNR) iNfin=iNR;
      for (G4int in=0; in<iNfin; in++) VZ2[z][in] = CalculateNuclMass(z,in+Nmin,s_value);
      iNin2[z]=true;
    }
    G4int dNn=n-Nmin;
    if(dNn<0)              // --> The minimum N must be reduced 
    {
#ifdef qdebug
      if(n<iNmin[z])
      {
        G4cout<<"-->>G4QPDGCode::GetNucM:Z="<<z<<", S=2 with N="<<n<<"<"<<iNmin[z]<<G4endl;
        iNmin[z]=n;
      }
#endif
      return CalculateNuclMass(z,n,s_value);
    }
    else if(dNn<iNL[z]) return VZ2[z][dNn]; // Found in DAM
    else                   // --> The maximum N must be increased
    {
#ifdef qdebug
      if(dNn>iNmax)
      {
        G4cout<<"**>>G4QPDGCode::GetNucM:Z="<<z<<", S=2 with dN="<<dNn<<">"<<iNmax<<G4endl;
        iNmax=dNn;
      }
      if(dNn>iNran[z])
      {
        G4cout<<">G4QPDGCode::GetNucM:Z="<<z<<", S=2 with dN="<<dNn<<">"<<iNran[z]<<G4endl;
        iNran[z]=dNn;
      }
#endif
      return CalculateNuclMass(z,n,s_value);
    }
  }
  else if(s_value==-2)     // ******************> S=-2 nucleus
  {
    G4int Nmin=iNF[z];     // Minimun N(Z) for the Dynamic Associative Memory (DAM)
    if(!iNin8[z])          // =--=> This element is already initialized
    {
#ifdef idebug
      G4cout<<"*>G4QPDGCode::GetNucM:Z="<<z<<", S=-2 is initialized. F="<<iNin8[z]<<G4endl;
#endif
      G4int iNfin=iNL[z];
      if(iNfin>iNR) iNfin=iNR;
      for (G4int in=0; in<iNfin; in++) VZ8[z][in] = CalculateNuclMass(z,in+Nmin,s_value);
      iNin8[z]=true;
    }
    G4int dNn=n-Nmin;
    if(dNn<0)              // --> The minimum N must be reduced 
    {
#ifdef qdebug
      if(n<iNmin[z])
      {
        G4cout<<"->>G4QPDGCode::GetNucM:Z="<<z<<", S=-2 with N="<<n<<"<"<<iNmin[z]<<G4endl;
        iNmin[z]=n;
      }
#endif
      return CalculateNuclMass(z,n,s_value);
    }
    else if(dNn<iNL[z]) return VZ8[z][dNn]; // Found in DAM
    else                   // --> The maximum N must be increased
    {
#ifdef qdebug
      if(dNn>iNmax)
      {
        G4cout<<"**>G4QPDGCode::GetNucM:Z="<<z<<", S=-2 with dN="<<dNn<<">"<<iNmax<<G4endl;
        iNmax=dNn;
      }
      if(dNn>iNran[z])
      {
        G4cout<<"G4QPDGCode::GetNucM:Z="<<z<<", S=-2 with dN="<<dNn<<">"<<iNran[z]<<G4endl;
        iNran[z]=dNn;
      }
#endif
      return CalculateNuclMass(z,n,s_value);
    }
  }
  else if(s_value==-3)     // ******************> S=-3 nucleus
  {
    G4int Nmin=iNF[z];     // Minimun N(Z) for the Dynamic Associative Memory (DAM)
    if(!iNin7[z])          // =--=> This element is already initialized
    {
#ifdef idebug
      G4cout<<"*>G4QPDGCode::GetNucM:Z="<<z<<", S=-3 is initialized. F="<<iNin7[z]<<G4endl;
#endif
      G4int iNfin=iNL[z];
      if(iNfin>iNR) iNfin=iNR;
      for (G4int in=0; in<iNfin; in++) VZ7[z][in] = CalculateNuclMass(z,in+Nmin,s_value);
      iNin7[z]=true;
    }
    G4int dNn=n-Nmin;
    if(dNn<0)              // --> The minimum N must be reduced 
    {
#ifdef qdebug
      if(n<iNmin[z])
      {
        G4cout<<"->>G4QPDGCode::GetNucM:Z="<<z<<", S=-3 with N="<<n<<"<"<<iNmin[z]<<G4endl;
        iNmin[z]=n;
      }
#endif
      return CalculateNuclMass(z,n,s_value);
    }
    else if(dNn<iNL[z]) return VZ7[z][dNn]; // Found in DAM
    else                   // --> The maximum N must be increased
    {
#ifdef qdebug
      if(dNn>iNmax)
      {
        G4cout<<"**>G4QPDGCode::GetNucM:Z="<<z<<", S=-3 with dN="<<dNn<<">"<<iNmax<<G4endl;
        iNmax=dNn;
      }
      if(dNn>iNran[z])
      {
        G4cout<<"G4QPDGCode::GetNucM:Z="<<z<<", S=-3 with dN="<<dNn<<">"<<iNran[z]<<G4endl;
        iNran[z]=dNn;
      }
#endif
      return CalculateNuclMass(z,n,s_value);
    }
  }
  else if(s_value==3)      // ******************> S=3 nucleus
  {
    G4int Nmin=iNF[z];     // Minimun N(Z) for the Dynamic Associative Memory (DAM)
    if(!iNin3[z])          // =--=> This element is already initialized
    {
#ifdef idebug
      G4cout<<"**>G4QPDGCode::GetNucM:Z="<<z<<", S=3 is initialized. F="<<iNin3[z]<<G4endl;
#endif
      G4int iNfin=iNL[z];
      if(iNfin>iNR) iNfin=iNR;
      for (G4int in=0; in<iNfin; in++) VZ3[z][in] = CalculateNuclMass(z,in+Nmin,s_value);
      iNin3[z]=true;
    }
    G4int dNn=n-Nmin;
    if(dNn<0)              // --> The minimum N must be reduced 
    {
#ifdef qdebug
      if(n<iNmin[z])
      {
        G4cout<<"-->>G4QPDGCode::GetNucM:Z="<<z<<", S=3 with N="<<n<<"<"<<iNmin[z]<<G4endl;
        iNmin[z]=n;
      }
#endif
      return CalculateNuclMass(z,n,s_value);
    }
    else if(dNn<iNL[z]) return VZ3[z][dNn]; // Found in DAM
    else                   // --> The maximum N must be increased
    {
#ifdef qdebug
      if(dNn>iNmax)
      {
        G4cout<<"**>>G4QPDGCode::GetNucM:Z="<<z<<", S=3 with dN="<<dNn<<">"<<iNmax<<G4endl;
        iNmax=dNn;
      }
      if(dNn>iNran[z])
      {
        G4cout<<">G4QPDGCode::GetNucM:Z="<<z<<", S=3 with dN="<<dNn<<">"<<iNran[z]<<G4endl;
        iNran[z]=dNn;
      }
#endif
      return CalculateNuclMass(z,n,s_value);
    }
  }
  else if(s_value==-4)     // ******************> S=-4 nucleus
  {
    G4int Nmin=iNF[z];     // Minimun N(Z) for the Dynamic Associative Memory (DAM)
    if(!iNin6[z])          // =--=> This element is already initialized
    {
#ifdef idebug
      G4cout<<"*>G4QPDGCode::GetNucM:Z="<<z<<", S=-4 is initialized. F="<<iNin6[z]<<G4endl;
#endif
      G4int iNfin=iNL[z];
      if(iNfin>iNR) iNfin=iNR;
      for (G4int in=0; in<iNfin; in++) VZ6[z][in] = CalculateNuclMass(z,in+Nmin,s_value);
      iNin6[z]=true;
    }
    G4int dNn=n-Nmin;
    if(dNn<0)              // --> The minimum N must be reduced 
    {
#ifdef qdebug
      if(n<iNmin[z])
      {
        G4cout<<"->>G4QPDGCode::GetNucM:Z="<<z<<", S=-4 with N="<<n<<"<"<<iNmin[z]<<G4endl;
        iNmin[z]=n;
      }
#endif
      return CalculateNuclMass(z,n,s_value);
    }
    else if(dNn<iNL[z]) return VZ6[z][dNn]; // Found in DAM
    else                   // --> The maximum N must be increased
    {
#ifdef qdebug
      if(dNn>iNmax)
      {
        G4cout<<"**>G4QPDGCode::GetNucM:Z="<<z<<", S=-4 with dN="<<dNn<<">"<<iNmax<<G4endl;
        iNmax=dNn;
      }
      if(dNn>iNran[z])
      {
        G4cout<<"G4QPDGCode::GetNucM:Z="<<z<<", S=-4 with dN="<<dNn<<">"<<iNran[z]<<G4endl;
        iNran[z]=dNn;
      }
#endif
      return CalculateNuclMass(z,n,s_value);
    }
  }
  else if(s_value==4)      // ******************> S=4 nucleus
  {
    G4int Nmin=iNF[z];     // Minimun N(Z) for the Dynamic Associative Memory (DAM)
    if(!iNin4[z])          // =--=> This element is already initialized
    {
#ifdef idebug
      G4cout<<"*>G4QPDGCode::GetNucM:Z="<<z<<", S=4 is initialized. F="<<iNin4[z]<<G4endl;
#endif
      G4int iNfin=iNL[z];
      if(iNfin>iNR) iNfin=iNR;
      for (G4int in=0; in<iNfin; in++) VZ4[z][in] = CalculateNuclMass(z,in+Nmin,s_value);
      iNin4[z]=true;
    }
    G4int dNn=n-Nmin;
    if(dNn<0)              // --> The minimum N must be reduced 
    {
#ifdef qdebug
      if(n<iNmin[z])
      {
        G4cout<<"-->>G4QPDGCode::GetNucM:Z="<<z<<", S=4 with N="<<n<<"<"<<iNmin[z]<<G4endl;
        iNmin[z]=n;
      }
#endif
      return CalculateNuclMass(z,n,s_value);
    }
    else if(dNn<iNL[z]) return VZ4[z][dNn]; // Found in DAM
    else                   // --> The maximum N must be increased
    {
#ifdef qdebug
      if(dNn>iNmax)
      {
        G4cout<<"**>>G4QPDGCode::GetNucM:Z="<<z<<", S=4 with dN="<<dNn<<">"<<iNmax<<G4endl;
        iNmax=dNn;
      }
      if(dNn>iNran[z])
      {
        G4cout<<">G4QPDGCode::GetNucM:Z="<<z<<", S=4 with dN="<<dNn<<">"<<iNran[z]<<G4endl;
        iNran[z]=dNn;
      }
#endif
      return CalculateNuclMass(z,n,s_value);
    }
  }
  else
  {
#ifdef qdebug
    if(s_value<Smin || s_value>Smax)
    {
      if(s_value<Smin) Smin=s_value;
      if(s_value>Smax) Smax=s_value;
      G4cout<<">>G4QPDGCode::GetNucM:Z="<<z<<" with S="<<s_value<<",N="<<n<<" (Improve)"<<G4endl;
    }
#endif
    rm=CalculateNuclMass(z,n,s_value);
  }
#ifdef debug
  G4cout<<"G4QPDGCode::GetMass:GetNuclMass="<<rm<<",Z="<<z<<",N="<<n<<",S="<<s_value<<G4endl;
#endif
  return rm;
} // End of GetNuclearMass

// Calculate a nuclear mass for Z (a#of protons), N (a#of neutrons), & S (a#of lambdas) 
G4double G4QPDGCode::CalculateNuclMass(G4int z, G4int n, G4int s_value)
{
  static const G4double mP  = QHaM(21); // Proton
  static const G4double mN  = QHaM(20); // Neutron
  static const G4double mL  = QHaM(22); // Lambda
  static const G4double mD  = G4Deuteron::Deuteron()->GetPDGMass(); // Deuteron (H-2)
  static const G4double mT  = G4Triton::Triton()->GetPDGMass();     // Triton (H-3)
  static const G4double mHe3= G4He3::He3()->GetPDGMass();           // Hetrium (He-3)
  static const G4double mAl = G4Alpha::Alpha()->GetPDGMass();       // Alpha (He-4)
  static const G4double dmP = mP+mP;    // DiProton
  static const G4double dmN = mN+mN;    // DiNeutron
  static const G4double dmL = mL+mL;    // DiLambda
  static const G4double dLN = mL+mN;    // LambdaNeutron
  static const G4double dLP = mL+mP;    // LambdaProton
  static const G4double mSm = QHaM(23); // Sigma-
  static const G4double mSp = QHaM(25); // Sigma+
  static const G4double dSP = mSp+mP;   // ProtonSigma+
  static const G4double dSN = mSm+mN;   // NeutronSigma-
  static const G4double dnS = dSN+mN;   // 2NeutronsSigma-
  static const G4double dpS = dSP+mP;   // 2ProtonsSigma+
  static const G4double mXm = QHaM(26); // Ksi-
  static const G4double mXz = QHaM(27); // Ksi0
  static const G4double mOm = QHaM(44); // Omega-
  static const G4double dXN = mXm+mN;   // NeutronKsi-
  static const G4double dXP = mXz+mP;   // ProtonKsi0
  static const G4double dOP = mOm+mP;   // ProtonOmega-
  static const G4double dON = mOm+mN;   // NeutronOmega-
  static const G4double mK  = QHaM(18); // Charged Kaon
  static const G4double mK0 = QHaM(17); // Neutral Kaon
  static const G4double mPi = QHaM(15); // Charged Pion
  //////////static const G4double mPi0= QHaM(14);
  //static const G4int    nSh = 164;
  //static G4double sh[nSh] = {0.,                        // Fake element for C++ N=Z=0
  //                             -4.315548,   2.435504,  -1.170501,   3.950887,   5.425238,
  //                             13.342524,  15.547586,  22.583129,  23.983480,  30.561036,
  //                             33.761971,  41.471027,  45.532156,  53.835880,  58.495514,
  //                             65.693445,  69.903344,  76.899581,  81.329361,  88.979438,
  //                             92.908703, 100.316636, 105.013393, 113.081686, 118.622601,
  //                            126.979113, 132.714435, 141.413182, 146.433488, 153.746754,
  //                            158.665225, 165.988967, 170.952395, 178.473011, 183.471531,
  //                            191.231310, 196.504414, 204.617158, 210.251108, 218.373984,
  //                            223.969281, 232.168660, 237.925619, 246.400505, 252.392471,
  //                            260.938644, 267.191321, 276.107788, 282.722682, 291.881502,
  //                            296.998590, 304.236025, 309.562296, 316.928655, 322.240263,
  //                            329.927236, 335.480630, 343.233705, 348.923475, 356.911659,
  //                            362.785757, 370.920926, 376.929998, 385.130316, 391.197741,
  //                            399.451554, 405.679971, 414.101869, 420.346260, 428.832412,
  //                            435.067445, 443.526983, 449.880034, 458.348602, 464.822352,
  //                            473.313779, 479.744524, 488.320887, 495.025069, 503.841579,
  //                            510.716379, 519.451976, 525.036156, 532.388151, 537.899017,
  //                            545.252264, 550.802469, 558.402181, 564.101100, 571.963429,
  //                            577.980340, 586.063802, 592.451334, 600.518525, 606.832027,
  //                            614.831626, 621.205330, 629.237413, 635.489106, 643.434167,
  //                            649.691284, 657.516479, 663.812101, 671.715021, 678.061128,
  //                            686.002970, 692.343712, 700.360477, 706.624091, 714.617848,
  //                            721.100390, 729.294717, 735.887170, 744.216084, 751.017094,
  //                            759.551764, 766.377807, 775.080204, 781.965673, 790.552795,
  //                            797.572494, 806.088030, 813.158751, 821.655631, 828.867137,
  //                            836.860955, 842.183292, 849.195302, 854.731798, 861.898839,
  //                            867.783606, 875.313342, 881.443441, 889.189065, 895.680189,
  //                            903.679729, 910.368085, 918.579876, 925.543547, 933.790028,
  //                            940.811396, 949.122548, 956.170201, 964.466810, 971.516490,
  //                            979.766905, 986.844659, 995.113552,1002.212760,1010.418770,
  //                           1018.302560,1025.781870,1033.263560,1040.747880,1048.234460,
  //                           1055.723430,1063.214780,1070.708750,1078.204870,1085.703370,
  //                           1093.204260,1100.707530,1108.213070};
  //static const G4double b1=8.09748564; // MeV
  //static const G4double b2=-0.76277387;
  //static const G4double b3=83.487332;  // MeV
  //static const G4double b4=0.090578206;// 2*b4
  //static const G4double b5=0.676377211;// MeV
  //static const G4double b6=5.55231981; // MeV
  static const G4double b7=25.;        // MeV (Lambda binding energy predexponent)
  // --- even-odd difference is 3.7(MeV)/X
  // --- S(X>151)=-57.56-5.542*X**1.05
  static const G4double b8=10.5;       // (Lambda binding energy exponent)
  //static const G4double b9=-1./3.;
  static const G4double a2=0.13;       // LambdaBindingEnergy for deutron+LambdaSystem(MeV)
  static const G4double a3=2.2;        // LambdaBindingEnergy for (t/He3)+LambdaSystem(MeV)
  static const G4double um_value=931.49432;  // Umified atomic mass unit (MeV)
  //static const G4double me =0.511;     // electron mass (MeV) :: n:8.071, p:7.289
  static const G4double eps =0.0001;   // security addition for multybaryons
  //static G4double c[9][9]={// z=1     =2     =3     =4     =5     =6     =7     =8     =9
  //                 {13.136,14.931,25.320,38.000,45.000,55.000,65.000,75.000,85.000},//n=1
  //     {14.950, 2.425,11.680,18.374,27.870,35.094,48.000,60.000,72.000},  //n=2
  //     {25.930,11.390,14.086,15.770,22.921,28.914,39.700,49.000,60.000},  //n=3
  //     {36.830,17.594,14.908, 4.942,12.416,15.699,24.960,32.048,45.000},  //n=4
  //     {41.860,26.110,20.946,11.348,12.051,10.650,17.338,23.111,33.610},  //n=5
  //     {45.000,31.598,24.954,12.607, 8.668, 0.000, 5.345, 8.006,16.780},  //n=6
  //     {50.000,40.820,33.050,20.174,13.369, 3.125, 2.863, 2.855,10.680},  //n=7
  //     {55.000,48.810,40.796,25.076,16.562, 3.020, 0.101,-4.737,1.9520},  //n=8
  //     {60.000,55.000,50.100,33.660,23.664, 9.873, 5.683,-0.809,0.8730}}; //n=9
  if(z>107)
  {
#ifdef debug
    G4cout<<"***G4QPDGCode::CalcNuclMass: Z="<<z<<">107, N="<<n<<", S="<<s_value<<G4endl;
#endif
    return 256000.;
  }
  G4int Z=z;
  G4int N=n;
  G4int S=s_value;
#ifdef debug
  G4cout<<"G4QPDGCode::CalcNuclMass called with Z="<<Z<<",N="<<N<<", S="<<S<<G4endl;
#endif
  if(!N&&!Z&&!S)
  {
#ifdef debug
    //G4cout<<"G4QPDGCode::GetNuclMass(0,0,0)="<<mPi0<<"#0"<<G4endl;
#endif
    //return mPi0; // Pi0 mass (dangerose)
    return 0.;     // Photon mass
  }
  else if(!N&&!Z&&S==1) return mL;
  else if(!N&&Z==1&&!S) return mP;
  else if(N==1&&!Z&&!S) return mN;
  else if(!N&&!Z&&S>1) return mL*S+eps;
  else if(!N&&Z>1&&!S) return mP*Z+eps;
  else if(N>1&&!Z&&!S) return mN*N+eps;
  G4int A=Z+N;
  G4int Bn=A+S;
  //if((Z<0||N<0)&&!Bn)
  //{
  //  if     (N<0) return Bn*mL-Z*mK - N*mK0+eps*   S ;
  //  else              return Bn*mL+N*mPi-A*mK +eps*(N+S);
  //}
  //if(A<0&&Bn>=0)                      // Bn*LAMBDA's + (-(N+Z))*antiKaons
  //{
  //  if     (N<0&&Z<0) return Bn*mL-Z*mK -N*mK0+eps*   S ;
  //  else if(N<0)      return Bn*mL+Z*mPi-A*mK0+eps*(Z+S);
  //  else              return Bn*mL+N*mPi-A*mK +eps*(N+S);
  //}
  if(!Bn)                     // => "GS Mesons - octet" case (without eta&pi0)
  {
    if     (!S&&Z<0) return mPi*N;
    else if(!S&&N<0) return mPi*Z;
    else if ( (N == 1 && S == -1) || (N == -1 && S == 1) ) 
      return mK0; // Simple decision
    else if ( (S == 1 && Z == -1) || (S == -1 && Z == 1) ) 
      return mK;  // Simple decision
    else if(S>0)                                  // General decision
    {
      if     (-Z>S) return S*mK-(S+Z)*mPi+eps;
      else if(Z>=0) return S*mK0+Z*mPi+eps;
      else          return (S+Z)*mK0-Z*mK+eps;
    }
    else if(S<0)                                  // General decision
    {
      if     (Z>-S) return -S*mK+(S+Z)*mPi+eps;
      else if(Z<=0) return -S*mK0-Z*mPi+eps;
      else          return -(S+Z)*mK0+Z*mK+eps;
    }
  }
  else if(Bn==1)                    // => "GS Baryons - octet" case (withoutSigma0)
  {
    if     (Z== 1 && N== 0 && S== 0) return mP;
    else if(Z== 0 && N== 1 && S== 0) return mN;
    else if(Z== 0 && N== 0 && S== 1) return mL;
    else if(Z== 1 && N==-1 && S== 1) return mSp;  // Lower than Sigma+ (Simp.Decis)
    else if(Z==-1 && N== 1 && S== 1) return mSm;  // Lower than Sigma- (Simp.Decis)
    else if(Z== 0 && N==-1 && S== 2) return mXz;  // Lower than Xi0    (Simp.Decis)
    else if(Z==-1 && N== 0 && S== 2) return mXm;  // Lower than Xi-    (Simp.Decis)
    else if(Z==-1 && N==-1 && S== 3) return mOm;  // Lower than Omega- (Simp.Decis)
    else if(!S&&Z<0) return mN-mPi*Z+eps;         // Negative Isonuclei
    else if(!S&&N<0) return mP-mPi*N+eps;         // Positive Isonuclei
    else if(S==1)                                 // --> General decision
    {
      if     (N>1)   return mSm+(N-1)*mPi+eps;    // (Sigma-)+(n*PI-)
      else if(Z>1)   return mSp+(Z-1)*mPi+eps;    // (Sigma+)+(n*PI+)
    }
    else if(S==2)                                 // --> General decision
    {
      if     (N>0)   return mXm+N*mPi+eps;        // (Xi-)+(n*PI-)
      else if(Z>0)   return mXz+Z*mPi+eps;        // (Xi0)+(n*PI+)
    }
    else if(S==3)                                 // --> General decision
    {
      if     (N>-1)  return mOm+(N+1)*mPi+eps;    // (Omega-)+(n*PI-)
      else if(Z>-1)  return mOm+(Z+1)*mPi+eps;    // (Omega-)+(n*PI+)
    }
    else if(S>3)                                  // --> General Omega- decision
    {
      if   (-Z>S-2)  return mOm+(S-3)*mK +(2-Z-S)*mPi+eps;
      else if(Z>-1)  return mOm+(S-3)*mK0+(Z+1)+mPi+eps;
      else           return mOm+(S+Z-2)*mK0-(Z+1)*mK+eps;
    }
  }
  else if(Bn==2) // => "GS Baryons - decuplet" case (NP,LP, and LN are made below)
  {
    if     (Z== 2 && N== 0 && S== 0) return dmP;
    //else if(Z== 1 && N== 1 && S== 0) return 1875.61282; // Exact deuteron PDG Mass
    else if(Z== 1 && N== 1 && S== 0) return mD;   // Exact deuteron PDG Mass
    else if(Z== 0 && N== 2 && S== 0) return dmN;
    else if(Z== 2 && N==-1 && S== 1) return dSP;
    else if(Z== 1 && N== 0 && S== 1) return dLP;
    else if(Z== 0 && N== 1 && S== 1) return dLN;
    else if(Z==-1 && N== 2 && S== 1) return dSN;
    else if(Z== 1 && N==-1 && S== 2) return dXP;
    else if(Z== 0 && N== 0 && S== 2) return dmL;
    else if(Z==-1 && N== 1 && S== 2) return dXN;
    else if(Z== 0 && N==-1 && S== 3) return dOP;
    else if(Z==-1 && N== 0 && S== 3) return dON;
    else if(!S&&Z<0) return dmN-mPi*Z+eps;        // Negative Isonuclei
    else if(!S&&N<0) return dmP-mPi*N+eps;        // Positive Isonuclei
    else if(S==1)                                 // --> General decision
    {
      if     (N>2)   return dSP+(N-2)*mPi+eps;    // (nSigma-)+(n*PI-)
      else if(Z>2)   return dSN+(Z-1)*mPi+eps;    // (pSigma+)+(n*PI+)
    }
    else if(S==2)                                 // --> General decision
    {
      if     (N>1)   return dXN+(N-1)*mPi+eps;    // (nXi-)+(n*PI-)
      else if(Z>1)   return dXP+(Z-1)*mPi+eps;    // (pXi0)+(n*PI+)
    }
    else if(S==3)                                 // --> General decision
    {
      if     (N>0)   return dON+N*mPi+eps;        // (nOmega-)+(n*PI-)
      else if(Z>0)   return dOP+Z*mPi+eps;        // (pOmega-)+(n*PI+)
    }
    else if(S>3)                                  // --> General Omega- decision
    {
      if   (-Z>S-2)  return dON+(S-3)*mK +(2-Z-S)*mPi+eps;
      else if(Z>0)   return dOP+(S-3)*mK0+Z+mPi+eps;
      else           return dOP+(S+Z-3)*mK0-Z*mK+eps;
    }
    //else if(S>0)                                // @@ Implement General Decision
    //{
    //  //#ifdef debug
    //  G4cout<<"***G4QPDGCode::GetNuclMass:B=2, Z="<<Z<<",N="<<N<<",S="<<S<<G4endl;
    //  //#endif
    //  return bigM;                              // Exotic dibaryons (?)
    //}
  }
  else if(Bn==3)
  {
    if(!S)                                      // Bn=3
    {
      if     (Z==1 && N== 2) return mT;         // tritium
      else if(Z==2 && N== 1) return mHe3;       // hetrium
    }
    if(S== 1 && Z==-1 && N== 3)                 // nnSig-
    {
      if     (Z==-1 && N== 3) return dnS;       // nnSig-
      else if(Z== 3 && N==-1) return dpS;       // ppSig+
    }
  }
  else if(!S)
  {
    if(Z==2 && N==2) return mAl;                // Alpha
    else if(Z<0) return A*mN-Z*mPi+eps;         // Multybaryon Negative Isonuclei
    else if(Z>A) return A*mP+(Z-A)*mPi+eps;     // Multybaryon Positive Isonuclei
  }
  // === Start mesonic extraction ===
  G4double km_value=0.;               // Mass Sum of K mesons (G4E::DecayAntiStrang.)
  G4int Zm=Z;
  G4int Nm=N;
  G4int Sm=S;
  if(S<0&&Bn>0)                       // NEW: the free mass minimum
  {
    if(Zm>=-S)                        // Enough charge for K+'s
    {
      km_value=-S*mK;                 // Anti-Lambdas are compensated by protons
      Zm+=S;
    }
    else if(Zm>0)
    {
      km_value=Zm*mK-(S+Zm)*mK0;      // Anti-Lambdas are partially compensated by neutrons
      Zm=0;
      Nm+=S+Zm;
    }
  }
  else Sm=0;                          // No alternative calculations
  // Old algorithm
  G4double k=0.;                      // Mass Sum of K mesons
  if(S<0&&Bn>0)                       // OLD @@ Can be improved by K0/K+ search of minimum
  {
    G4int sH=(-S)/2;                  // SmallHalfS || The algorithm must be the same
    G4int bH=-S-sH;                   // BigHalhS   || as in G4QE::DecayAntiStrange
    if(Z>0 && Z>N)
    {
      if(Z>=bH)                       // => "Enough protons in nucleus" case
      {
        if(N>=sH)
        {
          k=bH*mK+sH*mK0;
          Z-=bH;
          N-=sH;
        }
        else
        {
          G4int dN=Z-N;
          if(dN>=-S)
          {
            k=-S*mK;
            Z+=S;
          }
          else
          {
            G4int sS=(-S-dN)/2;
            G4int bS=-S-dN-sS;
            sS+=dN;
            if(Z>=sS&&N>=bS)
            {
              k=sS*mK+bS*mK0;
              Z-=sS;
              N-=bS;
            }
            else if(Z<sS)
            {
              G4int dS=-S-Z;
              k=Z*mK+dS*mK0;
              N-=dS;
              Z=0;
            }
            else
            {
              G4int dS=-S-N;
              k=dS*mK+N*mK0;
              Z-=dS;
              N=0;
            }
          }
        }
      }
      else // Must not be here
      {
#ifdef debug
        G4cout<<"***G4QPDGC::CalcNuclMass:Antimatter? Z="<<Z<<",N="<<N<<",S="<<S<<G4endl;
#endif
        return 0.;                          // @@ Antiparticles aren't implemented @@
      }
    }
    else if(N>=bH)
    {
      if(Z>=sH)
      {
        k=sH*mK+bH*mK0;
        Z-=sH;
        N-=bH;
      }
      else
      {
        G4int dN=N-Z;
        if(dN>=-S)
        {
          k=-S*mK0;
          N+=S;
        }
        else
        {
          G4int sS=(-S-dN)/2;
          G4int bS=-S-dN-sS;
          bS+=dN;
          if(N>=bS&&Z>=sS)
          {
            k=sS*mK+bS*mK0;
            Z-=sS;
            N-=bS;
          }
          else if(N<bS)
          {
            G4int dS=-S-N;
            k=dS*mK+N*mK0;
            Z-=dS;
            N=0;
          }
          else
          {
            G4int dS=-S-Z;
            k=Z*mK+dS*mK0;
            N-=dS;
            Z=0;
          }
        }
      }
    }
    else // Must not be here
    {
      return 0.;                            // @@ Antiparticles aren't implemented @@
#ifdef debug
        G4cout<<"***G4QPDGC::CalcNuclMass:Antimatter? N="<<N<<",Z="<<Z<<",S="<<S<<G4endl;
#endif
    }
    S=0;
  }
  if(N<0)
  {
    k+=-N*mPi;
    Z+=N;
    N=0;
  }
  if(Z<0)
  {
    k+=-Z*mPi;
    N+=Z;
    Z=0;
  }
  A=Z+N;
  if (!A) return k+S*mL+S*eps;              // @@ multy LAMBDA states are not implemented
  G4double m_value=k+A*um_value;            // Expected mass in atomic units
  //G4double D=N-Z;                         // Isotopic shift of the nucleus
  if ( (A+S < 1 && k==0.) || Z < 0 || N < 0 ) 
  {  // @@ Can be generalized to anti-nuclei
#ifdef debug
    G4cout<<"**G4QPDGCode::CalcNuclMass:A="<<A<<"<1 || Z="<<Z<<"<0 || N="<<N<<"<0"<<G4endl;
    //@@throw G4QException("***G4QPDGCode::GetNuclMass: Impossible nucleus");
#endif
    return 0.;                              // @@ Temporary
  }
  if     (!Z) return k+N*(mN+.1)+S*(mL+.1); // @@ n+LAMBDA states are not implemented
  else if(!N) return k+Z*(mP+1.)+S*(mL+.1); // @@ p+LAMBDA states are not implemented
  //else if(N<=9&&Z<=9) m_value+=1.433e-5*pow(double(Z),2.39)-Z*me+c[N-1][Z-1];// Geant4 Comp.now
  else 
  {
    //G4double fA=A;
    // IsInTable is inside G4NucleiProperties::GetNuclearMass
    // G4ParticleTable::GetParticleTable()->FindIon(Zm,Am,0,Zm) creates new Ion!
    if(A==256 && Z==128) m_value=256000.;
    else                 m_value=k+G4NucleiProperties::GetNuclearMass(A,Z);
  }
  //@@//G4double maxM= k+Z*mP+N*mN+S*mL+eps;      // @@ eps -- Wings of the Mass parabola
  //@@//if(m_value>maxM) m_value=maxM;
  G4double mm_value=m_value;
  if(Sm<0)                                  // For the new algorithm of calculation 
  {
    if(Nm<0)
    {
      km_value+=-Nm*mPi;
      Zm+=Nm;
      Nm=0;
    }
    if(Zm<0)
    {
      km_value+=-Zm*mPi;
      Nm+=Zm;
      Zm=0;
    }
    G4int Am=Zm+Nm;
    if(!Am) return km_value+eps;
    mm_value=km_value+Am*um_value;          // Expected mass in atomic units
    //G4double Dm=Nm-Zm;                    // Isotopic shift of the nucleus
    if ( (Am < 1 && km_value==0.) || Zm < 0 || Nm < 0 ) 
    {   // @@ Can be generalized to anti-nuclei
#ifdef debug
      G4cerr<<"*G4QPDGCode::CalcNucM:A="<<Am<<"<1 || Z="<<Zm<<"<0 || N="<<Nm<<"<0"<<G4endl;
#endif
    }
    if     (!Zm) return km_value+Nm*(mN+.1);
    else if(!Nm) return km_value+Zm*(mP+1.);
    //else if(Nm<=9&&Zm<=9) mm_value+=1.433e-5*pow(double(Zm),2.39)-Zm*me+c[Nm-1][Zm-1];// Geant4
    else 
    {
      // IsInTable is inside G4NucleiProperties::GetNuclearMass
      // G4ParticleTable::GetParticleTable()->FindIon(Zm,Am,0,Zm) creates new Ion!
      mm_value=km_value+G4NucleiProperties::GetNuclearMass(Am,Zm);
    }
    //@@//G4double mM= km_value+Zm*mP+Nm*mN+eps;
    //@@//if(mm_value>mM) mm_value=mM;
  }
  if(m_value>mm_value) m_value=mm_value;
  if(S>0)
  {
    G4double bs=0.;
    if     (A==2) bs=a2;
    else if(A==3) bs=a3;
    else if(A>3)  bs=b7*exp(-b8/(A+1.));
    m_value+=S*(mL-bs);
  }  
#ifdef debug
  G4cout<<"G4QPDGCode::CalcNuclMass: >->-> OUT <-<-< m="<<m_value<<G4endl;
#endif
  //if(z==8&&n==9&&!s_value) G4cout<<"G4QPDGC::GetNucM:m="<<m_value<<",mm="<<mm_value<<G4endl;
  return m_value;
}

// Get a Quark Content {d,u,s,ad,au,as} for the PDG
G4QContent G4QPDGCode::GetQuarkContent() const
{
  G4int            a=0; // particle
  if(thePDGCode<0) a=1; // anti-particle
  G4int d=0;
  G4int u=0;
  G4int s_value=0;
  G4int ad=0;
  G4int au=0;
  G4int as=0;
  G4int ab=abs(thePDGCode);
  if     (!ab)
  {
    G4cerr<<"***G4QPDGCode::GetQuarkContent: PDG=0, return (0,0,0,0,0,0)"<<G4endl;
    return G4QContent(0,0,0,0,0,0);
  }
  else if(ab<4) // @@ for SU(6) : ab<7
  {
    if     (thePDGCode== 1) return G4QContent(1,0,0,0,0,0);
    else if(thePDGCode== 2) return G4QContent(0,1,0,0,0,0);
    else if(thePDGCode== 3) return G4QContent(0,0,1,0,0,0);
    else if(thePDGCode==-1) return G4QContent(0,0,0,1,0,0);
    else if(thePDGCode==-2) return G4QContent(0,0,0,0,1,0);
    else if(thePDGCode==-3) return G4QContent(0,0,0,0,0,1);
  }
  else if(ab<99)
  {
    if     (thePDGCode== 11) return G4QContent(1,0,0,0,1,0);
    else if(thePDGCode==-11) return G4QContent(0,1,0,1,0,0);
    else if(thePDGCode== 13) return G4QContent(1,0,0,0,1,0);
    else if(thePDGCode==-13) return G4QContent(0,1,0,1,0,0);
    else if(thePDGCode== 15) return G4QContent(1,0,0,0,1,0);
    else if(thePDGCode==-15) return G4QContent(0,1,0,1,0,0);
#ifdef debug
    if     (ab==22) G4cout<<"-W-G4QPDGC::GetQuarkCont: For the Photon? - Return 0"<<G4endl;
    else if(ab==10) G4cout<<"-W-G4QPDGC::GetQuarkCont: For Chipolino? - Return 0"<<G4endl;
    else G4cout<<"-W-G4QPDGCode::GetQuarkCont: For PDG="<<thePDGCode<<" Return 0"<<G4endl;
#endif
    return G4QContent(0,0,0,0,0,0); // Photon, bosons, leptons
  }
  else if(ab<80000000) // Baryons & Mesons
  {
    G4int c=ab/10;     // temporary (quarks)
    G4int f=c%10;      // (1) low quark
    G4int v=c/10;      // (2,3) temporary(B) or (2) final(M) (high quarks, high quark)
    G4int t=0;         // (3)prototype of highest quark (B)
#ifdef sdebug
    G4cout<<"G4QPDGCode::GetQuarkContent: a="<<ab<<", c="<<c<<", f="<<f<<", v="<<v<<G4endl;
#endif
    if(v>10)           // Baryons
    {
      t=v/10;          // (3) highest quark
      v%=10;           // (2) high quark
      if     (f==1)
      {
        if(a) ad++;
        else   d++;
      }
      else if(f==2)
      {
        if(a) au++;
        else   u++;
      }
      else if(f==3)
      {
        if(a) as++;
        else   s_value++;
      }
      else G4cerr<<"*G4QPDGC::GetQCont:1 PDG="<<thePDGCode<<","<<f<<","<<v<<","<<t<<G4endl;
      if     (v==1)
      {
        if(a) ad++;
        else   d++;
      }
      else if(v==2)
      {
        if(a) au++;
        else   u++;
      }
      else if(v==3)
      {
        if(a) as++;
        else   s_value++;
      }
      else G4cerr<<"*G4QPDGC::GetQCont:2 PDG="<<thePDGCode<<","<<f<<","<<v<<","<<t<<G4endl;
      if     (t==1)
      {
        if(a) ad++;
        else   d++;
      }
      else if(t==2)
      {
        if(a) au++;
        else   u++;
      }
      else if(t==3)
      {
        if(a) as++;
        else   s_value++;
      }
      else G4cerr<<"*G4QPDGC::GetQCont:3 PDG="<<thePDGCode<<","<<f<<","<<v<<","<<t<<G4endl;
      return G4QContent(d,u,s_value,ad,au,as);
    }
    else        // Mesons
    {
      if(f==v)
      {
        if     (f==1) return G4QContent(1,0,0,1,0,0);
        else if(f==2) return G4QContent(0,1,0,0,1,0);
        else if(f==3) return G4QContent(0,0,1,0,0,1);
        else G4cerr<<"*G4QPDGC::GetQC:4 PDG="<<thePDGCode<<","<<f<<","<<v<<","<<t<<G4endl;
      }
      else
      {
        if     (f==1 && v==2)
        {
          if(a)return G4QContent(1,0,0,0,1,0);
          else return G4QContent(0,1,0,1,0,0);
        }
        else if(f==1 && v==3)
        {
          if(a)return G4QContent(0,0,1,1,0,0);
          else return G4QContent(1,0,0,0,0,1);
        }
        else if(f==2 && v==3)
        {
          if(a)return G4QContent(0,0,1,0,1,0);
          else return G4QContent(0,1,0,0,0,1);
        }
        else G4cerr<<"*G4QPDGC::GetQC:5 PDG="<<thePDGCode<<","<<f<<","<<v<<","<<t<<G4endl;
      }
    }
  }
  else                    
  {
    G4int szn=ab-90000000;
    G4int ds=0;
    G4int dz=0;
    G4int dn=0;
    if(szn<-100000)
    {
      G4int ns_value=(-szn)/1000000+1;
      szn+=ns_value*1000000;
      ds+=ns_value;
    }
    else if(szn<-100)
    {
      G4int nz=(-szn)/1000+1;
      szn+=nz*1000;
      dz+=nz;
    }
    else if(szn<0)
    {
      G4int nn=-szn;
      szn=0;
      dn+=nn;
    }
    G4int sz =szn/1000;
    G4int n  =szn%1000;
    if(n>700)
    {
      n-=1000;
      dz--;
    }
    G4int z  =sz%1000-dz;
    if(z>700)
    {
      z-=1000;
      ds--;
    }
    s_value  =sz/1000-ds;
    G4int b=z+n+s_value;
    d=n+b;
    u=z+b;
    if      (d<0&&u<0&&s_value<0) return G4QContent(0,0,0,-d,-u,-s_value);
    else if (u<0&&s_value<0)      return G4QContent(d,0,0,0,-u,-s_value);
    else if (d<0&&s_value<0)      return G4QContent(0,u,0,-d,0,-s_value);
    else if (d<0&&u<0)            return G4QContent(0,0,s_value,-d,-u,0);
    else if (u<0)                 return G4QContent(d,0,s_value,0,-u,0);
    else if (s_value<0)           return G4QContent(d,u,0,0,0,-s_value);
    else if (d<0)                 return G4QContent(0,u,s_value,-d,0,0);
    else                          return G4QContent(d,u,s_value,0,0,0);
  }
  return G4QContent(0,0,0,0,0,0);
}

// Quark Content of pseudo-meson for q_i->q_o Quark Exchange
G4QContent G4QPDGCode::GetExQContent(G4int i, G4int o)  const
{
  G4QContent cQC(0,0,0,0,0,0);
  if     (!i)   cQC.IncAD();
  else if(i==1) cQC.IncAU();
  else if(i==2) cQC.IncAS();
  else G4cerr<<"***G4QPDGCode::GetExQContent: strange entering quark i="<<i<<G4endl;
  if     (!o)   cQC.IncD();
  else if(o==1) cQC.IncU();
  else if(o==2) cQC.IncS();
  else G4cerr<<"***G4QPDGCode::GetExQContent: strange exiting quark o="<<o<<G4endl;
  return cQC;
}

// Relative Cross Index for q_i->q_o (q=0(d),1(u),2(s), 7 means there is no such cluster)
G4int G4QPDGCode::GetRelCrossIndex(G4int i, G4int o)  const
{
  //                          pD,nD,DD,DD,ppDnnDpDDnDD n, p, L,S-,S+,X-,X0,O-
  static const G4int b01[16]={ 7,15, 7, 7, 7, 7, 7, 7, 1, 7, 7,-1, 7, 7, 7, 7};
  static const G4int b02[16]={ 7, 7, 7, 7, 7, 7, 7, 7, 2, 7, 7, 7, 7, 7, 7, 7};
  static const G4int b10[16]={18, 7, 7, 7, 7, 7, 7, 7, 7,-1, 7, 7,-2, 7, 7, 7}; // Iso+Bary
  static const G4int b12[16]={ 7, 7, 7, 7, 7, 7, 7, 7, 7, 1, 7, 7, 7, 7, 7, 7};
  static const G4int b20[16]={ 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,-2, 7,-3, 7,-4, 7};
  static const G4int b21[16]={ 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,-1,-3, 7,-3, 7, 7};
  //                          nn,np,pp,Ln,Lp,LL,nS,pS,nnpnppLnnLpnLppLLnLLpnnS
  static const G4int d01[16]={ 1, 1, 7, 1, 7, 7,-3, 7, 1, 7, 1, 1, 7, 1, 7,-5};
  static const G4int d02[16]={ 3, 3, 7, 2, 7, 7, 7, 7, 3, 3, 3, 3, 7, 7, 7, 7};
  static const G4int d10[16]={ 7,-1,-1, 7,-1, 7, 7,-4, 7,-1, 7,-1,-1, 7,-1, 7}; //B=2,3
  static const G4int d12[16]={ 7, 2, 2, 7, 1, 7, 7, 7, 2, 2, 7, 2, 2, 7, 7, 7};
  static const G4int d20[16]={ 7, 7, 7,-3,-3,-2, 7,-5, 7, 7, 7,-3,-3,-3,-3, 7};
  static const G4int d21[16]={ 7, 7, 7,-2,-2,-1,-6, 7, 7, 7,-2,-2, 7,-2,-2, 7};
  //                          nn,np,pp,Ln,Lp,nn,np,pp! n, p,nn,np,pp,Ln,Lp
  static const G4int m01[15]={ 1, 1, 7, 1, 7, 1, 1, 7, 1, 7, 1, 1, 7, 1, 7};// Multibaryons
  static const G4int m02[15]={ 3, 3, 7, 3, 3, 7, 7, 7, 3, 3, 3, 3, 7, 7, 7};// @@ Regular !
  static const G4int m10[15]={ 7,-1,-1, 7,-1, 7,-1,-1, 7,-1, 7,-1,-1, 7,-1};// 01->1,10->-1
  static const G4int m12[15]={ 7, 2, 2, 2, 2, 7, 7, 7, 2, 2, 7, 2, 2, 7, 7};// 12->2,21->-2
  static const G4int m20[15]={ 7, 7, 7,-3,-3, 7,-3,-3, 7, 7, 7,-3,-3,-3,-3};// 02->3,20->-3
  static const G4int m21[15]={ 7, 7, 7,-2,-2,-2,-2, 7, 7, 7,-2,-2, 7,-2,-2};
  //static const G4int fragmStart = 82;    // "Isonuclei are added & Leptons are added"
  static const G4int fragmStart = 45;    // "Isonuclei & Leptons are added, Reduced"

  if(theQCode<fragmStart) return 7;
  G4int sub=theQCode-fragmStart;
  if ( (sub > 1 && sub < 8) || (sub > 12 && sub <16)) return 7; // SuperIso & SuperStrange 
  G4int rel=sub;                         // case of nuclear baryons and isonuclei
  if     (sub>31) rel =(sub-32)%15;      // case of heavy fragments (BaryNum>3)
  else if(sub>15) rel = sub-16;          // case of nuclear di-baryon & tri-baryons
#ifdef debug
  G4cout<<"G4QPDGCode::RelGetCrossIndex:i="<<i<<",o="<<o<<",su="<<sub<<",re="<<rel<<G4endl;
#endif
  //if ( (rel > 5 && rel < 9) || (rel > 13 && rel <16)) return 7; // SuperStrange're closed
  if     (!i)                            // ==> input quark = 0(d) (d=-1/3)
  {
    if     (!o)       return 0;          // -> output quark = 0(d) => 0 = the same cluster
    else if(o==1)                        // -> output quark = 1(u) (ch--)
    {
      if     (sub<16) return b01[rel];
      else if(sub<32) return d01[rel];
      else            return m01[rel];
    }
    else if(o==2)
    {
      if     (sub<16) return b02[rel];
      else if(sub<32) return d02[rel];
      else            return m02[rel];
    }
  }
  else if(i==1)                          // ==> input quark = 1(u) (u=2/3)
  {
    if     (!o)                          // -> output quark = 0(d) (ch++)
    {
      if     (sub<16) return b10[rel];
      else if(sub<32) return d10[rel];
      else            return m10[rel];
    }
    else if(o==1)     return 0;         // -> output quark = 1(u) => 0 = the same cluster
    else if(o==2)
    {
      if     (sub<16) return b12[rel];
      else if(sub<32) return d12[rel];
      else            return m12[rel];
    }
  }
  else if(i==2)
  {
    if     (!o)
    {
      if     (sub<16) return b20[rel];
      else if(sub<32) return d20[rel];
      else            return m20[rel];
    }
    else if(o==1)
    {
      if     (sub<16) return b21[rel];
      else if(sub<32) return d21[rel];
      else            return m21[rel];
    }
    else if(o==2)     return 0;
  }
  return 7;
}

// Get number of Combinations for q_i->q_o
G4int G4QPDGCode::GetNumOfComb(G4int i, G4int o) const
{
  if(i>-1&&i<3)
  {
    G4int shiftQ=GetRelCrossIndex(i, o);
    G4int sQCode=theQCode;                     // QCode of the parent cluster
    if     (shiftQ==7) return 0;               // A parent cluster doesn't exist
    else if(!shiftQ) sQCode+=shiftQ;           // Shift QCode of son to QCode of his parent
    G4QPDGCode parent;                         // Create a temporary (fake) parent cluster
    parent.InitByQCode(sQCode);                // Initialize it by Q-Code
    G4QContent parentQC=parent.GetQuarkContent();// Quark Content of the parent cluster
    if     (!o)   return parentQC.GetD();
    else if(o==1) return parentQC.GetU();
    else if(o==2) return parentQC.GetS();
    else G4cerr<<"***G4QPDGCode:::GetNumOfComb: strange exiting quark o="<<o<<G4endl;
  }
  else G4cerr<<"***G4QPDGCode:::GetNumOfComb: strange entering quark i="<<i<<G4endl;
  return 0;
}

// Get a total number of Combinations for q_i
G4int G4QPDGCode::GetTotNumOfComb(G4int i) const
{
  G4int tot=0;
  if(i>-1&&i<3) for(int j=0; j<3; j++) tot+=GetNumOfComb(i, j);
  else G4cerr<<"***G4QPDGCode:::GetTotNumOfComb: strange entering quark i="<<i<<G4endl;
  return tot;
}

// Converts nuclear PDGCode to Z(#of protons), N(#of neutrons), S(#of lambdas) values
void G4QPDGCode::ConvertPDGToZNS(G4int nucPDG, G4int& z, G4int& n, G4int& s_value)
{
  if(nucPDG>80000000&&nucPDG<100000000)            // Condition of conversion
  {
    z=0;
    n=0;
    s_value=0;
    G4int r=nucPDG;
    if(r==90000000) return;
    G4int cn =r%1000;                              // candidate to #of neutrons
    if(cn)
    {
      if(cn>500) cn-=1000;                         // AntiNeutrons
      n=cn;                                        // Increment neutrons
      r-=cn;                                       // Subtract them from the residual
      if(r==90000000) return;
    }
    G4int cz =r%1000000;                           // candidate to #of neutrons
    if(cz)
    {
      if(cz>500000) cz-=1000000;                   // AntiProtons
      z=cz/1000;                                   // Number of protons
      r-=cz;                                       // Subtract them from the residual
      if(r==90000000) return;
    }
    G4int cs =r%10000000;                           // candidate to #of neutrons
    if(cs)
    {
      if(cs>5000000) cs-=10000000;                 // AntiLambda
      s_value=cs/1000000;                          // Number of Lambdas
    }
  }
  return;
}

// Only for irreducable DiQaDiQ! L1!=R1 && L1!=R2 && L2!=R1 && L2!=R2
std::pair<G4int,G4int> G4QPDGCode::MakeTwoBaryons(G4int L1, G4int L2, G4int R1, G4int R2)
{
  G4int dl=0;
  G4int ul=0;
  //G4int sl=0;
  if     (L1==1) ++dl;
  else if(L1==2) ++ul;
  //else if(L1==3) ++sl;
  if     (L2==1) ++dl;
  else if(L2==2) ++ul;
  //else if(L2==3) ++sl;
  if     (R1==1) ++dl;
  else if(R1==2) ++ul;
  //else if(R1==3) ++sl;
  if     (R2==1) ++dl;
  else if(R2==2) ++ul;
  //else if(R2==3) ++sl;
  if     (dl==2 && ul==2) return make_pair(1114,2212); // @@ can be (2112,2224)
  else if(dl==1 && ul==2) return make_pair(3112,2212);
  else if(dl==0 && ul==2) return make_pair(3212,3212); // @@ can be (3312,2212)
  else if(dl==2 && ul==1) return make_pair(3222,2112);
  else if(dl==1 && ul==1) return make_pair(3312,2112); // @@ can be (3322,2212)
  else if(dl==2 && ul==0) return make_pair(3112,3112); // @@ can be (3322,1122)
  //#ifdef debug
  else G4cout<<"-Warning-G4QPDGCode::MakeTwoBaryons: Irreduceble? L1="<<L1<<",L2="<<L2
             <<",R1="<<R1<<",R2="<<R2<<G4endl;
  //#endif
  return make_pair(2212,2112);                         // @@ Just theMinimum, makeException
}
