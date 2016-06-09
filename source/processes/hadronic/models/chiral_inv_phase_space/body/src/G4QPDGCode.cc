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
// $Id: G4QPDGCode.cc,v 1.44 2004/11/09 11:11:17 mkossov Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
//      ---------------- G4QPDGCode ----------------
//             by Mikhail Kossov, Sept 1999.
//      class for Hadron definitions in CHIPS Model
// -------------------------------------------------------------------

//#define debug
//#define pdebug
//#define sdebug

#include "G4QPDGCodeVector.hh"
#include <cmath>
#include <cstdlib>
using namespace std;

G4QPDGCode::G4QPDGCode(G4int PDGCode): thePDGCode(PDGCode)
{
#ifdef sdebug
  G4cout<<"G4QPDGCode:Constructer is called with PDGCode="<<PDGCode<<G4endl;  
#endif
  if(PDGCode) theQCode=MakeQCode(PDGCode);
  else        
  {
#ifdef sdebug
    G4cout<<"***G4QPDGCode: Constructed with PDGCode=0, QCode=-2"<<G4endl;  
#endif
    theQCode=-2;
  }
#ifdef sdebug
  G4cout<<"G4QPDGCode:Constructer(PDG) the QCode="<<theQCode<<G4endl;  
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
  thePDGCode =rhs.thePDGCode;
  theQCode   =rhs.theQCode;

  return *this;
}

G4QPDGCode::~G4QPDGCode() {}

// Standard output for QPDGCode
ostream& operator<<(ostream& lhs, G4QPDGCode& rhs)
//       =========================================
{
  lhs << "[ PDG=" << rhs.GetPDGCode() << ", Q=" << rhs.GetQCode() << "]";
  return lhs;
}

// Standard output for const QPDGCode
ostream& operator<<(ostream& lhs, const G4QPDGCode& rhs)
//       ===============================================
{
  lhs << "[ PDG=" << rhs.GetPDGCode() << ", Q=" << rhs.GetQCode() << "]";
  return lhs;
}

// Overloading of QPDGCode addition
G4int operator+(const G4QPDGCode& lhs, const G4QPDGCode& rhs)
//         =======================================================
{
  G4int  s  = lhs.GetPDGCode();
  return s += rhs.GetPDGCode();
}
G4int operator+(const G4QPDGCode& lhs, const G4int& rhs)
//         =======================================================
{
  G4int  s  = lhs.GetPDGCode();
  return s += rhs;
}
G4int operator+(const G4int& lhs, const G4QPDGCode& rhs)
//         =======================================================
{
  G4int  s  = lhs;
  return s += rhs.GetPDGCode();
}

// Overloading of QPDGCode subtraction
G4int operator-(const G4QPDGCode& lhs, const G4QPDGCode& rhs)
//         =======================================================
{
  G4int  s  = lhs.GetPDGCode();
  return s -= rhs.GetPDGCode();
}
G4int operator-(const G4QPDGCode& lhs, const G4int& rhs)
//         =======================================================
{
  G4int  s  = lhs.GetPDGCode();
  return s -= rhs;
}
G4int operator-(const G4int& lhs, const G4QPDGCode& rhs)
//         =======================================================
{
  G4int  s  = lhs;
  return s -= rhs.GetPDGCode();
}

// Overloading of QPDGCode multiplication
G4int operator*(const G4QPDGCode& lhs, const G4QPDGCode& rhs)
//         =======================================================
{
  G4int  s  = lhs.GetPDGCode();
  return s *= rhs.GetPDGCode();
}
G4int operator*(const G4QPDGCode& lhs, const G4int& rhs)
//         =======================================================
{
  G4int  s  = lhs.GetPDGCode();
  return s *= rhs;
}
G4int operator*(const G4int& lhs, const G4QPDGCode& rhs)
//         =======================================================
{
  G4int  s  = lhs;
  return s *= rhs.GetPDGCode();
}

// Overloading of QPDGCode division
G4int operator/(const G4QPDGCode& lhs, const G4QPDGCode& rhs)
{//   =======================================================
  G4int  s  = lhs.GetPDGCode();
  return s /= rhs.GetPDGCode();
}
G4int operator/(const G4QPDGCode& lhs, const G4int& rhs)
{//   ==================================================
  G4int  s  = lhs.GetPDGCode();
  return s /= rhs;
}
G4int operator/(const G4int& lhs, const G4QPDGCode& rhs)
{//   ==================================================
  G4int  s  = lhs;
  return s /= rhs.GetPDGCode();
}

// Overloading of QPDGCode residual
G4int operator%(const G4QPDGCode& lhs, const G4int& rhs)
{//   ==================================================
  G4int  s  = lhs.GetPDGCode();
  return s %= rhs;
}

// TRUE if it is not RealNeutral (111,221,331 etc), FALSE if it is.
G4bool G4QPDGCode::TestRealNeutral(const G4int& PDGCode)
{//    =================================================
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
{//   ===========================================
  //static const G4int modi = 81;  // Q Codes for more than di-baryon nuclei
  //static const G4int modi = 89;  // Q Codes for more than di-baryon nuclei "IsoNuclei"
  static const G4int modi = 122; // Q Codes for more than quarta-baryon nuclei "Lept/Hyper"
  static G4int qC[modi] = { 11,   12,   13,   14,   15,   16,   22,   23,   24,   25, // 10
                            37,  110,  220,  330,  111,  211,  221,  311,  321,  331, // 20
						  2112, 2212, 3122, 3112, 3212, 3222, 3312, 3322,  113,  213, // 30
						   223,  313,  323,  333, 1114, 2114, 2214, 2224, 3124, 3114, // 40
                          3214, 3224, 3314, 3324, 3334,  115,  215,  225,  315,  325, // 50
                           335, 2116, 2216, 3126, 3116, 3216, 3226, 3316, 3326,  117, // 60
                           217,  227,  317,  327,  337, 1118, 2118, 2218, 2228, 3128, // 70
                          3118, 3218, 3228, 3318, 3328, 3338,  119,  219,  229,  319, // 80
                           329,  339, 90002999 , 89999003 , 90003998 ,
                           89998004 , 90003999 , 89999004 , 90004998 , 89998005 ,     // 90
                           90000001 , 90001000 , 91000000 , 90999001 , 91000999 ,
                           91999000 , 91999999 , 92999000 , 90000002 , 90001001 ,     //100
                           90002000 , 91000001 , 91001000 , 92000000 , 90999002 ,
                           91001999 , 90001002 , 90002001 , 91000002 , 91001001 ,     //110
                           91002000 , 92000001 , 92001000 , 90999003 , 90001003 ,
						   90002002 , 90003001 , 91001002 , 91002001 , 92000002 ,     //120
                           92001001 , 92002000};
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

// Make a Q Code out of the PDG Code
G4int G4QPDGCode::MakeQCode(const G4int& PDGCode)
{//   ===========================================
  G4int PDGC=abs(PDGCode);        // Qcode is always not negative
  G4int s=0;
  G4int z=0;
  G4int n=0;
  if (PDGC>100000000)             // Not supported
  {
#ifdef debug
    G4cout<<"***G4QPDGCode::MakeQCode: Unknown in Q-System code: "<<PDGCode<<G4endl;
#endif
    return -2;
  }
  else if (PDGC>80000000&&PDGC<100000000) // Try to convert the NUCCoding to PDGCoding
  {
    if(PDGC==90000000) return 6;
    ConvertPDGToZNS(PDGC, z, n, s);
    G4int b=n+z+s;                                 // Baryon number
#ifdef debug
    G4cout<<"***G4QPDGCode::Z="<<z<<",N="<<n<<",S="<<s<<G4endl;
#endif
    if(b<0)                                        // ---> Baryons & Fragments
	{
	  b=-b;
      n=-n;
      z=-z;
      s=-s;
      PDGC=90000000+s*1000000+z*1000+n;            // New PDGC for anti-baryons
    }
    else if(!b)                                    // --> Mesons
	{
      //G4bool anti=false;                           // For the PDG conversion
      if(z<0)                                      // --> Mesons conversion
	  {
        n=-n;
        z=-z;
        s=-s;
        //anti=true;                                 // For the PDG conversion
      }
      if(!z)
	  {
        if(s>0)
	    {
          n=-n;
          s=-s;
          //anti=true;                               // For the PDG conversion
        }
        if     (s==-1) return 17;                  // K0
        else if(s==-2) return -1;                  // K0+K0 chipolino
        else           return -2;                  // Not supported by Q Code
      }
      else                                         // --> z>0
	  {
        if(z==1)
        {
          if   (s==-1) return 18;                  // K+
          else         return 15;                  // pi+
        }
        else if(z==2)  return -1;                  // Chipolino
        else           return -2;                  // Not supported by Q Code
      }
    } // End of meson case
    if(b>0)                                        // --> Baryon case
	{
      if(b==1)
	  {
        if(!s)                                     // --> Baryons
		{
          if(z==-1)    return 34;                  // Delta-
          else if(!z)  return 91;                  // neutron
          else if(z==1)return 91;                  // proton
          else if(z==2)return 37;                  // Delta++
          else if(z==3||z==-2)return -1;           // Delta+pi Chipolino
          else         return -2;                  // Not supported by Q Code
        }
        else if(s==1)                              // --> Hyperons
		{
          if(z==-1)    return 93;                  // Sigma-
          else if(!z)  return 92;                  // Lambda (@@ 24->Sigma0)
          else if(z==1)return 94;                  // Sigma+
          else if(z==2||z==-2) return -1;          // Sigma+pi Chipolino
          else         return -2;                  // Not supported by Q Code
        }
        else if(s==2)                              // --> Xi Hyperons
		{
          if(z==-1)    return 95;                  // Xi-
          else if(!z)  return 96;                  // Xi0
          else if(z==1||z==-2)return -1;           // Xi+pi Chipolino
          else         return -2;                  // Not supported by Q Code
        }
        else if(s==3)                              // --> Xi Hyperons
		{
          if(z==-1)    return 97;                  // Omega-
          else if(!z||z==-2)  return -1;           // Omega+pi Chipolino
          else         return -2;                  // Not supported by Q Code
        }
      }
      else
	  {
        if(b==2)
        {
          if     (PDGC==90002999) return 82;       // p DEL++
          else if(PDGC==89999003) return 83;       // n DEL-
          else if(PDGC==90003998) return 84;       // DEL++ DEL++
          else if(PDGC==89998004) return 85;       // DEL-  DEL-
          else if(PDGC==90999002) return 104;      // n Sigma-
          else if(PDGC==91001999) return 105;      // p Sigma+
        }
        if(b==3)
        {
          if     (PDGC==90003999) return 86;       // p p DEL++
          else if(PDGC==89999004) return 87;       // n n DEL-
          else if(PDGC==90004998) return 88;       // p DEL++ DEL++
          else if(PDGC==89998005) return 89;       // n DEL-  DEL-
          else if(PDGC==90999003) return 113;      // n n Sigma-
        }
      }
    }
  }
  if (PDGC<80000000)              // ----> Direct Baryons & Mesons
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
    else if(PDGC==110) return 11; // Low R-P: Sigma (pi,pi S-wave)
    else if(PDGC==220) return 12; // Midle Regeon-Pomeron
    else if(PDGC==330) return 13; // High Regeon-Pomeron
    G4int p=PDGC/10;              // Quark Content
    G4int r=PDGC%10;              // 2s+1
    G4int         Q= 0;
    if     (!r)
    {
#ifdef pdebug
      G4cout<<"***G4QPDGCode::MakeQCode: (0) Unknown in Q-System code: "<<PDGCode<<G4endl;
#endif
      return -2;
    }
    else if(r==1) Q=13;
    else if(r==2) Q=19;
    else if(r==3) Q=27;
    else if(r==4) Q=33;
    else if(r==5) Q=44;
    else if(r==6) Q=50;
    else if(r==7) Q=58;
    else if(r==8) Q=64;
    else if(r==9) Q=75;
    if(r%2)                 // Mesons are all the same
	{
      if     (p==11) return Q+=1;
      else if(p==21) return Q+=2;
      else if(p==22) return Q+=3;
      else if(p==31) return Q+=4;
      else if(p==32) return Q+=5;
      else if(p==33) return Q+=6;
      else
      {
#ifdef pdebug
        G4cout<<"***G4QPDGCode::MakeQCode:(1) Unknown in Q-System code: "<<PDGCode<<G4endl;
#endif
        return -2;
      }
	}
    else                    // Baryons
	{
      G4int s=r/2;
      if(s%2)               // N Family
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
#ifdef pdebug
          G4cout<<"**G4QPDGCode::MakeQCode:(2) Unknown in Q-System code:"<<PDGCode<<G4endl;
#endif
          return -2;
        }
	  }
	  else                  // Delta Family
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
#ifdef pdebug
          G4cout<<"**G4QPDGCode::MakeQCode:(3) Unknown in Q-System code:"<<PDGCode<<G4endl;
#endif
          return -2;
        }
	  }
	}
  }
  else                        // Nuclear Fragments
  {
    G4int d=n+n+z+s;          // a#of d quarks
    G4int u=n+z+z+s;          // a#of u quarks
    G4int t=d+u+s;            // tot#of quarks
    if(t%3 || abs(t)<3)       // b=0 are in mesons
    {
#ifdef pdebug
      G4cout<<"***G4QPDGCode::MakeQCode: Unknown PDGCode="<<PDGCode<<", t="<<t<<G4endl;
#endif
      return -2;
	}
    else
	{
      G4int b=t/3;            // baryon number
      if(b==1)                // baryons
      {
        if     (s==0&&u==1&&d==2) return 90;
        else if(s==0&&u==2&&d==1) return 91;
        else if(s==1&&u==1&&d==1) return 92;
        else
        {
#ifdef pdebug
          G4cout<<"**G4QPDGCode::MakeQCode:(5) Unknown in Q-System code:"<<PDGCode<<G4endl;
#endif
          return -2;
        }
      }
      else if(b==2)           // di-baryons
      {
        if     (s==0&&u==2&&d==4) return 98;
        else if(s==0&&u==3&&d==3) return 99;
        else if(s==0&&u==4&&d==2) return 100;
        else if(s==1&&u==2&&d==3) return 101;
        else if(s==1&&u==3&&d==2) return 102;
        else if(s==2&&u==2&&d==2) return 103;
        else
        {
#ifdef pdebug
          G4cout<<"**G4QPDGCode::MakeQCode:(6) Unknown in Q-System code:"<<PDGCode<<G4endl;
#endif
          return -2;
        }
      }
      else if(b==3)           // tri-baryons
      {
        if     (s==0&&u==4&&d==5) return 106; // pnn
        else if(s==0&&u==5&&d==4) return 107; // npp
        else if(s==1&&u==3&&d==5) return 108; // Lnn
        else if(s==1&&u==4&&d==4) return 109; // Lnp
        else if(s==1&&u==5&&d==3) return 110; // Lpp
        else if(s==2&&u==3&&d==4) return 111; // LLn
        else if(s==2&&u==4&&d==3) return 112; // LLp
        else if(s==1&&u==2&&d==6) return 113; // SIG-nn
        else
        {
#ifdef pdebug
          G4cout<<"**G4QPDGCode::MakeQCode:(7) Unknown in Q-System code:"<<PDGCode<<G4endl;
#endif
          return -2;
        }
      }
      G4int c=b/2;
      G4int g=b%2;
      G4int h=c*3;
      //G4int Q=57+c*15;
      //G4int Q=65+c*15;           // "IsoNuclei"
      G4int Q=83+c*15;           // "Leptons/Hyperons"
      u-=h;
      d-=h;
      if(g)
      {
        if     (s==0&&u==1&&d==2) return Q+= 9;
        else if(s==0&&u==2&&d==1) return Q+=10;
        else if(s==1&&u==0&&d==2) return Q+=11;
        else if(s==1&&u==1&&d==1) return Q+=12;
        else if(s==1&&u==2&&d==0) return Q+=13;
        else if(s==2&&u==0&&d==1) return Q+=14;
        else if(s==2&&u==1&&d==0) return Q+=15;
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
        if     (s==0&&u==-1&&d== 1) return Q+=1;
        else if(s==0&&u== 0&&d== 0) return Q+=2;
        else if(s==0&&u== 1&&d==-1) return Q+=3;
        else if(s==1&&u==-1&&d== 0) return Q+=4;
        else if(s==1&&u== 0&&d==-1) return Q+=5;
        else if(s==2&&u==-2&&d== 0) return Q+=6;
        else if(s==2&&u==-1&&d==-1) return Q+=7;
        else if(s==2&&u== 0&&d==-2) return Q+=8;
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
#ifdef pdebug
  G4cout<<"***G4QPDGCode::MakeQCode: () Unknown in Q-System code: "<<PDGCode<<G4endl;
#endif
// return -2; not reachable statement  
}

// Get the mean mass value for the PDG
G4double G4QPDGCode::GetMass()
{//      =====================
  //static const int nM = 72;
  //static const int nM = 80; // "Isobars"
  static const int nM = 90; // "Leptons/Hyperons"
  static G4double m[nM]={.511, 0., 105.658, 0., 1777.03, 0.,     0., 91.188, 80.423, 140.00
    ,120.000,    980.,     780.,    1500.,  134.98,  139.57,  547.3, 497.67, 493.68, 957.78
    ,939.566, 938.272, 1115.684,  1197.44, 1192.55, 1189.37, 1321.3, 1314.9,   770.,   770.
    , 781.94,   896.1,   891.66, 1019.413,   1232.,   1232.,  1232.,  1232., 1519.5, 1387.2
    , 1383.7,  1382.8,    1535.,   1531.8, 1672.45,  1318.1, 1318.1,  1275., 1432.4, 1425.6
    ,  1525.,   1680.,    1680.,    1820.,   1915.,   1915.,  1915.,  2025.,  2025.,  1691.
    ,  1691.,   1667.,    1776.,    1776.,   1854.,   1950.,  1950.,  1950.,  1950.,  2100.
    ,  2030.,   2030.,    2030.,    2127.,   2127.,   2252.,  2020.,  2020.,  2044.,  2045.
	, 2045., 2297., 2170.272, 2171.565, 2464., 2464., 3108.544, 3111.13, 3402.272, 3403.565
    };
  //old//static G4int dn[15]={-1, 0, 1,-1, 0,-2,-1, 0, 0, 1,-1, 0, 1,-1, 0}; 
  //old//static G4int dz[15]={ 1, 0,-1, 0,-1, 0,-1,-2, 1, 0, 1, 0,-1, 0,-1};
  //old//static G4int ds[15]={ 0, 0, 0, 1, 1, 2, 2, 2, 0, 0, 1, 1, 1, 2, 2};
  //
  //static G4int dn[15]={ 1, 0, 1, 0,-1, 0,-1, 2, 1, 0, 1, 0, 1, 0,-1}; 
  //static G4int dz[15]={ 0, 1,-1, 0, 1,-1, 0, 0, 1, 2, 0, 1,-1, 0, 1};
  //static G4int ds[15]={ 0, 0, 1, 1, 1, 2, 2, 0, 0, 0, 1, 1, 2, 2, 2};
  //...........................................................................
  G4int ab=theQCode;
  if(ab<0&&thePDGCode<80000000||!thePDGCode)
  {
#ifdef debug
    if(thePDGCode!=10)
      G4cout<<"**G4QPDGCode::GetMass:m=100000.,QC="<<theQCode<<",PDG="<<thePDGCode<<G4endl;
#endif
    return 100000.;
  }
  else if(ab>-1 && ab<nM)
  {
#ifdef debug
    G4cout<<"G4QPDGCode::GetMass:sm="<<m[ab]<<",Q="<<theQCode<<",PDG="<<thePDGCode<<G4endl;
#endif
    return m[ab];            // Get mass from the table
  }
  //if(szn==0) return m[15];                    // @@ mPi0   @@ MK ?
  if(thePDGCode==90000000)
  {
#ifdef debug
    G4cout<<"G4QPDGCode::GetMass:***m=0, QC="<<theQCode<<",PDG="<<thePDGCode<<G4endl;
#endif
    return 0.;
  }
  G4int s=0;
  G4int z=0;
  G4int n=0;
  ConvertPDGToZNS(thePDGCode, z, n, s);
  G4double rm=GetNuclMass(z,n,s);
#ifdef debug
  G4cout<<"G4QPDGCode::GetMass:GetNucMass="<<rm<<",Z="<<z<<",N="<<n<<",S="<<s<<G4endl;
#endif
  return rm;
}

// Get the width value for the PDG
G4double G4QPDGCode::GetWidth()
//      =====================
{
  //static const int nW = 72;
  //static const int nW = 80; // "Isobars"
  static const int nW = 90; // "Leptons/Hypernuclei"
  static G4double width[nW] = {0.,0.,0.,0.,0.,0.,0.,2.495,2.118,10.
    ,  10.,  70., 450., 112.,   0.,   0., .00118,  0.,   0., .203
    ,   0.,   0.,   0.,   0.,   0.,   0.,   0.,    0., 160., 160.
    , 8.41, 50.5, 50.8, 4.43, 120., 120., 120.,  120., 15.6,  39.
    ,  36., 35.8,  10.,   9.,   0., 107., 107., 185.5, 109., 98.5
    ,  76., 130., 130.,  80., 120., 120., 120.,   20.,  20., 160.
    , 160., 168., 159., 159.,  87., 300., 300.,  300., 300., 200.
    , 180., 180., 180.,  99.,  99.,  55., 387.,  387., 208., 198.
	, 198., 149., 120., 120., 170., 170., 120.,  120., 170., 170.};
  G4int ab=abs(theQCode);
  if(ab<nW) return width[ab];
  return 0.;             // @@ May be real width should be implemented for nuclei e.g. pp
}

// Get a nuclear mass for Z (a#of protons), N (a#of neutrons), & S (a#of lambdas) 
G4double G4QPDGCode::GetNuclMass(G4int Z, G4int N, G4int S)
//      ===================================================
{
  //static const G4double bigM= 1000000.;                   // Big Mass
  static const G4double mP  = G4QPDGCode(2212).GetMass(); // Proton
  static const G4double mN  = G4QPDGCode(2112).GetMass(); // Neutron
  static const G4double mL  = G4QPDGCode(3122).GetMass(); // Lambda
  static const G4double dmP = mP+mP;                      // DiProton
  static const G4double dmN = mN+mN;                      // DiNeutron
  static const G4double dmL = mL+mL;                      // DiLambda
  static const G4double dLN = mL+mN;                      // LambdaNeutron
  static const G4double dLP = mL+mP;                      // LambdaProton
  static const G4double mSm = G4QPDGCode(3112).GetMass(); // Sigma-
  static const G4double mSp = G4QPDGCode(3222).GetMass(); // Sigma+
  static const G4double dSP = mSp+mP;                     // ProtonSigma+
  static const G4double dSN = mSm+mN;                     // NeutronSigma-
  static const G4double dnS = dSN+mN;                     // 2NeutronsSigma-
  static const G4double mXm = G4QPDGCode(3312).GetMass(); // Ksi-
  static const G4double mXz = G4QPDGCode(3322).GetMass(); // Ksi0
  static const G4double mOm = G4QPDGCode(3334).GetMass(); // Omega-
  static const G4double dXN = mXm+mN;                     // NeutronKsi-
  static const G4double dXP = mXz+mP;                     // ProtonKsi0
  static const G4double dOP = mOm+mP;                     // ProtonOmega-
  static const G4double dON = mOm+mN;                     // NeutronOmega-
  static const G4double mK  = G4QPDGCode( 321).GetMass();
  static const G4double mK0 = G4QPDGCode( 311).GetMass();
  static const G4double mPi = G4QPDGCode( 211).GetMass();
  //////////static const G4double mPi0= G4QPDGCode( 111).GetMass();
  static const G4int    nSh = 164;
  static G4double sh[nSh] = {0.,                        // Fake element for C++ N=Z=0
                               -4.315548,   2.435504,  -1.170501,   3.950887,   5.425238,
                                13.342524,  15.547586,  22.583129,  23.983480,  30.561036,
                                33.761971,  41.471027,  45.532156,  53.835880,  58.495514,
                                65.693445,  69.903344,  76.899581,  81.329361,  88.979438,
                                92.908703, 100.316636, 105.013393, 113.081686, 118.622601,
                               126.979113, 132.714435, 141.413182, 146.433488, 153.746754,
                               158.665225, 165.988967, 170.952395, 178.473011, 183.471531,
                               191.231310, 196.504414, 204.617158, 210.251108, 218.373984,
                               223.969281, 232.168660, 237.925619, 246.400505, 252.392471,
                               260.938644, 267.191321, 276.107788, 282.722682, 291.881502,
                               296.998590, 304.236025, 309.562296, 316.928655, 322.240263,
                               329.927236, 335.480630, 343.233705, 348.923475, 356.911659,
                               362.785757, 370.920926, 376.929998, 385.130316, 391.197741,
                               399.451554, 405.679971, 414.101869, 420.346260, 428.832412,
                               435.067445, 443.526983, 449.880034, 458.348602, 464.822352,
                               473.313779, 479.744524, 488.320887, 495.025069, 503.841579,
                               510.716379, 519.451976, 525.036156, 532.388151, 537.899017,
                               545.252264, 550.802469, 558.402181, 564.101100, 571.963429,
                               577.980340, 586.063802, 592.451334, 600.518525, 606.832027,
                               614.831626, 621.205330, 629.237413, 635.489106, 643.434167,
                               649.691284, 657.516479, 663.812101, 671.715021, 678.061128,
                               686.002970, 692.343712, 700.360477, 706.624091, 714.617848,
                               721.100390, 729.294717, 735.887170, 744.216084, 751.017094,
                               759.551764, 766.377807, 775.080204, 781.965673, 790.552795,
                               797.572494, 806.088030, 813.158751, 821.655631, 828.867137,
                               836.860955, 842.183292, 849.195302, 854.731798, 861.898839,
                               867.783606, 875.313342, 881.443441, 889.189065, 895.680189,
                               903.679729, 910.368085, 918.579876, 925.543547, 933.790028,
                               940.811396, 949.122548, 956.170201, 964.466810, 971.516490,
                               979.766905, 986.844659, 995.113552,1002.212760,1010.418770,
                              1018.302560,1025.781870,1033.263560,1040.747880,1048.234460,
                              1055.723430,1063.214780,1070.708750,1078.204870,1085.703370,
                              1093.204260,1100.707530,1108.213070};
  static const G4double b1=8.09748564; // MeV
  static const G4double b2=-0.76277387;
  static const G4double b3=83.487332;  // MeV
  static const G4double b4=0.090578206;// 2*b4
  static const G4double b5=0.676377211;// MeV
  static const G4double b6=5.55231981; // MeV
  static const G4double b7=25.;        // MeV (Lambda binding energy predexponent)
  // even-odd difference is 3.7(MeV)/X
  // S(X>151)=-57.56-5.542*X**1.05
  static const G4double b8=10.5;       // (Lambda binding energy exponent)
  static const G4double b9=-1./3.;
  static const G4double a2=0.13;       // LambdaBindingEnergy for deutron+LambdaSystem(MeV)
  static const G4double a3=2.2;        // LambdaBindingEnergy for (t/He3)+LambdaSystem(MeV)
  static const G4double um=931.49432;  // Umified atomic mass unit (MeV)
  static const G4double me =0.511;     // electron mass (MeV) :: n:8.071, p:7.289
  static const G4double eps =0.0001;   // security addition for multybaryons
  static G4double c[9][9]={// z=1     =2     =3     =4     =5     =6     =7     =8     =9
                   {13.136,14.931,25.320,38.000,45.000,55.000,65.000,75.000,85.000},  //n=1
				   {14.950, 2.425,11.680,18.374,27.870,35.094,48.000,60.000,72.000},  //n=2
				   {25.930,11.390,14.086,15.770,22.921,28.914,39.700,49.000,60.000},  //n=3
				   {36.830,17.594,14.908, 4.942,12.416,15.699,24.960,32.048,45.000},  //n=4
				   {41.860,26.110,20.946,11.348,12.051,10.650,17.338,23.111,33.610},  //n=5
				   {45.000,31.598,24.954,12.607, 8.668, 0.000, 5.345, 8.006,16.780},  //n=6
				   {50.000,40.820,33.050,20.174,13.369, 3.125, 2.863, 2.855,10.680},  //n=7
				   {55.000,48.810,40.796,25.076,16.562, 3.020, 0.101,-4.737,1.9520},  //n=8
				   {60.000,55.000,50.100,33.660,23.664, 9.873, 5.683,-0.809,0.8730}}; //n=9
#ifdef debug
  G4cout<<"G4QPDGCode::GetNuclMass called with Z="<<Z<<",N="<<N<<", S="<<S<<G4endl;
#endif
  if(!N&&!Z&&!S)
  {
#ifdef debug
    //G4cout<<"G4QPDGCode::GetNuclMass(0,0,0)="<<mPi0<<"#0"<<G4endl;
#endif
    //return mPi0;
    return 0.;
  }
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
    else if(N==1&&S==-1||N==-1&&S==1) return mK0; // Simple decision
    else if(S==1&&Z==-1||S==-1&&Z==1) return mK;  // Simple decision
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
    else if(Z== 1 && N== 1 && S== 0) return 1875.6134; // Exact deuteron PDG Mass
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
  else if(Z==-1 && N== 3 && S== 1) return dnS;    // Bn=3
  else if(!S&&Z<0) return A*mN-Z*mPi+eps;         // Multybaryon Negative Isonuclei
  else if(!S&&Z>A) return A*mP+(Z-A)*mPi+eps;     // Multybaryon Positive Isonuclei
  // === Start mesonic extraction ===
  G4double km=0.;                     // Mass Sum of K mesons (G4E::DecayAntiStrang.)
  G4int Zm=Z;
  G4int Nm=N;
  G4int Sm=S;
  if(S<0&&Bn>0)                       // NEW: the free mass minimum
  {
    if(Zm>=-S)                        // Enough charge for K+'s
	{
      km=-S*mK;                       // Anti-Lambdas are compensated by protons
	  Zm+=S;
    }
    else if(Zm>0)
	{
      km=Zm*mK-(S+Zm)*mK0;            // Anti-Lambdas are partially compensated by neutrons
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
        G4cout<<"***G4QPDGC::GetNuclMass:Antimatter? Z="<<Z<<",N="<<N<<",S="<<S<<G4endl;
#endif
        return 0.;            // @@ Antiparticles aren't implemented @@
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
      return 0.;                   // @@ Antiparticles aren't implemented @@
#ifdef debug
        G4cout<<"***G4QPDGC::GetNuclMass:Antimatter? N="<<N<<",Z="<<Z<<",S="<<S<<G4endl;
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
  if (!A) return k+S*mL+S*eps;     // @@ multy LAMBDA states are not implemented
  G4double m=k+A*um;                 // Expected mass in atomic units
  G4double D=N-Z;                    // Isotopic shift of the nucleus
  if(A+S<1&&k==0.||Z<0||N<0)         // @@ Can be generalized to anti-nuclei
  {
#ifdef debug
    G4cout<<"***G4QPDGCode::GetNuclMass:A="<<A<<"<1 || Z="<<Z<<"<0 || N="<<N<<"<0"<<G4endl;
    //@@throw G4QException("***G4QPDGCode::GetNuclMass: Impossible nucleus");
#endif
    return 0.;                       // @@ Temporary
  }
  if     (!Z) return k+N*(mN+.1)+S*(mL+.1);  // @@ n+LAMBDA states are not implemented
  else if(!N) return k+Z*(mP+1.)+S*(mL+.1);  // @@ p+LAMBDA states are not implemented
  else if(N<=9&&Z<=9) m+=1.433e-5*pow(double(Z),2.39)-Z*me+c[N-1][Z-1];
  else 
  {
    G4double fA=A;
    if(G4NucleiPropertiesTable::IsInTable(Z,A))m=k+G4NucleiProperties::GetNuclearMass(A,Z);
    else
     m+=-sh[Z]-sh[N]+b1*D*D*pow(fA,b2)+b3*(1.-2./(1.+exp(b4*D)))+Z*Z*(b5*pow(fA,b9)+b6/fA);
  }
  G4double maxM= k+Z*mP+N*mN+S*mL+eps;       // @@ eps -- Wings of the Mass parabola
  if(m>maxM) m=maxM;
  G4double mm=m;
  if(Sm<0)                           // For the new algorithm of calculation 
  {
    if(Nm<0)
    {
      km+=-Nm*mPi;
      Zm+=Nm;
      Nm=0;
    }
    if(Zm<0)
    {
      km+=-Zm*mPi;
      Nm+=Zm;
      Zm=0;
    }
    G4int Am=Zm+Nm;
    if(!Am) return km+eps;
    mm=km+Am*um;                     // Expected mass in atomic units
    G4double Dm=Nm-Zm;               // Isotopic shift of the nucleus
    if(Am<1&&km==0.||Zm<0||Nm<0)     // @@ Can be generalized to anti-nuclei
    {
#ifdef debug
      G4cerr<<"**G4QPDGCode::GetNucM:A="<<Am<<"<1 || Z="<<Zm<<"<0 || N="<<Nm<<"<0"<<G4endl;
#endif
    }
    if     (!Zm) return km+Nm*(mN+.1);
    else if(!Nm) return km+Zm*(mP+1.);
    else if(Nm<=9&&Zm<=9) mm+=1.433e-5*pow(double(Zm),2.39)-Zm*me+c[Nm-1][Zm-1];
    else 
    {
      G4double fA=Am;
      if(G4NucleiPropertiesTable::IsInTable(Zm,Am)) mm=km+G4NucleiProperties::GetNuclearMass(Am,Zm);
      else mm+=-sh[Zm]-sh[Nm]+b1*Dm*Dm*pow(fA,b2)+b3*(1.-2./(1.+exp(b4*Dm)))+Zm*Zm*(b5*pow(fA,b9)+b6/Am);
    }
    G4double mM= km+Zm*mP+Nm*mN+eps;
    if(mm>mM) mm=mM;
  }
  if(m>mm) m=mm;
  if(S>0)
  {
    G4double bs=0.;
    if     (S==2) bs=a2;
    else if(S==3) bs=a3;
    else if(S>3)  bs=b7*exp(-b8/(A+1.));
    m+=S*(mL-bs);
  }  
#ifdef debug
  G4cout<<"G4QPDGCode::GetNuclMass: >>>OUT<<< m="<<m<<G4endl;
#endif
  return m;
}

// Get a Quark Content {d,u,s,ad,au,as} for the PDG
G4QContent G4QPDGCode::GetQuarkContent() const
//      ======================================
{
  G4int            a=0; // particle
  if(thePDGCode<0) a=1; // anti-particle
  G4int d=0;
  G4int u=0;
  G4int s=0;
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
#ifdef debug
    if     (ab==22) G4cout<<"-W-G4QPDGCode::GetQuarkContent: For the Photon? - Return 0"<<G4endl;
    else if(ab==10) G4cout<<"-W-G4QPDGCode::GetQuarkContent: For the Chipolino? - Return 0"<<G4endl;
	else G4cout<<"-W-G4QPDGCode::GetQuarkContent: For PDG="<<thePDGCode<<" Return 0"<<G4endl;
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
	G4cout<<"G4QPDGCode::GetQuarkContent: ab="<<ab<<", c="<<c<<", f="<<f<<", v="<<v<<G4endl;
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
        else   s++;
	  }
      else G4cerr<<"***G4QPDGCode::GetQContent:1 PDG="<<thePDGCode<<","<<f<<","<<v<<","<<t<<G4endl;
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
        else   s++;
	  }
      else G4cerr<<"***G4QPDGCode::GetQContent:2 PDG="<<thePDGCode<<","<<f<<","<<v<<","<<t<<G4endl;
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
        else   s++;
	  }
      else G4cerr<<"***G4QPDGCode::GetQCont:3 PDG="<<thePDGCode<<","<<f<<","<<v<<","<<t<<G4endl;
      return G4QContent(d,u,s,ad,au,as);
	}
    else        // Mesons
	{
      if(f==v)
	  {
        if     (f==1) return G4QContent(1,0,0,1,0,0);
        else if(f==2) return G4QContent(0,1,0,0,1,0);
        else if(f==3) return G4QContent(0,0,1,0,0,1);
        else G4cerr<<"***G4QPDGCode::GetQCont:4 PDG="<<thePDGCode<<",i="<<f<<","<<v<<","<<t<<G4endl;
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
        else G4cerr<<"***G4QPDGCode::GetQCont:5 PDG="<<thePDGCode<<","<<f<<","<<v<<","<<t<<G4endl;
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
      G4int ns=(-szn)/1000000+1;
      szn+=ns*1000000;
      ds+=ns;
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
    s  =sz/1000-ds;
    G4int b=z+n+s;
    d=n+b;
    u=z+b;
    if      (d<0&&u<0&&s<0) return G4QContent(0,0,0,-d,-u,-s);
    else if (u<0&&s<0)      return G4QContent(d,0,0,0,-u,-s);
    else if (d<0&&s<0)      return G4QContent(0,u,0,-d,0,-s);
    else if (d<0&&u<0)      return G4QContent(0,0,s,-d,-u,0);
    else if (u<0)           return G4QContent(d,0,s,0,-u,0);
    else if (s<0)           return G4QContent(d,u,0,0,0,-s);
    else if (d<0)           return G4QContent(0,u,s,-d,0,0);
    else                    return G4QContent(d,u,s,0,0,0);
  }
  return G4QContent(0,0,0,0,0,0);
}

// Quark Content of pseudo-meson for q_i->q_o Quark Exchange
G4QContent G4QPDGCode::GetExQContent(G4int i, G4int o)  const
{//        ==================================================
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
{//   =====================================================
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
  static const G4int fragmStart = 82;    // "Isonuclei are added && Leptons are added"

  if(theQCode<fragmStart) return 7;
  G4int sub=theQCode-fragmStart;
  if(sub>1&&sub<8||sub==15) return 7; //@@Why they are in clusters?-Residuals(?)
  G4int rel=sub;                         // case of nuclear baryons and isonuclei
  if     (sub>31)rel =(sub-32)%15;       // case of heavy fragments (BaryNum>3)
  else if(sub>15)rel = sub-16;           // case of nuclear di-baryon & tri-baryons
#ifdef debug
  G4cout<<"G4QPDGCode::RelGetCrossIndex:i="<<i<<",o="<<o<<",su="<<sub<<",re="<<rel<<G4endl;
#endif
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
{//   ================================================
  if(i>-1&&i<3)
  {
    G4int shiftQ=GetRelCrossIndex(i, o);
    G4int sQCode=theQCode;                        // QCode of the parent cluster
    if     (shiftQ==7) return 0;                  // A parent cluster doesn't exist
    else if(!shiftQ) sQCode+=shiftQ;              // Shift QCode of son to QCode of his parent
    G4QPDGCode parent;                            // Create a temporary (fake) parent cluster
    parent.InitByQCode(sQCode);                   // Initialize it by Q-Code
    G4QContent parentQC=parent.GetQuarkContent(); // Quark Content of the parent cluster
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
{//   ==========================================
  G4int tot=0;
  if(i>-1&&i<3) for(int j=0; j<3; j++) tot+=GetNumOfComb(i, j);
  else G4cerr<<"***G4QPDGCode:::GetTotNumOfComb: strange entering quark i="<<i<<G4endl;
  return tot;
}

// Converts nuclear PDGCode to Z(#of protons), N(#of neutrons), S(#of lambdas) values
void G4QPDGCode::ConvertPDGToZNS(G4int nucPDG, G4int& z, G4int& n, G4int& s)
{//  =======================================================================
  if(nucPDG>80000000&&nucPDG<100000000)            // Condition of conversion
  {
    G4int r=nucPDG-90000000;
    if(!r)
	{
      z=0;
      n=0;
      s=0;
      return;
	}
    // Antinucleus extraction
    if(r<-200000)                                  // Negative -> anLambdas              
    {
      G4int ns=(-r-200000)/1000000+1;
      r+=ns*1000000;                               // Get out aL from PDG
      s=-ns;                                       // Remember aLambdas
    }
    if(r<-200)                                     // Negative -> aProtons
    {
      G4int nz=(-r-200)/1000+1;
      r+=nz*1000;                                  // Get out aP from PDG
      z=-nz;                                       // Remember aProtons
    }
    if(r<0)                                        // Negative -> aNeutrons
    {
      G4int nn=-r;
      r=0;                                         // Get out aN from PDG
      n=-nn;                                       // Remember aNeutrons
    }
    G4int sz =r/1000;                              // Residual to analize
    n+=r%1000;                                     // A#of Neutrons
    if(n>700)                                      // AntiNutrons
    {
      n-=1000;
      z++;
    }
    z+=sz%1000;                                    // A#of Protons
    if(z>700)                                      // AntiProtons
    {
      z-=1000;
      s++;
    }
    s+=sz/1000;                                    // A#of Lambdas
  }
  return;
}





