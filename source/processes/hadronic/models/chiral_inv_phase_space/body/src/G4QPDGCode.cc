// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4QPDGCode.cc,v 1.3 2000-09-10 13:58:58 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4QPDGCode ----------------
//             by Mikhail Kossov, Sept 1999.
//      class for Hadron definitions in CHIPS Model
// -------------------------------------------------------------------

//#define debug
//#define pdebug

#include "G4QPDGCodeVector.hh"

G4QPDGCode::G4QPDGCode(G4int PDGCode): thePDGCode(PDGCode)
{
  if(PDGCode) theQCode=MakeQCode(PDGCode);
  else        
  {
#ifdef debug
    cerr<<"***G4QPDGCode: Constructed with PDGCode=0, QCode=-2"<<endl;  
#endif
    theQCode=-2;
  }
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

G4QPDGCode::~G4QPDGCode() {};

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
  static const G4int modi = 81;      // Q Codes for more than di-baryon nuclei
  static G4int qC[modi] = { 22,  110,  220,  330,  111,  211,  221,  311,  321,  331,
                          2112, 2212, 3122, 3112, 3212, 3222, 3312, 3322,  113,  213,
                           223,  313,  323,  333, 1114, 2114, 2214, 2224, 3124, 3114,
                          3214, 3224, 3314, 3324, 3334,  115,  215,  225,  315,  325,
                           335, 2116, 2216, 3126, 3116, 3216, 3226, 3316, 3326,  117,
                           217,  227,  317,  327,  337, 1118, 2118, 2218, 2228, 3128,
                          3118, 3218, 3228, 3318, 3328, 3338,  119,  219,  229,  319,
                           329,  339,  90000001 ,  90001000 ,  91000000 ,  90000002 ,
                           90001001 ,  90002000 ,  91000001 ,  91001000 ,  92000000};
  static G4int aC[15] = {1,1000,999001,1000000,1000999,1999000,1999999,        // sum 1
                         2,1001,2000,1000001,1001000,1999001,2000000,2000999}; // sum 2
  if      (QCode<0)
  {
    cerr<<"***G4QPDGCode::MakePDGCode: negative Q Code ="<<QCode<<endl;
    return 0;
  }
  else if (QCode>=modi)
  {
    G4int q=QCode-modi;              // Starting BarNum=3
    G4int a=q/15+1;                  // BarNum/2
    G4int b=q%15;
    return 90000000+a*1001+aC[b];
  }
  return qC[theQCode];
}

// Make a Q Code out of the PDG Code
G4int G4QPDGCode::MakeQCode(const G4int& PDGCode)
{//   ===========================================
  static const G4int NUCPDG=90000000;
  G4int PDGC=abs(PDGCode);        // Qcode is always not negative
  if (PDGC<NUCPDG)                // Baryons & Mesons
  {
    if(PDGC>10000000)  return -2; // Impossible code
    if     (PDGC==10)  return -1; // Chipolino
    else if(PDGC==22)  return  0; // Photon
    else if(PDGC==110) return  1; // Low R-P: Sigma (pi,pi S-wave)
    else if(PDGC==220) return  2; // Midle Regeon-Pomeron
    else if(PDGC==330) return  3; // High Regeon-Pomeron
    G4int p=PDGC/10;              // Quark Content
    G4int r=PDGC%10;              // 2s+1
    G4int         Q= 0;
    if     (!r)
    {
      cerr<<"***G4QPDGCode::MakeQCode: (0) Unknown in Q-System code: "<<PDGCode<<endl;
      return -2;
    }
    else if(r==1) Q= 3;
    else if(r==2) Q= 9;
    else if(r==3) Q=17;
    else if(r==4) Q=23;
    else if(r==5) Q=34;
    else if(r==6) Q=40;
    else if(r==7) Q=48;
    else if(r==8) Q=54;
    else if(r==9) Q=66;     // @@ Not yet implemented family of mesons with spin 4
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
        cerr<<"***G4QPDGCode::MakeQCode: (1) Unknown in Q-System code: "<<PDGCode<<endl;
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
          cerr<<"***G4QPDGCode::MakeQCode: (2) Unknown in Q-System code: "<<PDGCode<<endl;
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
          cerr<<"***G4QPDGCode::MakeQCode: (3) Unknown in Q-System code: "<<PDGCode<<endl;
          return -2;
        }
	  }
	}
  }
  else                        // Nuclear Fragments
  {
    G4int r=PDGC-NUCPDG;      // cut the fake 90000000
    if(!r) return -2;
    G4int n=r%1000;           // a#of neutrons
    G4int a=r/1000;           // a#of protons + 1000*Lambdas
    G4int z=a%1000;           // a#of protons
    G4int s=a/1000;           // a#of s quarks
    G4int d=n+n+z+s;          // a#of d quarks
    G4int u=n+z+z+s;          // a#of u quarks
    G4int t=d+u+s;            // tot#of quarks
    if(t%3 || t<3)            // b=0 are in mesons
    {
      cerr<<"***G4QPDGCode::MakeQCode: Unknown PDGCode="<<PDGCode<<", t="<<t<<endl;
      return -2;
	}
    else
	{
      G4int b=t/3;            // baryon number
      if(b==1)                // baryons
      {
        if     (s==0&&u==1&&d==2) return 72;
        else if(s==0&&u==2&&d==1) return 73;
        else if(s==1&&u==1&&d==1) return 74;
        else
        {
          cerr<<"***MakeQCode: (5) Unknown in Q-System code: "<<PDGCode<<endl;
          return -2;
        }
      }
      else if(b==2)           // di-baryons
      {
        if     (s==0&&u==2&&d==4) return 75;
        else if(s==0&&u==3&&d==3) return 76;
        else if(s==0&&u==4&&d==2) return 77;
        else if(s==1&&u==2&&d==3) return 78;
        else if(s==1&&u==3&&d==2) return 79;
        else if(s==2&&u==2&&d==2) return 80;
        else
        {
          cerr<<"***G4QPDGCode::MakeQCode: (6) Unknown in Q-System code: "<<PDGCode<<endl;
          return -2;
        }
      }
      G4int c=b/2;
      G4int g=b%2;
      G4int h=c*3;
      G4int Q=57+c*15;
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
          cerr<<"***G4QPDGCode::MakeQCode: (7) Unknown in Q-System code: "<<PDGCode<<endl;
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
          cerr<<"***G4QPDGCode::MakeQCode: (8) Unknown in Q-System code: "<<PDGCode<<endl;
#endif
          return -2;
        }
	  }
	}
  }
  cerr<<"***G4QPDGCode::MakeQCode: () Unknown in Q-System code: "<<PDGCode<<endl;
  return -2;
}

// Get the mean mass value for the PDG
G4double G4QPDGCode::GetMass()
{//      =====================
  static const int nM = 72;
  static G4double m[nM] = {0., 980., 780., 1500., 134.98, 139.57, 547.3, 497.67, 493.68, 957.78
    , 939.5656, 938.2723, 1115.684,  1197.44, 1192.55, 1189.37, 1321.3, 1314.9,   770.,   770.
    ,   781.94,    896.1,   891.66, 1019.413,   1232.,   1232.,  1232.,  1232., 1519.5, 1387.2
    ,   1383.7,   1382.8,    1535.,   1531.8, 1672.45,  1318.1, 1318.1,  1275., 1432.4, 1425.6
    ,    1525.,    1680.,    1680.,    1820.,   1915.,   1915.,  1915.,  2025.,  2025.,  1691.
    ,    1691.,    1667.,    1776.,    1776.,   1854.,   1950.,  1950.,  1950.,  1950.,  2100.
    ,    2030.,    2030.,    2030.,    2127.,   2127.,   2252.,  2020.,  2020.,  2044.,  2045.
    ,    2045.,    2297.
    };
  //static G4int n[15]={-1, 0, 1,-1, 0,-2,-1, 0, 0, 1,-1, 0, 1,-1, 0}; 
  //static G4int z[15]={ 1, 0,-1, 0,-1, 0,-1,-2, 1, 0, 1, 0,-1, 0,-1};
  //static G4int s[15]={ 0, 0, 0, 1, 1, 2, 2, 2, 0, 0, 1, 1, 1, 2, 2};
  //
  static G4int n[15]={ 1, 0, 1, 0,-1, 0,-1, 2, 1, 0, 1, 0, 1, 0,-1}; 
  static G4int z[15]={ 0, 1,-1, 0, 1,-1, 0, 0, 1, 2, 0, 1,-1, 0, 1};
  static G4int s[15]={ 0, 0, 1, 1, 1, 2, 2, 0, 0, 0, 1, 1, 2, 2, 2};
  //...........................................................................
  G4int ab=theQCode;
  G4int szn=thePDGCode-90000000;
  if(szn==0) return 0.;
  if(ab<0 && szn>0)
  {
    G4int s = szn/1000000;
    G4int zn= szn%1000000;
    G4int z = zn /1000;
    G4int n = zn %1000;
    return GetNuclMass(z,n,s);
  }
  else if(ab<0)
  {
#ifdef pdebug
    cerr<<"***G4QPDGCode::GetMass: m=100000., QCode="<<theQCode<<",PDGCode="<<thePDGCode<<endl;
#endif
    return 100000.;
  }
  if(ab<nM) return m[ab];
  G4int a=ab-66;
  G4int c=a%15;
  G4int b=a/15;      // b_base/2
  if (ab<81)
  {
    if     (ab==72) return GetNuclMass(0,1,0); //  n
    else if(ab==73) return GetNuclMass(1,0,0); //  p
    else if(ab==74) return GetNuclMass(0,0,1); //  L
    else if(ab==75) return GetNuclMass(0,2,0); // nn
    else if(ab==76) return GetNuclMass(1,1,0); //  d
    else if(ab==77) return GetNuclMass(2,0,0); // pp
    else if(ab==78) return GetNuclMass(0,1,1); // nL
    else if(ab==79) return GetNuclMass(1,0,1); // pL
    else if(ab==80) return GetNuclMass(0,0,2); // LL
  }
  else return GetNuclMass(b+z[c],b+n[c],s[c]);
}

// Get the width value for the PDG
G4double G4QPDGCode::GetWidth()
//      =====================
{
  static const int nW = 72;
  static G4double width[nW] = {   0.,  70., 450., 112.,   0.,   0., .00118,  0.,   0., .203
                              ,   0.,   0.,   0.,   0.,   0.,   0.,   0.,    0., 160., 160.
                              , 8.41, 50.5, 50.8, 4.43, 120., 120., 120.,  120., 15.6,  39.
                              ,  36., 35.8,  10.,   9.,   0., 107., 107., 185.5, 109., 98.5
                              ,  76., 130., 130.,  80., 120., 120., 120.,   20.,  20., 160.
                              , 160., 168., 159., 159.,  87., 300., 300.,  300., 300., 200.
                              , 180., 180., 180.,  99.,  99.,  55., 387.,  387., 208., 198.
                              , 198., 149.
                              };
  G4int ab=abs(theQCode);
  if(ab<nW) return width[ab];
  return 0.;                  // @@ May be real width should be implemented for nuclei e.g. pp
}

// Get a nuclear mass for Z (a#of protons), N (a#of neutrons), & S (a#of lambdas) 
G4double G4QPDGCode::GetNuclMass(G4int Z, G4int N, G4int S)
//      ===================================================
{
  static const G4double delMi = G4QPDGCode(1114).GetMass();
  static const G4double delPP = G4QPDGCode(2224).GetMass();
  static const G4double mP    = G4QPDGCode(2212).GetMass();
  static const G4double mN    = G4QPDGCode(2112).GetMass();
  static const G4double mL    = G4QPDGCode(3122).GetMass();
  static const G4double mK    = G4QPDGCode( 321).GetMass();
  static const G4double mK0   = G4QPDGCode( 311).GetMass();
  static const G4int nSh = 164;
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
  static G4double b1=8.09748564; // MeV
  static G4double b2=-0.76277387;
  static G4double b3=83.487332;  // MeV
  static G4double b4=0.090578206;// 2*b4
  static G4double b5=0.676377211;// MeV
  static G4double b6=5.55231981; // MeV
  static G4double b7=25.;        // MeV (Lambda binding energy predexponent)
  // even-odd difference is 3.7(MeV)/X
  // S(X>151)=-57.56-5.542*X**1.05
  static G4double b8=10.5;       // (Lambda binding energy exponent)
  static G4double b9=-1./3.;
  static G4double a2=0.13;       // Lambda binding energy for the deutron+Lambda system (MeV)
  static G4double a3=2.2;        // Lambda binding energy for the (t or He3)+Lambda system (MeV)
  static G4double ml=1115.684;   // Lambda mass (MeV)
  static G4double um=931.49432;  // Umified atomic mass unit (MeV)
  static G4double me =0.511;     // electron mass (MeV)
  static G4double c[9][9]={
                   {13.136,14.931,25.320,38.000,45.000,55.000,65.000,75.000,85.000},
				   {14.950, 2.425,11.680,18.374,27.870,35.094,48.000,60.000,72.000},
				   {25.930,11.390,14.086,15.770,22.921,28.914,39.700,49.000,60.000},
				   {36.830,17.594,14.908, 4.942,12.416,15.699,24.960,32.048,45.000},
				   {41.860,26.110,20.946,11.348,12.051,10.650,17.338,23.111,33.610},
				   {45.000,31.598,24.954,12.607, 8.668, 0.000, 5.345, 8.006,16.780},
				   {50.000,40.820,33.050,20.174,13.369, 3.125, 2.863, 2.855,10.680},
				   {55.000,48.810,40.796,25.076,16.562, 3.020, 0.101,-4.737,1.9520},
				   {60.000,55.000,50.100,33.660,23.664, 9.873, 5.683,-0.809,0.8730}};
  G4double k=0.;                      // Mass Sum of K+ mesons
  if(S<0)                             // @@ Can be improved by K0/K+ search of minimum
  {
    if(Z+N+S<1) return 0.;
    if(Z>0)
	{
      if(Z>=-S)                         // => "Enough protons in nucleus" case
	  {
        k=-S*mK;
        Z+=S;
	  }
      else
	  {
        k=Z*mK-(S+Z)*mK0;
        N+=S+Z;
        if(N<0) return 0.;
        Z=0;
	  }
	}
    else if(N>0)
	{
      if(N>=-S)                         // => "Enough neutrons in nucleus" case
	  {
        k=-S*mK0;
        N+=S;
	  }
      else
	  {
        k=N*mK0-(S+N)*mK;
        Z+=S+N;
        if(Z<0) return 0.;
        N=0;
	  }
	}
    else return 0.;
    S=0;
  }
  G4double A=Z+N;                     // Baryon Number of the nucleus
  G4double m=k+A*um;                  // Expected mass in atomic units
  G4double D=N-Z;                     // Isotopic shift of the nucleus
  if (Z==-1||N==-1)                   // Iso-Nuclei
  {
    if(Z==-1) return k+delMi+(N-1)*mN+S*mL;
    else      return k+delPP+(Z-1)*mP+S*mL;
  }
  else if(A<1||Z<0||N<0)              // @@ Can be generalized to anti-nuclei
  {
    //G4cerr<<"***G4QPDGCode::GetNuclMass: A="<<A<<"<1 || Z="<<Z<<"<0 || N="<<N<<"<0"<<G4endl;
    //@@G4Exception("***G4QPDGCode::GetNuclMass: Impossible nucleus");
    return 0.;                        // @@ Temporary
  }
  if     (!Z) return k+N*mN+S*mL+(N+S)*.001;  // @@ n+LAMBDA states are not implemented
  else if(!N) return k+Z*mP+S*mL+(Z+S)*.001;  // @@ p+LAMBDA states are not implemented
  else if(N<=9&&Z<=9) m+=1.433e-5*pow(double(Z),2.39)-Z*me+c[N-1][Z-1];
  else m+=-sh[Z]-sh[N]+b1*D*D*pow(A,b2)+b3*(1.-2./(1.+exp(b4*D)))+Z*Z*(b5*pow(A,b9)+b6/A);
  //else m=G4NucleiProperties::GetNuclearMass(Z+N,Z);
  G4double maxM= k+Z*mP+N*mN+S*mL+.001;       // @@ .001 ?? Wings of the Mass parabola
  if(m>maxM) m=maxM;
  if(S>0)
  {
    G4double bs=0.;
    if     (S==2) bs=a2;
    else if(S==3) bs=a3;
    else if(S>3)  bs=b7*exp(-b8/(A+1.));
    m+=S*(mL-bs);
  }  

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
  if     (ab==22) return G4QContent(0,0,0,0,0,0); // Photon
  else if(ab<90000000) // Baryons & Mesons
  {
    G4int c=ab/10;     // temporary (quarks)
    G4int f=c%10;      // (1) low quark
    G4int v=c/10;      // (2,3) temporary(B) or (2) final(M) (high quarks, high quark)
    G4int t=0;         // (3)prototype of highest quark (B)
#ifdef debug
	cout<<"G4QPDGCode::GetQuarkContent: ab="<<ab<<", c="<<c<<", f="<<f<<", v="<<v<<endl;
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
      else cerr<<"*G4PDGCode::GetQuarkContent:1 PDG="<<thePDGCode<<","<<f<<","<<v<<","<<t<<endl;
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
      else cerr<<"*G4PDGCode::GetQuarkContent:2 PDG="<<thePDGCode<<","<<f<<","<<v<<","<<t<<endl;
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
      else cerr<<"*G4PDGCode::GetQuarkContent:3 PDG="<<thePDGCode<<","<<f<<","<<v<<","<<t<<endl;
      return G4QContent(d,u,s,ad,au,as);
	}
    else        // Mesons
	{
      if(f==v)
	  {
        if     (f==1) return G4QContent(1,0,0,1,0,0);
        else if(f==2) return G4QContent(0,1,0,0,1,0);
        else if(f==3) return G4QContent(0,0,1,0,0,1);
        else cerr<<"*G4PDGCode::GetQuarkContent:4 P="<<thePDGCode<<","<<f<<","<<v<<","<<t<<endl;
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
        else cerr<<"*G4PDGCode::GetQuarkContent:5 P="<<thePDGCode<<","<<f<<","<<v<<","<<t<<endl;
	  }
	}
  }
  else
  {
    G4int b=ab-90000000;
    G4int n=b%1000;
    G4int c=b/1000;
    G4int p=c%1000;
    s      =c/1000;
    b=s+p+n;
    d=b+n;
    u=b+p;
    if(a)return G4QContent(0,0,0,d,u,s);
    else return G4QContent(d,u,s,0,0,0);
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
  else cerr<<"***G4QPDGCode::GetExQContent: strange entering quark i="<<i<<endl;
  if     (!o)   cQC.IncD();
  else if(o==1) cQC.IncU();
  else if(o==2) cQC.IncS();
  else cerr<<"***G4QPDGCode::GetExQContent: strange exiting quark o="<<o<<endl;
  return cQC;
}


// Relative Cross Index for q_i->q_o (7 means there is no such cluster)
G4int G4QPDGCode::GetRelCrossIndex(G4int i, G4int o)  const
{//   =====================================================
  static const G4int b01[ 3]={ 1, 7, 7};
  static const G4int b02[ 3]={ 2, 7, 7};
  static const G4int b10[ 3]={ 7,-1, 7};                                     // Baryons
  static const G4int b12[ 3]={ 7, 1, 7};
  static const G4int b20[ 3]={ 7, 7,-2};
  static const G4int b21[ 3]={ 7, 7,-1};
  static const G4int d01[ 6]={ 1, 1, 7, 1, 7, 7};
  static const G4int d02[ 6]={ 3, 3, 7, 2, 7, 7};
  static const G4int d10[ 6]={ 7,-1,-1, 7,-1, 7};                            // Dibaryons
  static const G4int d12[ 6]={ 7, 2, 2, 7, 1, 7};
  static const G4int d20[ 6]={ 7, 7, 7,-3,-3,-2};
  static const G4int d21[ 6]={ 7, 7, 7,-2,-2,-1};
  static const G4int m01[15]={ 1, 7, 1, 1, 7, 1, 7, 1, 1, 7, 1, 7, 1, 1, 7}; // Multibaryons
  static const G4int m02[15]={ 3, 3, 3, 3, 7, 7, 7, 3, 3, 7, 3, 3, 7, 7, 7}; // @@ Regular !!!!
  static const G4int m10[15]={ 7,-1, 7,-1,-1, 7,-1, 7,-1,-1, 7,-1, 7,-1,-1}; // 01-> 1, 10-> -1
  static const G4int m12[15]={ 2, 2, 7, 2, 2, 7, 7, 7, 2, 2, 2, 2, 7, 7, 7}; // 12-> 2, 21-> -2
  static const G4int m20[15]={ 7, 7, 7,-3,-3,-3,-3, 7, 7, 7,-3,-3, 7,-3,-3}; // 02-> 3, 20-> -3 
  static const G4int m21[15]={ 7, 7,-2,-2, 7,-2,-2, 7, 7, 7,-2,-2,-2,-2, 7};
  // === I<->O Reversed version ===
  //static const G4int b01[ 3]={ 7,-1, 7};
  //static const G4int b02[ 3]={ 7, 7,-2};
  //static const G4int b10[ 3]={ 1, 7, 7};                                     // Baryons
  //static const G4int b12[ 3]={ 7, 7,-1};
  //static const G4int b20[ 3]={ 2, 7, 7};
  //static const G4int b21[ 3]={ 7, 1, 7};
  //static const G4int d01[ 6]={ 7,-1,-1, 7,-1, 7};
  //static const G4int d02[ 6]={ 7, 7, 7,-3,-3,-2};
  //static const G4int d10[ 6]={ 1, 1, 7, 1, 7, 7};                            // Dibaryons
  //static const G4int d12[ 6]={ 7, 7, 7,-2,-2,-1};
  //static const G4int d20[ 6]={ 3, 3, 7, 2, 7, 7};
  //static const G4int d21[ 6]={ 7, 2, 2, 7, 1, 7};
  //static const G4int m01[15]={ 7,-1, 7,-1,-1, 7,-1, 7,-1,-1, 7,-1, 7,-1,-1}; // Multibaryons
  //static const G4int m02[15]={ 7, 7, 7,-3,-3,-3,-3, 7, 7, 7,-3,-3, 7,-3,-3}; // @@ Regular !!!!
  //static const G4int m10[15]={ 1, 7, 1, 1, 7, 1, 7, 1, 1, 7, 1, 7, 1, 1, 7}; // 01-> 1, 10-> -1
  //static const G4int m12[15]={ 7, 7,-2,-2, 7,-2,-2, 7, 7, 7,-2,-2,-2,-2, 7}; // 12-> 2, 21-> -2
  //static const G4int m20[15]={ 3, 3, 3, 3, 7, 7, 7, 3, 3, 7, 3, 3, 7, 7, 7}; // 02-> 3, 20-> -3 
  //static const G4int m21[15]={ 2, 2, 7, 2, 2, 7, 7, 7, 2, 2, 2, 2, 7, 7, 7};
  static const G4int fragmStart = 72;

  if(theQCode<fragmStart)
  {
#ifdef debug
    cerr<<"***G4QPDGCode::RelGetCrossIndex: not nucleus QCode="<<theQCode<<endl;
#endif
    return 7;
  }
  G4int sub=theQCode-fragmStart;
  G4int rel=sub;                         // case of nuclear baryons
  if     (sub>8)rel =(sub-9)%15;         // case of heavy fragments (BaryNum>2)
  else if(sub>2)rel = sub-3;             // case of nuclear di-baryon
#ifdef debug
	cout<<"G4QPDGCode::RelGetCrossIndex:i="<<i<<",o="<<o<<",sub="<<sub<<",rel="<<rel<<endl;
#endif
  if     (!i)
  {
    if     (!o)      return 0;
    else if(o==1)
    {
      if     (sub<3) return b01[rel];
      else if(sub<9) return d01[rel];
      else           return m01[rel];
    }
    else if(o==2)
    {
      if     (sub<3) return b02[rel];
      else if(sub<9) return d02[rel];
      else           return m02[rel];
    }
    else cerr<<"***G4QPDGCode::RelGetCrossIndex: strange exiting quark (i=0) o="<<o<<endl;
  }
  else if(i==1)
  {
    if     (!o)
    {
      if     (sub<3) return b10[rel];
      else if(sub<9) return d10[rel];
      else           return m10[rel];
    }
    else if(o==1) return 0;
    else if(o==2)
    {
      if     (sub<3) return b12[rel];
      else if(sub<9) return d12[rel];
      else           return m12[rel];
    }
    else cerr<<"***G4QPDGCode::RelGetCrossIndex: strange exiting quark (i=0) o="<<o<<endl;
  }
  else if(i==2)
  {
    if     (!o)
    {
      if     (sub<3) return b20[rel];
      else if(sub<9) return d20[rel];
      else           return m20[rel];
    }
    else if(o==1)
    {
      if     (sub<3) return b21[rel];
      else if(sub<9) return d21[rel];
      else           return m21[rel];
    }
    else if(o==2) return 0;
    else cerr<<"***G4QPDGCode::RelGetCrossIndex: strange exiting quark (i=0) o="<<o<<endl;
  }
  else cerr<<"***G4QPDGCode::RelGetCrossIndex: strange entering quark i="<<i<<endl;
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
    else cerr<<"***G4QPDGCode:::GetNumOfComb: strange exiting quark o="<<o<<endl;
  }
  else cerr<<"***G4QPDGCode:::GetNumOfComb: strange entering quark i="<<i<<endl;
  return 0;
}

// Get a total number of Combinations for q_i
G4int G4QPDGCode::GetTotNumOfComb(G4int i) const
{//   ==========================================
  G4int tot=0;
  if(i>-1&&i<3) for(int j=0; j<3; j++) tot+=GetNumOfComb(i, j);
  else cerr<<"***G4QPDGCode:::GetTotNumOfComb: strange entering quark i="<<i<<endl;
  return tot;
}





