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
// $Id: G4QBesIKJY.cc,v 1.4 2009-11-10 17:13:46 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//      ---------------- G4QBesIKJY ----------------
//             by Mikhail Kossov, Sept 1999.
//      class  for Bessel I0/I1 and K0/K1 functions in CHIPS Model
// -------------------------------------------------------------------
// Short description: Bessel functions class (can be substituted)
// -------------------------------------------------------------------

//#define debug
//#define pdebug

#include "G4QBesIKJY.hh"

// Constructor
G4QBesIKJY::G4QBesIKJY(G4QBIType type)
{
  ex=false;
  switch (type)
  {
    case BessI0:
      nu=0;
      ik=true;
      ij=true;
      break;
    case BessI1:
      nu=1;
      ik=true;
      ij=true;
      break;
    case EBessI0:
      nu=0;
      ex=true;
      ik=true;
      ij=true;
      break;
    case EBessI1:
      nu=1;
      ex=true;
      ik=true;
      ij=true;
      break;
    case BessJ0:
      nu=0;
      ik=true;
      ij=false;
      break;
    case BessJ1:
      nu=1;
      ik=true;
      ij=false;
      break;
    case BessK0:
      nu=0;
      ik=false;
      ij=true;
      break;
    case BessK1:
      nu=1;
      ik=false;
      ij=true;
      break;
    case EBessK0:
      nu=0;
      ex=true;
      ik=false;
      ij=true;
      break;
    case EBessK1:
      nu=1;
      ex=true;
      ik=false;
      ij=true;
      break;
    case BessY0:
      nu=0;
      ik=false;
      ij=false;
      break;
    case BessY1:
      nu=1;
      ik=false;
      ij=false;
      break;
  }
}

G4QBesIKJY::~G4QBesIKJY(){;}

G4double G4QBesIKJY::operator() (G4double X) const
{
  static const G4int nat1 = 15;            // a # of attempts to reach the X<1 accuracy
  static const G4int nat2 = nat1+nat1;     // a # of attempts to reach the X<5 accuracy
  static const G4int npi = 25;
  static const G4int npil = npi-1;
  static const G4int npk = 17;
  static const G4int npkl = npk-1;
  static const G4int npj = 20;
  static const G4int npjl = npj-1;
  static const G4complex CI(0,1);
  static const G4double EPS = 1.E-15;
  static const G4double Z1  = 1.;
  static const G4double HF  = Z1/2;
  static const G4double R8  = HF/4;
  static const G4double R32 = R8/4;
  static const G4double PI  = 3.14159265358979324;
  static const G4double CE  = 0.57721566490153286;
  static const G4double PIH = PI/2;
  static const G4double PI4 = PIH/2;   // PI/4
  static const G4double PI3 = PIH+PI4; // 3*PI/4
  static const G4double RPIH = 2./PI;
  static const G4double RPI2 = RPIH/4;

  static const G4double CI0[npi]={+1.00829205458740032E0,
             +.00845122624920943E0,+.00012700630777567E0,+.00007247591099959E0,
             +.00000513587726878E0,+.00000056816965808E0,+.00000008513091223E0,
             +.00000001238425364E0,+.00000000029801672E0,-.00000000078956698E0,
             -.00000000033127128E0,-.00000000004497339E0,+.00000000001799790E0,
             +.00000000000965748E0,+.00000000000038604E0,-.00000000000104039E0,
             -.00000000000023950E0,+.00000000000009554E0,+.00000000000004443E0,
             -.00000000000000859E0,-.00000000000000709E0,+.00000000000000087E0,
             +.00000000000000112E0,-.00000000000000012E0,-.00000000000000018E0};

  static const G4double CI1[npi]={+.975800602326285926E0,
           -.024467442963276385E0,-.000277205360763829E0,-.000009732146728020E0,
           -.000000629724238640E0,-.000000065961142154E0,-.000000009613872919E0,
           -.000000001401140901E0,-.000000000047563167E0,+.000000000081530681E0,
           +.000000000035408148E0,+.000000000005102564E0,-.000000000001804409E0,
           -.000000000001023594E0,-.000000000000052678E0,+.000000000000107094E0,
           +.000000000000026120E0,-.000000000000009561E0,-.000000000000004713E0,
           +.000000000000000829E0,+.000000000000000743E0,-.000000000000000080E0,
           -.000000000000000117E0,+.000000000000000011E0,+.000000000000000019E0};

  static const G4double CK0[npk]={+.988408174230825800E0,-.011310504646928281E0,
           +.000269532612762724E0,-.000011106685196665E0,+.000000632575108500E0,
           -.000000045047337641E0,+.000000003792996456E0,-.000000000364547179E0,
           +.000000000039043756E0,-.000000000004579936E0,+.000000000000580811E0, 
           -.000000000000078832E0,+.000000000000011360E0,-.000000000000001727E0,
           +.000000000000000275E0,-.000000000000000046E0,+.000000000000000008E0};

  static const G4double CK1[npk]={+1.035950858772358331E0,+.035465291243331114E0,
            -.000468475028166889E0,+.000016185063810053E0,-.000000845172048124E0,
            +.000000057132218103E0,-.000000004645554607E0,+.000000000435417339E0,
            -.000000000045757297E0,+.000000000005288133E0,-.000000000000662613E0,
            +.000000000000089048E0,-.000000000000012726E0,+.000000000000001921E0,
            -.000000000000000305E0,+.000000000000000050E0,-.000000000000000009E0};

  static const G4double CA[npk]={+.157727971474890120E0,-.008723442352852221E0,
          +.265178613203336810E0,-.370094993872649779E0,+.158067102332097261E0,
          -.034893769411408885E0,+.004819180069467605E0,-.000460626166206275E0,
          +.000032460328821005E0,-.000001761946907762E0,+.000000076081635924E0,
          -.000000002679253531E0,+.000000000078486963E0,-.000000000001943835E0,
          +.000000000000041253E0,-.000000000000000759E0,+.000000000000000012E0};

  static const G4double CB[npk]={-.021505111449657551E0,-.275118133043518791E0,
          +.198605634702554156E0,+.234252746109021802E0,-.165635981713650413E0,
          +.044621379540669282E0,-.006932286291523188E0,+.000719117403752303E0,
          -.000053925079722939E0,+.000003076493288108E0,-.000000138457181230E0,
          +.000000005051054369E0,-.000000000152582850E0,+.000000000003882867E0,
          -.000000000000084429E0,+.000000000000001587E0,-.000000000000000026E0};

  static const G4complex CC[npj]={
    G4complex(+.998988089858965153E0,-.012331520578544144E0),
    G4complex(-.001338428549971856E0,-.012249496281259475E0),
    G4complex(-.000318789878061893E0,+.000096494184993423E0),
    G4complex(+.000008511232210657E0,+.000013655570490357E0),
    G4complex(+.000000691542349139E0,-.000000851806644426E0),
    G4complex(-.000000090770101537E0,-.000000027244053414E0),
    G4complex(+.000000001454928079E0,+.000000009646421338E0),
    G4complex(+.000000000926762487E0,-.000000000683347518E0),
    G4complex(-.000000000139166198E0,-.000000000060627380E0),
    G4complex(+.000000000003237975E0,+.000000000021695716E0),
    G4complex(+.000000000002535357E0,-.000000000002304899E0),
    G4complex(-.000000000000559090E0,-.000000000000122554E0),
    G4complex(+.000000000000041919E0,+.000000000000092314E0),
    G4complex(+.000000000000008733E0,-.000000000000016778E0),
    G4complex(-.000000000000003619E0,+.000000000000000754E0),
    G4complex(+.000000000000000594E0,+.000000000000000462E0),
    G4complex(-.000000000000000010E0,-.000000000000000159E0),
    G4complex(-.000000000000000024E0,+.000000000000000025E0),
    G4complex(+.000000000000000008E0,+.000000000000000000E0),
    G4complex(-.000000000000000001E0,-.000000000000000001E0)};

  static const G4double CD[npk]={+0.648358770605264921E0,-1.191801160541216873E0,
         +1.287994098857677620E0,-0.661443934134543253E0,+0.177709117239728283E0,
         -0.029175524806154208E0,+0.003240270182683857E0,-0.000260444389348581E0,
         +0.000015887019239932E0,-0.000000761758780540E0,+0.000000029497070073E0,
         -0.000000000942421298E0,+0.000000000025281237E0,-0.000000000000577740E0,
         +0.000000000000011386E0,-0.000000000000000196E0,+0.000000000000000003E0};

  static const G4double EE[npk]={-.040172946544414076E0,-.444447147630558063E0,
          -.022719244428417736E0,+.206644541017490520E0,-.086671697056948524E0,
          +.017636703003163134E0,-.002235619294485095E0,+.000197062302701541E0,
          -.000012885853299241E0,+.000000652847952359E0,-.000000026450737175E0,
          +.000000000878030117E0,-.000000000024343279E0,+.000000000000572612E0,
          -.000000000000011578E0,+.000000000000000203E0,-.000000000000000003E0};

  static const G4complex CF[npj]={
    G4complex(+1.001702234853820996E0,+.037261715000537654E0),
     G4complex(+.002255572846561180E0,+.037145322479807690E0),
     G4complex(+.000543216487508013E0,-.000137263238201907E0),
     G4complex(-.000011179461895408E0,-.000019851294687597E0),
     G4complex(-.000000946901382392E0,+.000001070014057386E0),
     G4complex(+.000000111032677121E0,+.000000038305261714E0),
     G4complex(-.000000001294398927E0,-.000000011628723277E0),
     G4complex(-.000000001114905944E0,+.000000000759733092E0),
     G4complex(+.000000000157637232E0,+.000000000075476075E0),
     G4complex(-.000000000002830457E0,-.000000000024752781E0),
     G4complex(-.000000000002932169E0,+.000000000002493893E0),
     G4complex(+.000000000000617809E0,+.000000000000156198E0),
     G4complex(-.000000000000043162E0,-.000000000000103385E0),
     G4complex(-.000000000000010133E0,+.000000000000018129E0),
     G4complex(+.000000000000003973E0,-.000000000000000709E0),
     G4complex(-.000000000000000632E0,-.000000000000000520E0),
     G4complex(+.000000000000000006E0,+.000000000000000172E0),
     G4complex(+.000000000000000027E0,-.000000000000000026E0),
     G4complex(-.000000000000000008E0,-.000000000000000000E0),
     G4complex(+.000000000000000001E0,+.000000000000000001E0)};
  // -------------------------------------------------------------------------------------
  G4double H=0.;                     // Prototype of the result
  if (ij)                            // I/K Bessel functions
  {
    if (ik)                          // I0/I1/EI0/EI1 Bessel functions (symmetric)
    {
      G4double V=std::abs(X);
      G4double CJ=0.;                // Prototype of the element of the CI0/CI1 matrix
      if (V < 8.)
      {
        G4double HFV=HF*V;
        G4double Y=HFV*HFV;
        G4int    V3=nu+1;
        G4int    XL=V3+1;
        G4int    XLI=XL+1;
        G4int    XLD=XLI+1;
        G4int    W1=XLD+1;
        G4double A0=1.;
        G4double DY=Y+Y;
        G4double A1=1.+DY/(XLI*V3);
        G4double A2=1.+Y*(4.+(DY+Y)/(XLD*XL))/(W1*V3);
        G4double B0=1.;
        G4double B1=1.-Y/XLI;
        G4double B2=1.-Y*(1.-Y/(XLD+XLD))/W1;
        G4int    V1=3-XL;
        G4double V2=V3+V3;
        G4double C=0.;
        for (G4int N=3; N<=30; N++)
        {
          G4double C0=C;
          G4double FN=N;
                   W1=W1+2;
          G4int    W2=W1-1;
          G4int    W3=W2-1;
          G4int    W4=W3-1;
          G4int    W5=W4-1;
          G4int    W6=W5-1;
                   V1=V1+1;
                   V2=V2+1;
                   V3=V3+1;
          G4double U1=FN*W4;
          G4double E=V3/(U1*W3);
          G4double U2=E*Y;
          G4double F1=1.+Y*V1/(U1*W1);
          G4double F2=(1.+Y*V2/(V3*W2*W5))*U2;
          G4double F3=-Y*Y*U2/(W4*W5*W5*W6);
          G4double A=F1*A2+F2*A1+F3*A0;
          G4double B=F1*B2+F2*B1+F3*B0;
                   C=A/B;
          if (std::abs(C0-C) < EPS*std::abs(C)) break;
          A0=A1; A1=A2; A2=A; B0=B1; B1=B2; B2=B;
        }
        H=C;
        if (nu==1) H*=HF*X;
        if (ex) H*=std::exp(-V);
      }
      else
      {
        G4double P=16./V-1.;
        G4double ALFA=P+P;
        G4double B1=0.;
        G4double B2=0.;
        for (G4int I=npil; I>=0; I--)
        {
          if (!nu) CJ=CI0[I];
          else     CJ=CI1[I];
          G4double B0=CJ+ALFA*B1-B2;
                   B2=B1;
                   B1=B0;
        }
        H=std::sqrt(RPI2/V)*(B1-P*B2);
        if (nu && X < 0.) H=-H;
        if (!ex) H*=std::exp(V);
      }
    }
    else                             // K0/K1/EK0/EK1 Bessel functions
    {
#ifdef debug
      G4cout<<"G4BesIKJY: >>>>>>>>>>>>>> K is called, X="<<X<<",n="<<nu<<",E="<<ex<<G4endl;
#endif
      G4double CK=0.;                // Prototype of the element of the CI0/CI1 matrix
      if (X < 0.)
      {
        G4cout<<"G4BesIKJY::NegativeArg in K-BessFun X="<<X<<", n="<<nu<<",E="<<ex<<G4endl;
        return H;
      }
      else if (X < 1.)
      {
#ifdef debug
        G4cout<<"G4BesIKJY: >>>> [ X < 1 ] is called, X="<<X<<",n="<<nu<<",E="<<ex<<G4endl;
#endif
        G4double B=HF*X;
        G4double BK=-(std::log(B)+CE);
        G4double F=BK;
        G4double P=HF;
        G4double Q=HF;
        G4double C=1.;
        G4double D=B*B;
        G4double BK1=P;
        for (G4int N=1; N<=nat1; N++)  // @@ "nat" can be increased
        {
          G4double FN=N;
                   P/=FN;
                   Q/=FN;
                   F=(F+P+Q)/FN;
                   C*=D/FN;
          G4double G=C*(P-FN*F);
          G4double R=C*F;
                   BK=BK+R;
                   BK1=BK1+G;
          if (BK1*R+std::abs(G)*BK < EPS*BK*BK1) break;
        }
        if (nu==1) H=BK1/B;
        else       H=BK;
        if (ex) H*=std::exp(X);
      }
      else if (X < 5.)
      {
#ifdef debug
        G4cout<<"G4BesIKJY: >>>> [ X < 5 ] is called, X="<<X<<",n="<<nu<<",E="<<ex<<G4endl;
#endif
        G4int NUS=0;              // @@ NU**2 for future NU>1 applications
        if (nu==1) NUS=1;
        G4double DNUS=NUS+NUS;
        G4double XN=DNUS+DNUS;
        G4double A=9.-XN;
        G4double B=25.-XN;
        G4double C=768*X*X;
        G4double HX=16*X;
        G4double C0=HX+HX+HX;;
        G4double A0=1.;
        G4double A1=(HX+7.+XN)/A;
        G4double A2=(C+C0*(XN+23.)+XN*(XN+62.)+129.)/(A*B);
        G4double B0=1.;
        G4double B1=(HX+9.-XN)/A;
        G4double B2=(C+C0*B)/(A*B)+1.;
                 C=0.;
        for (G4int N=3; N<=nat2; N++)
        {
          C0=C;
          G4double FN=N;
          G4double FN2=FN+FN;
          G4double FNP=FN2+1.;
          G4double FN1=FN2-1.;
          G4double FNM=FN1-4.;
          G4double FN3=FN1/(FN2-3.);
          G4double FN4=12*FN*FN-(1.-XN);
          G4double FN5=16*FN1*X;
          G4double RAN=1./(FNP*FNP-XN);
          G4double F1=FN3*(FN4-20*FN)+FN5;
          G4double F2=28*FN-FN4-8.+FN5;
          G4double F3=FN3*(FNM*FNM-XN);
                   A=(F1*A2+F2*A1+F3*A0)*RAN;
                   B=(F1*B2+F2*B1+F3*B0)*RAN;
                   C=A/B;
          if (std::abs(C0-C) < EPS*std::abs(C)) break;
          A0=A1; A1=A2; A2=A; B0=B1; B1=B2; B2=B;
        }
        H=C/std::sqrt(RPIH*X);
        if (!ex) H*=std::exp(-X);
      }
      else
      {
#ifdef debug
        G4cout<<"G4BesIKJY: >>> [ X >= 5 ] is called, X="<<X<<",n="<<nu<<",E="<<ex<<G4endl;
#endif
        G4double P=10./X-1.;
        G4double ALFA=P+P;
        G4double B1=0.;
        G4double B2=0.;
#ifdef debug
        G4cout<<"G4BesIKJY: >>> [ X >= 5 ] is called, X="<<X<<",n="<<nu<<",E="<<ex<<G4endl;
#endif
        for (G4int I=npkl; I>=0; I--)
        {
          if (!nu) CK=CK0[I];
          else     CK=CK1[I];
          G4double B0=CK+ALFA*B1-B2;
                   B2=B1;
                   B1=B0;
        }
        H=std::sqrt(PIH/X)*(B1-P*B2);
        if (!ex) H*=std::exp(-X);
      }
    }
  }
  else
  {
    if (!ik && X < 0.)
    {
      G4cout<<"G4BesIKJY::NegativeArgument in Y BesselFunction X="<<X<<", nu="<<nu<<G4endl;
      return H;
    }
    G4double V=std::abs(X);
    if (!nu)                          // J0/Y0 Bessel functions
    {
      if (V < 8.)
      {
        G4double P=R32*V*V-1.;
        G4double ALFA=P+P;
        G4double B1=0.;
        G4double B2=0.;
        for (G4int IT=npkl; IT>=0; IT--)
        {
          G4double B0=CA[IT]+ALFA*B1-B2;
                   B2=B1;
                   B1=B0;
        }
        H=B1-P*B2;
        if (!ik)
        {
          B1=0.;
          B2=0.;
          for (G4int JT=npkl; JT>=0; JT--)
          {
            G4double B0=CB[JT]+ALFA*B1-B2;
                     B2=B1;
                     B1=B0;
          }
          H=RPIH*(CE+std::log(HF*X))*H+B1-P*B2;
        }
      }
      else
      {
        G4double P=10./V-1.;
        G4double ALFA=P+P;
        G4complex CB1(0.,0.);
        G4complex CB2(0.,0.);
        for (G4int IT=npjl; IT>=0; IT--)
        {
          G4complex CB0=CC[IT]+ALFA*CB1-CB2;
                    CB2=CB1;
                    CB1=CB0;
        }
        CB1=std::sqrt(RPIH/V)*std::exp(CI*(V-PI4))*(CB1-P*CB2);
        if (ik) H=real(CB1);
        else    H=real(-CI*CB1);
      }
    }
    else                          // J1/Y1 Bessel functions
    {
      if (V < 8.)
      {
        G4double Y=R8*V;
        G4double Y2=Y*Y;
        G4double P=Y2+Y2-1.;
        G4double ALFA=P+P;
        G4double B1=0.;
        G4double B2=0.;
        for (G4int IT=npkl; IT>=0; IT--)
        {
          G4double B0=CD[IT]+ALFA*B1-B2;
                   B2=B1;
                   B1=B0;
        }
        H=Y*(B1-P*B2);
        if (!ik)
        {
          B1=0.;
          B2=0.;
          for (G4int JT=npkl; JT>=0; JT--)
          {
            G4double B0=EE[JT]+ALFA*B1-B2;
                     B2=B1;
                     B1=B0;
          }
          H=RPIH*((CE+std::log(HF*X))*H-1./X)+Y*(B1-B2);
        }
      }
      else
      {
        G4double P=10./V-1.;
        G4double ALFA=P+P;
        G4complex CB1(0.,0.);
        G4complex CB2(0.,0.);
        for (G4int IT=npjl; IT>=0; IT--)
        {
          G4complex CB0=CF[IT]+ALFA*CB1-CB2;
                    CB2=CB1;
                    CB1=CB0;
        }
        CB1=std::sqrt(RPIH/V)*std::exp(CI*(V-PI3))*(CB1-P*CB2);
        if (ik) H=real(CB1);
        else    H=real(-CI*CB1);
      }
      if (X < 0.) H=-H;
    }
  }
  return H;
}
