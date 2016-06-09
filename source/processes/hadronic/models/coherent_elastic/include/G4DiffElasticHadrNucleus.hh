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
// $Id: G4DiffElasticHadrNucleus.hh,v 1.15 2006/06/29 20:08:59 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//

//  High energy hadron-nucleus elastic scattering
//  Kinetic energy T > 1 GeV
//  N. Starkov 2003.
//
//  Modifications:
//  14.11.05 The HE elastic scattering on proton is added (N.Starkov)
//  30.05.06 New version, which not uses elastic data (N.Starkov)
//  30.05.06 cleanup (V.Ivanchenko)
//

#ifndef G4DiffElasticHadrNucleus_h
#define G4DiffElasticHadrNucleus_h 1

#include "G4Nucleus.hh"
#include "G4DynamicParticle.hh"
#include "G4HadronValues.hh"
#include "G4IntegrHadrNucleus.hh"

//  ###############################################
class NucleusParameters
{

public:

  NucleusParameters() {;}
  virtual ~NucleusParameters(){;}    

  NucleusParameters & operator= (const NucleusParameters &t) 
  {
    if(this!=&t)
      {
	R1    = t.R1;
	R2    = t.R2;
	Aeff  = t.Aeff;
	Pnucl = t.Pnucl;
      }
    return *this;
  }

  void GetNucleusParameters(G4int AtWeight);
   
  G4double R1, R2, Aeff, Pnucl;
};

//  ################################################ 
class G4DiffElasticHadrNucleus: public //G4HadronValues
                                      G4IntegrHadrNucleus
{
    
public:
  G4DiffElasticHadrNucleus() : // G4HadronValues(),
    G4IntegrHadrNucleus() 
  {
    Factors();
    Binom();
    r0       = 1.1;   //  The WS's parameters
    r01      = 1.16;
    rAmax    = 2.5;
    NpointsB = 500;
    FmToGeV  = 5.068;
    Fm2ToGeV2= FmToGeV*FmToGeV;
  }

  virtual ~G4DiffElasticHadrNucleus()   {;}

  G4double HadrNuclDifferCrSec(
                          const  G4DynamicParticle *  aHadron,
                                 G4Nucleus         *  aNucleus, 
                                 G4double            Q2);
    
  G4double HadrNuclDifferCrSecT(
                           const  G4DynamicParticle *  aHadron,
                                  G4Nucleus         *  aNucleus, 
                                  G4double            Q2,
                                  G4int                Kind);
    
  G4double Thickness(G4int A, G4double b);

  void     GetIntegrandB(G4int  Anucleus);
 
  G4double GetIntegrandS(G4int  Anucleus, G4double b, G4int Kind);

  void     GetNucleusParameters(G4Nucleus   * aNucleus);
  
  G4double HadronProtonDiffCrSec(G4double Q2);

  void     GetKinematics(const G4DynamicParticle * aHadron);

  void     GetParametersHP(const G4DynamicParticle * aHadron);

  G4double LineInterpol(G4double p0, G4double p1,
			G4double c0, G4double c1, 
			G4double p);
//  ====================================================
  G4double MyJ0(G4double x)
  {
    G4double x1, x2, x3, x4, x5, x6, x7, f0, th0, res;

    x2 = x*x/9.0;

    if(std::fabs(x)<=3.0)
      {
         x3 = x2*x2;
         x4 = x3*x2;
         x5 = x4*x2;
         x6 = x5*x2;
         x7 = x6*x2;
        res = 1.0-2.2499997*x2+1.2656208*x3-
                  0.3163866*x4+0.0444479*x5-
                  0.0039444*x6+0.0002100*x7;
      }
    else
      {
         x1 = x/3.0;
         x3 = x2*x1;
         x4 = x3*x1;
         x5 = x4*x1;
         x6 = x5*x1;
         f0 =    0.7978456    -0.00000077/x-
                 0.00552740/x2-0.00009512/x3+
                 0.00137237/x4-0.00072805/x5+
	         0.00014476/x6;
         th0 = x-0.78539816   -0.04166397/x1-
                 0.00003954/x2+0.00262573/x3-
                 0.00054125/x4-0.00029333/x5+
                 0.00013558/x6;
	 res = f0*std::cos(th0)/std::sqrt(x);
      }
    return res;
  }
//  +++++++++++++++++++++++++++++++++++++++++++++
  G4double MyI0(G4double x)
  {
     G4double p1=1.0,       p2=3.5156229, p3=3.0899424,
              p4=1.2067492, p5=0.2659732, p6=3.360768e-2,
              p7=4.5813e-3;
     G4double q1=0.39894228,  q2=1.328592e-2,  q3=2.25319e-3,
              q4=-1.57565e-3, q5=9.16281e-3,   q6=-2.057706e-2,
              q7=2.635537e-2, q8=-1.647633e-2, q9=3.922377e-3;

     G4double k1=3.75, ax=std::fabs(x), y=0.0, result=0.0;
 
     if(ax<k1)
       {
	 G4double xx = x/k1;
	 y = xx*xx;
	 result = p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))));
       } 
     else 
       {
	 y = k1/ax;
	 result =  std::exp(ax)/std::sqrt(ax)*(q1+y*(q2+y*(q3+y*(q4+
		   y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))));
       }
     return result;
  }

//  private:

// ++++++++++++++++++++++++++++++++++++++++++++++++++
  void Binom()
  {       
    G4int N, M;
    G4double  Fact1, J;

    for(N=0; N<240; N++)
      {
	J = 1;
	for(M=0; M<=N; M++)
	  {
            Fact1 = 1;
            if ((N>0) && (N>M) && M>0) 
	      {
		J = J*(N-M+1)/M;
		Fact1 = J;
	      }
            SetBinom[N][M] = Fact1;
// G4cout<<" N M "<<N<<"  "<<M<<" F  "<<Fact1<<G4endl;
          }
      }
    return;
  }
// +++++++++++++++++++++++++++++++++++++++++++++++++++
  G4double binom(G4int N, G4int M)
  {
    G4double Fact1 = 1;

    if((N>1) & (N>=M))
      Fact1 = Factorials[N]/Factorials[M]/Factorials[N-M];
    return Fact1;
  }
// +++++++++++++++++++++++++++++++++++++++++++++++++++
  G4double Factorial(G4int N)
  {
    G4double  Res=1;

    if(N == 0) return Res;

    if(N < 100) for(G4int M = 1; M<=N; M++)  Res = Res*M;         

    if(N >= 100)  Res = 2.50662827*
		    std::exp(static_cast<double>(-N-1))*
		    std::pow(static_cast<double>(N+1),N+0.5)*
		    (1.+1./12./(N+1)+1./288./(N+1)/(N+1)-
		     139./51840./(N+1)/(N+1)/(N+1)-
		     571./2488320./(N+1)/(N+1)/(N+1)/(N+1));
    return Res;
  }
// +++++++++++++++++++++++++++++++++++++++++++++++++++
  void Factors()
  {
    G4int ii, ll, mm;
    G4double Sum1, Fac1, Fac3, Sum3;

    Factorials[0] = 1;

    for( ii = 1; ii<140; ii++)
      {
	Factorials[ii] = Factorial(ii);
	if(ii >= 100) Mnoj[ii] = 3.03;  //  there is the saturation
	else
	  {
	    Sum1 = 0;
	    Fac1 = 1;
	    for( ll = 0; ll<=ii; ll++)
	      {
		Fac1 = binom(ii,ll);
		Fac3 = 1;
		Sum3 = 1;
		for( mm = 1; mm<=ii-ll; mm++)
		  {
		    Fac3 = binom(ii-ll,mm);
		    Sum3 = Sum3 + 1/Fac3;
		  }  //  mm
		Sum1 = Sum1 + Sum3/Fac1;
	      }      //  ll
	    Mnoj[ii] = Sum1;
	  }       //  else
      }           //  ii
    Mnoj[0] = 1;
  }   //   Factors
//  --------------------------------------------------------

public:

  G4double Fm2ToGeV2, FmToGeV;

  G4double  SigTot, Slope, ReOnIm, HdrE, ConstU, Q2;
  G4double  ProtonM, HadronM, HadronE, Sh, PM2, HM2;

  G4double  EcmP, EcmH, Kcm, CosCM, SqrtS, MaxT, FmaxT, DsigDt;
  G4double  Slope1, Slope2, Coeff1, Coeff2, IntConst, MaxTR;
  G4double  Slope0, Coeff0;
  G4int     NumbPointsB, NumbPointsT, HadrCode;
  G4String  HadronName;

  G4double  Mnoj[250], Factorials[250];
  G4double  ReIntegrand[1000], ImIntegrand[1000], Thick[1000];
  G4double  rAfm, rAGeV, stepB, MomentumCMN, MomentumH;

  G4double  R1, R2, Pnucl, Aeff, AIm, ARe, Dem, DIm, InCoh, InCohI;
  G4double  r0, r01, rAmax;
  G4double  BoundaryP[7], BoundaryTL[7], BoundaryTG[7];
  G4int     NpointsB, HadrCodes[7], Intern;
  G4double  SetBinom[240][240];
};

#endif
