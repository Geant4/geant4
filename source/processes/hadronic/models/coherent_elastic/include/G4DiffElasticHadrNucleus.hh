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
// $Id: G4DiffElasticHadrNucleus.hh,v 1.13 2005/11/25 19:17:32 dennis Exp $
// GEANT4 tag $Name: geant4-08-00 $
//

//
//  G4DiffElasticHadrNucleus header file
//
//  High energy hadron-nucleus elastic scattering
//  Kinetic energy T > 1 GeV
//  N.  Starkov 2003.
//
//  Modifications:
//  14.11.05 The HE elastic scattering on proton is added (N.Starkov)
//  23.11.05 int -> G4int, fabs -> abs (V.Ivanchenko)
//

#ifndef G4DiffElasticHadrNucleus_h
#define G4DiffElasticHadrNucleus_h 1

#include "G4Nucleus.hh"
#include "G4DynamicParticle.hh"
#include "G4HadronValues.hh"
#include "G4IntegrHadrNucleus.hh"

class G4DiffElasticHadrNucleus: public G4IntegrHadrNucleus
{
    
public:
  
  G4DiffElasticHadrNucleus();

  ~G4DiffElasticHadrNucleus();

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

     if(std::abs(x)<=3.0)
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

     G4double k1=3.75, ax=std::abs(x), y=0.0, result=0.0;
 
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

private:

  void Factors();

// ++++++++++++++++++++++++++++++++++++++++++++++++++
  G4double binom(G4int N, G4int M)
       {       
          G4double  Fact1 = 1;
       if ((N>1) & (N>=M)) 
          {
            Fact1 = Factorials[N]/Factorials[M]/
                      Factorials[N-M];
          }
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
                         (1+1/12/(N+1)+1/288/(N+1)/(N+1)-
                         139/51840/(N+1)/(N+1)/(N+1)-
                         571/2488320/(N+1)/(N+1)/(N+1)/(N+1));

      return Res;
    }

//  --------------------------------------------------------
public:

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
};

#endif
