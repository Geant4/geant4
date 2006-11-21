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
// $Id: G4ElasticHadrNucleusHE.hh,v 1.30 2006-11-21 19:38:53 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4ElasticHadrNucleusHe.hh

//  The generator of high energy hadron-nucleus elastic scattering
//  The hadron kinetic energy T > 1 GeV
//  N.  Starkov 2003.

//  19.05.04 Variant for G4 6.1: The 'ApplyYourself' was changed

//  November 2005 - The HE elastic scattering on proton is added
//  N. Starkov

#ifndef G4ElasticHadrNucleusHE_h
#define G4ElasticHadrNucleusHE_h 1

#include <vector>

#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleChange.hh"
#include "G4Nucleus.hh"

#include "G4HadronicInteraction.hh"

static const G4int  ONQ0     = 5;   //  The initial number of steps on Q2
static const G4int  ONQ2     = 150; //  The total number of steps on Q2
static const G4int  ONE      = 5;   //  The number of steps on E
static const G4int  AreaNumb = 6;   //  The number of order steps on E
static const G4int  ONQ2XE   = ONQ2*ONE; //  The dimension of a distr. func. array
static const G4int  MaxN     = 10;  
static const G4int  NENERGY  = 30;  
static const G4int  NQTABLE  = NENERGY*ONQ2;  

class ElasticData
{
public:

  ElasticData(const G4ParticleDefinition* h, 
	      G4int AtomWeight,
	      G4double* en);

  ~ElasticData(){;}

  void GetNucleusParameters(G4int Nucleus);
  const G4ParticleDefinition* Hadron() {return hadr;}

  //  ============================================
  void   fillQ2limit();
  G4int  GetNumberE(G4double E);

  const G4ParticleDefinition*  hadr;
  G4int     AtomicWeight;
  G4int     dnkE[NENERGY];
  G4double  TableE[NENERGY];
  G4double  TableQ2[ONQ2];
  G4double  TableCrossSec[NQTABLE];
  G4double  maxQ2, dQ2;
  G4double  R1, R2, Pnucl, Aeff;

};

//  ############################################################
class G4ElasticHadrNucleusHE : public G4HadronicInteraction
{
public:

  G4ElasticHadrNucleusHE();

  virtual ~G4ElasticHadrNucleusHE();

  G4HadFinalState * ApplyYourself(const G4HadProjectile &aTrack,
				  G4Nucleus       &G4Nucleus);

  G4double SampleT(const G4ParticleDefinition* p,
		   G4double pTotLabMomentum, G4int Z, G4int N);

private:

  G4double InterPol(G4double X1, G4double X2, G4double X3,
		    G4double Y1, G4double Y2, G4double Y3, 
		    G4double X);

  G4double HadronNucleusQ2_2(const G4ParticleDefinition * aHadron,
			     G4int                  AWeight,
			     G4double   q, G4double aLabMom,
			     G4int kk, ElasticData * pElD);
  G4double GetLightFq2(G4int N, G4double Q, G4int Step);

  G4int    GetBinom(G4int m, G4int n);


//  ======================================================
  G4DynamicParticle  aHad;

  G4float  GetFt(G4double T);
  G4float  GetDistrFun(G4double Q2);
  G4double GetQ2(G4double Ran);
  G4double GetQ2_2(G4int  N, G4double * Q, 
		   G4double * F, G4double R);
  G4double HadronProtonQ2(const G4ParticleDefinition * aHadron,
			  G4double TotMom);

  void     GetKinematics(const G4ParticleDefinition * aHadron,
			 G4double TotMom);

  void     GetParametersHP(const G4ParticleDefinition * aHadron,
			   G4double TotMom);

  G4double LineInterpol(G4double p0, G4double p2,
			G4double c1, G4double c2,
			G4double p);

  G4double SigTot, Slope, ReOnIm, HdrE, ConstU, Q2;
  G4double ProtonM, HadronM, HadronE, Sh, PM2, HM2;

  G4double EcmP, EcmH, Kcm, CosCM, SqrtS, MaxT, FmaxT, DsigDt;
  G4double Slope1, Slope2, Coeff1, Coeff2, IntConst, MaxTR;
  G4double Slope0, Coeff0;
  G4double BoundaryP[7], BoundaryTL[7], BoundaryTG[7];
  G4double MomentumH;
  G4int    HadrCodes[7];
  G4int    verboselevel;
//  ======================================================
  G4double MbToGeV2;
  G4double sqMbToGeV;
  G4double Fm2ToGeV2;
  G4double emin, emax, deltae;

  G4int     HadrCode;
  G4String  HadronName;
  G4double  R1, R2, Pnucl, Aeff;
  G4double  Energy[NENERGY];
// +++++++++++++++++++++++++++++++++++++++++++++++++++
  std::vector<ElasticData*> SetOfElasticData;

  G4int     Nstep,         // The number of steps on Q2
            iKindWork,     // 
            iContr,        //
            iPoE;          // The number of steps on E
  G4int     iTypeWork, CurrentN;
  G4double  aNucleon;

  G4double  dEbeg1, dEend1, dQ2, maxQ2;
  G4double  SetBinom[240][240];

//  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
          }
        }
      return;
    }

  void GetHadronValues(const G4ParticleDefinition * aHadron,
                                  G4double TotMom);

  G4double  HadrTot, HadrSlope, HadrReIm,  DDSect2, DDSect3, MomentumCM;

  G4double  Q2res;

};     //   The end of the class description

#endif
