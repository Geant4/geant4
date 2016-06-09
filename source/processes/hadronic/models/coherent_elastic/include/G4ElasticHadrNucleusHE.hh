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
// $Id: G4ElasticHadrNucleusHE.hh,v 1.33 2006/12/13 15:45:22 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
// G4ElasticHadrNucleusHe.hh

//  The generator of high energy hadron-nucleus elastic scattering
//  The hadron kinetic energy T > 1 GeV
//  N.  Starkov 2003.
//
//  19.05.04 Variant for G4 6.1: The 'ApplyYourself' was changed
//  19.11.05 The HE elastic scattering on proton is added (N.Starkov)
//  16.11.06 General redesign (N.Starkov)
//  23.11.06 General cleanup, ONQ0=3 (V.Ivanchenko)
//

#ifndef G4ElasticHadrNucleusHE_h
#define G4ElasticHadrNucleusHE_h 1

#include <vector>

#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleChange.hh"
#include "G4Nucleus.hh"

#include "G4HadronicInteraction.hh"

static const G4int  ONQ0     = 3;   //  The initial number of steps on Q2
static const G4int  ONQ2     = 150; //  The total number of steps on Q2
static const G4int  NENERGY  = 30;  
static const G4int  NQTABLE  = NENERGY*ONQ2;  

class ElasticData
{
public:

  ElasticData(const G4ParticleDefinition* h, 
	      G4int AtomWeight);

  ~ElasticData(){;}

  void GetNucleusParameters(G4int Nucleus);
  const G4ParticleDefinition* Hadron() {return hadr;}

  //  ============================================
  void   fillQ2limit();

  const G4ParticleDefinition*  hadr;
  G4int     AtomicWeight;
  G4int     dnkE[NENERGY];
  G4double  TableQ2[ONQ2];
  G4double  TableCrossSec[NQTABLE];
  G4double  maxQ2, dQ2;
  G4double  R1, R2, Pnucl, Aeff;
  G4double  massGeV;
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
  G4double GeV2;
  G4double emin, emax, deltae;

  G4int     HadrCode;
  //  G4String  HadronName;
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
