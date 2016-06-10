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
// $Id: G4ElasticHadrNucleusHE.hh 90756 2015-06-09 07:43:33Z gcosmo $
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
//  15.05.07 Redesign and cleanup (V.Ivanchenko)
//  18.05.07 Cleanup (V.Grichine)
//  19.04.12 Fixed reproducibility violation (A.Ribon)
//  12.06.12 Fixed warnings of shadowed variables (A.Ribon)
//

#ifndef G4ElasticHadrNucleusHE_h
#define G4ElasticHadrNucleusHE_h 1

#include <vector>

#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleChange.hh"
#include "G4Nucleus.hh"
#include "G4HadronElastic.hh"

class G4NistManager;

static const G4int  NHADRONS = 26;  //  Number of hadrons for which model is applied
static const G4int  ONQ0     = 5;   //  The initial number of steps on Q2
static const G4int  ONQ2     = 100; //  The total number of steps on Q2
static const G4int  NENERGY  = 30;  
static const G4int  ZMAX     = 93;  
static const G4int  NQTABLE  = NENERGY*ONQ2;  


///////////////////////////////////////////////////////////////////////
//
//

class G4ElasticData
{
public:

  G4ElasticData(const G4ParticleDefinition* h, 
		G4int Z, G4double A, G4double* eGeV);

  ~G4ElasticData(){}

  const G4ParticleDefinition* Hadron() {return hadr;}

private:
  void DefineNucleusParameters(G4double A);
  const G4ParticleDefinition*  hadr;

  // hide assignment operator
  G4ElasticData & operator=(const G4ElasticData &right);
  G4ElasticData(const G4ElasticData&);

public:
  G4int     AtomicWeight;
  G4double  R1, R2, Pnucl, Aeff;
  G4double  limitQ2;
  G4double  massGeV;
  G4double  mass2GeV2;
  G4double  massA;
  G4double  massA2;
  G4int     dnkE[NENERGY];
  G4double  maxQ2[NENERGY];

  G4double  TableQ2[ONQ2];
  G4double  TableCrossSec[NQTABLE];
};

/////////////////////////////////////////////////////////////////////
//
//

class G4ElasticHadrNucleusHE : public G4HadronElastic
{
public:

  G4ElasticHadrNucleusHE(const G4String& name = "hElasticGlauber");

  virtual ~G4ElasticHadrNucleusHE();

  virtual G4double SampleInvariantT(const G4ParticleDefinition* p, 
				    G4double plab, 
				    G4int Z, G4int A);

  virtual void ModelDescription(std::ostream&) const;

  G4double SampleT(const G4ParticleDefinition* p, 
		   G4double plab, 
		   G4int Z, G4int A);

  G4double HadronNucleusQ2_2(G4ElasticData * pElD, G4int Z, 
			     G4double plabGeV, G4double tmax);

  void DefineHadronValues(G4int Z);

  G4double GetLightFq2(G4int Z, G4int A, G4double Q);
  G4double GetHeavyFq2(G4int Z, G4int Nucleus, G4double *LineFq2); 

  G4double GetQ2_2(G4int  N, G4double * Q, 
		   G4double * F, G4double R);

  G4double LineInterpol(G4double p0, G4double p2,
			G4double c1, G4double c2,
			G4double p);

  G4double HadrNucDifferCrSec(G4int Z, G4int Nucleus, G4double Q2); 

  void InterpolateHN(G4int n, const G4double EnP[], 
		     const G4double C0P[], const G4double C1P[], 
		     const G4double B0P[], const G4double B1P[]);

  // hide assignment operator
  G4ElasticHadrNucleusHE & operator=(const G4ElasticHadrNucleusHE &right);
  G4ElasticHadrNucleusHE(const G4ElasticHadrNucleusHE&);

  G4double GetBinomCof( G4int n, G4int m );

  G4double GetFt(G4double Q2);

  G4double GetDistrFun(G4double Q2);

  G4double GetQ2(G4double Ran);

  G4double HadronProtonQ2(const G4ParticleDefinition * aHadron,
                          G4double inLabMom);

  void     GetKinematics(const G4ParticleDefinition * aHadron,
		          G4double MomentumH);
private:

  void     Binom();

  //  fields

  G4int    iHadrCode;
  G4int    iHadron;
  G4int    HadronCode[NHADRONS];
  G4int    HadronType[NHADRONS];
  G4int    HadronType1[NHADRONS];

  // protection energy and momemtum

  G4double lowestEnergyLimit;  
  G4double plabLowLimit;
  G4double dQ2;  

  // transition between internal and CLHEP units

  G4double MbToGeV2;
  G4double sqMbToGeV;
  G4double Fm2ToGeV2;
  G4double GeV2;
  G4double protonM;     // GeV
  G4double protonM2;    // GeV^2

  // projectile kinematics in GeV

  G4double hMass;
  G4double hMass2;
  G4double hLabMomentum;
  G4double hLabMomentum2;
  G4double MomentumCM;
  G4double HadrEnergy;

  // nucleaus parameters

  G4double  R1, R2, Pnucl, Aeff;
  G4int     NumbN;

  // elastic parameters

  G4double  HadrTot, HadrSlope, HadrReIm, TotP, 
            DDSect2, DDSect3, ConstU, FmaxT;

  // momentum limits for different models of hadron/nucleon scatetring
  G4double BoundaryP[7], BoundaryTL[7], BoundaryTG[7];

  // parameterisation of scattering

  G4double Slope1, Slope2, Coeff1, Coeff2, MaxTR;
  G4double Slope0, Coeff0;

  G4double aAIm, aDIm, Dtot11;

  G4double        Energy[NENERGY];
  G4double        LowEdgeEnergy[NENERGY];

  G4double        SetBinom[240][240];

  static G4ElasticData*  SetOfElasticData[NHADRONS][ZMAX];
  G4NistManager*  nistManager;
  static G4Mutex eldata_m[NHADRONS][ZMAX];

};     //   The end of the class description

////////////////////////////////////////////////////////////////

inline
G4double G4ElasticHadrNucleusHE::LineInterpol(G4double p1, G4double p2,
					      G4double c1, G4double c2,
					      G4double p)
{
//  G4cout<<"  LineInterpol: p1 p2 c1 c2 "<<p1<<"  "<<p2<<"  "
//        <<c1<<"  "<<c2<<"  c  "<<c1+(p-p1)*(c2-c1)/(p2-p1)<<G4endl;


  return c1+(p-p1)*(c2-c1)/(p2-p1);
}

////////////////////////////////////////////////////////////////

inline
void G4ElasticHadrNucleusHE::InterpolateHN(G4int n, const G4double EnP[], 
                   const G4double C0P[], const G4double C1P[], 
                   const G4double B0P[], const G4double B1P[])
{
  G4int i; 

  for(i=1; i<n; i++) if(hLabMomentum <= EnP[i]) break;
 
  if(i == n) i = n - 1;

  Coeff0 = LineInterpol(EnP[i], EnP[i-1], C0P[i], C0P[i-1], hLabMomentum);
  Coeff1 = LineInterpol(EnP[i], EnP[i-1], C1P[i], C1P[i-1], hLabMomentum);
  Slope0 = LineInterpol(EnP[i], EnP[i-1], B0P[i], B0P[i-1], hLabMomentum);
  Slope1 = LineInterpol(EnP[i], EnP[i-1], B1P[i], B1P[i-1], hLabMomentum);

//  G4cout<<"  InterpolHN:  n i "<<n<<"  "<<i<<"  Mom "
//        <<hLabMomentum<<G4endl;
}

////////////////////////////////////////////////////////////////

inline
G4double G4ElasticHadrNucleusHE::GetBinomCof( G4int numN, G4int numM )
{
  if ( numN >= numM && numN <= 240) return SetBinom[numN][numM];
  else                     return 0.;
}

////////////////////////////////////////////////////////////////

inline
G4double G4ElasticHadrNucleusHE::GetDistrFun(G4double Q2)
{
  return GetFt(Q2)/FmaxT;
}

#endif
