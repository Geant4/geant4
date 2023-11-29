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
// G4ElasticHadrNucleusHe.hh

//  The generator of high energy hadron-nucleus elastic scattering
//  The hadron kinetic energy T > 1 GeV
//  N.Starkov 2003.
//
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
#include <iostream>
#include <fstream>

#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleChange.hh"
#include "G4Nucleus.hh"
#include "G4HadronElastic.hh"
#include "G4Threading.hh"

class G4NistManager;

static const G4int  NHADRONS = 26;     // Number of allowed hadrons 
static const G4int  ONQ2     = 102;    // Number of points on Q2
static const G4int  NENERGY  = 24;
static const G4int  ZMAX     = 93;

///////////////////////////////////////////////////////////////////////

class G4ElasticData
{

friend class G4ElasticHadrNucleusHE;

public:

  G4ElasticData(const G4ParticleDefinition* h, G4int Z, G4int A, 
                const G4double* e);

  ~G4ElasticData() {}

private:

  void DefineNucleusParameters(G4int A);

  // hide assignment operator
  G4ElasticData & operator=(const G4ElasticData &right);
  G4ElasticData(const G4ElasticData&);

  G4double  R1, R2, Pnucl, Aeff;
  G4double  dQ2;
  G4double  massA;
  G4double  massA2;
  G4double  maxQ2[NENERGY];
  std::vector<G4double> fCumProb[NENERGY];
};

/////////////////////////////////////////////////////////////////////

class G4ElasticHadrNucleusHE : public G4HadronElastic
{
public:

  explicit G4ElasticHadrNucleusHE(const G4String& name = "hElasticGlauber");

  ~G4ElasticHadrNucleusHE() override;

  G4double SampleInvariantT(const G4ParticleDefinition* p, G4double plab, 
			    G4int Z, G4int A) override;

  void InitialiseModel() override;

  void ModelDescription(std::ostream&) const override;

private:

  G4double HadronNucleusQ2_2(const G4ElasticData *pElD, G4double plabGeV, 
                             G4double tmax);

  void DefineHadronValues(G4int Z);
  G4int FillFq2(G4int A); 

  G4double GetLightFq2(G4int Z, G4int A, G4double Q);

  G4double GetQ2_2(G4int N, G4int Nmax, 
		   const std::vector<G4double>& F, G4double rand);

  G4double HadrNucDifferCrSec(G4int A, G4double Q2); 

  void InterpolateHN(G4int n, const G4double EnP[], 
		     const G4double C0P[], const G4double C1P[], 
		     const G4double B0P[], const G4double B1P[]);

  G4double GetFt(G4double Q2);

  G4double HadronProtonQ2(G4double plab, G4double tmax);

  void Binom();

  void FillData(const G4ParticleDefinition* p, G4int idx, G4int Z);

  void InFileName(std::ostringstream&, const G4ParticleDefinition* p, G4int Z);

  void OutFileName(std::ostringstream&, const G4ParticleDefinition* p, G4int Z);

  G4bool ReadLine(std::ifstream&, std::vector<G4double>&);

  void WriteLine(std::ofstream&, std::vector<G4double>&);

  inline G4double LineInterpol(G4double p0, G4double p2,
                               G4double c1, G4double c2, G4double p);

  inline G4double GetBinomCof( G4int n, G4int m );

  // hide assignment operator
  G4ElasticHadrNucleusHE & operator=(const G4ElasticHadrNucleusHE &right);
  G4ElasticHadrNucleusHE(const G4ElasticHadrNucleusHE&);

  //  fields
  G4int    iHadrCode;
  G4int    iHadron;
  G4int    iHadron1;
  static const G4int fHadronCode[NHADRONS];
  static const G4int fHadronType[NHADRONS];
  static const G4int fHadronType1[NHADRONS];

  static G4bool fStoreToFile;
  static G4bool fRetrieveFromFile;

  // momemtum limits
  G4double ekinLowLimit;
  G4double dQ2;  

  // projectile kinematics in GeV
  G4double hMass;
  G4double hMass2;
  G4double hLabMomentum;
  G4double hLabMomentum2;
  G4double HadrEnergy;

  // elastic parameters
  G4double HadrTot, HadrSlope, HadrReIm, TotP; 
  G4double DDSect2, DDSect3, ConstU;

  // momentum limits for different models of hadron/nucleon scatetring
  G4double BoundaryP[7], BoundaryTL[7], BoundaryTG[7];

  // parameterisation of scattering
  G4double Slope1, Slope2, Coeff1, Coeff2;
  G4double Slope0, Coeff0;

  G4double aAIm, aDIm, Dtot11;

  // nucleaus parameters
  G4double  R1, R2, Pnucl, Aeff, Q2max;

  static G4double fLineF[ONQ2];
  static G4double fEnergy[NENERGY];
  static G4double fLowEdgeEnergy[NENERGY];
  static G4double fBinom[240][240];

  static G4ElasticData* fElasticData[NHADRONS][ZMAX];
  G4NistManager*  nistManager;
  const char* fDirectory;

  G4bool isMaster;

#ifdef G4MULTITHREADED
  static G4Mutex elasticMutex;
#endif

};

////////////////////////////////////////////////////////////////

inline
G4double G4ElasticHadrNucleusHE::LineInterpol(G4double p1, G4double p2,
					      G4double c1, G4double c2,
					      G4double p)
{
  return c1+(p-p1)*(c2-c1)/(p2-p1);
}

////////////////////////////////////////////////////////////////

inline
G4double G4ElasticHadrNucleusHE::GetBinomCof(G4int numN, G4int numM)
{
  return (numN >= numM && numN < 240) ? fBinom[numN][numM] : 0.0;
}

////////////////////////////////////////////////////////////////

#endif
