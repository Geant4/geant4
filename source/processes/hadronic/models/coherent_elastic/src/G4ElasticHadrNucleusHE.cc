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
//  The generator of high energy hadron-nucleus elastic scattering
//  The hadron kinetic energy T > 1 GeV
//  N.Starkov 2003.
//
//  19.11.05 The HE elastic scattering on proton is added (N.Starkov)
//  16.11.06 The low energy boundary is shifted to T = 400 MeV (N.Starkov)
//  23.11.06 General cleanup, ONQ0=3, use pointer instead of particle name (VI)
//  02.05.07 Scale sampled t as p^2 (VI)
//  15.05.07 Redesign and cleanup (V.Ivanchenko)
//  17.05.07 cleanup (V.Grichine)
//  19.04.12 Fixed reproducibility violation (A.Ribon)
//  12.06.12 Fixed warnings of shadowed variables (A.Ribon)
//

#include  "G4ElasticHadrNucleusHE.hh"
#include  "G4PhysicalConstants.hh"
#include  "G4SystemOfUnits.hh"
#include  "Randomize.hh"
#include  "G4ios.hh"
#include  "G4ParticleTable.hh"
#include  "G4NucleiProperties.hh"
#include  "G4IonTable.hh"
#include  "G4Proton.hh"
#include  "G4PionPlus.hh"
#include  "G4PionMinus.hh"
#include  "G4NistManager.hh"
#include  "G4ProductionCutsTable.hh"
#include  "G4MaterialCutsCouple.hh"
#include  "G4Material.hh"
#include  "G4Element.hh"
#include  "G4Log.hh"
#include  "G4Exp.hh"

using namespace std;

const G4int G4ElasticHadrNucleusHE::fHadronCode[] = 
{211,-211,2112,2212,321,-321,130,310,311,-311,
 3122,3222,3112,3212,3312,3322,3334,
 -2212,-2112,-3122,-3222,-3112,-3212,-3312,-3322,-3334};

const G4int G4ElasticHadrNucleusHE::fHadronType[] = 
{2,3,6,0,4,5,4,4,4,5,
 0,0,0,0,0,0,0,
 1,7,1,1,1,1,1,1,1};

const G4int G4ElasticHadrNucleusHE::fHadronType1[] = 
{3,4,1,0,5,6,5,5,5,6,
 0,0,0,0,0,0,0,
 2,2,2,2,2,2,2,2,2};

G4double G4ElasticHadrNucleusHE::fLineF[]  = {0.0};
G4double G4ElasticHadrNucleusHE::fEnergy[] = {0.0};
G4double G4ElasticHadrNucleusHE::fLowEdgeEnergy[] = {0.0};
G4double G4ElasticHadrNucleusHE::fBinom[240][240] = {{0.0}};

G4ElasticData* 
G4ElasticHadrNucleusHE::fElasticData[NHADRONS][ZMAX] = {{nullptr}};

#ifdef G4MULTITHREADED
  G4Mutex G4ElasticHadrNucleusHE::elasticMutex = G4MUTEX_INITIALIZER;
#endif

G4bool G4ElasticHadrNucleusHE::fStoreToFile = false;
G4bool G4ElasticHadrNucleusHE::fRetrieveFromFile = false;

const G4double invGeV    =  1.0/CLHEP::GeV;
const G4double MbToGeV2  =  2.568;
const G4double GeV2      =  CLHEP::GeV*CLHEP::GeV;
const G4double invGeV2   =  1.0/GeV2;
const G4double protonM   =  CLHEP::proton_mass_c2*invGeV;
const G4double protonM2  =  protonM*protonM;

///////////////////////////////////////////////////////////////

G4ElasticData::G4ElasticData(const G4ParticleDefinition* p, 
			     G4int Z, G4int A, const G4double* e) 
{ 
  G4double massGeV  = p->GetPDGMass()*invGeV;
  G4double mass2GeV2= massGeV*massGeV;

  DefineNucleusParameters(A);
  G4double limitQ2 = 35./(R1*R1);     //  (GeV/c)^2

  massA  = G4NucleiProperties::GetNuclearMass(A, Z)*invGeV;
  massA2 = massA*massA; 
  /*
  G4cout << " G4ElasticData for " << p->GetParticleName()
	 << " Z= " << Z << " A= " << A << " R1= " << R1 
	 << " R2= " << R2 << G4endl;  
  */
  for(G4int kk = 0; kk<NENERGY; ++kk) 
  {
    G4double elab = e[kk] + massGeV;
    G4double plab2= e[kk]*(e[kk] + 2.0*massGeV);
    G4double Q2m  = 4.0*plab2*massA2/(mass2GeV2 + massA2 + 2.*massA*elab);

    if(Z == 1 && p == G4Proton::Proton()) { Q2m *= 0.5; }

    maxQ2[kk] = Q2m;
    /*
    G4cout << " Ekin= " << e[kk] << " Q2m= " << Q2m 
	   << " limitQ2= " << limitQ2 << G4endl;
    */
  }

  dQ2 = limitQ2/(G4double)(ONQ2-2);
}

/////////////////////////////////////////////////////////////////////////

void G4ElasticData::DefineNucleusParameters(G4int A)
{
  switch (A) {
    case 207:
    case 208:
      R1       = 20.5;
      R2       = 15.74;
      Pnucl    = 0.4;
      Aeff     = 0.7;
      break;
    case 237:
    case 238:
      R1       = 21.7;    
      R2       = 16.5;
      Pnucl    = 0.4;
      Aeff     = 0.7;
      break;
    case 90:
    case 91:
      R1    = 16.5;
      R2    = 11.62;
      Pnucl = 0.4;
      Aeff  = 0.7;
      break;
    case 58:
    case 59:
      R1    = 15.75;
      R2    = 9.9;
      Pnucl = 0.45;
      Aeff  = 0.85;
      break;
    case 48:
    case 47:
      R1    = 14.0;
      R2    = 9.26;
      Pnucl = 0.31;
      Aeff  = 0.75;
      break;
    case 40:
    case 41:
      R1    = 13.3;
      R2    = 9.26;
      Pnucl = 0.31;
      Aeff  = 0.75;
      break;
    case 28:
    case 29:
      R1    = 12.0;
      R2    = 7.64;
      Pnucl = 0.253;
      Aeff  = 0.8;
      break;
    case 16:
      R1    = 10.50;
      R2    = 5.5;
      Pnucl = 0.7;
      Aeff  = 0.98;
      break;
    case 12:
      R1    = 9.3936;
      R2    = 4.63;
      Pnucl = 0.7;
      Aeff  = 1.0;
      break;
    case 11:
      R1    = 9.0;
      R2    = 5.42;
      Pnucl = 0.19;
      Aeff  = 0.9;
      break;
    case 9:
      R1    = 9.9;
      R2    = 6.5;
      Pnucl = 0.690;
      Aeff  = 0.95;
      break;
    case 4:
      R1    = 5.3;   
      R2    = 3.7;
      Pnucl = 0.4;
      Aeff  = 0.75;
      break;
    case 1:
      R1    = 4.5;   
      R2    = 2.3;
      Pnucl = 0.177;
      Aeff  = 0.9;
      break;
    default:
      R1    = 4.45*G4Exp(G4Log((G4double)(A - 1))*0.309)*0.9;
      R2    = 2.3 *G4Exp(G4Log((G4double)A)* 0.36);

      if(A < 100 && A > 3) { Pnucl = 0.176 + 0.00275*A; }
      else                 { Pnucl = 0.4; }
      //G4cout<<" Deault: A= "<<A<<"  R1 R2 Aeff Pnucl "<<R1<<"  "<<R2<<"  "
      //      <<Aeff<<"  "<<Pnucl<<G4endl;

      if(A >= 100)               { Aeff = 0.7; }
      else if(A < 100 && A > 75) { Aeff = 1.5 - 0.008*A; }
      else                       { Aeff = 0.9; }
      break;
  }
  //G4cout<<" Result: A= "<<A<<"  R1 R2 Aeff Pnucl "<<R1<<"  "<<R2<<"  "
  //      <<Aeff<<"  "<<Pnucl<<G4endl;
}

////////////////////////////////////////////////////////////////////

G4ElasticHadrNucleusHE::G4ElasticHadrNucleusHE(const G4String& name)
  : G4HadronElastic(name), fDirectory(nullptr), isMaster(false)
{
  dQ2 = hMass = hMass2 = hLabMomentum = hLabMomentum2 = HadrEnergy 
    = R1 = R2 = Pnucl = Aeff = HadrTot = HadrSlope = HadrReIm = TotP = DDSect2
    = DDSect3 = ConstU = Slope1 = Slope2 = Coeff1 = Coeff2
    = Slope0 = Coeff0 = aAIm = aDIm = Dtot11 = Q2max = 0.0;
  iHadrCode = iHadron = iHadron1 = 0;

  verboseLevel = 0;
  ekinLowLimit = 400.0*CLHEP::MeV;

  BoundaryP[0]=9.0; BoundaryTG[0]=5.0;BoundaryTL[0]=0.;
  BoundaryP[1]=20.0;BoundaryTG[1]=1.5;BoundaryTL[1]=0.;
  BoundaryP[2]=5.0; BoundaryTG[2]=1.0;BoundaryTL[2]=1.5;
  BoundaryP[3]=8.0; BoundaryTG[3]=3.0;BoundaryTL[3]=0.;
  BoundaryP[4]=7.0; BoundaryTG[4]=3.0;BoundaryTL[4]=0.;
  BoundaryP[5]=5.0; BoundaryTG[5]=2.0;BoundaryTL[5]=0.;
  BoundaryP[6]=5.0; BoundaryTG[6]=1.5;BoundaryTL[6]=3.0;

  nistManager = G4NistManager::Instance();

  if(fEnergy[0] == 0.0) {
#ifdef G4MULTITHREADED
    G4MUTEXLOCK(&elasticMutex);
    if(fEnergy[0] == 0.0) {
#endif
      isMaster = true;
      Binom();
      // energy in GeV
      fEnergy[0] = 0.4;
      fEnergy[1] = 0.6;
      fEnergy[2] = 0.8;
      fEnergy[3] = 1.0;
      fLowEdgeEnergy[0] = 0.0;
      fLowEdgeEnergy[1] = 0.5;
      fLowEdgeEnergy[2] = 0.7;
      fLowEdgeEnergy[3] = 0.9;
      G4double f = G4Exp(G4Log(10.)*0.1);
      G4double e = f*f;
      for(G4int i=4; i<NENERGY; ++i) {
	fEnergy[i] = e;
	fLowEdgeEnergy[i] = e/f;
	e *= f*f;
      }
      if(verboseLevel > 0) {
	G4cout << "### G4ElasticHadrNucleusHE: energy points in GeV" << G4endl;
	for(G4int i=0; i<NENERGY; ++i) {
	  G4cout << "  " << i << "   " << fLowEdgeEnergy[i] 
		 << "  " << fEnergy[i] << G4endl;
	}
      }
#ifdef G4MULTITHREADED
    }
    G4MUTEXUNLOCK(&elasticMutex);
#endif
  }
}

///////////////////////////////////////////////////////////////////

void G4ElasticHadrNucleusHE::ModelDescription(std::ostream& outFile) const
{
  outFile << "G4ElasticHadrNucleusHE is a hadron-nucleus elastic scattering\n"
	  << "model developed by N. Starkov which uses a Glauber model\n"
	  << "parameterization to calculate the final state.  It is valid\n"
	  << "for all hadrons with incident momentum above 0.4 GeV/c.\n";
}

///////////////////////////////////////////////////////////////////

G4ElasticHadrNucleusHE::~G4ElasticHadrNucleusHE()
{
  if(isMaster) {
    for(G4int j = 0; j < NHADRONS; ++j) {
      for(G4int k = 0; k < ZMAX; ++k) {
	G4ElasticData* ptr = fElasticData[j][k]; 
	if(ptr) { 
	  delete ptr;
	  fElasticData[j][k] = nullptr;
	  for(G4int l = j+1; l < NHADRONS; ++l) {
	    if(ptr == fElasticData[l][k]) { fElasticData[l][k] = nullptr; }
	  }
	}
      }
    }
    delete fDirectory;
    fDirectory = nullptr;
  }
}

///////////////////////////////////////////////////////////////////

void G4ElasticHadrNucleusHE::InitialiseModel()
{
  if(!isMaster) { return; }
  G4ProductionCutsTable* theCoupleTable=
    G4ProductionCutsTable::GetProductionCutsTable();
  G4int numOfCouples = (G4int)theCoupleTable->GetTableSize();
  
  for(G4int i=0; i<2; ++i) {
    const G4ParticleDefinition* p = G4PionPlus::PionPlus();
    if(1 == i) { p = G4PionMinus::PionMinus(); } 
    iHadrCode = fHadronCode[i]; 
    iHadron   = fHadronType[i];
    iHadron1  = fHadronType1[i];
    hMass     = p->GetPDGMass()*invGeV;
    hMass2    = hMass*hMass;
    for(G4int j=0; j<numOfCouples; ++j) {
      auto mat = theCoupleTable->GetMaterialCutsCouple(j)->GetMaterial();
      auto elmVec = mat->GetElementVector();
      std::size_t numOfElem = mat->GetNumberOfElements();
      for(std::size_t k=0; k<numOfElem; ++k) {
        G4int Z = std::min((*elmVec)[k]->GetZasInt(), ZMAX-1);
        if(!fElasticData[i][Z]) { 
          if(1 == i && Z > 1) { 
            fElasticData[1][Z] = fElasticData[0][Z]; 
          } else {
            FillData(p, i, Z);
          } 
        }
      }
    }
  }
}

////////////////////////////////////////////////////////////////////

G4double 
G4ElasticHadrNucleusHE::SampleInvariantT(const G4ParticleDefinition* p,
					 G4double inLabMom, 
					 G4int iZ, G4int A)
{
  G4double mass = p->GetPDGMass();
  G4double kine = sqrt(inLabMom*inLabMom + mass*mass) - mass;
  if(kine <= ekinLowLimit) {
    return G4HadronElastic::SampleInvariantT(p,inLabMom,iZ,A);
  }
  G4int Z = std::min(iZ,ZMAX-1);
  G4double Q2 = 0.0;
  iHadrCode = p->GetPDGEncoding();

  // below computations in GeV/c
  hMass  = mass*invGeV;
  hMass2 = hMass*hMass;
  G4double plab = inLabMom*invGeV;
  G4double tmax = pLocalTmax*invGeV2;

  if(verboseLevel > 1) {
    G4cout<< "G4ElasticHadrNucleusHE::SampleT: " 
	  << " for " << p->GetParticleName() 
	  << " at Z= " << Z << " A= " << A
	  << " plab(GeV)= " << plab
	  << " hadrCode= " << iHadrCode
	  << G4endl;
  }
  iHadron = -1;
  G4int idx;
  for(idx=0; idx<NHADRONS; ++idx) {
    if(iHadrCode == fHadronCode[idx]) { 
      iHadron = fHadronType[idx];
      iHadron1 = fHadronType1[idx];
      break; 
    }
  }
  // Hadron is not in the list
  if(0 > iHadron) { return 0.0; }

  if(Z==1) {
    Q2 = HadronProtonQ2(plab, tmax);

    if (verboseLevel>1) {
      G4cout<<"  Proton : Q2  "<<Q2<<G4endl;
    }
  } else {
    const G4ElasticData* ElD1 = fElasticData[idx][Z];

    // Construct elastic data
    if(!ElD1) { 
      FillData(p, idx, Z); 
      ElD1 = fElasticData[idx][Z];
      if(!ElD1) { return 0.0; }
    }

    // sample scattering
    Q2 = HadronNucleusQ2_2(ElD1, plab, tmax);

    if(verboseLevel > 1) {
      G4cout<<" SampleT: Q2(GeV^2)= "<<Q2<< "  t/tmax= " 
            << Q2/tmax <<G4endl;
    }
  }
  return Q2*GeV2;
}

////////////////////////////////////////////////////////////////

void G4ElasticHadrNucleusHE::FillData(const G4ParticleDefinition* p, 
                                      G4int idx, G4int Z)
{
#ifdef G4MULTITHREADED
  G4MUTEXLOCK(&elasticMutex);
  if(!fElasticData[idx][Z]) {
#endif
    G4int A = G4lrint(nistManager->GetAtomicMassAmu(Z));
    G4ElasticData* pElD = new G4ElasticData(p, Z, A, fEnergy);
    if(fRetrieveFromFile) { 
      std::ostringstream ss;
      InFileName(ss, p, Z); 
      std::ifstream infile(ss.str(), std::ios::in);
      for(G4int i=0; i<NENERGY; ++i) {
	if(ReadLine(infile, pElD->fCumProb[i])) {
	  continue;
	} else {
	  fRetrieveFromFile = false;
          break;
	}
      }
      infile.close();
    }
    R1     = pElD->R1;
    R2     = pElD->R2;
    Aeff   = pElD->Aeff;
    Pnucl  = pElD->Pnucl;
    dQ2    = pElD->dQ2;
    if(verboseLevel > 0) {
      G4cout<<"### FillData for " << p->GetParticleName() 
	    << " Z= " << Z << " idx= " << idx << " iHadron= " << iHadron 
	    <<" iHadron1= " << iHadron1 << " iHadrCode= " << iHadrCode
            <<"\n   R1= " << R1 << " R2= " << R2 << " Aeff= " << Aeff 
	    <<" Pnucl= " << Pnucl << G4endl;
    }

    if(!fRetrieveFromFile) {  
      for(G4int i=0; i<NENERGY; ++i) {
	G4double T = fEnergy[i];
	hLabMomentum2 = T*(T + 2.*hMass);
	hLabMomentum  = std::sqrt(hLabMomentum2);
	HadrEnergy = hMass + T;
	DefineHadronValues(Z);
	Q2max = pElD->maxQ2[i];

	G4int length  = FillFq2(A); 
	(pElD->fCumProb[i]).reserve(length);
	G4double norm = 1.0/fLineF[length-1];

	if(verboseLevel > 0) {
	  G4cout << "### i= " << i << " Z= " << Z << " A= " << A 
		 << " length= " << length << " Q2max= " << Q2max << G4endl;
	}

	(pElD->fCumProb[i]).push_back(0.0);
	for(G4int ii=1; ii<length-1; ++ii) {
	  (pElD->fCumProb[i]).push_back(fLineF[ii]*norm);
	  if(verboseLevel > 2) {
	    G4cout << "    ii= " << ii << " val= " 
		   << (pElD->fCumProb[i])[ii] << G4endl;
	  }
	}
	(pElD->fCumProb[i]).push_back(1.0);
      }
    }

    if(fStoreToFile) {
      std::ostringstream ss;
      OutFileName(ss, p, Z); 
      std::ofstream fileout(ss.str());
      for(G4int i=0; i<NENERGY; ++i) {
	WriteLine(fileout, pElD->fCumProb[i]);
      }
      fileout.close();
    }
    
    if(verboseLevel > 0) {
      G4cout << " G4ElasticHadrNucleusHE::FillData done for idx= " << idx
	     << " for " << p->GetParticleName() << " Z= " << Z
	     << " A= " << A << G4endl;
    }
    fElasticData[idx][Z] = pElD;

#ifdef G4MULTITHREADED
  }
  G4MUTEXUNLOCK(&elasticMutex);
#endif  
}

////////////////////////////////////////////////////////////////

void G4ElasticHadrNucleusHE::InterpolateHN(G4int n, const G4double EnP[], 
                   const G4double C0P[], const G4double C1P[], 
                   const G4double B0P[], const G4double B1P[])
{
  G4int i; 

  for(i=1; i<n; ++i) { if(hLabMomentum <= EnP[i]) { break; } }
  if(i == n) { i = n - 1; }

  Coeff0 = LineInterpol(EnP[i], EnP[i-1], C0P[i], C0P[i-1], hLabMomentum);
  Coeff1 = LineInterpol(EnP[i], EnP[i-1], C1P[i], C1P[i-1], hLabMomentum);
  Slope0 = LineInterpol(EnP[i], EnP[i-1], B0P[i], B0P[i-1], hLabMomentum);
  Slope1 = LineInterpol(EnP[i], EnP[i-1], B1P[i], B1P[i-1], hLabMomentum);

//  G4cout<<"  InterpolHN:  n i "<<n<<"  "<<i<<"  Mom "
//        <<hLabMomentum<<G4endl;
}

//////////////////////////////////////////////////////////////////////////

G4double 
G4ElasticHadrNucleusHE::HadronNucleusQ2_2(const G4ElasticData* pElD,
                                          G4double plab, G4double tmax)
{
  G4double ekin  = std::sqrt(hMass2 + plab*plab) - hMass;

  if(verboseLevel > 1) {
    G4cout<<"Q2_2: ekin(GeV)= " << ekin << "  plab(GeV/c)= " << plab 
          <<"  tmax(GeV2)= " << tmax <<G4endl;
  }
  // Find closest energy bin
  G4int idx; 
  for(idx=0; idx<NENERGY-1; ++idx) {
    if(ekin <= fLowEdgeEnergy[idx+1]) { break; }
  }
  //G4cout << "   idx= " << idx << G4endl;

  // Select kinematics for node energy
  R1    = pElD->R1;
  dQ2   = pElD->dQ2;
  Q2max = pElD->maxQ2[idx];
  G4int length = (G4int)(pElD->fCumProb[idx]).size();

  G4double Rand = G4UniformRand();

  G4int iNumbQ2 = 0;
  for(iNumbQ2=1; iNumbQ2<length; ++iNumbQ2) {
    if(Rand <= (pElD->fCumProb[idx])[iNumbQ2]) { break; }
  }
  iNumbQ2 = std::min(iNumbQ2, length - 1);
  G4double Q2 = GetQ2_2(iNumbQ2, length, pElD->fCumProb[idx], Rand);
  Q2 = std::min(Q2, Q2max);
  Q2 *= tmax/Q2max;

  if(verboseLevel > 1) {
    G4cout<<" HadrNucleusQ2_2(2): Q2= "<<Q2<<" iNumbQ2= " << iNumbQ2 
	  << " rand= " << Rand << " Q2max= " << Q2max 
          << " tmax= " << tmax << G4endl;
  }
  return Q2;
}       

///////////////////////////////////////////////////////////////////////
//
//  The randomization of one dimensional array 
//

G4double G4ElasticHadrNucleusHE::GetQ2_2(G4int kk, G4int kmax,
					 const std::vector<G4double>& F, 
                                         G4double ranUni)
{
  //G4cout << "GetQ2_2 kk= " << kk << " kmax= " << kmax << "  size= " 
  //	 << F.size() << "  rand= " << ranUni << G4endl;
  if(kk == kmax-1) {
    G4double X1 = dQ2*kk;
    G4double F1 = F[kk-1];
    G4double X2 = Q2max;
    G4double xx = R1*(X2 - X1);
    xx = (xx > 20.) ? 0.0 : G4Exp(-xx);
    G4double Y = X1 - G4Log(1.0 - (ranUni - F1)*(1.0 - xx)/(1.0 - F1))/R1;
    return Y;
  } 
  G4double F1, F2, F3, X1, X2, X3;

  if(kk == 1 || kk == 0) {
    F1 = F[0]; 
    F2 = F[1];
    F3 = F[2];
    X1 = 0.0;
    X2 = dQ2;
    X3 = dQ2*2;
  } else {
    F1 = F[kk-2];
    F2 = F[kk-1];
    F3 = F[kk];
    X1 = dQ2*(kk-2);
    X2 = dQ2*(kk-1);
    X3 = dQ2*kk;
  }
  if(verboseLevel > 1) {
    G4cout << "GetQ2_2 kk= " << kk << " X2= " << X2 << " X3= " << X3 
	   << " F2= " << F2 << " F3= " << F3 << " Rndm= " << ranUni << G4endl;
  }

  G4double F12 = F1*F1;
  G4double F22 = F2*F2;
  G4double F32 = F3*F3;

  G4double D0  = F12*F2+F1*F32+F3*F22-F32*F2-F22*F1-F12*F3;

  if(verboseLevel > 2) {
    G4cout << "       X1= " << X1 << " F1= " << F1 << "  D0= " 
           << D0 << G4endl; 
  }
  G4double Y;
  if(std::abs(D0) < 1.e-9) { 
    Y = X2 + (ranUni - F2)*(X3 - X2)/(F3 - F2);
  } else {
    G4double DA = X1*F2+X3*F1+X2*F3-X3*F2-X1*F3-X2*F1;
    G4double DB = X2*F12+X1*F32+X3*F22-X2*F32-X3*F12-X1*F22;
    G4double DC = X3*F2*F12+X2*F1*F32+X1*F3*F22
	           -X1*F2*F32-X2*F3*F12-X3*F1*F22;
    Y = (DA*ranUni*ranUni + DB*ranUni + DC)/D0;
  }
  return Y;
}

////////////////////////////////////////////////////////////////////////

G4int G4ElasticHadrNucleusHE::FillFq2(G4int A) 
{
  G4double curQ2, curSec;
  G4double curSum = 0.0;
  G4double totSum = 0.0;

  G4double ddQ2 = dQ2*0.1;
  G4double Q2l  = 0.0;

  G4int ii = 0;
  for(ii=1; ii<ONQ2-1; ++ii) {
    curSum = curSec = 0.0;

    for(G4int jj=0; jj<10; ++jj) {
      curQ2 = Q2l+(jj + 0.5)*ddQ2;
      if(curQ2 >= Q2max) { break; }
      curSec = HadrNucDifferCrSec(A, curQ2);
      curSum += curSec;
    }
    G4double del = (curQ2 >= Q2max) ? Q2max - Q2l : dQ2;
    Q2l    += del;
    curSum *= del*0.1;
    totSum += curSum;
    fLineF[ii] = totSum;
    if (verboseLevel>2) {
      G4cout<<ii << ". FillFq2: A= " << A << " Q2= "<<Q2l<<" dQ2= "
	    <<dQ2<<" Tot= "<<totSum << " dTot " <<curSum
	    <<" curSec= " <<curSec<<G4endl;
    }
    if(totSum*1.e-4 > curSum || Q2l >= Q2max) { break; }
  }
  ii = std::min(ii, ONQ2-2);
  curQ2 = Q2l;
  G4double xx = R1*(Q2max - curQ2);
  if(xx > 0.0) {
    xx = (xx > 20.) ? 0.0 : G4Exp(-xx);
    curSec = HadrNucDifferCrSec(A, curQ2);
    totSum += curSec*(1.0 - xx)/R1;
  }
  fLineF[ii + 1] = totSum;
  if (verboseLevel>1) {
    G4cout << "### FillFq2 done curQ2= " << curQ2 << " Q2max= "<< Q2max 
           << " sumG= " << fLineF[ONQ2-2] << "  totSum= " << totSum 
	   << " Nbins= " << ii + 1 << G4endl;
  }
  return ii + 2;
}

////////////////////////////////////////////////////////////////////////

G4double G4ElasticHadrNucleusHE::GetLightFq2(G4int Z, G4int A, G4double Q2)
{
  // Scattering off proton
  if(Z == 1) 
  {
    G4double SqrQ2  = std::sqrt(Q2);
    G4double valueConstU = 2.*(hMass2 + protonM2) - Q2;

    G4double y = (1.-Coeff1-Coeff0)/HadrSlope*(1.-G4Exp(-HadrSlope*Q2))
      + Coeff0*(1.-G4Exp(-Slope0*Q2))
      + Coeff2/Slope2*G4Exp(Slope2*valueConstU)*(G4Exp(Slope2*Q2)-1.)
      + 2.*Coeff1/Slope1*(1./Slope1-(1./Slope1+SqrQ2)*G4Exp(-Slope1*SqrQ2));

    return y;
  }

  // The preparing of probability function  

  G4double prec = A > 208  ?  1.0e-7 : 1.0e-6;

  G4double    Stot     = HadrTot*MbToGeV2;     //  Gev^-2
  G4double    Bhad     = HadrSlope;         //  GeV^-2
  G4double    Asq      = 1+HadrReIm*HadrReIm;
  G4double    Rho2     = std::sqrt(Asq);

  if(verboseLevel >1) {
    G4cout<<" Fq2 Before for i Tot B Im "<<HadrTot<<"  "<<HadrSlope<<"  "
      <<HadrReIm<<G4endl;
  }
  if(verboseLevel > 1) {
    G4cout << "GetFq2: Stot= " << Stot << " Bhad= " << Bhad 
           <<"  Im "<<HadrReIm 
           << " Asq= " << Asq << G4endl;
    G4cout << "R1= " << R1 << " R2= " << R2 << " Pnucl= " << Pnucl <<G4endl;
  }
  G4double    R12      = R1*R1;
  G4double    R22      = R2*R2;
  G4double    R12B     = R12+2*Bhad;
  G4double    R22B     = R22+2*Bhad;

  G4double    Norm     = (R12*R1-Pnucl*R22*R2); // HP->Aeff;

  G4double    R13      = R12*R1/R12B;
  G4double    R23      = Pnucl*R22*R2/R22B;
  G4double    Unucl    = Stot/twopi*R13/Norm;
  G4double    UnucRho2 = -Unucl*Rho2;

  G4double    FiH      = std::asin(HadrReIm/Rho2);
  G4double    NN2      = R23/R13;

  if(verboseLevel > 2) {
    G4cout << "UnucRho2= " << UnucRho2 << " FiH= " << FiH << " NN2= " << NN2 
	   << " Norm= " << Norm << G4endl;
  }
  G4double    Prod0 = 0.;
  G4double    N1    = -1.0;

  for(G4int i1 = 1; i1<= A; ++i1) ////++++++++++  i1
    {
      N1 *= (-Unucl*Rho2*(A-i1+1)/(G4double)i1);
      G4double Prod1 = 0.;
      G4double N2    = -1.;

      for(G4int i2 = 1; i2<=A; ++i2) ////+++++++++ i2
        {
          N2 *= (-Unucl*Rho2*(A-i2+1)/(G4double)i2);
          G4double Prod2 = 0; 
          G4double N5    = -1/NN2;
	  for(G4int j2=0; j2<= i2; ++j2) ////+++++++++ j2
            {
              G4double Prod3 = 0;
              G4double exp2  = 1./((G4double)j2/R22B+(G4double)(i2-j2)/R12B);
              N5 *= (-NN2);
              G4double N4 = -1./NN2;
	      for(G4int j1=0; j1<=i1; ++j1) ////++++++++ j1
		{
		  G4double exp1  = 1./((G4double)j1/R22B+(G4double)(i1-j1)/R12B);
		  G4double dddd  = 0.25*(exp1+exp2);
		  N4    *= (-NN2);
		  Prod3 += 
                    N4*exp1*exp2*(1.-G4Exp(-Q2*dddd))*GetBinomCof(i1,j1)/dddd;
		}                                   // j1
	      Prod2 += Prod3*N5*GetBinomCof(i2,j2);
	    }                                      // j2
	  Prod1 += Prod2*N2*std::cos(FiH*(i1-i2));

	  if (std::abs(Prod2*N2/Prod1)<prec) break;
        }                                         // i2
      Prod0 += Prod1*N1;
      if(std::abs(N1*Prod1/Prod0) < prec) break;
    }                                           // i1

  const G4double fact = 0.25*CLHEP::pi/MbToGeV2; 
  Prod0 *= fact;  //  This is in mb

  if(verboseLevel>1) {
    G4cout << "GetLightFq2 Z= " << Z << " A= " << A 
	   <<" Q2= " << Q2 << " Res= " << Prod0 << G4endl;
  }
  return Prod0;
}

///////////////////////////////////////////////////////////////////

G4double 
G4ElasticHadrNucleusHE::HadrNucDifferCrSec(G4int A, G4double aQ2)
{
  //   ------ All external kinematical variables are in MeV -------
  //            ------ but internal in GeV !!!!  ------

  // Scattering of proton
  if(A == 1) 
  {
    G4double SqrQ2  = std::sqrt(aQ2);
    G4double valueConstU = hMass2 + protonM2-2*protonM*HadrEnergy - aQ2;
    
    BoundaryTL[0] = Q2max;
    BoundaryTL[1] = Q2max;
    BoundaryTL[3] = Q2max;
    BoundaryTL[4] = Q2max;
    BoundaryTL[5] = Q2max;

    G4double dSigPodT = HadrTot*HadrTot*(1+HadrReIm*HadrReIm)*
      ( Coeff1*G4Exp(-Slope1*SqrQ2)+
	Coeff2*G4Exp( Slope2*(valueConstU)+aQ2)+
	(1-Coeff1-Coeff0)*G4Exp(-HadrSlope*aQ2)+
	Coeff0*G4Exp(-Slope0*aQ2) )*2.568/(16*pi);

    return dSigPodT;
  }

  G4double    Stot     = HadrTot*MbToGeV2; 
  G4double    Bhad     = HadrSlope; 
  G4double    Asq      = 1+HadrReIm*HadrReIm;
  G4double    Rho2     = std::sqrt(Asq);
  G4double    R12      = R1*R1;
  G4double    R22      = R2*R2;
  G4double    R12B     = R12+2*Bhad;
  G4double    R22B     = R22+2*Bhad;
  G4double    R12Ap    = R12+20;
  G4double    R22Ap    = R22+20;
  G4double    R13Ap    = R12*R1/R12Ap;
  G4double    R23Ap    = R22*R2*Pnucl/R22Ap;
  G4double    R23dR13  = R23Ap/R13Ap;
  G4double    R12Apd   = 2/R12Ap;
  G4double    R22Apd   = 2/R22Ap;
  G4double R12ApdR22Ap = 0.5*(R12Apd+R22Apd);

  G4double DDSec1p  = (DDSect2+DDSect3*G4Log(0.53*HadrEnergy/R1));
  G4double DDSec2p  = (DDSect2+DDSect3*G4Log(0.53*HadrEnergy/
                             std::sqrt((R12+R22)*0.5)));
  G4double DDSec3p  = (DDSect2+DDSect3*G4Log(0.53*HadrEnergy/R2));

  G4double    Norm     = (R12*R1-Pnucl*R22*R2)*Aeff;
  G4double    R13      = R12*R1/R12B;
  G4double    R23      = Pnucl*R22*R2/R22B;
  G4double    Unucl    = Stot/(twopi*Norm)*R13;
  G4double    UnuclScr = Stot/(twopi*Norm)*R13Ap;
  G4double    SinFi    = HadrReIm/Rho2;
  G4double    FiH      = std::asin(SinFi);
  G4double    N        = -1;
  G4double    N2       = R23/R13;

  G4double ImElasticAmpl0 = 0;
  G4double ReElasticAmpl0 = 0;
  G4double exp1;

  for(G4int i=1; i<=A; ++i) {
    N  *= (-Unucl*Rho2*(A-i+1)/(G4double)i);
    G4double N4 = 1;
    G4double medTot = R12B/(G4double)i;
    G4double Prod1  = G4Exp(-aQ2*R12B/(G4double)(4*i))*medTot;

    for(G4int l=1; l<=i; ++l) {
      exp1 = l/R22B+(i-l)/R12B;
      N4 *= (-N2*(i-l+1)/(G4double)l);
      G4double expn4 = N4/exp1;
      Prod1  += expn4*G4Exp(-aQ2/(exp1*4));
      medTot += expn4;
    }  // end l

    G4double dcos = N*std::cos(FiH*i);
    ReElasticAmpl0  += Prod1*N*std::sin(FiH*i);
    ImElasticAmpl0  += Prod1*dcos;
    if(std::abs(Prod1*N/ImElasticAmpl0) < 0.000001) break;
  }      // i

  static const G4double pi25 = CLHEP::pi/2.568;
  ImElasticAmpl0 *= pi25;   // The amplitude in mB
  ReElasticAmpl0 *= pi25;   // The amplitude in mB

  G4double C1 = R13Ap*R13Ap*0.5*DDSec1p;
  G4double C2 = 2*R23Ap*R13Ap*0.5*DDSec2p;
  G4double C3 = R23Ap*R23Ap*0.5*DDSec3p;
  
  G4double N1p  = 1;
  G4double Din1 = 0.5*(C1*G4Exp(-aQ2/8*R12Ap)/2*R12Ap-
		       C2/R12ApdR22Ap*G4Exp(-aQ2/(4*R12ApdR22Ap))+
		       C3*R22Ap/2*G4Exp(-aQ2/8*R22Ap));

  G4double DTot1 = 0.5*(C1*0.5*R12Ap-C2/R12ApdR22Ap+C3*R22Ap*0.5);

  for(G4int i=1; i<= A-2; ++i) {
    N1p *= (-UnuclScr*Rho2*(A-i-1)/(G4double)i);
    G4double N2p  = 1;
    G4double Din2 = 0;
    G4double DmedTot = 0;
    G4double BinCoeff = 1.0;
    for(G4int l=0; l<=i; ++l) {
      if(l > 0) { BinCoeff *= (i-l+1)/(G4double)l; }

      exp1  = l/R22B+(i-l)/R12B;
      G4double exp1p = exp1+R12Apd;
      G4double exp2p = exp1+R12ApdR22Ap;
      G4double exp3p = exp1+R22Apd;

      Din2 += N2p*BinCoeff*(C1/exp1p*G4Exp(-aQ2/(4*exp1p))-
			    C2/exp2p*G4Exp(-aQ2/(4*exp2p))+
			    C3/exp3p*G4Exp(-aQ2/(4*exp3p)));

      DmedTot += N2p*BinCoeff*(C1/exp1p-C2/exp2p+C3/exp3p);

      N2p *= -R23dR13;
    }     // l

    G4double dcos = N1p*std::cos(FiH*i)/(G4double)((i+2)*(i+1));
    Din1  += Din2*dcos;
    DTot1 += DmedTot*dcos;
    
    if(std::abs(Din2*N1p/Din1) < 0.000001) break;
  }           //  i
  G4double gg = (G4double)(A*(A-1)*4)/(Norm*Norm); 

  Din1  *= (-gg);
  DTot1 *= 5*gg;

  //  ----------------  dSigma/d|-t|,  mb/(GeV/c)^-2  -----------------

  G4double DiffCrSec2 = (ReElasticAmpl0*ReElasticAmpl0+
			 (ImElasticAmpl0+Din1)*
			 (ImElasticAmpl0+Din1))/twopi;

  Dtot11 = DTot1;
  aAIm   = ImElasticAmpl0;
  aDIm   = Din1;

  return DiffCrSec2;  //  dSig/d|-t|,  mb/(GeV/c)^-2
} 

////////////////////////////////////////////////////////////////

void G4ElasticHadrNucleusHE::DefineHadronValues(G4int Z)
{
  G4double sHadr = 2.*HadrEnergy*protonM+protonM2+hMass2;
  G4double sqrS  = std::sqrt(sHadr);
  
  if(verboseLevel>2) {
    G4cout << "GetHadrValues: Z= " << Z << " iHadr= " << iHadron 
	   << " E(GeV)= " << HadrEnergy << " sqrS= " << sqrS
	   << " plab= " << hLabMomentum   
	   <<"  E - m  "<<HadrEnergy - hMass<< G4endl;
  }
  G4double TotN = 0.0;
  G4double logE = G4Log(HadrEnergy);
  G4double logS = G4Log(sHadr);
           TotP = 0.0;

  switch (iHadron) {
  case 0:                  //  proton, neutron
  case 6:

    if(hLabMomentum > 10) {
      TotP = TotN = 7.5*logE - 40.12525 + 103*G4Exp(-logS*0.165);// mb

    } else {
      // ==================  neutron  ================

      if( hLabMomentum > 1.4 ) {
	TotN = 33.3+15.2*(hLabMomentum2-1.35)/
	  (G4Exp(G4Log(hLabMomentum)*2.37)+0.95);
		
      } else if(hLabMomentum > 0.8) {
	G4double A0 = logE + 0.0513;
	TotN = 33.0 + 25.5*A0*A0;  
      } else {
	G4double A0 = logE - 0.2634;  // log(1.3)
	TotN = 33.0 + 30.*A0*A0*A0*A0;
      }
      //  =================  proton  ===============

      if(hLabMomentum >= 1.05) {
	TotP = 39.0+75.*(hLabMomentum-1.2)/(hLabMomentum2*hLabMomentum+0.15);
      } else if(hLabMomentum >= 0.7) {
	G4double A0 = logE + 0.3147;
	TotP = 23.0 + 40.*A0*A0;
      } else {
	TotP = 23.+50.*G4Exp(G4Log(G4Log(0.73/hLabMomentum))*3.5);
      }
    }
    HadrTot = 0.5*(TotP+TotN);
    //  ...................................................
    //  Proton slope
    if(hLabMomentum >= 2.)       { HadrSlope = 5.44 + 0.88*logS; }
    else if(hLabMomentum >= 0.5) { HadrSlope = 3.73*hLabMomentum-0.37; }
    else                         { HadrSlope = 1.5; }

    //  ...................................................
    if(hLabMomentum >= 1.2) {
      HadrReIm  = 0.13*(logS - 5.8579332)*G4Exp(-logS*0.18);
    } else if(hLabMomentum >= 0.6) { 
      HadrReIm = -75.5*(G4Exp(G4Log(hLabMomentum)*0.25)-0.95)/
	(G4Exp(G4Log(3*hLabMomentum)*2.2)+1);     
    } else {
      HadrReIm = 15.5*hLabMomentum/(27*hLabMomentum2*hLabMomentum+2);
    }
    //  ...................................................
    DDSect2   = 2.2;                              //mb*GeV-2
    DDSect3   = 0.6;                               //mb*GeV-2
    //  ================== lambda  ==================
    if( iHadrCode == 3122) {
      HadrTot   *= 0.88;
      HadrSlope *=0.85;
      //  ================== sigma +  ==================
    } else if( iHadrCode == 3222) {
      HadrTot   *=0.81;
      HadrSlope *=0.85;
      //  ================== sigma 0,-  ==================
    } else if(iHadrCode == 3112 || iHadrCode == 3212 ) {
      HadrTot   *=0.88;
      HadrSlope *=0.85;
      //  ===================  xi  =================
    } else if( iHadrCode == 3312 || iHadrCode == 3322 ) {
      HadrTot   *=0.77;
      HadrSlope *=0.75;
      //  =================  omega  =================
    } else if( iHadrCode == 3334) {
      HadrTot   *=0.78;
      HadrSlope *=0.7;
    }
    break;
    //  ===========================================================
  case 1:              //   antiproton
  case 7:              //   antineutron

    HadrTot   = 5.2+5.2*logE + 123.2/sqrS;     //  mb
    HadrSlope = 8.32+0.57*logS;                //(GeV/c)^-2

    if( HadrEnergy < 1000 ) {
      HadrReIm  = 0.06*(sqrS-2.236)*(sqrS-14.14)*G4Exp(-logS*0.8);
    } else {
      HadrReIm  = 0.6*(logS - 5.8579332)*G4Exp(-logS*0.25);
    }
    DDSect2   = 11;                            //mb*(GeV/c)^-2
    DDSect3   = 3;                             //mb*(GeV/c)^-2
    //  ================== lambda  ==================
    if( iHadrCode == -3122) {
      HadrTot   *= 0.88;
      HadrSlope *=0.85;
      //  ================== sigma +  ==================
    } else if( iHadrCode == -3222) {
      HadrTot   *=0.81;
      HadrSlope *=0.85;
      //  ================== sigma 0,-  ==================
    } else if(iHadrCode == -3112 || iHadrCode == -3212 ) {
      HadrTot   *=0.88;
      HadrSlope *=0.85;
    //  ===================  xi  =================
    } else if( iHadrCode == -3312 || iHadrCode == -3322 ) {
      HadrTot   *=0.77;
      HadrSlope *=0.75;
      //  =================  omega  =================
    } else if( iHadrCode == -3334) {
      HadrTot   *=0.78;
      HadrSlope *=0.7;
    }
    break;
    //  -------------------------------------------
  case 2:             //   pi plus, pi minus
  case 3:

    if(hLabMomentum >= 3.5) {
      TotP = 10.6+2.*logE + 25.*G4Exp(-logE*0.43); // mb
      //  =========================================
    } else if(hLabMomentum >= 1.15) {
      G4double x = (hLabMomentum - 2.55)/0.55; 
      G4double y = (hLabMomentum - 1.47)/0.225;
      TotP = 3.2*G4Exp(-x*x) + 12.*G4Exp(-y*y) + 27.5;
      //  =========================================
    } else if(hLabMomentum >= 0.4) {
      TotP  = 88*(logE+0.2877)*(logE+0.2877)+14.0;
    //  =========================================
    } else {
      G4double x = (hLabMomentum - 0.29)/0.085;
      TotP = 20. + 180.*G4Exp(-x*x);
    }
    //  -------------------------------------------

    if(hLabMomentum >= 3.0 ) {
      TotN = 10.6 + 2.*logE + 30.*G4Exp(-logE*0.43); // mb
    } else if(hLabMomentum >= 1.3) {
      G4double x = (hLabMomentum - 2.1)/0.4;
      G4double y = (hLabMomentum - 1.4)/0.12;
      TotN = 36.1+0.079 - 4.313*logE + 3.*G4Exp(-x*x) + 1.5*G4Exp(-y*y);
    } else if(hLabMomentum >= 0.65) {
      G4double x = (hLabMomentum - 0.72)/0.06;
      G4double y = (hLabMomentum - 1.015)/0.075;
      TotN = 36.1 + 10.*G4Exp(-x*x) + 24*G4Exp(-y*y);
    } else if(hLabMomentum >= 0.37) {
      G4double x = G4Log(hLabMomentum/0.48);
      TotN = 26. + 110.*x*x;
    } else {
      G4double x = (hLabMomentum - 0.29)/0.07;
      TotN = 28.0 + 40.*G4Exp(-x*x);
    }
    HadrTot = (TotP+TotN)*0.5;
    //  ........................................
    HadrSlope = 7.28+0.245*logS;        // GeV-2
    HadrReIm  = 0.2*(logS - 4.6051702)*G4Exp(-logS*0.15);

    DDSect2   = 0.7;                               //mb*GeV-2
    DDSect3   = 0.27;                              //mb*GeV-2

    break;
    //  ==========================================================
  case 4:            //  K plus

    HadrTot   = 10.6+1.8*logE + 9.0*G4Exp(-logE*0.55);  // mb
    if(HadrEnergy>100) { HadrSlope = 15.0; }
    else { HadrSlope = 1.0+1.76*logS - 2.84/sqrS; }   // GeV-2

    HadrReIm  = 0.4*(sHadr-20)*(sHadr-150)*G4Exp(-G4Log(sHadr+50)*2.1);
    DDSect2   = 0.7;                             //mb*GeV-2
    DDSect3   = 0.21;                            //mb*GeV-2
    break;
    //  =========================================================
  case 5:              //   K minus

    HadrTot   = 10+1.8*logE + 25./sqrS; // mb
    HadrSlope = 6.98+0.127*logS;         // GeV-2
    HadrReIm  = 0.4*(sHadr-20)*(sHadr-20)*G4Exp(-G4Log(sHadr+50)*2.1);
    DDSect2   = 0.7;                             //mb*GeV-2
    DDSect3   = 0.27;                            //mb*GeV-2
    break;
  }   
  //  =========================================================
  if(verboseLevel>2) {
    G4cout << "HadrTot= " << HadrTot << " HadrSlope= " << HadrSlope
	   << " HadrReIm= " << HadrReIm << " DDSect2= " << DDSect2
	   << " DDSect3= " << DDSect3 << G4endl;
  }
  if(Z != 1) return;

  // Scattering of protons

  Coeff0 = Coeff1 = Coeff2 = 0.0;
  Slope0 = Slope1 = 1.0;
  Slope2 = 5.0;

  // data for iHadron=0
  static const G4double EnP0[6]={1.5,3.0,5.0,9.0,14.0,19.0};
  static const G4double C0P0[6]={0.15,0.02,0.06,0.08,0.0003,0.0002};
  static const G4double C1P0[6]={0.05,0.02,0.03,0.025,0.0,0.0};
  static const G4double B0P0[6]={1.5,2.5,3.0,4.5,1.4,1.25};
  static const G4double B1P0[6]={5.0,1.0,3.5,4.0,4.8,4.8};
      
  // data for iHadron=6,7
  static const G4double EnN[5]={1.5,5.0,10.0,14.0,20.0};
  static const G4double C0N[5]={0.0,0.0,0.02,0.02,0.01};
  static const G4double C1N[5]={0.06,0.008,0.0015,0.001,0.0003};
  static const G4double B0N[5]={1.5,2.5,3.8,3.8,3.5};
  static const G4double B1N[5]={1.5,2.2,3.6,4.5,4.8};

  // data for iHadron=1
  static const G4double EnP[2]={1.5,4.0};
  static const G4double C0P[2]={0.001,0.0005};
  static const G4double C1P[2]={0.003,0.001};
  static const G4double B0P[2]={2.5,4.5};
  static const G4double B1P[2]={1.0,4.0};

  // data for iHadron=2
  static const G4double EnPP[4]={1.0,2.0,3.0,4.0};
  static const G4double C0PP[4]={0.0,0.0,0.0,0.0};
  static const G4double C1PP[4]={0.15,0.08,0.02,0.01};
  static const G4double B0PP[4]={1.5,2.8,3.8,3.8};
  static const G4double B1PP[4]={0.8,1.6,3.6,4.6};

  // data for iHadron=3
  static const G4double EnPPN[4]={1.0,2.0,3.0,4.0};
  static const G4double C0PPN[4]={0.0,0.0,0.0,0.0};
  static const G4double C1PPN[4]={0.0,0.0,0.0,0.0};
  static const G4double B0PPN[4]={1.5,2.8,3.8,3.8};
  static const G4double B1PPN[4]={0.8,1.6,3.6,4.6};

  // data for iHadron=4
  static const G4double EnK[4]={1.4,2.33,3.0,5.0};
  static const G4double C0K[4]={0.0,0.0,0.0,0.0};
  static const G4double C1K[4]={0.01,0.007,0.005,0.003};
  static const G4double B0K[4]={1.5,2.0,3.8,3.8};
  static const G4double B1K[4]={1.6,1.6,1.6,1.6};

  // data for iHadron=5
  static const G4double EnKM[2]={1.4,4.0};
  static const G4double C0KM[2]={0.006,0.002};
  static const G4double C1KM[2]={0.00,0.00};
  static const G4double B0KM[2]={2.5,3.5};
  static const G4double B1KM[2]={1.6,1.6};

  switch(iHadron) {
  case 0:

    if(hLabMomentum <BoundaryP[0]) {
      InterpolateHN(6,EnP0,C0P0,C1P0,B0P0,B1P0);
    }
    Coeff2 = 0.8/hLabMomentum2;
    break; 

  case 6:

    if(hLabMomentum < BoundaryP[1]) {
      InterpolateHN(5,EnN,C0N,C1N,B0N,B1N);
    }
    Coeff2 = 0.8/hLabMomentum2;
    break; 

  case 1:
  case 7:
    if(hLabMomentum <  BoundaryP[2]) {
      InterpolateHN(2,EnP,C0P,C1P,B0P,B1P);
    }
    break; 

  case 2:

    if(hLabMomentum < BoundaryP[3]) {
      InterpolateHN(4,EnPP,C0PP,C1PP,B0PP,B1PP);
    }
    Coeff2 = 0.02/hLabMomentum;
    break; 

  case 3:

    if(hLabMomentum < BoundaryP[4]) {
      InterpolateHN(4,EnPPN,C0PPN,C1PPN,B0PPN,B1PPN);
    }
    Coeff2 = 0.02/hLabMomentum;
    break;
 
  case 4:

    if(hLabMomentum < BoundaryP[5]) {
      InterpolateHN(4,EnK,C0K,C1K,B0K,B1K);
    }
    if(hLabMomentum < 1) { Coeff2 = 0.34; }
    else  { Coeff2 = 0.34/(hLabMomentum2*hLabMomentum); }
    break; 

  case 5:
    if(hLabMomentum < BoundaryP[6]) {
      InterpolateHN(2,EnKM,C0KM,C1KM,B0KM,B1KM);
    }
    if(hLabMomentum < 1) { Coeff2 = 0.01; }
    else  { Coeff2 = 0.01/(hLabMomentum2*hLabMomentum); }
    break; 
  }

  if(verboseLevel > 2) {
    G4cout<<"  HadrVal : Plasb  "<<hLabMomentum
	  <<"  iHadron  "<<iHadron<<"  HadrTot  "<<HadrTot<<G4endl;
  }
}

///////////////////////////////////////////////////////////////////

G4double G4ElasticHadrNucleusHE::GetFt(G4double Q2)
{
  G4double Fdistr=0;
  G4double SqrQ2 = std::sqrt(Q2);
 
  Fdistr = (1-Coeff1-Coeff0) / HadrSlope*(1-G4Exp(-HadrSlope*Q2))
    + Coeff0*(1-G4Exp(-Slope0*Q2))
    + Coeff2/Slope2*G4Exp(Slope2*ConstU)*(G4Exp(Slope2*Q2)-1)
    + 2*Coeff1/Slope1*(1/Slope1-(1/Slope1+SqrQ2)*G4Exp(-Slope1*SqrQ2));

  if (verboseLevel>1) {
    G4cout<<"Old:  Coeff0 Coeff1 Coeff2 "<<Coeff0<<"  "
          <<Coeff1<<"  "<<Coeff2<<"  Slope Slope0 Slope1 Slope2 "
          <<HadrSlope<<"  "<<Slope0<<"  "<<Slope1<<"  "<<Slope2
          <<"  Fdistr "<<Fdistr<<G4endl;
  }
  return Fdistr;
}

///////////////////////////////////////////////////////////////////

G4double 
G4ElasticHadrNucleusHE::HadronProtonQ2(G4double plab, G4double tmax)
{
  hLabMomentum  = plab;
  hLabMomentum2 = hLabMomentum*hLabMomentum;
  HadrEnergy    = std::sqrt(hMass2 + hLabMomentum2);
  DefineHadronValues(1);

  G4double Sh = 2.0*protonM*HadrEnergy+protonM2+hMass2; // GeV
  ConstU = 2*protonM2+2*hMass2-Sh;

  BoundaryTL[0] = tmax;
  BoundaryTL[1] = tmax;
  BoundaryTL[3] = tmax;
  BoundaryTL[4] = tmax;
  BoundaryTL[5] = tmax;

  G4double MaxTR = (plab < BoundaryP[iHadron1]) ? 
    BoundaryTL[iHadron1] : BoundaryTG[iHadron1]; 

  if (verboseLevel>1) {
    G4cout<<"3  GetKin. : iHadron1  "<<iHadron1
	  <<"  Bound.P[iHadron1] "<<BoundaryP[iHadron1]
	  <<"  Bound.TL[iHadron1] "<<BoundaryTL[iHadron1]
	  <<"  Bound.TG[iHadron1] "<<BoundaryTG[iHadron1]
	  <<"  MaxT MaxTR "<<tmax<<"  "<<MaxTR<<G4endl;
  }

  G4double rand = G4UniformRand();

  G4double DDD0=MaxTR*0.5, DDD1=0.0, DDD2=MaxTR;

  G4double norm  = 1.0/GetFt(MaxTR);
  G4double delta = GetFt(DDD0)*norm - rand;

  static const G4int maxNumberOfLoops = 10000;
  G4int loopCounter = -1;
  while ( (std::abs(delta) > 0.0001) && 
          ++loopCounter < maxNumberOfLoops )  /* Loop checking, 10.08.2015, A.Ribon */
    {
      if(delta>0) 
      {
        DDD2 = DDD0;
        DDD0 = (DDD0+DDD1)*0.5;
      }
      else if(delta<0.0)
      {
        DDD1 = DDD0; 
        DDD0 = (DDD0+DDD2)*0.5;
      }
      delta = GetFt(DDD0)*norm - rand;
    }
  return (loopCounter >= maxNumberOfLoops) ? 0.0 : DDD0;
}

///////////////////////////////////////////////////////////////////

void G4ElasticHadrNucleusHE::Binom()
{
  for(G4int N = 0; N < 240; ++N) {
    G4double J = 1.0;
    for(G4int M = 0; M <= N; ++M) {
      G4double Fact1 = 1.0;
      if (N > 0 && N > M && M > 0 ) {
	J *= (G4double)(N-M+1)/(G4double)M;
	Fact1 = J;
      }
      fBinom[N][M] = Fact1;
    }
  }
}

///////////////////////////////////////////////////////////

void 
G4ElasticHadrNucleusHE::InFileName(std::ostringstream& ss, 
				   const G4ParticleDefinition* p, G4int Z)
{
  if(!fDirectory) {
    fDirectory = G4FindDataDir("G4LEDATA");
    if (fDirectory) { 
      ss << fDirectory << "/";
    }
  }
  OutFileName(ss, p, Z);
}

///////////////////////////////////////////////////////////

void 
G4ElasticHadrNucleusHE::OutFileName(std::ostringstream& ss, 
				    const G4ParticleDefinition* p, G4int Z)
{
  ss << "hedata/" << p->GetParticleName() << Z << ".dat";
}

///////////////////////////////////////////////////////////

G4bool G4ElasticHadrNucleusHE::ReadLine(std::ifstream& infile, 
					std::vector<G4double>& v)
{
  G4int n(0);
  infile >> n;
  if (infile.fail()) { return false; }
  if(n > 0) {
    v.reserve(n);
    G4double x(0.0);
    for(G4int i=0; i<n; ++i) {
      infile >> x;
      if (infile.fail()) { return false; }
      v.emplace_back(x);
    }
  }
  return true;
}

///////////////////////////////////////////////////////////

void G4ElasticHadrNucleusHE::WriteLine(std::ofstream& outfile, 
				       std::vector<G4double>& v)
{
  std::size_t n = v.size();
  outfile << n << G4endl;
  if(n > 0) {
    for(std::size_t i=0; i<n; ++i) {
      outfile << v[i] << " ";
    }
    outfile << G4endl;
  }
}

///////////////////////////////////////////////////////////
