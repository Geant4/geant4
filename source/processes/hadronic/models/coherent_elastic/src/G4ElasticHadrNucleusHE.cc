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
// $Id: G4ElasticHadrNucleusHE.cc 94236 2015-11-09 11:00:13Z gcosmo $
//
//
//  The generator of high energy hadron-nucleus elastic scattering
//  The hadron kinetic energy T > 1 GeV
//  N.  Starkov 2003.
//
//  19.05.04 Variant for G4 6.1: The 'ApplyYourself' was changed
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
#include  "G4IonTable.hh"
#include  "G4Proton.hh"
#include  "G4NistManager.hh"
#include  "G4Log.hh"
#include  "G4Exp.hh"

#include <cmath>

using namespace std;

//ANDREA-> MT Fix
#include "G4AutoLock.hh"
G4ElasticData* G4ElasticHadrNucleusHE::SetOfElasticData[NHADRONS][ZMAX];
G4Mutex G4ElasticHadrNucleusHE::eldata_m[NHADRONS][ZMAX];
namespace {
  G4Mutex aMutex = G4MUTEX_INITIALIZER;
  G4bool onlyOnceInit  = true;
  G4bool onlyOnceDestroy = true;
}
//ANDREA<-

///////////////////////////////////////////////////////////////
//
//

G4ElasticData::G4ElasticData(const G4ParticleDefinition* p, 
			     G4int Z, G4double AWeight, G4double* eGeV)
{ 
  hadr     = p;
  massGeV  = p->GetPDGMass()/GeV;
  mass2GeV2= massGeV*massGeV;
  AtomicWeight = G4lrint(AWeight);

  DefineNucleusParameters(AWeight);

  limitQ2 = 35./(R1*R1);     //  (GeV/c)^2

  G4double dQ2 = limitQ2/(ONQ2 - 1.);

  TableQ2[0] = 0.0;

  for(G4int ii = 1; ii < ONQ2; ii++) 
  {
    TableQ2[ii] = TableQ2[ii-1]+dQ2;
  }

  massA  = AWeight*amu_c2/GeV;
  massA2 = massA*massA; 

  for(G4int kk = 0; kk < NENERGY; kk++) 
  {
    dnkE[kk] = 0;
    G4double elab = eGeV[kk] + massGeV;
    G4double plab2= eGeV[kk]*(eGeV[kk] + 2.0*massGeV);
    G4double Q2m  = 4.0*plab2*massA2/(mass2GeV2 + massA2 + 2.*massA2*elab);

    if(Z == 1 && p == G4Proton::Proton()) Q2m *= 0.5;

    maxQ2[kk] = std::min(limitQ2, Q2m);
    TableCrossSec[ONQ2*kk] = 0.0;

//    G4cout<<" kk  eGeV[kk] "<<kk<<"  "<<eGeV[kk]<<G4endl;
  }
}

/////////////////////////////////////////////////////////////////////////
//
//

void G4ElasticData::DefineNucleusParameters(G4double A)
{
  switch (AtomicWeight)
    {
    case 207:
    case 208:
      R1       = 20.5;
//      R1       = 17.5;
//      R1       = 21.3;    
      R2       = 15.74;
//      R2       = 10.74;

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
      R1    = 16.5*1.0;
      R2    = 11.62;
      Pnucl = 0.4;
      Aeff  = 0.7;
      break;
    case 58:
    case 59:
      R1    = 15.0*1.05;
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
//      Pnucl = 0.5397;
      Aeff  = 1.0;
      break;
    case 11:
      R1    = 9.0;
      R2    = 5.42;
      Pnucl = 0.19;
//      Pnucl = 0.5397;
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
      R1    = 4.45*G4Exp(G4Log(A - 1.)*0.309)*0.9;
      R2    = 2.3 *G4Exp(G4Log(A)* 0.36);

      if(A < 100 && A > 3) Pnucl = 0.176 + 0.00275*A;
      else                 Pnucl = 0.4;

//G4cout<<" Deault: A= "<<A<<"  R1 R2 Aeff Pnucl "<<R1<<"  "<<R2<<"  "
//      <<Aeff<<"  "<<Pnucl<<G4endl;

      if(A >= 100)               Aeff = 0.7;
      else if(A < 100 && A > 75) Aeff = 1.5 - 0.008*A;
      else                       Aeff = 0.9;
      break;
    }
//G4cout<<" Result: A= "<<A<<"  R1 R2 Aeff Pnucl "<<R1<<"  "<<R2<<"  "
//      <<Aeff<<"  "<<Pnucl<<G4endl;
}

////////////////////////////////////////////////////////////////////
//
//  The constructor for the generating of events
//

G4ElasticHadrNucleusHE::G4ElasticHadrNucleusHE(const G4String& name)
  : G4HadronElastic(name)
{
  //ANDREA->
  G4AutoLock l(&aMutex);
  if ( onlyOnceInit ) {
      for ( int i = 0 ; i< NHADRONS ; ++i) {
          for (int j = 0 ; j<ZMAX ; ++j ) {
              SetOfElasticData[i][j]=0;
              G4MUTEXINIT(eldata_m[i][j]);
          }
      }
      onlyOnceInit = false;
  }
  l.unlock();
  //ANDREA<-

  dQ2 = hMass = hMass2 =  hLabMomentum = hLabMomentum2 = MomentumCM = HadrEnergy 
    = R1 = R2 = Pnucl = Aeff = HadrTot = HadrSlope = HadrReIm = TotP = DDSect2
    = DDSect3 = ConstU = FmaxT = Slope1 = Slope2 = Coeff1 = Coeff2 = MaxTR 
    = Slope0 = Coeff0 = aAIm = aDIm = Dtot11 = 0.0;
  NumbN = iHadrCode = iHadron = 0;

  verboseLevel = 0;
  plabLowLimit = 20.0*MeV;
  lowestEnergyLimit = 0.0;
  //Description();

  MbToGeV2  =  2.568;
  sqMbToGeV =  1.602;
  Fm2ToGeV2 =  25.68;
  GeV2      =  GeV*GeV;
  protonM   =  proton_mass_c2/GeV;
  protonM2  =  protonM*protonM;

  BoundaryP[0]=9.0;BoundaryTG[0]=5.0;BoundaryTL[0]=0.;
  BoundaryP[1]=20.0;BoundaryTG[1]=1.5;BoundaryTL[1]=0.;
  BoundaryP[2]=5.0; BoundaryTG[2]=1.0;BoundaryTL[2]=1.5;
  BoundaryP[3]=8.0; BoundaryTG[3]=3.0;BoundaryTL[3]=0.;
  BoundaryP[4]=7.0; BoundaryTG[4]=3.0;BoundaryTL[4]=0.;
  BoundaryP[5]=5.0; BoundaryTG[5]=2.0;BoundaryTL[5]=0.;
  BoundaryP[6]=5.0; BoundaryTG[6]=1.5;BoundaryTL[6]=3.0;

  Binom();
  // energy in GeV
  Energy[0] = 0.4;
  Energy[1] = 0.6;
  Energy[2] = 0.8;
  LowEdgeEnergy[0] = 0.0;
  LowEdgeEnergy[1] = 0.5;
  LowEdgeEnergy[2] = 0.7;
  G4double e = 1.0;
  G4double f = G4Exp(G4Log(10.)*0.1);
  for(G4int i=3; i<NENERGY; i++) {
    Energy[i] = e;
    LowEdgeEnergy[i] = e/f;
    e *= f*f;
  }
  nistManager = G4NistManager::Instance();

  // PDG code for hadrons
  G4int ic[NHADRONS] = {211,-211,2112,2212,321,-321,130,310,311,-311,
			3122,3222,3112,3212,3312,3322,3334,
			-2212,-2112,-3122,-3222,-3112,-3212,-3312,-3322,-3334};
  // internal index 
  G4int id[NHADRONS] = {2,3,6,0,4,5,4,4,4,5,
			0,0,0,0,0,0,0,
			1,7,1,1,1,1,1,1,1};

  G4int id1[NHADRONS] = {3,4,1,0,5,6,5,5,5,6,
                        0,0,0,0,0,0,0,
                        2,2,2,2,2,2,2,2,2};

  for(G4int j=0; j<NHADRONS; j++) 
  {
    HadronCode[j]  = ic[j];
    HadronType[j]  = id[j];
    HadronType1[j] = id1[j];

    for(G4int k = 0; k < ZMAX; k++) { SetOfElasticData[j][k] = 0; }
  } 
}


void G4ElasticHadrNucleusHE::ModelDescription(std::ostream& outFile) const
{

    outFile << "G4ElasticHadrNucleusHE is a hadron-nucleus elastic scattering\n"
            << "model developed by N. Starkov which uses a Glauber model\n"
            << "parameterization to calculate the final state.  It is valid\n"
            << "for all hadrons with incident energies above 1 GeV.\n";

}


///////////////////////////////////////////////////////////////////
//
//

G4ElasticHadrNucleusHE::~G4ElasticHadrNucleusHE()
{
    //ANDREA->
  G4AutoLock l(&aMutex);
  if ( onlyOnceDestroy ) {
      for(G4int j = 0; j < NHADRONS; j++)
        {
          for(G4int k = 0; k < ZMAX; k++)
            {
              if ( SetOfElasticData[j][k] ) {
                  delete SetOfElasticData[j][k];
                  SetOfElasticData[j][k]=0;
                  G4MUTEXDESTROY(eldata_m[j][k]);
              }
            }
        }
      onlyOnceDestroy = false;
  }
}

////////////////////////////////////////////////////////////////////
//
//

G4double 
G4ElasticHadrNucleusHE::SampleInvariantT(const G4ParticleDefinition* p,
					 G4double inLabMom, 
					 G4int iZ, G4int N)
{
  G4int Z = iZ;
  if(Z >= ZMAX) { Z = ZMAX-1; }
  G4double plab  = inLabMom/GeV;   // (GeV/c)
  G4double Q2 = 0;

  iHadrCode = p->GetPDGEncoding();

  NumbN = N;

  if(verboseLevel > 1)
  {
    G4cout<< " G4ElasticHadrNucleusHE::SampleT: " 
	  << " for " << p->GetParticleName() 
	  << " at Z= " << Z << " A= " << N
	  << " plab(GeV)= " << plab
	  << G4endl;
  }
  G4int idx;

  for( idx = 0 ; idx < NHADRONS; idx++) 
  {
    if(iHadrCode == HadronCode[idx]) break;
  }

  // Hadron is not in the list

  if( idx >= NHADRONS ) return Q2;

  iHadron = HadronType[idx];
  iHadrCode = HadronCode[idx];

  if(Z==1)
    {
      hMass  = p->GetPDGMass()/GeV;
      hMass2 = hMass*hMass;

      G4double T = sqrt(plab*plab+hMass2)-hMass;

      if(T > 0.4) Q2 = HadronProtonQ2(p, plab);

      if (verboseLevel>1)
	G4cout<<"  Proton : Q2  "<<Q2<<G4endl;
    }
  else
    {
      G4AutoLock l(&(eldata_m[idx][Z]));//ANDREA
      G4ElasticData* ElD1 = SetOfElasticData[idx][Z];

      // Construct elastic data
      if(!ElD1) 
	{
	  G4double AWeight = nistManager->GetAtomicMassAmu(Z);
	  ElD1 = new  G4ElasticData(p, Z, AWeight, Energy);
	  SetOfElasticData[idx][Z] = ElD1;
    
	  if(verboseLevel > 1)
	    {
	      G4cout<< " G4ElasticHadrNucleusHE::SampleT:  new record " << idx
		    << " for " << p->GetParticleName() << " Z= " << Z
		    << G4endl;
	    }
	}  
      hMass          = ElD1->massGeV;
      hMass2         = ElD1->mass2GeV2;
      G4double M     = ElD1->massA;
      G4double M2    = ElD1->massA2;
      G4double plab2 = plab*plab;
      G4double Q2max = 4.*plab2*M2/
	(hMass2 + M2 + 2.*M*std::sqrt(plab2 + hMass2));

      // sample scattering
      G4double T = sqrt(plab2+hMass2)-hMass;

      if(T > 0.4) Q2 = HadronNucleusQ2_2(ElD1, Z, plab, Q2max);

      if(verboseLevel > 1)
	G4cout<<" SampleT: Q2(GeV^2)= "<<Q2<< "  t/tmax= " << Q2/Q2max <<G4endl;
    }
  return  Q2*GeV2;
}

//////////////////////////////////////////////////////////////////////////
//
//

G4double 
G4ElasticHadrNucleusHE::SampleT(const G4ParticleDefinition* p,
				G4double inLabMom, 
				G4int Z, G4int N)
{
  return SampleInvariantT(p, inLabMom, Z, N);
}

//////////////////////////////////////////////////////////////////////////
//
//

G4double G4ElasticHadrNucleusHE::
                          HadronNucleusQ2_2(G4ElasticData* pElD, G4int Z, 
                                            G4double plab, G4double tmax)
{
  //ANDREA: Important notice on this function
  //        For MT we are sharing among threads the G4ElasticData classes
  //        The *call* to this function is proteced
  //        with a mutex in the calling function
  G4double LineFq2[ONQ2];

  G4double Rand = G4UniformRand();

  G4int      iNumbQ2 = 0;
  G4double   Q2 = 0.0;

  G4double ptot2 = plab*plab;
  G4double ekin  = std::sqrt(hMass2 + ptot2) - hMass;

  if(verboseLevel > 1)
    G4cout<<"Q2_2: ekin  plab  "<<ekin<<"    "<<plab<<"  tmax "<<tmax<<G4endl;

  // Find closest energy bin
  G4int NumbOnE; 
  for( NumbOnE = 0; NumbOnE < NENERGY-1; NumbOnE++ ) 
  {
    if( ekin <= LowEdgeEnergy[NumbOnE+1] ) break;
  }
  G4double* dNumbQ2 = pElD->TableQ2;

  if(NumbOnE >= NENERGY-1) { NumbOnE = NENERGY-2; }
  G4int index = NumbOnE*ONQ2;

  // Select kinematics for node energy
  G4double T     = Energy[NumbOnE];
  hLabMomentum2  = T*(T + 2.*hMass);
  G4double Q2max = pElD->maxQ2[NumbOnE];
  G4int length   = pElD->dnkE[NumbOnE];

  // Build vector
  if(length == 0) 
    {
      R1    = pElD->R1;
      R2    = pElD->R2;
      Aeff  = pElD->Aeff;
      Pnucl = pElD->Pnucl;
      hLabMomentum = std::sqrt(hLabMomentum2);
 
      DefineHadronValues(Z);

      if(verboseLevel >0)
	{
	  G4cout<<"1  plab  T  "<<plab<<"  "<<T<<"  sigTot  B  ReIm  "
		<<HadrTot<<"  "<<HadrSlope<<"  "<<HadrReIm<<G4endl;
	  G4cout<<"  R1  R2  Aeff  p  "<<R1<<"  "<<R2<<"  "<<Aeff<<"  "
		<<Pnucl<<G4endl;
	}

      //pElD->CrossSecMaxQ2[NumbOnE] = 1.0;

      if(verboseLevel > 1)
	G4cout<<" HadrNucleusQ2_2: NumbOnE= " << NumbOnE 
	      << " length= " << length 
	      << " Q2max= " << Q2max 
	      << " ekin= " << ekin <<G4endl;
    
      pElD->TableCrossSec[index] = 0;


      dQ2 = pElD->TableQ2[1]-pElD->TableQ2[0];

      GetHeavyFq2(Z, NumbN, LineFq2);  //  %%%%%%%%%%%%%%%%%%%%%%%%%

      for(G4int ii=0; ii<ONQ2; ++ii)
	{
	  //if(verboseLevel > 2)
	  //  G4cout<<"  ii LineFq2  "<<ii<<"  "<<LineFq2[ii]/LineFq2[ONQ2-1]
	  //	<<"  dF(q2) "<<LineFq2[ii]-LineFq2[ii-1]<<G4endl;

	  pElD->TableCrossSec[index+ii] = LineFq2[ii]/LineFq2[ONQ2-1];
	}
    
      pElD->dnkE[NumbOnE] = ONQ2;
      length = ONQ2;
    } 

  G4double* dNumbFQ2 = &(pElD->TableCrossSec[index]);

  for( iNumbQ2 = 1; iNumbQ2<length; ++iNumbQ2) 
    {
      if(Rand <= pElD->TableCrossSec[index+iNumbQ2]) break;
    }
  if(iNumbQ2 >= ONQ2) { iNumbQ2 = ONQ2 - 1; }
  Q2 = GetQ2_2(iNumbQ2, dNumbQ2, dNumbFQ2, Rand);

  if(tmax < Q2max) Q2 *= tmax/Q2max;

  if(verboseLevel > 1)
    G4cout<<" HadrNucleusQ2_2(2): Q2= "<<Q2<<" iNumbQ2= " << iNumbQ2 
	  << " rand= " << Rand << G4endl;
  
  return Q2;
}       

///////////////////////////////////////////////////////////////////////
//
//  The randomization of one dimensional array 
//
G4double G4ElasticHadrNucleusHE::GetQ2_2(G4int kk, G4double * Q,
					 G4double * F, G4double ranUni)
{
  G4double ranQ2;

  G4double F1  = F[kk-2];
  G4double F2  = F[kk-1];
  G4double F3  = F[kk];
  G4double X1  = Q[kk-2];
  G4double X2  = Q[kk-1];
  G4double X3  = Q[kk];

  if(verboseLevel > 2) 
    G4cout << "GetQ2_2 kk= " << kk << " X2= " << X2 << " X3= " << X3 
	   << " F2= " << F2 << " F3= " << F3 << " Rndm= " << ranUni << G4endl;

  if(kk == 1 || kk == 0)
  {
     F1  = F[0]; 
     F2  = F[1];
     F3  = F[2];
     X1  = Q[0];
     X2  = Q[1];
     X3  = Q[2];
  }

  G4double F12 = F1*F1;
  G4double F22 = F2*F2;
  G4double F32 = F3*F3;

  G4double D0  = F12*F2+F1*F32+F3*F22-F32*F2-F22*F1-F12*F3;

  if(verboseLevel > 2) 
    G4cout << "       X1= " << X1 << " F1= " << F1 << "  D0= " 
           << D0 << G4endl; 

  if(std::abs(D0) < 0.00000001)
    { 
      ranQ2 = X2 + (ranUni - F2)*(X3 - X2)/(F3 - F2);
    }
  else    
    {
      G4double DA = X1*F2+X3*F1+X2*F3-X3*F2-X1*F3-X2*F1;
      G4double DB = X2*F12+X1*F32+X3*F22-X2*F32-X3*F12-X1*F22;
      G4double DC = X3*F2*F12+X2*F1*F32+X1*F3*F22
	           -X1*F2*F32-X2*F3*F12-X3*F1*F22;
      ranQ2 = (DA*ranUni*ranUni + DB*ranUni + DC)/D0;
    }
  return ranQ2;         //  MeV^2
}

////////////////////////////////////////////////////////////////////////
//
//
G4double G4ElasticHadrNucleusHE::GetHeavyFq2(G4int Z, G4int Nucleus, G4double* LineF) 
{
  G4int ii, jj, aSimp;
  G4double curQ2, curSec;
  G4double curSum = 0.0;
  G4double totSum = 0.0;

  G4double ddQ2 = dQ2/20;
  G4double Q2l  = 0;

  LineF[0] = 0;
  for(ii = 1; ii<ONQ2; ii++)
    {
      curSum = 0;
      aSimp  = 4;   

      for(jj = 0; jj<20; jj++)
	{
	  curQ2 = Q2l+jj*ddQ2;

	  curSec  = HadrNucDifferCrSec(Z, Nucleus, curQ2);
	  curSum += curSec*aSimp;

	  if(aSimp > 3) aSimp = 2;
	  else          aSimp = 4;

	  if(jj == 0 && verboseLevel>2)
	    G4cout<<"  Q2  "<<curQ2<<"  AIm  "<<aAIm<<"  DIm  "<<aDIm
		  <<"  Diff  "<<curSec<<"  totSum  "<<totSum<<G4endl;
	}

      Q2l    += dQ2;
      curSum *= ddQ2/2.3;   //  $$$$$$$$$$$$$$$$$$$$$$$
      totSum += curSum;

      LineF[ii] = totSum;
	
      if (verboseLevel>2)
	G4cout<<"  GetHeavy: Q2  dQ2  totSum  "<<Q2l<<"  "<<dQ2<<"  "<<totSum
	      <<"  curSec  "
	      <<curSec<<"  totSum  "<< totSum<<"  DTot "
	      <<curSum<<G4endl;
    }      
  return totSum;
}

////////////////////////////////////////////////////////////////////////
//
//

G4double G4ElasticHadrNucleusHE::GetLightFq2(G4int Z, G4int Nucleus, 
                                             G4double Q2)
{
  // Scattering of proton
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

  G4double prec = Nucleus > 208  ?  1.0e-7 : 1.0e-6;

  G4double    Stot     = HadrTot*MbToGeV2;     //  Gev^-2
  G4double    Bhad     = HadrSlope;         //  GeV^-2
  G4double    Asq      = 1+HadrReIm*HadrReIm;
  G4double    Rho2     = std::sqrt(Asq);

//  if(verboseLevel >1)
    G4cout<<" Fq2 Before for i Tot B Im "<<HadrTot<<"  "<<HadrSlope<<"  "
      <<HadrReIm<<G4endl;

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
  G4double    Unucl    = Stot/twopi/Norm*R13;
  G4double    UnucRho2 = -Unucl*Rho2;

  G4double    FiH      = std::asin(HadrReIm/Rho2);
  G4double    NN2      = R23/R13;

  if(verboseLevel > 2) 
    G4cout << "UnucRho2= " << UnucRho2 << " FiH= " << FiH << " NN2= " << NN2 
	   << " Norm= " << Norm << G4endl;

  G4double    dddd;
 
  G4double    Prod0    = 0;
  G4double    N1       = -1.0;
  //G4double    Tot0     = 0;
  G4double    exp1;

  G4double    Prod3 ;
  G4double    exp2  ;
  G4double    N4, N5, N2, Prod1, Prod2;
  G4int    i1, i2, j1, j2;

  for(i1 = 1; i1<= Nucleus; i1++) ////++++++++++  i1
    {
      N1    = -N1*Unucl*(Nucleus-i1+1)/i1*Rho2;
      Prod1 = 0;
      //Tot0  = 0;
      N2    = -1;

      for(i2 = 1; i2<=Nucleus; i2++) ////+++++++++ i2
        {
          N2    = -N2*Unucl*(Nucleus-i2+1)/i2*Rho2;
          Prod2 = 0; 
          N5    = -1/NN2;
	  for(j2=0; j2<= i2; j2++) ////+++++++++ j2
            {
              Prod3 = 0;
              exp2  = 1/(j2/R22B+(i2-j2)/R12B);
              N5    = -N5*NN2;
              N4    = -1/NN2;
	      for(j1=0; j1<=i1; j1++) ////++++++++ j1
		{
		  exp1  = 1/(j1/R22B+(i1-j1)/R12B);
		  dddd  = exp1+exp2;
		  N4    = -N4*NN2;
		  Prod3 = Prod3+N4*exp1*exp2*
		    (1-G4Exp(-Q2*dddd/4))/dddd*4*SetBinom[i1][j1];
               }                                   // j1
	      Prod2 = Prod2 +Prod3*N5*SetBinom[i2][j2];
	    }                                      // j2
	  Prod1 = Prod1 + Prod2*N2*std::cos(FiH*(i1-i2));

	  if (std::fabs(Prod2*N2/Prod1)<prec) break;
        }                                         // i2
      Prod0   = Prod0 + Prod1*N1;
      if(std::fabs(N1*Prod1/Prod0) < prec) break;

    }                                           // i1

  Prod0 *= 0.25*pi/MbToGeV2;  //  This is in mb

  if(verboseLevel>1) 
    G4cout << "GetLightFq2 Z= " << Z << " A= " << Nucleus 
	   <<" Q2= " << Q2 << " Res= " << Prod0 << G4endl;
  return Prod0;
}
//  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
G4double G4ElasticHadrNucleusHE::
   HadrNucDifferCrSec(G4int Z, G4int , G4double aQ2)
{
//   ------ All external kinematical variables are in MeV -------
//            ------ but internal in GeV !!!!  ------

  G4int NWeight = G4lrint(nistManager->GetAtomicMassAmu(Z)); 

  G4double    theQ2 = aQ2;   ///GeV/GeV;  

  // Scattering of proton
  if(NWeight == 1) 
  {
    G4double SqrQ2  = std::sqrt(aQ2);
    G4double valueConstU = hMass2 + protonM2-2*protonM*HadrEnergy - aQ2;

    G4double MaxT = 4*MomentumCM*MomentumCM;

     BoundaryTL[0] = MaxT;
     BoundaryTL[1] = MaxT;
     BoundaryTL[3] = MaxT;
     BoundaryTL[4] = MaxT;
     BoundaryTL[5] = MaxT;

    G4double dSigPodT;

    dSigPodT = HadrTot*HadrTot*(1+HadrReIm*HadrReIm)*
                 (
                  Coeff1*G4Exp(-Slope1*SqrQ2)+
                  Coeff2*G4Exp( Slope2*(valueConstU)+aQ2)+
                  (1-Coeff1-Coeff0)*G4Exp(-HadrSlope*aQ2)+
                 +Coeff0*G4Exp(-Slope0*aQ2)
//                +0.1*(1-std::fabs(CosTh))
                  )*2.568/(16*pi);

    return dSigPodT;
  }

    G4double    Stot     = HadrTot*MbToGeV2; 
    G4double    Bhad     = HadrSlope; 
    G4double    Asq      = 1+HadrReIm*HadrReIm;
    G4double    Rho2     = std::sqrt(Asq);
    G4double    Pnuclp   = 0.001;
                Pnuclp   = Pnucl;
    G4double    R12      = R1*R1;
    G4double    R22      = R2*R2;
    G4double    R12B     = R12+2*Bhad;
    G4double    R22B     = R22+2*Bhad;
    G4double    R12Ap    = R12+20;
    G4double    R22Ap    = R22+20;
    G4double    R13Ap    = R12*R1/R12Ap;
    G4double    R23Ap    = R22*R2/R22Ap*Pnuclp;
    G4double    R23dR13  = R23Ap/R13Ap;
    G4double    R12Apd   = 2/R12Ap;
    G4double    R22Apd   = 2/R22Ap;
    G4double R12ApdR22Ap = 0.5*(R12Apd+R22Apd);

    G4double DDSec1p  = (DDSect2+DDSect3*G4Log(0.53*HadrEnergy/R1));
    G4double DDSec2p  = (DDSect2+DDSect3*G4Log(0.53*HadrEnergy/
                             std::sqrt((R12+R22)/2)));
    G4double DDSec3p  = (DDSect2+DDSect3*G4Log(0.53*HadrEnergy/R2));

    G4double    Norm     = (R12*R1-Pnucl*R22*R2)*Aeff;
    G4double    Normp    = (R12*R1-Pnuclp*R22*R2)*Aeff;
    G4double    R13      = R12*R1/R12B;
    G4double    R23      = Pnucl*R22*R2/R22B;
    G4double    Unucl    = Stot/(twopi*Norm)*R13;
    G4double    UnuclScr = Stot/(twopi*Normp)*R13Ap;
    G4double    SinFi    = HadrReIm/Rho2;
    G4double    FiH      = std::asin(SinFi);
    G4double    N        = -1;
    G4double    N2       = R23/R13;

    G4double    ImElasticAmpl0  = 0;
    G4double    ReElasticAmpl0  = 0;

    G4double    exp1;
    G4double    N4;
    G4double    Prod1, Tot1=0, medTot, DTot1, DmedTot;
    G4int       i;

    for( i=1; i<=NWeight; i++)
    {
      N       = -N*Unucl*(NWeight-i+1)/i*Rho2;
      N4      = 1;
      Prod1   = G4Exp(-theQ2/i*R12B/4)/i*R12B;
      medTot  = R12B/i;

       for(G4int l=1; l<=i; l++)
       {
         exp1    = l/R22B+(i-l)/R12B;
         N4      = -N4*(i-l+1)/l*N2;
         Prod1   = Prod1+N4/exp1*G4Exp(-theQ2/(exp1*4));
         medTot  = medTot+N4/exp1;
       }  // end l

      ReElasticAmpl0  = ReElasticAmpl0+Prod1*N*std::sin(FiH*i);
      ImElasticAmpl0  = ImElasticAmpl0+Prod1*N*std::cos(FiH*i);
      Tot1            = Tot1+medTot*N*std::cos(FiH*i);
      if(std::fabs(Prod1*N/ImElasticAmpl0) < 0.000001) break;
    }      // i

    ImElasticAmpl0 = ImElasticAmpl0*pi/2.568;   // The amplitude in mB
    ReElasticAmpl0 = ReElasticAmpl0*pi/2.568;   // The amplitude in mB
    Tot1           = Tot1*twopi/2.568;

    G4double C1 = R13Ap*R13Ap*0.5*DDSec1p;
    G4double C2 = 2*R23Ap*R13Ap*0.5*DDSec2p;
    G4double C3 = R23Ap*R23Ap*0.5*DDSec3p;

    G4double N1p  = 1;

    G4double Din1 = 0.5;     

    Din1  = 0.5*(C1*G4Exp(-theQ2/8*R12Ap)/2*R12Ap-
                 C2/R12ApdR22Ap*G4Exp(-theQ2/(4*R12ApdR22Ap))+
                 C3*R22Ap/2*G4Exp(-theQ2/8*R22Ap));

    DTot1 = 0.5*(C1/2*R12Ap-C2/R12ApdR22Ap+C3*R22Ap/2);

    G4double exp1p;
    G4double exp2p;
    G4double exp3p;
    G4double N2p;
    G4double Din2, BinCoeff;

    BinCoeff = 1;

    for( i = 1; i<= NWeight-2; i++)
    {
      N1p     = -N1p*UnuclScr*(NWeight-i-1)/i*Rho2;
      N2p     = 1;
      Din2    = 0;
      DmedTot = 0;
        for(G4int l = 0; l<=i; l++) 
        {
          if(l == 0)      BinCoeff = 1;
          else if(l !=0 ) BinCoeff = BinCoeff*(i-l+1)/l;

          exp1  = l/R22B+(i-l)/R12B;
          exp1p = exp1+R12Apd;
          exp2p = exp1+R12ApdR22Ap;
          exp3p = exp1+R22Apd;

          Din2  = Din2 + N2p*BinCoeff*
	    (C1/exp1p*G4Exp(-theQ2/(4*exp1p))-
	     C2/exp2p*G4Exp(-theQ2/(4*exp2p))+
	     C3/exp3p*G4Exp(-theQ2/(4*exp3p)));

	  DmedTot = DmedTot + N2p*BinCoeff*
	    (C1/exp1p-C2/exp2p+C3/exp3p);

	  N2p   = -N2p*R23dR13;
	}     // l

	Din1  = Din1+Din2*N1p/*Mnoj[i]*//((i+2)*(i+1))*std::cos(FiH*i);
	DTot1 = DTot1+DmedTot*N1p/*Mnoj[i]*//((i+2)*(i+1))*std::cos(FiH*i);
 
	if(std::fabs(Din2*N1p/Din1) < 0.000001) break;
    }           //  i

    Din1 = -Din1*NWeight*(NWeight-1)*4/(Normp*Normp);

    DTot1 = DTot1*NWeight*(NWeight-1)*4/(Normp*Normp);

    DTot1 *= 5;   //  $$$$$$$$$$$$$$$$$$$$$$$$ 
//     Din1  *= 0.2;    //   %%%%%%%%%%%%%%%%%%%%%%%   proton
//     Din1 *= 0.05;    //   %%%%%%%%%%%%%%%%%%%%%%%  pi+
//  ----------------  dSigma/d|-t|,  mb/(GeV/c)^-2  -----------------

    G4double DiffCrSec2 = (ReElasticAmpl0*ReElasticAmpl0+
                           (ImElasticAmpl0+Din1)*
                           (ImElasticAmpl0+Din1))/twopi;

    Tot1   = Tot1-DTot1;
     //  Tott1  = Tot1*1.0;
    Dtot11 = DTot1;
    aAIm   = ImElasticAmpl0;
    aDIm   = Din1;

    return DiffCrSec2;  //  dSig/d|-t|,  mb/(GeV/c)^-2
}   // function
//  ##############################################

////////////////////////////////////////////////////////////////
//
//

void  G4ElasticHadrNucleusHE::DefineHadronValues(G4int Z)
{
  HadrEnergy = std::sqrt(hMass2 + hLabMomentum2);

  G4double sHadr = 2.*HadrEnergy*protonM+protonM2+hMass2;
  G4double sqrS  = std::sqrt(sHadr);
  G4double Ecm   = 0.5*(sHadr-hMass2+protonM2)/sqrS;
  MomentumCM     = std::sqrt(Ecm*Ecm-protonM2);
  
  if(verboseLevel>2)
    G4cout << "GetHadrVall.: Z= " << Z << " iHadr= " << iHadron 
	   << " E(GeV)= " << HadrEnergy << " sqrS= " << sqrS
	   << " plab= " << hLabMomentum   
	   <<"  E - m  "<<HadrEnergy - hMass<< G4endl;

  G4double TotN = 0.0;
  G4double logE = G4Log(HadrEnergy);
  G4double logS = G4Log(sHadr);
           TotP = 0.0;

  switch (iHadron)
    {
    case 0:                  //  proton, neutron
    case 6:

      if(hLabMomentum > 10)
	TotP = TotN = 7.5*logE - 40.12525 + 103*G4Exp(-G4Log(sHadr)*0.165);// mb

      else
	{
// ==================  neutron  ================

////	  if(iHadrCode == 2112) 


	  if( hLabMomentum > 1.4 )
	    TotN = 33.3+15.2*(hLabMomentum2-1.35)/
	      (G4Exp(G4Log(hLabMomentum)*2.37)+0.95);
		
	  else if(hLabMomentum > 0.8)
	    {
	      G4double A0 = logE + 0.0513;
	      TotN = 33.0 + 25.5*A0*A0;  
	    }
	  else 
	    {
	      G4double A0 = logE - 0.2634;  // log(1.3)
	      TotN = 33.0 + 30.*A0*A0*A0*A0;
	    }
//  =================  proton  ===============
//       else if(iHadrCode == 2212) 
	  {
	    if(hLabMomentum >= 1.05)
              {
		TotP = 39.0+75.*(hLabMomentum-1.2)/
		  (hLabMomentum2*hLabMomentum+0.15);
              }

	    else if(hLabMomentum >= 0.7)
	      {
		 G4double A0 = logE + 0.3147;
		 TotP = 23.0 + 40.*A0*A0;
	      }
	    else 
              {
		TotP = 23.+50.*G4Exp(G4Log(G4Log(0.73/hLabMomentum))*3.5);
	      }
	  }
	}

//      HadrTot = 0.5*(82*TotP+126*TotN)/104;  //  $$$$$$$$$$$$$$$$$$
      HadrTot = 0.5*(TotP+TotN);
//  ...................................................
      //  Proton slope
      if(hLabMomentum >= 2.)       HadrSlope = 5.44 + 0.88*logS;

      else if(hLabMomentum >= 0.5) HadrSlope = 3.73*hLabMomentum-0.37;

      else                         HadrSlope = 1.5;

//  ...................................................
      if(hLabMomentum >= 1.2)
	HadrReIm  = 0.13*(logS - 5.8579332)*G4Exp(-G4Log(sHadr)*0.18);
       
      else if(hLabMomentum >= 0.6) 
	HadrReIm = -75.5*(G4Exp(G4Log(hLabMomentum)*0.25)-0.95)/
	  (G4Exp(G4Log(3*hLabMomentum)*2.2)+1);     

      else 
	HadrReIm = 15.5*hLabMomentum/(27*hLabMomentum2*hLabMomentum+2);
//  ...................................................
      DDSect2   = 2.2;                              //mb*GeV-2
      DDSect3   = 0.6;                               //mb*GeV-2
      //  ================== lambda  ==================
      if( iHadrCode == 3122)
	{
	  HadrTot   *= 0.88;
	  HadrSlope *=0.85;
	}
      //  ================== sigma +  ==================
      else if( iHadrCode == 3222)
	{
	  HadrTot   *=0.81;
	  HadrSlope *=0.85;
	}
      //  ================== sigma 0,-  ==================
      else if(iHadrCode == 3112 || iHadrCode == 3212 )
	{
	  HadrTot   *=0.88;
	  HadrSlope *=0.85;
	}
      //  ===================  xi  =================
      else if( iHadrCode == 3312 || iHadrCode == 3322 )
	{
	  HadrTot   *=0.77;
	  HadrSlope *=0.75;
	}
      //  =================  omega  =================
      else if( iHadrCode == 3334)
	{
	  HadrTot   *=0.78;
	  HadrSlope *=0.7;
	}

      break;
//  ===========================================================
    case 1:              //   antiproton
    case 7:              //   antineutron

      HadrTot   = 5.2+5.2*logE + 123.2/sqrS;     //  mb
      HadrSlope = 8.32+0.57*logS;                //(GeV/c)^-2

      if( HadrEnergy < 1000 )
	HadrReIm  = 0.06*(sqrS-2.236)*(sqrS-14.14)*G4Exp(-G4Log(sHadr)*0.8);
      else
	HadrReIm  = 0.6*(logS - 5.8579332)*G4Exp(-G4Log(sHadr)*0.25);

      DDSect2   = 11;                            //mb*(GeV/c)^-2
      DDSect3   = 3;                             //mb*(GeV/c)^-2
      //  ================== lambda  ==================
      if( iHadrCode == -3122)
	{
	  HadrTot   *= 0.88;
	  HadrSlope *=0.85;
	}
      //  ================== sigma +  ==================
      else if( iHadrCode == -3222)
	{
	  HadrTot   *=0.81;
	  HadrSlope *=0.85;
	}
      //  ================== sigma 0,-  ==================
      else if(iHadrCode == -3112 || iHadrCode == -3212 )
	{
	  HadrTot   *=0.88;
	  HadrSlope *=0.85;
	}
      //  ===================  xi  =================
      else if( iHadrCode == -3312 || iHadrCode == -3322 )
	{
	  HadrTot   *=0.77;
	  HadrSlope *=0.75;
	}
      //  =================  omega  =================
      else if( iHadrCode == -3334)
	{
	  HadrTot   *=0.78;
          HadrSlope *=0.7;
	}

      break;
//  -------------------------------------------
    case 2:             //   pi plus, pi minus
    case 3:

      if(hLabMomentum >= 3.5)
	TotP = 10.6+2.*logE + 25.*G4Exp(-G4Log(HadrEnergy)*0.43); // mb
//  =========================================
      else if(hLabMomentum >= 1.15)
	{
          G4double x = (hLabMomentum - 2.55)/0.55; 
	  G4double y = (hLabMomentum - 1.47)/0.225;
	  TotP = 3.2*G4Exp(-x*x) + 12.*G4Exp(-y*y) + 27.5;
	}
//  =========================================
      else if(hLabMomentum >= 0.4)
	{
	TotP  = 88*(logE+0.2877)*(logE+0.2877)+14.0;
        }
//  =========================================
      else 
	{
	  G4double x = (hLabMomentum - 0.29)/0.085;
	  TotP = 20. + 180.*G4Exp(-x*x);
	}
//  -------------------------------------------

      if(hLabMomentum >= 3.0 )
	TotN = 10.6 + 2.*logE + 30.*G4Exp(-G4Log(HadrEnergy)*0.43); // mb

      else if(hLabMomentum >= 1.3) 
	{
          G4double x = (hLabMomentum - 2.1)/0.4;
          G4double y = (hLabMomentum - 1.4)/0.12;
	  TotN = 36.1+0.079 - 4.313*logE + 3.*G4Exp(-x*x) + 
                                              1.5*G4Exp(-y*y);
	}
      else if(hLabMomentum >= 0.65)
	{
          G4double x = (hLabMomentum - 0.72)/0.06;
          G4double y = (hLabMomentum - 1.015)/0.075;
	  TotN = 36.1 + 10.*G4Exp(-x*x) + 24*G4Exp(-y*y);
	}
      else if(hLabMomentum >= 0.37)
	{
	  G4double x = G4Log(hLabMomentum/0.48);
	  TotN = 26. + 110.*x*x;
	}
      else 
	{
          G4double x = (hLabMomentum - 0.29)/0.07;
	  TotN = 28.0 + 40.*G4Exp(-x*x);
	}
      HadrTot = (TotP+TotN)/2;
//  ........................................
      HadrSlope = 7.28+0.245*logS;        // GeV-2
      HadrReIm  = 0.2*(logS - 4.6051702)*G4Exp(-G4Log(sHadr)*0.15);

      DDSect2   = 0.7;                               //mb*GeV-2
      DDSect3   = 0.27;                              //mb*GeV-2

      break;
//  ==========================================================
    case 4:            //  K plus

      HadrTot   = 10.6+1.8*logE + 9.0*G4Exp(-G4Log(HadrEnergy)*0.55);  // mb
      if(HadrEnergy>100) HadrSlope = 15.0;
      else HadrSlope = 1.0+1.76*logS - 2.84/sqrS;   // GeV-2

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
  if(verboseLevel>2)
    G4cout << "HadrTot= " << HadrTot << " HadrSlope= " << HadrSlope
	   << " HadrReIm= " << HadrReIm << " DDSect2= " << DDSect2
	   << " DDSect3= " << DDSect3 << G4endl;

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

  switch(iHadron)
    {
    case 0 :

      if(hLabMomentum <BoundaryP[0])
	InterpolateHN(6,EnP0,C0P0,C1P0,B0P0,B1P0);

      Coeff2 = 0.8/hLabMomentum/hLabMomentum;
      break; 

    case  6 :

      if(hLabMomentum < BoundaryP[1])
	InterpolateHN(5,EnN,C0N,C1N,B0N,B1N);

      Coeff2 = 0.8/hLabMomentum/hLabMomentum;
      break; 

    case 1 :
    case 7 :
      if(hLabMomentum <  BoundaryP[2])
	InterpolateHN(2,EnP,C0P,C1P,B0P,B1P);
      break; 

    case 2 :

      if(hLabMomentum < BoundaryP[3])
	InterpolateHN(4,EnPP,C0PP,C1PP,B0PP,B1PP);

      Coeff2 = 0.02/hLabMomentum;
      break; 

    case 3 :

      if(hLabMomentum < BoundaryP[4])
	InterpolateHN(4,EnPPN,C0PPN,C1PPN,B0PPN,B1PPN);

      Coeff2 = 0.02/hLabMomentum;
      break;
 
    case 4 :

      if(hLabMomentum < BoundaryP[5])
	InterpolateHN(4,EnK,C0K,C1K,B0K,B1K);

      if(hLabMomentum < 1) Coeff2 = 0.34;
      else  Coeff2 = 0.34/hLabMomentum2/hLabMomentum;
      break; 

    case 5 :
      if(hLabMomentum < BoundaryP[6])
	InterpolateHN(2,EnKM,C0KM,C1KM,B0KM,B1KM);

      if(hLabMomentum < 1) Coeff2 = 0.01;
      else  Coeff2 = 0.01/hLabMomentum2/hLabMomentum;
      break; 
    }

  if(verboseLevel > 2) 
    G4cout<<"  HadrVal : Plasb  "<<hLabMomentum
	  <<"  iHadron  "<<iHadron<<"  HadrTot  "<<HadrTot<<G4endl;
}

//  =====================================================
void  G4ElasticHadrNucleusHE::
       GetKinematics(const G4ParticleDefinition * aHadron,
                           G4double MomentumH)
{
  if (verboseLevel>1)
    G4cout<<"1  GetKin.: HadronName MomentumH "
	  <<aHadron->GetParticleName()<<"  "<<MomentumH<<G4endl;

  DefineHadronValues(1);

  G4double Sh     = 2.0*protonM*HadrEnergy+protonM2+hMass2; // GeV

  ConstU = 2*protonM2+2*hMass2-Sh;

  G4double MaxT = 4*MomentumCM*MomentumCM;

  BoundaryTL[0] = MaxT; //2.0;
  BoundaryTL[1] = MaxT;
  BoundaryTL[3] = MaxT;
  BoundaryTL[4] = MaxT;
  BoundaryTL[5] = MaxT;

  G4int NumberH=0;

  while(iHadrCode!=HadronCode[NumberH]) NumberH++;  /* Loop checking, 10.08.2015, A.Ribon */

  NumberH = HadronType1[NumberH];   

  if(MomentumH<BoundaryP[NumberH]) MaxTR = BoundaryTL[NumberH];
  else MaxTR = BoundaryTG[NumberH];

  if (verboseLevel>1)
    G4cout<<"3  GetKin. : NumberH  "<<NumberH
	  <<"  Bound.P[NumberH] "<<BoundaryP[NumberH]
	  <<"  Bound.TL[NumberH] "<<BoundaryTL[NumberH]
	  <<"  Bound.TG[NumberH] "<<BoundaryTG[NumberH]
	  <<"  MaxT MaxTR "<<MaxT<<"  "<<MaxTR<<G4endl;

//     GetParametersHP(aHadron, MomentumH);
}
//  ============================================================
G4double G4ElasticHadrNucleusHE::GetFt(G4double Q2)
{
  G4double Fdistr=0;
  G4double SqrQ2 = std::sqrt(Q2);
 
  Fdistr = (1-Coeff1-Coeff0) //-0.0*Coeff2*G4Exp(ConstU))
    /HadrSlope*(1-G4Exp(-HadrSlope*Q2))
    + Coeff0*(1-G4Exp(-Slope0*Q2))
    + Coeff2/Slope2*G4Exp(Slope2*ConstU)*(G4Exp(Slope2*Q2)-1)
    + 2*Coeff1/Slope1*(1/Slope1-(1/Slope1+SqrQ2)*G4Exp(-Slope1*SqrQ2));

  if (verboseLevel>1)
    G4cout<<"Old:  Coeff0 Coeff1 Coeff2 "<<Coeff0<<"  "
          <<Coeff1<<"  "<<Coeff2<<"  Slope Slope0 Slope1 Slope2 "
          <<HadrSlope<<"  "<<Slope0<<"  "<<Slope1<<"  "<<Slope2
          <<"  Fdistr "<<Fdistr<<G4endl;
  return Fdistr;
}
//  +++++++++++++++++++++++++++++++++++++++
G4double G4ElasticHadrNucleusHE::GetQ2(G4double Ran)
{
  G4double DDD0=MaxTR*0.5, DDD1=0.0, DDD2=MaxTR, delta;
  G4double Q2=0;

  FmaxT = GetFt(MaxTR);
  delta = GetDistrFun(DDD0)-Ran;

  const G4int maxNumberOfLoops = 10000;
  G4int loopCounter = -1;
  while ( (std::fabs(delta) > 0.0001) && 
          ++loopCounter < maxNumberOfLoops )  /* Loop checking, 10.08.2015, A.Ribon */
    {
      if(delta>0) 
      {
        DDD2 = DDD0;
        DDD0 = (DDD0+DDD1)*0.5;
      }
      else if(delta<0)
      {
        DDD1 = DDD0; 
        DDD0 = (DDD0+DDD2)*0.5;
      }
      delta = GetDistrFun(DDD0)-Ran;
    }
  if ( loopCounter >= maxNumberOfLoops ) {
    return 0.0;
  }

  Q2 = DDD0;

  return Q2;
}
//  ++++++++++++++++++++++++++++++++++++++++++
G4double G4ElasticHadrNucleusHE::
       HadronProtonQ2(const G4ParticleDefinition * p,
                            G4double inLabMom)
{

  hMass         = p->GetPDGMass()/GeV;
  hMass2        = hMass*hMass;
  hLabMomentum  = inLabMom;
  hLabMomentum2 = hLabMomentum*hLabMomentum;
  HadrEnergy    = sqrt(hLabMomentum2+hMass2); 

  G4double Rand = G4UniformRand();

  GetKinematics(p, inLabMom);

  G4double Q2 = GetQ2(Rand);

  return Q2;
}

//  ===========================================

///////////////////////////////////////////////////////////////////
//
//  

void G4ElasticHadrNucleusHE::Binom()
{
  G4int N, M;
  G4double  Fact1, J;

  for(N = 0; N < 240; N++)
  {
    J = 1;

      for( M = 0; M <= N; M++ )
      {
	  Fact1 = 1;

	  if ( ( N > 0 ) && ( N > M ) && M > 0 )
          {
              J *= G4double(N-M+1)/G4double(M);
              Fact1 = J;
          }
	  SetBinom[N][M] = Fact1;
      }
  }
  return;
}


//
//
///////////////////////////////////////////////////////////

