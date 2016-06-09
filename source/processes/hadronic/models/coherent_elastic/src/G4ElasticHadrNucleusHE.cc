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
// $Id: G4ElasticHadrNucleusHE.cc,v 1.73 2007/06/14 10:00:45 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-00 $
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
//

#include  "G4ElasticHadrNucleusHE.hh"
#include  "Randomize.hh"
#include  "G4ios.hh"
#include  "G4ParticleTable.hh"
#include  "G4IonTable.hh"
#include  "G4Proton.hh"
#include  "G4NistManager.hh"

using namespace std;

///////////////////////////////////////////////////////////////
//
//

G4ElasticData::G4ElasticData(const G4ParticleDefinition* p, 
			     G4int Z, G4double AWeight, G4double* eGeV)
{ 
  hadr     = p;
  massGeV  = p->GetPDGMass()/GeV;
  mass2GeV2= massGeV*massGeV;
  AtomicWeight = G4int(AWeight + 0.5);
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
      R1    = 16.5*1.1;
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
    case 16:
      R1    = 10.50;
      R2    = 5.5;
      Pnucl = 0.7;
      Aeff  = 0.98;
      break;
    case 9:
      R1    = 9.0;
      R2    = 7.0;
      Pnucl = 0.190;
      Aeff  = 0.9;
      break;
    case 4:
      R1    = 6.0;   
      R2    = 3.7;
      Pnucl = 0.4;
      Aeff  = 0.87;
      break;
    case 1:
      R1    = 4.5;   
      R2    = 2.3;
      Pnucl = 0.177;
      Aeff  = 0.9;
      break;
    default:
      R1    = 4.45*std::pow(A - 1.,0.309)*0.9;
      R2    = 2.3 *std::pow(A, 0.36);

      if(A < 100 && A > 3) Pnucl = 0.176 + 0.00275*A;
      else                 Pnucl = 0.4;

      if(A >= 100)               Aeff = 0.7;
      else if(A < 100 && A > 75) Aeff = 1.5 - 0.008*A;
      else                       Aeff = 0.9;
      break;
    }
}

////////////////////////////////////////////////////////////////////
//
//  The constructor for the generating of events
//

G4ElasticHadrNucleusHE::G4ElasticHadrNucleusHE()
  :G4HadronicInteraction("G4ElasticHadrNucleusHE")
{
  verboseLevel = 0;
  plabLowLimit = 20.0*MeV;
  lowestEnergyLimit = 0.0;

  MbToGeV2  =  2.568;
  sqMbToGeV =  1.602;
  Fm2ToGeV2 =  25.68;
  GeV2      =  GeV*GeV;
  protonM   =  proton_mass_c2/GeV;
  protonM2  =  protonM*protonM;

  //  BoundaryP[0]=19.0;BoundaryTG[0]=5.0;BoundaryTL[0]=0.;
  // BoundaryP[1]=20.0;BoundaryTG[1]=1.5;BoundaryTL[1]=0.;
  // BoundaryP[2]=5.0; BoundaryTG[2]=1.0;BoundaryTL[2]=1.5;
  // BoundaryP[3]=8.0; BoundaryTG[3]=3.0;BoundaryTL[3]=0.;
  // BoundaryP[4]=7.0; BoundaryTG[4]=3.0;BoundaryTL[4]=0.;
  // BoundaryP[5]=5.0; BoundaryTG[5]=2.0;BoundaryTL[5]=0.;
  // BoundaryP[6]=5.0; BoundaryTG[6]=1.5;BoundaryTL[6]=3.0;

  Binom();
  // energy in GeV
  Energy[0] = 0.4;
  Energy[1] = 0.6;
  Energy[2] = 0.8;
  LowEdgeEnergy[0] = 0.0;
  LowEdgeEnergy[1] = 0.5;
  LowEdgeEnergy[2] = 0.7;
  G4double e = 1.0;
  G4double f = std::pow(10.,0.1);
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

  for(G4int j=0; j<NHADRONS; j++) 
  {
    HadronCode[j] = ic[j];
    HadronType[j] = id[j];

    for(G4int k = 0; k < 93; k++) SetOfElasticData[j][k] = 0;
  } 
}

///////////////////////////////////////////////////////////////////
//
//

G4ElasticHadrNucleusHE::~G4ElasticHadrNucleusHE()
{
  for(G4int j = 0; j < NHADRONS; j++) 
  {
    for(G4int k = 0; k < 93; k++) 
    {
      if(SetOfElasticData[j][k]) delete SetOfElasticData[j][k];
    }
  } 
}

////////////////////////////////////////////////////////////////////
//
//

G4HadFinalState * G4ElasticHadrNucleusHE::ApplyYourself(
                          const  G4HadProjectile  &aTrack,
                                 G4Nucleus        &targetNucleus)
{
  theParticleChange.Clear();

  const G4HadProjectile* aParticle = &aTrack;
  G4double ekin = aParticle->GetKineticEnergy();

  if( ekin <= lowestEnergyLimit ) 
  {
    theParticleChange.SetEnergyChange(ekin);
    theParticleChange.SetMomentumChange(aTrack.Get4Momentum().vect().unit());
    return &theParticleChange;
  }

  G4double aTarget = targetNucleus.GetN();
  G4double zTarget = targetNucleus.GetZ();

  G4double plab = aParticle->GetTotalMomentum();

  if (verboseLevel >1)
  { 
    G4cout << "G4ElasticHadrNucleusHE: Incident particle plab=" 
	   << plab/GeV << " GeV/c " 
	   << " ekin(MeV) = " << ekin/MeV << "  " 
	   << aParticle->GetDefinition()->GetParticleName() << G4endl;
  }
  // Scattered particle referred to axis of incident particle

  const G4ParticleDefinition* theParticle = aParticle->GetDefinition();
  G4double m1 = theParticle->GetPDGMass();

  G4int Z = static_cast<G4int>(zTarget+0.5);
  G4int A = static_cast<G4int>(aTarget+0.5);
  G4int projPDG = theParticle->GetPDGEncoding();

  if (verboseLevel>1) 
  {
    G4cout << "G4ElasticHadrNucleusHE for " << theParticle->GetParticleName()
	   << " PDGcode= " << projPDG << " on nucleus Z= " << Z 
	   << " A= " << A 
	   << G4endl;
  }
  G4ParticleDefinition * theDef = 0;

  if      (Z == 1 && A == 1) theDef = G4Proton::Proton();
  else if (Z == 1 && A == 2) theDef = G4Deuteron::Deuteron();
  else if (Z == 1 && A == 3) theDef = G4Triton::Triton();
  else if (Z == 2 && A == 3) theDef = G4He3::He3();
  else if (Z == 2 && A == 4) theDef = G4Alpha::Alpha();
  else                       theDef = G4ParticleTable::GetParticleTable()->FindIon(Z,A,0,Z);
 
  G4double m2 = theDef->GetPDGMass();
  G4LorentzVector lv1 = aParticle->Get4Momentum();
  G4LorentzVector lv(0.0,0.0,0.0,m2);   
  lv += lv1;

  G4ThreeVector bst = lv.boostVector();
  lv1.boost(-bst);

  G4ThreeVector p1 = lv1.vect();
  G4double ptot = p1.mag();
  G4double tmax = 4.0*ptot*ptot;
  G4double t = 0.0;

  // Choose generator
  G4bool swave = false;

  // S-wave for very low energy
  if(plab < plabLowLimit) swave = true;

  // normal sampling
  if(!swave) {
    t = SampleT(theParticle,plab,Z,A);
    if(t > tmax) swave = true;
  }

  if(swave) t = G4UniformRand()*tmax;

  // Sampling in CM system
  G4double phi  = G4UniformRand()*twopi;
  G4double cost = 1. - 2.0*t/tmax;
  G4double sint;

  if( cost >= 1.0 ) 
  {
    cost = 1.0;
    sint = 0.0;
  }
  else if( cost <= -1.0) 
  {
    cost = -1.0;
    sint =  0.0;
  }
  else  
  {
    sint = std::sqrt((1.0-cost)*(1.0+cost));
  }  
  if (verboseLevel>1)
  { 
    G4cout << "cos(t)=" << cost << " std::sin(t)=" << sint << G4endl;
  }
  G4ThreeVector v1(sint*std::cos(phi),sint*std::sin(phi),cost);
  v1 *= ptot;
  G4LorentzVector nlv1( v1.x(), v1.y(), v1.z(), std::sqrt(ptot*ptot + m1*m1));

  nlv1.boost(bst); 

  G4double eFinal = nlv1.e() - m1;

  if (verboseLevel > 1)
  { 
    G4cout << "Scattered: "
	   << nlv1<<" m= " << m1 << " ekin(MeV)= " << eFinal 
	   << " Proj: 4-mom " << lv1 
	   <<G4endl;
  }
  if( eFinal < 0.0 ) 
  {
    G4cout << "G4ElasticHadrNucleusHE WARNING ekin= " << eFinal
	   << " after scattering of " 
	   << aParticle->GetDefinition()->GetParticleName()
	   << " p(GeV/c)= " << plab
	   << " on " << theDef->GetParticleName()
	   << G4endl;
    eFinal = 0.0;
    nlv1.setE(m1);
  }

  theParticleChange.SetMomentumChange( nlv1.vect().unit() );
  theParticleChange.SetEnergyChange(eFinal);
  
  G4LorentzVector nlv0 = lv - nlv1;
  G4double erec =  nlv0.e() - m2;

  if (verboseLevel > 1)
  { 
    G4cout << "Recoil: "
	   << nlv0<<" m= " << m2 << " ekin(MeV)= " << erec 
	   <<G4endl;
  }
  if(erec < 0.0) 
  {
    G4cout << "G4ElasticHadrNucleusHE WARNING Erecoil(MeV)= " << erec
	   << " after scattering of " 
	   << aParticle->GetDefinition()->GetParticleName()
	   << " p(GeV/c)= " << plab
	   << " on " << theDef->GetParticleName()
	   << G4endl;
    nlv0.setE(m2);
  }
  G4DynamicParticle * aSec = new G4DynamicParticle(theDef, nlv0);
  theParticleChange.AddSecondary(aSec);

  return &theParticleChange;
}

//////////////////////////////////////////////////////////////////////////
//
//

G4double G4ElasticHadrNucleusHE::SampleT( const G4ParticleDefinition* p,
                                                 G4double inLabMom, G4int Z, G4int N)
{
  G4double plab  = inLabMom/GeV;   // (GeV/c)
  G4double Q2 = 0;

  iHadrCode = p->GetPDGEncoding();

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
  G4double Q2max = 4.*plab2*M2/(hMass2 + M2 + 2.*M*std::sqrt(plab2 + hMass2));

  // sample scattering

  Q2 = HadronNucleusQ2_2(ElD1, Z, plab, Q2max);

  if(verboseLevel > 1)
    G4cout<<" SampleT: Q2(GeV^2)= "<<Q2<< "  t/tmax= " << Q2/Q2max <<G4endl;
  
  return  Q2*GeV2;
}

//////////////////////////////////////////////////////////////////////////
//
//

G4double G4ElasticHadrNucleusHE::HadronNucleusQ2_2(G4ElasticData* pElD, G4int Z, 
						   G4double plab, G4double tmax)
{

  G4double Rand = G4UniformRand();

  G4int      iNumbQ2 = 0;
  G4double   Q2 = 0.0, Buf = 0.0;

  G4double ptot2 = plab*plab;
  G4double ekin  = std::sqrt(hMass2 + ptot2) - hMass;

  // Find closest energy bin
  G4int NumbOnE; 

  for( NumbOnE = 0; NumbOnE < NENERGY-1; NumbOnE++ ) 
  {
    if( ekin <= LowEdgeEnergy[NumbOnE+1] ) break;
  }
  G4double* dNumbQ2 = pElD->TableQ2;

  G4int index = NumbOnE*ONQ2;

  G4double Weight= 1.0;
  G4double rmax  = 1.0;

  // Select kinematics for node energy
  G4double T     = Energy[NumbOnE];
  hLabMomentum2  = T*(T + 2.*hMass);
  G4double Q2max = pElD->maxQ2[NumbOnE];
  G4int length   = pElD->dnkE[NumbOnE];
  G4bool isIni   = false;

  // Build first part of the vector

  if(length == 0) 
  {
    isIni = true;
    R1    = pElD->R1;
    R2    = pElD->R2;
    Aeff  = pElD->Aeff;
    Pnucl = pElD->Pnucl;
    hLabMomentum = std::sqrt(hLabMomentum2);
    DefineHadronValues(Z);
    G4int AWeight = pElD->AtomicWeight;
    Weight = GetLightFq2(Z, AWeight, Q2max);
    pElD->CrossSecMaxQ2[NumbOnE] = Weight;

    if(verboseLevel > 1)
      G4cout<<" HadrNucleusQ2_2: NumbOnE= " << NumbOnE 
	    << " length= " << length 
	    << " Weight "<<Weight
	    << " Q2max= " << Q2max 
	    << " ekin= " << ekin <<G4endl;
    
    pElD->TableCrossSec[index] = 0;

    for(G4int ii=1; ii<ONQ0; ii++)
    {
	Q2 = pElD->TableQ2[ii];

	if(Q2 < Q2max) Buf = GetLightFq2(Z, AWeight, Q2)/Weight;
	else           Buf = 1.0;

	pElD->TableCrossSec[index+ii] = Buf;

	if(verboseLevel > 1)
	  G4cout<<" HadrNucleusQ2_2: ii= " << ii << " Q2= "
		<<Q2 <<" p= " <<Buf<<" B*W "<<Buf*Weight<<G4endl;
    }   // for ii

    rmax    = Buf;
    length  = ONQ0;
    pElD->dnkE[NumbOnE] = ONQ0;

  } 
  else 
  {
    rmax = pElD->TableCrossSec[index+length-1];
  }

  G4double* dNumbFQ2 = &(pElD->TableCrossSec[index]);

  // No more vector needed

  if(rmax >= Rand) 
  {

    for( iNumbQ2 = 1; iNumbQ2<length; iNumbQ2++ ) 
    {
      if(Rand <= pElD->TableCrossSec[index+iNumbQ2]) break;
    }

  
  } 
  else // Build second part of the vector
  {    
    if(!isIni) 
    {
      R1    = pElD->R1;
      R2    = pElD->R2;
      Aeff  = pElD->Aeff;
      Pnucl = pElD->Pnucl;
      hLabMomentum = std::sqrt(hLabMomentum2);
      DefineHadronValues(Z);
      Weight = pElD->CrossSecMaxQ2[NumbOnE];
    }
    G4int AWeight = pElD->AtomicWeight;
     
    // Stop building when find out the node

    for(iNumbQ2 = length; iNumbQ2<ONQ2; iNumbQ2++) 
    {

      Q2 = pElD->TableQ2[iNumbQ2];

      if(Q2 < Q2max) Buf = GetLightFq2(Z, AWeight, Q2)/Weight;
      else           Buf = 1.0;

      pElD->TableCrossSec[index+iNumbQ2] = Buf;

      if(verboseLevel > 1)
	G4cout<<" HadrNucleusQ2_2: NumbOnE= " << NumbOnE 
	      << " iNumbQ2= " << iNumbQ2 << " Q2= "
	      <<Q2 <<" Buf= " <<Buf<<" B*W "<<Buf*Weight<<G4endl;
      
      if(Rand <= Buf) 
      {
	pElD->dnkE[NumbOnE] = iNumbQ2+1;
	break;
      }
    }   
  }
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
  G4double F2  = F[kk-1];
  G4double F3  = F[kk];
  G4double X2  = Q[kk-1];
  G4double X3  = Q[kk];

  if(verboseLevel > 2) 
    G4cout << "GetQ2_2 kk= " << kk << " X2= " << X2 << " X3= " << X3 
	   << " F2= " << F2 << " F3= " << F3 << " Rndm= " << ranUni << G4endl;

  if(kk <= 2)
  {
      ranQ2 = X2 + (ranUni - F2)*(X3 - X2)/(F3 - F2);
      return ranQ2;
  }
  G4double F1  = F[kk-2];

  G4double F12 = F1*F1;
  G4double F22 = F2*F2;
  G4double F32 = F3*F3;

  G4double X1  = Q[kk-2];          //  MeV^2

  G4double D0  = F12*F2+F1*F32+F3*F22-F32*F2-F22*F1-F12*F3;

  if(verboseLevel > 2) 
    G4cout << "       X1= " << X1 << " F1= " << F1 << "  D0= " << D0 << G4endl; 

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

G4double G4ElasticHadrNucleusHE::GetLightFq2(G4int Z, G4int Nucleus, G4double Q2)
{
  // Scattering of proton
  if(Z == 1) 
  {
    G4double SqrQ2  = std::sqrt(Q2);
    G4double ConstU = 2.*(hMass2 + protonM2) - Q2;

    G4double y = (1.-Coeff1-Coeff0)/HadrSlope*(1.-std::exp(-HadrSlope*Q2))
      + Coeff0*(1.-std::exp(-Slope0*Q2))
      + Coeff2/Slope2*std::exp(Slope2*ConstU)*(std::exp(Slope2*Q2)-1.)
      + 2.*Coeff1/Slope1*(1./Slope1-(1./Slope1+SqrQ2)*std::exp(-Slope1*SqrQ2));

    return y;
  }

  // The preparing of probability function  

  G4double prec = Nucleus > 208  ?  1.0e-7 : 1.0e-6;

  G4double    Stot     = HadrTot*MbToGeV2;     //  Gev^-2
  G4double    Bhad     = HadrSlope;         //  GeV^-2
  G4double    Asq      = 1+HadrReIm*HadrReIm;
  G4double    Rho2     = std::sqrt(Asq);

  if(verboseLevel > 1) {
    G4cout << "Stot= " << Stot << " Bhad= " << Bhad << " Asq= " << Asq << G4endl;
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
  G4double    Tot0     = 0;
  G4double    exp1;

  G4double    Prod3 ;
  G4double    exp2  ;
  G4double    N4, N5, N2, Prod1, Prod2;
  G4int    i1, i2, m1, m2;

  for(i1 = 1; i1<= Nucleus; i1++) ////++++++++++  i1
    {
      N1    = -N1*Unucl*(Nucleus-i1+1)/i1*Rho2;
      Prod1 = 0;
      Tot0  = 0;
      N2    = -1;

      for(i2 = 1; i2<=Nucleus; i2++) ////+++++++++ i2
        {
          N2    = -N2*Unucl*(Nucleus-i2+1)/i2*Rho2;
          Prod2 = 0; 
          N5    = -1/NN2;
	  for(m2=0; m2<= i2; m2++) ////+++++++++ m2
            {
              Prod3 = 0;
              exp2  = 1/(m2/R22B+(i2-m2)/R12B);
              N5    = -N5*NN2;
              N4    = -1/NN2;
	      for(m1=0; m1<=i1; m1++) ////++++++++ m1
		{
		  exp1  = 1/(m1/R22B+(i1-m1)/R12B);
		  dddd  = exp1+exp2;
		  N4    = -N4*NN2;
		  Prod3 = Prod3+N4*exp1*exp2*
		    (1-std::exp(-Q2*dddd/4))/dddd*4*SetBinom[i1][m1];
               }                                   // m1
	      Prod2 = Prod2 +Prod3*N5*SetBinom[i2][m2];
	    }                                      // m2
	  Prod1 = Prod1 + Prod2*N2*std::cos(FiH*(i1-i2));

	  if (std::fabs(Prod2*N2/Prod1)<prec) break;
        }                                         // i2
      Prod0   = Prod0 + Prod1*N1;
      if(std::fabs(N1*Prod1/Prod0) < prec) break;
    }                                           // i1
      /*
  for(G4int i1 = 1; i1<= Nucleus; i1++) 
  {
    N1 *= UnucRho2*G4double(Nucleus-i1+1)/G4double(i1);
    Prod1 = 0;
    Tot0  = 0;
    N2    = -1;

    for(G4int i2 = 1; i2<=Nucleus; i2++) 
    {
      N2 *= UnucRho2*G4double(Nucleus-i2+1)/G4double(i2);
      Prod2 = 0; 
      N5    = -1.0/NN2;

      for(G4int m2=0; m2<= i2; m2++) 
      {
	Prod3 = 0;
	exp2  = 1.0/(m2/R22B+(i2-m2)/R12B);
	N5   *= (-NN2);
	N4    = -1.0/NN2;

	for(G4int m1=0; m1<=i1; m1++) 
        {
	  exp1   = 1.0/(m1/R22B+(i1-m1)/R12B);
	  dddd   = 0.25*(exp1 + exp2);
	  N4    *= (-NN2);
	  Prod3 += N4*exp1*exp2*SetBinom[i1][m1]*(1-std::exp(-Q2*dddd))/dddd;
	}                                   // m1
	Prod2 += Prod3*N5*SetBinom[i2][m2];
      }                                      // m2
      Prod1 += Prod2*N2*std::cos(FiH*(i1-i2));

      if (std::abs(Prod2*N2/Prod1)<prec) break;
    }                                         // i2
    Prod0  += Prod1*N1;

    if(std::abs(N1*Prod1/Prod0) < prec) break;
  }     
      */                                      // i1
  Prod0 *= 0.25*pi/MbToGeV2;  //  This is in mb
  if(verboseLevel>1) 
    G4cout << "GetLightFq2 Z= " << Z << " A= " << Nucleus 
	   <<" Q2= " << Q2 << " Res= " << Prod0 << G4endl;
  return Prod0;
}

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
    G4cout << "Z= " << Z << " iHadr= " << iHadron 
	   << " E(GeV)= " << HadrEnergy << " sqrS= " << sqrS
	   << " plab= " << hLabMomentum   
	   << G4endl;
  
  if( HadrEnergy - hMass <0.39 )
  {
      G4cout<<"ElasticHE(GetHadronValues): The energy T = "
	    <<(HadrEnergy-hMass)
	    <<" GeV is very low for this method! T = 0.400 GeV is set."
	    <<G4endl;
  }
  G4double TotP = 0.0, TotN = 0.0;
  G4double logE = std::log(HadrEnergy);
  G4double logS = std::log(sHadr);

  switch (iHadron)
    {
    case 0:                  //  proton, neutron
    case 6:

      TotP = TotN = 7.5*logE - 40.12525 + 103*std::pow(sHadr,-0.165); //  mb

      if(hLabMomentum<10.)
	{
	  // ==================  neutron  ================

	  if(iHadrCode == 2112) 
	    {    
	      if( hLabMomentum > 1.4 )
		TotN = 33.3+15.2*(hLabMomentum2-1.35)/(std::pow(hLabMomentum,2.37)+0.95);
		
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
	    }
	  else if(iHadrCode == 2212) 
           {
	     if(hLabMomentum>1.05)
	       TotP = 39.0+75.*(hLabMomentum-1.2)/(hLabMomentum2*hLabMomentum+0.15);
	    
	     else if(hLabMomentum > 0.7)
	       {
		 G4double A0 = logE + 0.3147;
		 TotP = 23.0 + 40.*A0*A0;
	       }
	     else 
	       TotP = 23.+50.*std::pow(std::log(0.73/hLabMomentum),3.5);
	   }
	 }
      HadrTot = 0.5*(TotP+TotN);
	
      //  Proton slope
      HadrSlope = 5.44 + 0.88*logS;     //  GeV-2

      if( hLabMomentum < 2.) HadrSlope = 1.5;

      if(hLabMomentum>1.2)
	 HadrReIm  = 0.13*(logS - 5.8579332)*std::pow(sHadr,-0.18);
       
      else if(hLabMomentum > 0.6)
	HadrReIm = -75.5*(std::pow(hLabMomentum,0.25)-0.95)/
	  (std::pow(3*hLabMomentum,2.2)+1);     
                                                                                    
      else 
	HadrReIm = 15.5*hLabMomentum/(27*hLabMomentum2*hLabMomentum+2);
      
      DDSect2   = 11;                              //mb*GeV-2
      DDSect3   = 3;                               //mb*GeV-2
      //  ================== lambda  ==================
      if( iHadrCode == 3122)
	{
	  HadrTot *= 0.88;
	  HadrSlope *=0.85;
	}
      //  ================== sigma +  ==================
      if( iHadrCode == 3222)
	{
	  HadrTot   *=0.81;
	  HadrSlope *=0.85;
	}
      //  ================== sigma 0,-  ==================
      if(iHadrCode == 3112 || iHadrCode == 3212 )
	{
	  HadrTot   *=0.88;
	  HadrSlope *=0.85;
	}
      //  ===================  xi  =================
      if( iHadrCode == 3312 || iHadrCode == 3322 )
	{
	  HadrTot   *=0.77;
	  HadrSlope *=0.75;
	}
      //  =================  omega  =================
      if( iHadrCode == 3334)
	{
	  HadrTot   *=0.78;
	  HadrSlope *=0.7;
	}

      break;
      
    case 1:              
    case 7:              //   antiproton

      HadrTot   = 5.2+5.2*logE + 123.2/sqrS;     //  mb
      HadrSlope = 8.32+0.57*logS; //GeV-2

      if( HadrEnergy < 1000 )
	HadrReIm  = 0.06*(sqrS-2.236)*(sqrS-14.14)*std::pow(sHadr,-0.8);
      else
	HadrReIm  = 0.6*(logS - 5.8579332)*std::pow(sHadr,-0.25);

      DDSect2   = 11;                                //mb*GeV-2
      DDSect3   = 3;                                 //mb*GeV-2
      //  ================== lambda  ==================
      if( iHadrCode == -3122)
	{
	  HadrTot *= 0.88;
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

      if(hLabMomentum>3.5)
	TotP = 10.6+2.*logE + 25.*std::pow(HadrEnergy,-0.43); // mb
      //  =========================================
      else if(hLabMomentum >1.15)
	{
          G4double x = (hLabMomentum - 2.55)/0.55; 
	  G4double y = (hLabMomentum - 1.47)/0.225;
	  TotP = 3.2*std::exp(-x*x) + 12.*std::exp(-y*y) + 27.5;
	}
      //  =========================================
      else if(hLabMomentum >0.4)
	TotP  = 88*(logE+0.2877)*(logE+0.2877)+14.0;
      
      //  =========================================
      else 
	{
	  G4double x = (hLabMomentum - 0.29)/0.085;
	  TotP = 20. + 180.*std::exp(-x*x);
	}
      //  -------------------------------------------

      if(hLabMomentum > 3.0 )
	TotN = 10.6 + 2.*logE + 30.*std::pow(HadrEnergy,-0.43); // mb

      else if(hLabMomentum > 1.3) 
	{
          G4double x = (hLabMomentum - 2.1)/0.4;
          G4double y = (hLabMomentum - 1.4)/0.12;
	  TotN = 36.1+0.079 - 4.313*logE + 3.*std::exp(-x*x) + 1.5*std::exp(-y*y);
	}
      else if(hLabMomentum > 0.65)
	{
          G4double x = (hLabMomentum - 0.72)/0.06;
          G4double y = (hLabMomentum - 1.015)/0.075;
	  TotN = 36.1 + 10.*std::exp(-x*x) + 24*std::exp(-y*y);
	}
      else if(hLabMomentum > 0.37)
	{
	  G4double x = std::log(hLabMomentum/0.48);
	  TotN = 26. + 110.*x*x;
	}
      else 
	{
          G4double x = (hLabMomentum - 0.29)/0.07;
	  TotN = 28.0 + 40.*std::exp(-x*x);
	}
      HadrTot = (TotP+TotN)/2;

      HadrSlope = 7.28+0.245*logS;        // GeV-2
      HadrReIm  = 0.2*(logS - 4.6051702)*std::pow(sHadr,-0.15);

      DDSect2   = 4.6;                               //mb*GeV-2
      DDSect3   = 1.33;                              //mb*GeV-2

      break;
      
    case 4:            //  K plus

      HadrTot   = 10.6+1.8*logE + 9.0*std::pow(HadrEnergy,-0.55);  // mb
      if(HadrEnergy>100) HadrSlope = 15.0;
      else HadrSlope = 1.0+1.76*logS - 2.84/sqrS;   // GeV-2

      HadrReIm  = 0.4*(sHadr-20)*(sHadr-150)*std::pow(sHadr+50,-2.1);
      DDSect2   = 3.5;                             //mb*GeV-2
      DDSect3   = 1.03;                            //mb*GeV-2
      break;
      
     case 5:              //   K minus

       HadrTot   = 10+1.8*logE + 25./sqrS; // mb
       HadrSlope = 6.98+0.127*logS;         // GeV-2
       HadrReIm  = 0.4*(sHadr-20)*(sHadr-20)*std::pow(sHadr+50,-2.1);
       DDSect2   = 3.5;                             //mb*GeV-2
       DDSect3   = 1.03;                            //mb*GeV-2
       break;
  }   
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

      if(hLabMomentum < EnP0[5])
	InterpolateHN(6,EnP0,C0P0,C1P0,B0P0,B1P0);
      break; 

    case  6 :
    case  7 :

      if(hLabMomentum < EnN[4])
	InterpolateHN(5,EnN,C0N,C1N,B0N,B1N);
      break; 

    case 1 :

      if(hLabMomentum < EnP[1])
	InterpolateHN(2,EnP,C0P,C1P,B0P,B1P);
      break; 

    case 2 :

      if(hLabMomentum < EnPP[3])
	InterpolateHN(4,EnPP,C0PP,C1PP,B0PP,B1PP);
      Coeff2 = 0.02/hLabMomentum;
      break; 

    case 3 :

      if(hLabMomentum < EnPPN[3])
	InterpolateHN(4,EnPPN,C0PPN,C1PPN,B0PPN,B1PPN);
      Coeff2 = 0.02/hLabMomentum;
      break;
 
    case 4 :

      if(hLabMomentum < EnK[3])
	InterpolateHN(4,EnK,C0K,C1K,B0K,B1K);
      Coeff2 = 0.34/hLabMomentum2/hLabMomentum;
      break; 

    case 5 :
      if(hLabMomentum < EnKM[1])
	InterpolateHN(2,EnKM,C0KM,C1KM,B0KM,B1KM);
      Coeff2 = 0.01/hLabMomentum2/hLabMomentum;
      break; 
  }
}

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

