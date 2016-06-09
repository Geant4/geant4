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
// $Id: G4ElasticHadrNucleusHE.cc,v 1.34 2006/07/06 18:11:11 vnivanch Exp $
// GEANT4 tag $Name: geant4-08-01-patch-02 $
//G4ElasticHadrNucleusHE.cc
//
//
//  G4ElasticHadrNucleusHE class
//
//  The generator of high energy hadron-nucleus elastic scattering
//  The hadron kinetic energy T > 1 GeV
//  N.  Starkov 2003.

//  Modifications:
//  19.05.04 Variant for G4 6.1: The 'ApplyYourself' was changed
//  14.11.05 The HE elastic scattering on proton is added (N.Starkov)
//  23.11.05 int -> G4int, fabs -> abs (V.Ivanchenko)
//  09.05.06 N.Starkov propose new method SampleT (V.Ivanchenko)
//  30.05.06 Variant without use of elastic data files (N.Starkov)
//  07.06.06 V.Ivanchenko fix problem of rotation of final state
//

#include  "G4ElasticHadrNucleusHE.hh"
#include  "Randomize.hh"
#include  <sstream>
#include  <iostream>
#include  <fstream>
#include  "G4ios.hh"
#include "G4ParticleTable.hh"

using namespace std;
//  ############################################################
ElasticData:: ElasticData(G4String HadrName, G4int AtomWeight)
   { 
//     G4Nucleus *pNucleus, InterNucleus;
     G4int      ii, kk;
     G4double   Expo=1, Kinetic[5] = {1.0, 1.5, 2.5, 4.0, 6.0};

//     InterNucleus.SetParameters(AtomWeight, AtomWeight/2);
//     pNucleus = &InterNucleus;
     hadrName = HadrName;
     nuclAtomicNumber = AtomWeight;


     NucPar.GetNucleusParameters(AtomWeight);

     theR1 = NucPar.R1;
     R2    = NucPar.R2;
     Aeff  = NucPar.Aeff;
     Pnucl = NucPar.Pnucl;
     theMaxQ2 = GetQ2limit(theR1);

     for(kk=0; kk<AreaNumb; kk++)
     {
       for(ii=0; ii<ONE;ii++)
       {
         TableE[ii+ONE*kk] = Kinetic[ii]*Expo;
       }
         Expo      *= 10;
       for(ii=0; ii<ONQ2XE; ii++)
         TableCrossSec[ii+kk*ONQ2XE] = 0;
     }
//     G4DiffElasticHadrNucleus::GetNucleusParameters(pNucleus);
//     R1 = G4ElasticHadrNucleusHE::R1;
   }
//  ###########################################################
   G4double ElasticData::GetQ2limit(G4double aR1)
   {
     G4double maxQ2_0 = 35/*1000*1000*//aR1/aR1;

     Nstep = ONQ2;
     dQ2   = maxQ2_0/(Nstep/*3./2.*/-1.);

     TableQ2[0] = 0.00001;

     for(G4int ii=1; ii<Nstep; ii++)
     {
          TableQ2[ii] = TableQ2[ii-1]+dQ2;
     }
//ii<Nstep/2 ? iQ2[ii-1]+dQ2 : iQ2[ii-1]+2*dQ2;

     return maxQ2_0;
   }
//  ############################################################
   G4int  ElasticData::GetNumberE(G4double E)
   {
     G4int ii;

     for(ii=0; ii<ONE*AreaNumb; ii++)
       if(E > TableE[ii]) break;
     return ii;
   }

//  ####### The constructor for the generating of events #######
   G4ElasticHadrNucleusHE::G4ElasticHadrNucleusHE()
       :  G4DiffElasticHadrNucleus()
   {
      ;
   }
//  ####### The constructor for the preparation of data ########
   G4ElasticHadrNucleusHE::
        G4ElasticHadrNucleusHE(const G4ParticleDefinition * aHadron1,
                                     G4Nucleus            * aNucleus) 
       :  G4DiffElasticHadrNucleus(), G4HadronicInteraction()
   {  
//  --------------  here all variables are in MeV  ----------------
    G4int ii, kk, ik;
           
    iKindWork = 3;
    iPoE      = 5;
    Nstep     = ONQ2;

    G4ParticleDefinition * aHadron = 
           const_cast<G4ParticleDefinition *>(aHadron1); 

    G4String sNameHdr    = aHadron->GetParticleName();
    G4int    iNnucl      = (int) aNucleus->GetN();
    G4double      Energies[5]; // = {1.95, 2.0, 4.0, 6.0, 10.0};
    G4double      Kinetic[5]  = {1.0, 1.5, 2.5, 4.0, 6.0};

    if(iNnucl==1) return;

    if(iNnucl==2 || iNnucl==3)
      G4Exception(" This nucleus is very light for this model !!!");

    if(iNnucl>238)
      G4Exception(" This nucleus is very heavy for this model !!!");


    G4DynamicParticle * innerHadron =  new G4DynamicParticle();
    innerHadron->SetDefinition(aHadron);

    G4ThreeVector  inVector(0.,0.,0.);

    G4double    mHadr   = aHadron->GetPDGMass();     //   MeV
    G4double    curE, curP;         

    for(ii=0; ii<iPoE; ii++) Energies[ii] = Kinetic[ii]+mHadr/1000;
  
    innerHadron->SetDefinition(aHadron);

    G4double    dPower=1000.0; 

    G4cout<<G4endl<<" The preparing of elastic Data array "<<G4endl;
    G4cout<<" for '"<<sNameHdr<<"' on nucleus "<<iNnucl<<G4endl<<G4endl;

    for( ik = 0; ik<AreaNumb; ik++) 
    {
      for( ii = 0; ii<iPoE; ii++)
      {
        curE         = Energies[ii]*dPower;
        curP         = std::sqrt(curE*curE-mHadr*mHadr);  //  MeV
        pTableE[ii+ik*iPoE] = curE;
        inVector.setZ(curP);                //  MeV
        innerHadron->SetMomentum(inVector);
//  ------ The creating of the distribution function array ------
        CreationArray(innerHadron, aNucleus);

        for( kk = 0; kk<ONQ2; kk++)
          pTableCrSec[ii*ONQ2+kk+ik*ONQ2XE] = iIntgr[kk];

	G4cout<<" Energy of "<<sNameHdr<<" =  "<<curE<<G4endl;
      }         //   ii   step Ei
      dPower *= 10.0; 
    }             //   ik   step power

  G4cout<<" The array for elastic scattering of the "<<sNameHdr<<G4endl;
  G4cout<<" on the Nuclues "<<iNnucl<<" has been prepared !"
        <<G4endl<<G4endl;

//  ---------   The writing of array into sName File  ---------

 G4String   sPath;

  if(getenv("G4ELASTICDATA")) 
   {
    sPath =  getenv("G4ELASTICDATA");
    sPath =  sPath+"/Elastic/";
   }
  else  sPath =       "./Elastic/";

  //  G4cout<<" Path : "<<sPath<<G4endl;
  std::ostringstream sNameFile;
  sNameFile<<sPath<<sNameHdr<<"_"<<iNnucl<<".dat"<<std::ends;

  //       sNameFile<<"Elastic//"<<sNameHdr<<"_"<<iNnucl<<".dat"<<ends;
  G4String str = sNameFile.str();
  std::ofstream TestFile(str, std::ios::out);
  TestFile.setf(std::ios::scientific, std::ios::floatfield);
  TestFile.precision(9);

  TestFile <<iPoE<<G4endl;

  G4int Mult;
  dPower = 1000.0;
  for( ik=0; ik<AreaNumb; ik++)
  {
    for(kk=0; kk<iPoE; kk++) 
    { 
      curE = Energies[kk]*dPower;
      Mult = (int) curE;
      TestFile<<"  "<<Mult<<" ";
    }
    TestFile<<G4endl;

    for(ii=0; ii<ONQ2; ii++)
    {
      for(kk=0; kk<iPoE; kk++) 
         TestFile << pTableCrSec[kk*ONQ2+ii+ik*ONQ2XE]<<" ";
      TestFile << G4endl;
    }        //   ii
    TestFile<<G4endl;
    dPower *= 10.0;
   }           //   ik

   TestFile.close();

  if(getenv("HEElastic_debug"))
  {
    G4cout<<" The array for elastic scattering of "<<sNameHdr<<G4endl;
    G4cout<<" on the Nuclues "<<iNnucl<<" has been written !"<<G4endl;
    G4cout<<" The Name of File is "<<sNameFile.str()<<G4endl<<G4endl;
  }

  iContr = 0;
 }   //  Constructor 

//  ######################################################
G4int G4ElasticHadrNucleusHE::         
             ReadOfData(G4ParticleDefinition * aHadron,
                        G4Nucleus            * aNucleus)
{
  G4int    ii, kk, ik;
  G4String sNameHdr  = aHadron->GetParticleName();
  G4int    iNnucl    = (int) aNucleus->GetN();

  if(iNnucl==2 || iNnucl==3)
    G4Exception(" This nucleus is very light for this model !!!");

  if(iNnucl>238)
    G4Exception(" This nucleus is very heavy for this model !!!");

  G4String sPath;

  if(getenv("G4ELASTICDATA")) 
    {
     sPath =  getenv("G4ELASTICDATA");
     sPath =  sPath+"/Elastic/";
    }
  else  sPath =       "./Elastic/";

  std::ostringstream sNameFile;
  sNameFile<<sPath<<sNameHdr<<"_"<<iNnucl<<".dat"<<std::ends;
  
  if(getenv("HEElastic_debug"))
    {
      G4cout <<" Reading file for: Hadron - "<<sNameHdr
	     <<".  Nucleus - "<<iNnucl<<G4endl;
      G4cout <<" The Name of File is "<<sNameFile.str()<<G4endl;
    }

  GetNucleusParameters(aNucleus);

  iContr         = iPoE;
  Nstep          = ONQ2;
  iKindWork      = 4;
  //  theMaxQ2       = GetQ2limit(theR1);   //   MeV^2

  G4String stri = sNameFile.str();
  std::ifstream TestFile(stri, std::ios::in);
  TestFile.setf(std::ios::scientific);

  if(TestFile)
    {
      ElasticData  ElD;
      if(getenv("HEElastic_debug")) G4cout<<" The file exists "<<G4endl;
      
      ElD.hadrName         = sNameHdr;
      ElD.nuclAtomicNumber = iNnucl;
     
      TestFile  >> iPoE;

      G4double  dPower;

      for(ik=0; ik<AreaNumb; ik++)
	{
	  dPower = std::pow(10.0,ik);
	  
	  for(kk=0; kk<iPoE; kk++)  TestFile >> ElD.TableE[kk+ik*iPoE];

	  for(ii=0; ii<ONQ2; ii++)
	    {

	      ElD.TableQ2[ii] = iQ2[ii];
	      for(kk=0; kk<iPoE; kk++) 
		TestFile >> ElD.TableCrossSec[kk*ONQ2+ii+ik*ONQ2XE];
	    }       //  ii
	}               //  ik
      
      TestFile.close();

      SetOfElasticData.push_back(ElD);

      if(getenv("HEElastic_debug"))
	{
	  G4cout<<" The array for elastic scattering of "<<sNameHdr;
	  G4cout<<" on the Nuclues "<<iNnucl<<G4endl;
	  G4cout<<" from the file "<<sNameFile.str()<<" has been readout !"
		<<G4endl<<G4endl;
	}
      return 1;
    }  //  if check

  else return -1;
}  //  ReadOfData

//  ++++++++++++++++++++   GetQ2limit  +++++++++++++++++++
     G4double G4ElasticHadrNucleusHE::
                   GetQ2limit(G4double aR1)
   {
     G4double mQ2 = 35*1000*1000/aR1/aR1;

     dQ2   = mQ2/(Nstep/*3./2.*/-1.);

     iQ2[0] = 1;

     for(G4int ii=1; ii<Nstep; ii++)
          iQ2[ii] = ii<Nstep/2 ? iQ2[ii-1]+dQ2 : iQ2[ii-1]+2*dQ2;

     return mQ2;
   }
//  ========================================================
G4bool G4ElasticHadrNucleusHE::
       GetHadronNucleusData(G4DynamicParticle * aParticle,
                            G4Nucleus         * aNucleus,
                            ElasticData       & ElD)
{
  ElD.Clean();

  G4int     nN       = (G4int) aNucleus->GetN();
  G4String  hadrName = aParticle->GetDefinition()->GetParticleName();

  G4ParticleDefinition * PartDef = aParticle->GetDefinition();

  size_t SizeData = SetOfElasticData.size();

  G4int NumberOfRecord = -1;

  if( SizeData!= 0)
    {
      for(size_t kk = 0; kk<SizeData; kk++)
	{
	  if( SetOfElasticData[kk].nuclAtomicNumber == nN &&
	      SetOfElasticData[kk].hadrName == hadrName)
	    {
	      NumberOfRecord = kk;
	      break;
	    }
	}
    }

  if(NumberOfRecord >= 0)
    ElD = SetOfElasticData[NumberOfRecord];
  else
    {
      G4int TestFile = ReadOfData(PartDef, aNucleus);
      SizeData       = SetOfElasticData.size();
      NumberOfRecord = SizeData-1;

      if(TestFile < 0)
	{
	  G4cout<<" File for elastic scattering of hadron  \""
		<<hadrName<<"\""<<G4endl;
	  G4cout<<"  on nucleus \""
		<<nN<<"\"  does not exist!"<<G4endl;
	  G4cout<<"  You must prepare it by costructor: "<<G4endl;
	  G4cout<<
	    "  G4ElasticHadrNucleusHE(G4DynamicParticle *, G4Nucleus *) "
		<<G4endl;
	  return false;
	}
      else ElD = SetOfElasticData[NumberOfRecord];
    }
  return true;
}

//  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
G4double G4ElasticHadrNucleusHE::SampleT1(
                          const G4ParticleDefinition * p,
                          G4double pTotLabMomentum, 
                          G4int Z,  G4int A)
{
  G4double aQ2 = 0.0;

  G4DynamicParticle    * dHadron = &aHad;
  G4ThreeVector HadrMomentum(0.0, 0.0, pTotLabMomentum);

  dHadron->SetDefinition(const_cast<G4ParticleDefinition*>(p));
  dHadron->SetMomentum(HadrMomentum);
    
  if(A > 1)
    {
      G4Nucleus   InterNucleus;
      ElasticData ElD;      

      InterNucleus.SetParameters(A, Z);
      if(GetHadronNucleusData(dHadron, &InterNucleus, ElD))
	aQ2 = RandomElastic1(dHadron, & ElD);
    }
  else  aQ2 = HadronProtonQ2(dHadron);

  return aQ2;
}

//  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
G4double G4ElasticHadrNucleusHE::SampleT(
                          const G4ParticleDefinition * p,
                          G4double pTotLabMomentum, 
                          G4int,  G4int A)
{
  G4double aQ2 = 0.0;

  G4DynamicParticle    * dHadron = &aHad;
  G4ThreeVector HadrMomentum(0.0, 0.0, pTotLabMomentum);
  G4ParticleDefinition* pd = const_cast<G4ParticleDefinition*>(p); 
  dHadron->SetDefinition(pd);
  dHadron->SetMomentum(HadrMomentum);
 
  if(A > 1)
    {          
      G4String hName = p->GetParticleName();
      G4int  Step = 1;
      ElasticData * ElD = 0;

      G4int NumberOfRecord = -1;
      G4int SizeData = SetOfElasticData.size();

      if( SizeData!= 0)
	{
	  for(G4int kk = 0; kk<SizeData; kk++)
	    {
	      if( SetOfElasticData[kk].nuclAtomicNumber == A &&
                  SetOfElasticData[kk].hadrName == hName)
		{
		  NumberOfRecord = kk;
		  break;
		}

	      if(NumberOfRecord >= 0)
		ElD = &SetOfElasticData[NumberOfRecord];
	    }
	}
      if(!ElD) 
	{
	  ElasticData      ElD2(hName, A);
	  SetOfElasticData.push_back(ElD2);
	  G4int n       = SetOfElasticData.size();
	  ElD = &SetOfElasticData[n-1];
	}
 
      if(ElD) aQ2 = HadronNucleusQ2_2(dHadron, A, Step, ElD);

   }   //  if Amass

  else  aQ2 = HadronProtonQ2(dHadron);

  return  aQ2;
}

//  ++++++++++++++++++  ApplayYourself  ++++++++++++++++++++
G4HadFinalState * G4ElasticHadrNucleusHE::ApplyYourself(
                          const  G4HadProjectile  &aHadron,
                                 G4Nucleus        &aNucleus)
{
  iContr = 133;

  if(iContr==137)
    G4cout<<G4endl<<"----------- Begin Applay ------------"<<G4endl;

  G4Nucleus * aNucl      = &aNucleus;

  theParticleChange.Clear();

//  ---------------  Hadron Definition  ---------------

  G4DynamicParticle    * dParticle = &aHad;
  G4ParticleDefinition * hadrDef = 
    const_cast<G4ParticleDefinition *>(aHadron.GetDefinition());

  G4ThreeVector  hadrMomentum = aHadron.Get4Momentum().vect();
  dParticle->SetDefinition(hadrDef);
  dParticle->SetMomentum(hadrMomentum);

//  ---------------  Nucleus Definition  ---------------
  G4double A = aNucl->GetN();
  G4int nA   = (G4int) A;
  G4int nZ   = (G4int) aNucl->GetZ();

  G4ParticleDefinition * secNuclDef = 0;

  if(nZ == 1 && nA == 1)      secNuclDef = G4Proton::Proton();
  else if(nZ == 1 && nA == 2) secNuclDef = G4Deuteron::Deuteron();
  else if(nZ == 1 && nA == 3) secNuclDef = G4Triton::Triton();
  else if(nZ == 2 && nA == 3) secNuclDef = G4He3::He3();
  else if(nZ == 2 && nA == 4) secNuclDef = G4Alpha::Alpha();
  else secNuclDef = G4ParticleTable::GetParticleTable()->FindIon(nZ,nA,0,nZ);

//  ----------------------------------------------------
// ----------------_ Randomization of Q2 as dSig/dt --------------

  G4double ranQ2 = SampleT(hadrDef,dParticle->GetTotalMomentum(),nZ,nA);

//  ----------------  Hadron kinematics ----------------

  G4double m1 = hadrDef->GetPDGMass();
  G4double m2 = secNuclDef->GetPDGMass();
  G4LorentzVector lv1 = dParticle->Get4Momentum();
  G4LorentzVector lv0(0.0,0.0,0.0,m2);
  //  G4cout << "lv0= " << lv0 << "  lv1= " << lv1 << "   m1= " 
  //	 << m1 << " ranQ2= " << ranQ2<<G4endl;
  G4LorentzVector lv  = lv0 + lv1;
  G4ThreeVector bst = lv.boostVector();
  lv1.boost(-bst);
  lv0.boost(-bst);
  G4ThreeVector p1 = lv1.vect();
  G4double ptot = p1.mag();

  // Sampling in CM system
  G4double phi  = G4UniformRand()*twopi;
  G4double cost = 1. - 0.5*ranQ2/(ptot*ptot);

  if(std::abs(cost) > 1.0) cost = 2.0*G4UniformRand() - 1.0;
  G4double sint = std::sqrt((1.0-cost)*(1.0+cost));

  // G4cout << "Entering elastic scattering 3"<<G4endl;
  if (verboseLevel > 1) 
    G4cout << "cos(t)=" << cost << " std::sin(t)=" << sint << G4endl;

  G4ThreeVector v1(sint*std::cos(phi),sint*std::sin(phi),cost);
  v1 *= ptot;
  G4LorentzVector nlv1(v1.x(),v1.y(),v1.z(),std::sqrt(ptot*ptot + m1*m1));
  G4LorentzVector nlv0 = lv0 + lv1 - nlv1;
  nlv0.boost(bst);
  nlv1.boost(bst); 

  G4double eFinal = nlv1.e() - m1;
  if (verboseLevel > 1) 
    G4cout << " P0= "<< nlv0 << "   P1= "
	   << nlv1<<" m= " << m1 << " ekin0= " << eFinal 
	   << " ekin1= " << nlv0.e() - m2 
	   <<G4endl;
  if(eFinal < 0.0) {
    G4cout << "G4HadronElastic WARNING ekin= " << eFinal
	   << " after scattering of " 
	   << hadrDef->GetParticleName()
	   << " p(GeV/c)= " << dParticle->GetTotalMomentum()/GeV
	   << " on " << secNuclDef->GetParticleName()
	   << G4endl;
    eFinal = 0.0;
  }

  theParticleChange.SetMomentumChange(nlv1.vect().unit());
  theParticleChange.SetEnergyChange(eFinal);
  
  G4DynamicParticle * aSec = new G4DynamicParticle(secNuclDef, nlv0);
  theParticleChange.AddSecondary(aSec);

  return &theParticleChange;
}

//  ########################################################
  G4double G4ElasticHadrNucleusHE::
               HadronNucleusQ2_2(G4DynamicParticle * dHadron,
                                 G4int                  nucN,
                                 G4int                  Step,
                                 ElasticData          * pElD)
{
  G4int  Step1=0, ii;
  G4double Rand = G4UniformRand();

  G4int      kk, NumbOnE, iNumbQ2;
  G4double * dNumbQ2, * dNumbFQ2;
  G4double aQ2 = 0.0;
  G4double Buf;
  Nstep = ONQ2;

  iContr = 2;

  G4String  hadrName = dHadron->GetDefinition()->GetParticleName();
  G4double  Energy = dHadron->GetTotalEnergy();

  NumbOnE = pElD->GetNumberE(Energy);

  dNumbQ2 = pElD->TableQ2;
  if(pElD->TableCrossSec[NumbOnE*ONQ2+2]<0.00001) Step=0;

  if(Step==0)
    {
      G4HadronValues::GetHadronValues(dHadron);
      NucPar = pElD->NucPar;
      aQ2 = pElD->theMaxQ2;

      Weight = GetLightFq2(nucN, aQ2, Step1, &pElD->NucPar);

      iIntgr[0] = 0;
      pElD->TableCrossSec[NumbOnE*ONQ2] = 0;
// ------------------------------------------
      for(ii=1; ii<Nstep/*3*/; ii++)
	{
	  aQ2      = pElD->TableQ2[ii];

	  Buf = GetLightFq2(nucN, aQ2, Step1, &pElD->NucPar);
	  pElD->TableCrossSec[NumbOnE*ONQ2+ii] = Buf/Weight;
	  RandMax  = Buf/Weight;

	}   // for ii
      CurrentN = 3;

    }  //  if Step

  ii = 0;
  while(pElD->TableCrossSec[NumbOnE*ONQ2+ii]>0.0001) 
    {
      RandMax  = pElD->TableCrossSec[NumbOnE*ONQ2+ii];
      CurrentN = ii;
      ii++;
    }

  dNumbFQ2 = &pElD->TableCrossSec[NumbOnE*ONQ2];

  if(RandMax > Rand)
    {
      kk = 0;
      while(Rand>pElD->TableCrossSec[NumbOnE*ONQ2+kk]) kk++;
      iNumbQ2 = kk;
      aQ2 = GetQ2_2(iNumbQ2, dNumbQ2, dNumbFQ2, Rand);
      return aQ2;
    }
  else
    {
      Step = 1;
      for(ii=CurrentN; ii<Nstep; ii++)
	{ 

	  aQ2  = pElD->TableQ2[ii];  //(ii-0)*dQ2;
	  Buf = GetLightFq2(nucN, aQ2, Step, &pElD->NucPar);
	  pElD->TableCrossSec[NumbOnE*ONQ2+ii] = Buf/Weight;
	  if(Buf/Weight>Rand) 
	    {
	      iNumbQ2 = ii;
	      aQ2 = GetQ2_2(iNumbQ2, dNumbQ2, dNumbFQ2, Rand);
	      CurrentN = ii;
	      RandMax  = Buf/Weight;
	      break;
		//          return Q2;    
	    }  //  if Rand
	}    //  for ii
    }      //  else    
  return  aQ2;
}        //  function

//  =========================================================
//  +++++++  The randomization of one dimensional array ++++++
  G4double G4ElasticHadrNucleusHE::GetQ2_2(G4int kk, G4double * Q,
                                       G4double * F, G4double ranUni)
  {
   G4double ranQ2;

   if(kk == 1)
   {
     ranQ2 = (*(Q+kk-1)+(ranUni-*(F+kk-1)))*
	      (*(Q+kk)-(*(Q+kk-1)))/(*(F+kk)-(*(F+kk-1)));
    return ranQ2;
   }

   G4double F1  = *(F+kk-2);
   G4double F2  = *(F+kk-1);
   G4double F3  = *(F+kk);

   G4double F12 = F1*F1;
   G4double F22 = F2*F2;
   G4double F32 = F3*F3;

   G4double X1  = *(Q+kk-2);          //  MeV^2
   G4double X2  = *(Q+kk-1);
   G4double X3  = *(Q+kk);

   G4double D0  = F12*F2+F1*F32+F3*F22-F32*F2-F22*F1-F12*F3;

   if(std::abs(D0) < 0.00000001 || kk == 2)
   { 
    ranQ2 = (*(Q+kk-1)+(ranUni-*(F+kk-1))*
                 (*(Q+kk)-*(Q+kk-1))
                 /(*(F+kk)-*(F+kk-1)));

    }
   else    
    {
      G4double DA = X1*F2+X3*F1+X2*F3-X3*F2-X1*F3-X2*F1;
      G4double DB = X2*F12+X1*F32+X3*F22-X2*F32-X3*F12-X1*F22;
      G4double DC = X3*F2*F12+X2*F1*F32+X1*F3*F22
                          -X1*F2*F32-X2*F3*F12-X3*F1*F22;
      ranQ2 = (DA*ranUni*ranUni+DB*ranUni+DC)/D0;
    }
  if(ranUni>0.99995)
   if(ranQ2<-10000)
  {
 G4cout<<" X1 X2 X3 "<<X1<<"  "<<X2<<"  "<<X3
       <<" F1 F2 F3 "<<F1<<"  "<<F2<<"  "<<F3<<G4endl;
 G4cout<<"  ranUni  "<<ranUni<<"  ranQ2 "<<ranQ2
       <<" D0 "<<D0<<G4endl;
  }
   return ranQ2;         //  MeV^2
  }
//  ==========================================================
G4double G4ElasticHadrNucleusHE::  
        GetLightFq2(G4int Nucleus, G4double aQ2, G4int,
           NucleusParameters *HP)
{

// ----------------  The preparing of probability function  ------------

//    G4double    MbToB    = 2.568;             //  from mb to GeV^-2
//    G4double    Pi1      = 3.1416;
    G4double    Stot     = HadrTot*MbToGeV2;     //  Gev^-2
    G4double    Bhad     = HadrSlope;         //  GeV^-2
    G4double    Asq      = 1+HadrReIm*HadrReIm;
    G4double    Rho2     = std::sqrt(Asq);

    G4double    R12      = HP->R1*HP->R1;
    G4double    R22      = HP->R2*HP->R2;
    G4double    R12B     = R12+2*Bhad;
    G4double    R22B     = R22+2*Bhad;
///      G4double    R12Bp    = R12+20;
///      G4double    R22Bp    = R22+20;
///      G4double    R13Bp    = R12*R1/R12Bp;
///      G4double    R23Bp    = R22*R2/R22Bp;
///      G4double    R12Ap    = R12+20;
///      G4double    R22Ap    = R22+20;
///      G4double    R13Ap    = R12*R1/R12Ap;
///      G4double    R23Ap    = R22*R2/R22Ap*PnuclP;
///      G4double    R23dR13  = R23Ap/R13Ap;
///      G4double    R12Apd   = 2/R12Ap;
///      G4double    R22Apd   = 2/R22Ap;

    G4double    Norm     = R12*HP->R1-HP->Pnucl*R22*HP->R2; //*HP->Aeff;
///      G4double    NormP    = R12*R1-PnuclP*R22*R2;
    G4double    R13      = R12*HP->R1/R12B;
    G4double    R23      = HP->Pnucl*R22*HP->R2/R22B;
    G4double    Unucl    = Stot/2/MyPi/Norm*R13;
////      G4double    Unclprod = Stot/2/Pi1/NormP*R13Ap;
    G4double    FiH      = std::asin(HadrReIm/Rho2);
    G4double    NN2      = R23/R13;

    G4double    dddd;
///      G4double    DDSec1p  = (DDSect2+
///                    DDSect3*std::log(1.06*2*Ehad/R1/4));

///      G4double    DDSec2p  = (DDSect2+
///                    DDSect3*std::log(1.06*2*Ehad/
//                     std::sqrt((R12+R22)/2)/4));

///      G4double    DDSec3p  = (DDSect2+
///                    DDSect3*std::log(1.06*2*Ehad/R2/4));

///      G4double    R12ApdR22Ap = 0.5*(R12Apd+R22Apd);

//                 iIntgr[0] = 0;

//    G4int ii;

//    for(ii=0; ii<Nstep-1; ii++) 
//    {

//    G4double    Q2       = iQ2[ii+1]/1000/1000;// internal in GeV^2

    G4double Prod0    = 0;
    G4double N1       = -1;
    G4double Tot0     = 0;
    G4double exp1;

    G4double Prod3 ;
    G4double exp2  ;
    G4double N4, N5, N2, Prod1, Prod2;
    G4int    i1, i2, m1, m2;


//    G4cout<<" Light 1 Q2  "<<Q2<<G4endl;

      for(i1 = 1; i1<= Nucleus; i1++) ////++++++++++  i1
      {
        N1    = -N1*Unucl*(Nucleus-i1+1)/i1*Rho2;
        Prod1 = 0;
        Tot0  = 0;
        N2    = -1;

        for(i2 = 1; i2<=Nucleus; i2++) ////+++++++++ i2
        {
          N2    = -N2*Unucl*(Nucleus-i2+1)/i2*Rho2;
          Prod2 = 0; //std::exp(-Q2/i2*R12B/4)/i2*R12B;
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
                            (1-std::exp(-aQ2*dddd/*(1/exp1+1/exp2)*//4))/
		            dddd/*(1/exp1+1/exp2)*/*4*SetBinom[i1][m1];
//                             Factorials1[i1]/
//                             Factorials1[m1]/
//                             Factorials1[i1-m1];
               }                                   // m1
               Prod2 = Prod2 +Prod3*N5*SetBinom[i2][m2];
//                 Factorials1[i2]/
//                 Factorials1[m2]/
//                Factorials1[i2-m2];

           }                                      // m2
         Prod1 = Prod1 + Prod2*N2*std::cos(FiH*(i1-i2));
//         Tot0  = Tot0  + Prod2*N2*std::sin(FiH*(i1-i2));

         if (std::fabs(Prod2*N2/Prod1)<1e-6) break;
        }                                         // i2
//        ImDistr = Tot0  + Tot0*N1;
         Prod0   = Prod0 + Prod1*N1;
         if(std::fabs(N1*Prod1/Prod0) < 1e-6) break;
      }                                           // i1
        Prod0        = Prod0*MyPi/MbToGeV2/4;  //  This is in mb
//        iIntgr[ii+1] = Prod0;  //iIntgr[ii-1]+Prod0;
//    }                                            //   ii (Q2)
   return Prod0;
  }

//  ##########################################################    
  G4double G4ElasticHadrNucleusHE::
               HadronNucleusQ2(G4DynamicParticle * aHadron,
                               G4Nucleus         & aNucleus)
  {
    G4int       nN;
    G4double    ranQ2 = 0.0;
    ElasticData ElD;

    nN = (G4int) aNucleus.GetN();

    if(nN<2) ranQ2 = HadronProtonQ2(aHadron);

    else
    {
      if(GetHadronNucleusData(aHadron, & aNucleus, ElD))
	ranQ2 = RandomElastic1(aHadron, & ElD);
    }

    return ranQ2;
  }
//  +++++++  The randomization of one dimensional array ++++++
  G4double G4ElasticHadrNucleusHE::RandomElastic0()
  {
   G4double ranQ2;
   G4double ranUni = G4UniformRand();

   G4int kk = 0;

   while( ranUni > iIntgr[kk])  kk ++;
         
   if(kk == 1 ) kk = 2;

   G4double F1  = iIntgr[kk-2];
   G4double F2  = iIntgr[kk-1];
   G4double F3  = iIntgr[kk];
   G4double F12 = F1*F1;
   G4double F22 = F2*F2;
   G4double F32 = F3*F3;
   G4double X1  = iQ2[kk-2];          //  MeV^2
   G4double X2  = iQ2[kk-1];
   G4double X3  = iQ2[kk];
   G4double D0  = F12*F2+F1*F32+F3*F22-F32*F2-F22*F1-F12*F3;

   if(D0 == 0) 
    ranQ2 = (iQ2[kk-1]+(ranUni-iIntgr[kk-1])*
                 (iQ2[kk]-iQ2[kk-1])
                 /(iIntgr[kk]-iIntgr[kk-1]));
   else    
    {
      G4double DA = X1*F2+X3*F1+X2*F3-X3*F2-X1*F3-X2*F1;
      G4double DB = X2*F12+X1*F32+X3*F22-X2*F32-X3*F12-X1*F22;
      G4double DC = X3*F2*F12+X2*F1*F32+X1*F3*F22
                          -X1*F2*F32-X2*F3*F12-X3*F1*F22;
      ranQ2 = (DA*ranUni*ranUni+DB*ranUni+DC)/D0;
    }
   return ranQ2;         //  MeV^2
  }     
// +++++++++ The randomization of two dimensoinal array   ++++++++
    G4double G4ElasticHadrNucleusHE::RandomElastic1(
                           const G4DynamicParticle *  aHadron,
                           const ElasticData  * ElD)
  {
    G4double forQ2[3];
    G4double curE = aHadron->GetTotalEnergy();   //  in  MeV
    G4double X1;
    G4double X2;
    G4double X3;
    G4double Y1;
    G4double Y2;
    G4double Y3;
    G4int    iiE = 0;
    G4int    k[3];

    while(curE > ElD->TableE[iiE] ) iiE++;

    G4double ranQ2;
    G4double ranUni = G4UniformRand();

    G4int kk;

    for(G4int mm=0; mm<3; mm++)
    {
       kk = 0;

       while( ranUni > 
            ElD->TableCrossSec[iiE*ONQ2+mm*ONQ2+kk-ONQ2])  
      {
        kk ++;
        k[mm]=kk;
      }     //   while   

//           if(kk == 1 ) kk = 2;
      X1 = ElD->TableCrossSec[iiE*ONQ2+mm*ONQ2+kk-ONQ2-1];
      X2 = ElD->TableCrossSec[iiE*ONQ2+mm*ONQ2+kk-ONQ2];
      X3 = ElD->TableCrossSec[iiE*ONQ2+mm*ONQ2+kk-ONQ2+1];
      Y1 = ElD->TableQ2[kk-1];            //  MeV^2
      Y2 = ElD->TableQ2[kk];
      Y3 = ElD->TableQ2[kk+1];

      forQ2[mm] = InterPol(X1, X2, X3, Y1, Y2, Y3, ranUni);
    }     //  mm

    X1 = ElD->TableE[iiE-1];
    X2 = ElD->TableE[iiE];
    X3 = ElD->TableE[iiE+1];
    Y1 = forQ2[0];             //  MeV^2
    Y2 = forQ2[1];
    Y3 = forQ2[2];

    ranQ2 = InterPol(X1, X2, X3, Y1, Y2, Y3, curE);

    return ranQ2;          //   MeV^2
 }
//  ++++++++++++ The filling of the one-dimensional array +++++++++++++
void  G4ElasticHadrNucleusHE::CreationArray(
              const G4DynamicParticle * aHadron,
                    G4Nucleus         * aNucleus)
  {
    GetNucleusParameters(aNucleus);

    //    Nstep  = ONQ2;
    // theMaxQ2  = GetQ2limit(theR1);             //   Mev^2

    if(aNucleus->GetN() > MaxN )
           ArrayForHeavy( aHadron, aNucleus);
    else
           ArrayForLight( aHadron, aNucleus);
  }
//  +++++++++++++ The  method for heavy nucleus  ++++++++++++++++
void G4ElasticHadrNucleusHE::ArrayForHeavy( 
            const G4DynamicParticle * aHadron,
                  G4Nucleus         * aNucleus)
  {
   GetNucleusParameters(aNucleus);

//////       G4double aNuc  = aNucleus->GetN();
   G4double intgrS=0., intgStep;
   G4int    ii, Kind=0;

   iIntgr[0] = 0;

   for( ii=0; ii<Nstep-1; ii++) 
   {

// ----------------  The integration  -------------

     G4double  curQ2  = iQ2[ii];             //  MeV^2
     G4double  ddQ2   =(iQ2[ii+1]-curQ2)/20; //  MeV^2
     G4double  curSum = 0; 
     G4double  curSec = 0;
     G4double aSimp = 4;
 
        for(G4int ll=0; ll<20; ll++)                     
        {                                              
//         curSec = aDiffElHadNcls.
//            HadrNuclDifferCrSec(aHadron, aNucleus, curQ2); 
           curSec = aDiffElHadNcls.
                     HadrNuclDifferCrSecT(aHadron, aNucleus, 
					  curQ2, Kind);

           curQ2  = curQ2 + ddQ2;                        
           curSum = curSum + curSec*aSimp;        

           if(aSimp > 3 ) aSimp = 2;
           else  aSimp = 4;
//           if(ll == 0) aSimp = 4;
           Kind = 1;
         }      //  ll                           
     intgStep = curSum*ddQ2/3;      
//  --------------------------------------------------------------------
     iIntgr[ii+1] = iIntgr[ii]+intgStep; 
     intgrS       = intgrS + intgStep;
   }   //  ii
//  -------------------  The normalization  ----------------
      for( ii=0; ii<Nstep; ii++)
               iIntgr[ii] = iIntgr[ii]/iIntgr[Nstep-1];
 }             //  end  subroutine
//  +++++++++++++++++  The  method for light nucleus  +++++++++++++++++++
void G4ElasticHadrNucleusHE::ArrayForLight( 
                          const G4DynamicParticle * aHadron,
                                G4Nucleus         * aNucleus)
{
// -------------- !!! Here the main unit is GeV !!!  ---------------
    G4HadronValues::GetHadronValues(aHadron);
    G4double aNuc     = aNucleus->GetN();
    G4int    Nucleus  = (G4int) aNuc;

    GetNucleusParameters(aNucleus);

    //    theMaxQ2    = GetQ2limit(theR1)/GeV/GeV;        //  GeV^2
// ----------------  The preparing of probability function  ------------

    G4double    MbToB    = 2.568;             //  from mb to GeV^-2
    G4double    Pi1      = pi;
    G4double    Stot     = HadrTot*MbToB;     //  Gev^-2
    G4double    Bhad     = HadrSlope;         //  GeV^-2
    G4double    Asq      = 1+HadrReIm*HadrReIm;
    G4double    Rho2     = std::sqrt(Asq);

    G4double    R12      = R1*R1;
    G4double    R22      = R2*R2;
    G4double    R12B     = R12+2*Bhad;
    G4double    R22B     = R22+2*Bhad;
///      G4double    R12Bp    = R12+20;
///      G4double    R22Bp    = R22+20;
///      G4double    R13Bp    = R12*R1/R12Bp;
///      G4double    R23Bp    = R22*R2/R22Bp;
///      G4double    R12Ap    = R12+20;
///      G4double    R22Ap    = R22+20;
///      G4double    R13Ap    = R12*R1/R12Ap;
///      G4double    R23Ap    = R22*R2/R22Ap*PnuclP;
///      G4double    R23dR13  = R23Ap/R13Ap;
///      G4double    R12Apd   = 2/R12Ap;
///      G4double    R22Apd   = 2/R22Ap;

    G4double    Norm     = R12*R1-Pnucl*R22*R2;
///      G4double    NormP    = R12*R1-PnuclP*R22*R2;
    G4double    R13      = R12*R1/R12B;
    G4double    R23      = Pnucl*R22*R2/R22B;
    G4double    Unucl    = Stot/2/Pi1/Norm*R13;
////      G4double    Unclprod = Stot/2/Pi1/NormP*R13Ap;
    G4double    FiH      = std::asin(HadrReIm/Rho2);
    G4double    NN2      = R23/R13;

///      G4double    DDSec1p  = (DDSect2+
///                    DDSect3*std::log(1.06*2*Ehad/R1/4));

///      G4double    DDSec2p  = (DDSect2+
///                    DDSect3*std::log(1.06*2*Ehad/
//                     std::sqrt((R12+R22)/2)/4));

///      G4double    DDSec3p  = (DDSect2+
///                    DDSect3*std::log(1.06*2*Ehad/R2/4));

///      G4double    R12ApdR22Ap = 0.5*(R12Apd+R22Apd);

                 iIntgr[0] = 0;

    G4int ii;

    for(ii=0; ii<Nstep-1; ii++) 
    {

    G4double    aQ2       = iQ2[ii+1]/1000/1000;// internal in GeV^2

    G4double    Prod0    = 0;
    G4double    N1       = -1;
    G4double    Tot0     = 0;
    G4double    exp1;

      for(G4int i1 = 1; i1<= Nucleus; i1++) ////++++++++++  i1
      {
        N1              = -N1*Unucl*(Nucleus-i1+1)/i1*Rho2;
        G4double Prod1  = 0;
        Tot0            = 0;
        G4double  N2    = -1;

        for(G4int i2 = 1; i2<=Nucleus; i2++) ////+++++++++ i2
        {
          N2    = -N2*Unucl*(Nucleus-i2+1)/i2*Rho2;
          G4double Prod2  = 0; //std::exp(-Q2/i2*R12B/4)/i2*R12B;
          G4double N5     = -1/NN2;

            for(G4int m2=0; m2<= i2; m2++) ////+++++++++ m2
            {
             G4double Prod3 = 0;
             G4double exp2  = m2/R22B+(i2-m2)/R12B;
             N5             = -N5*NN2;
             G4double N4    = -1/NN2;

               for(G4int m1=0; m1<=i1; m1++) ////++++++++ m1
               {
                 exp1  = m1/R22B+(i1-m1)/R12B;
                 N4    = -N4*NN2;
                 Prod3 = Prod3+N4/exp1/exp2*
                            (1-std::exp(-aQ2*(1/exp1+1/exp2)/4))/
                            (1/exp1+1/exp2)*4*
                             Factorials1[i1]/
                             Factorials1[m1]/
                             Factorials1[i1-m1];

               }                                   // m1

             Prod2 = Prod2 +Prod3*N5*
                 Factorials1[i2]/
                 Factorials1[m2]/
                 Factorials1[i2-m2];

           }                                      // m2
         Prod1 = Prod1 + Prod2*N2*std::cos(FiH*(i1-i2));
         Tot0  = Tot0  + Prod2*N2*std::sin(FiH*(i1-i2));

         if (std::fabs(Prod2*N2/Prod1)<1e-6) break;
        }                                         // i2
//        ImDistr = Tot0  + Tot0*N1;
         Prod0   = Prod0 + Prod1*N1;

         if(std::fabs(N1*Prod1/Prod0) < 1e-6) break;
      }                                           // i1
        Prod0        = Prod0*Pi1/2.568/4;  //  This is in mb
        iIntgr[ii+1] = Prod0;  //iIntgr[ii-1]+Prod0;
    }                                            //   ii (Q2)

      for(ii = 0; ii<Nstep; ii++)
                  iIntgr[ii]=iIntgr[ii]/iIntgr[Nstep-1];
   }           //  subroutine
//  ++++++++++++++++++++  The interpolation +++++++++++++++++++++
G4double G4ElasticHadrNucleusHE::InterPol(
                  G4double X1, G4double X2, G4double X3,
                  G4double Y1, G4double Y2, G4double Y3, 
                  G4double X)
   {
      G4double ranQ2;
      G4double F12 = X1*X1;
      G4double F22 = X2*X2;
      G4double F32 = X3*X3;
      G4double D0  = F12*X2+X1*F32+X3*F22-F32*X2-F22*X1-F12*X3;

      if(std::fabs(D0) < 1e-8 || D0 == 0) 
                  ranQ2 = (Y2+(X-X2)*(Y3-Y2) /(X3-X2));   //   MeV^2

      else    
      {
       G4double DA =  Y1*X2    +Y3*X1    +Y2*X3 -Y3*X2 -Y1*X3 -Y2*X1;
       G4double DB =  Y2*F12   +Y1*F32   +Y3*F22-Y2*F32-Y3*F12-Y1*F22;
       G4double DC =  Y3*X2*F12+Y2*X1*F32+Y1*X3*F22
                           -Y1*X2*F32-Y2*X3*F12-Y3*X1*F22;
             ranQ2 = (DA*X*X+DB*X+DC)/D0;           //   MeV^2
      }
      return  ranQ2;
    }
//  ============================================================
  G4float G4ElasticHadrNucleusHE::GetFt(G4double aQ2)
  {
    G4float Fdistr, SqrQ2 = std::sqrt(aQ2);

    Fdistr = (1-Coeff1-Coeff0/*-0.0*Coeff2*std::exp(ConstU)*/)
                   /Slope*(1-std::exp(-Slope*aQ2))

                  +Coeff0*(1-std::exp(-Slope0*aQ2))

                  +Coeff2/Slope2*std::exp(Slope2*ConstU)*
                      (std::exp(Slope2*aQ2)-1)

                  +2*Coeff1/Slope1*(1/Slope1-(1/Slope1+SqrQ2)*
                  std::exp(-Slope1*SqrQ2))
                   ;
    return Fdistr;
  }
//  +++++++++++++++++++++++++++++++++++++++
  G4float G4ElasticHadrNucleusHE::GetDistrFun(G4double aQ2)
  {
    G4float Fq2;

    Fq2   = GetFt(aQ2);

    return Fq2/FmaxT;
  }
//  +++++++++++++++++++++++++++++++++++++++
  G4double G4ElasticHadrNucleusHE::
              GetQ2(G4double Ran)
  {
    G4double DDD0=MaxTR*0.5, DDD1=0.0, DDD2=MaxTR, delta;
    G4double aQ2;

    FmaxT = GetFt(MaxTR);
    delta = GetDistrFun(DDD0)-Ran;

    while(std::fabs(delta) > 0.0001)
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

    aQ2 = DDD0;

    return aQ2;
  }
//  ++++++++++++++++++++++++++++++++++++++++++
  G4double G4ElasticHadrNucleusHE::
       HadronProtonQ2(G4DynamicParticle * aHadron)
  {
    G4double Ran = G4UniformRand(), aQ2;

    GetKinematics(aHadron);

//    G4cout<<"  HRQ2HE Ran "<<Ran<<G4endl;

    aQ2 = GetQ2(Ran);

    return aQ2*GeV*GeV;
  }

//  ===========================================

/*  End of file  */

