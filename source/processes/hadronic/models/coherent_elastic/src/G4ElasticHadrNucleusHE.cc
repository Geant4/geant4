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
// $Id: G4ElasticHadrNucleusHE.cc,v 1.36 2006-10-20 16:43:44 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//G4ElasticHadrNucleusHE.cc
//
//

//  The generator of high energy hadron-nucleus elastic scattering
//  The hadron kinetic energy T > 1 GeV
//  N.  Starkov 2003.

//  19.05.04 Variant for G4 6.1: The 'ApplyYourself' was changed

//  November 2005 - The HE elastic scattering on proton is added
//  N. Starkov

//#include  <strstream>
#include  <sstream>
#include  <iostream>
#include  <ostream>
#include  <fstream>

#include  "G4ElasticHadrNucleusHE.hh"
#include  "Randomize.hh"
#include  "G4ios.hh"
#include  "G4ParticleTable.hh"
#include  "G4IonTable.hh"

  using namespace std;
//  ############################################################
   ElasticData:: ElasticData(G4String HadrName, G4int AtomWeight)
   { 
     G4int      ii, kk;
     G4double   Expo=1, Kinetic[5] = {1.0, 1.5, 2.5, 4.0, 6.0};

     hadrName = HadrName;
     AtomicWeight = AtomWeight;

     GetNucleusParameters(AtomWeight);

     GetQ2limit(R1);

     verboselevel = 0;
  if(verboselevel == 1)
  G4cout<<" ElastData(constr). R1 R2 Aeff Pnucl "<<R1<<"  "<<R2<<"  "
        <<Aeff<<"  "<<Pnucl<<G4endl;

     for(kk=0; kk<AreaNumb; kk++)
     {
       for(ii=0; ii<ONE;ii++)
       {
         TableE[ii+ONE*kk] = Kinetic[ii]*Expo;

  if(verboselevel == 1)
  G4cout<<" Constr ElD: Po E ii  Expo  E  "<<ii<<"  "
        <<Expo<<"  "<< TableE[ii+ONE*kk]<<G4endl;
       }
         Expo      *= 10;
       for(ii=0; ii<ONQ2XE; ii++)
         TableCrossSec[ii+kk*ONQ2XE] = 0;

  if(verboselevel == 1)
  G4cout<<" ConstrElD: E "<<TableE[ii+ONE*kk-1]<<G4endl;
     }

  if(verboselevel == 1)
  G4cout<<" ConstrElD: maxQ2 "<<maxQ2<<"  Qi "
        <<TableQ2[1]<<"  "<<TableQ2[2]<<"  "
        <<TableQ2[3]<<G4endl;
   }
//  ###########################################################
   G4double ElasticData::GetQ2limit(G4double R1)
   {
     G4double maxQ2_0 = 35/*1000*1000*//R1/R1;     //  (GeV/c)^2

     Nstep = ONQ2;
     dQ2   = maxQ2_0/(Nstep/*3./2.*/-1.);

     maxQ2 = maxQ2_0;
     TableQ2[0] = 0.00001;

     for(G4int ii=1; ii<Nstep; ii++) TableQ2[ii] = TableQ2[ii-1]+dQ2;

     return maxQ2_0;
   }
//  ############################################################
   G4int  ElasticData::GetNumberE(G4double E)
   {
     G4int ii=0;

  if(verboselevel == 1)
  G4cout<<"  GetNumbE: before cycle ii = "<<ii<<"  TableE "
        <<TableE[ii]<<"  E  "<<E<<G4endl;

    while(E > TableE[ii])
    {
      ii++;

  if(verboselevel == 1)
  G4cout<<"  GetNumbE: cycle ii = "<<ii<<"  TableE "
        <<TableE[ii]<<"  E  "<<E<<G4endl;
    }

  if(verboselevel == 1)
  G4cout<<" GetNumbE: after cycle ii = "<<ii
        <<"  Table[ii] "<<TableE[ii]<<G4endl;

     return ii;
   }
//  ####### The constructor for the generating of events #######
   G4ElasticHadrNucleusHE::G4ElasticHadrNucleusHE()
   {
      MbToGeV2  =  2.568;
      sqMbToGeV =  1.602;
      Fm2ToGeV2 =  25.68;

      Binom();

      verboselevel = -2;
   }
//  ####### The constructor for the preparation of data ########
//  ++++++++++++++++++++   GetQ2limit  +++++++++++++++++++
     G4double G4ElasticHadrNucleusHE::
                   GetQ2limit(G4double R1)
   {
     G4double maxQ2 = 35*1000*1000/R1/R1;    //  (MeV/c)^2

     dQ2   = maxQ2/(Nstep/*3./2.*/-1.);

// SetOfElasticData.size();
//     for(G4int ii=1; ii<Nstep; ii++)
//          iQ2[ii] = ii<Nstep/2 ? iQ2[ii-1]+dQ2 : iQ2[ii-1]+2*dQ2;

     return maxQ2;
   } 
//  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   G4double G4ElasticHadrNucleusHE::SampleT(
                          const G4ParticleDefinition * p,
                          G4double inLabMom, 
                          G4int Z,  G4int N)
  {
    G4String               hName = p->GetParticleName();
    G4double               Q2, Ran, pTotLabMomentum = inLabMom/1000; // (GeV/c)
    G4int                  Amass=Z+N;
    G4ThreeVector HadrMomentum(0.0, 0.0, pTotLabMomentum);

    HadrCode = p->GetPDGEncoding();

 if(verboselevel == 1) 
 G4cout<<" SampleT: AMass "<<Amass<<"  HadrCode "<<HadrCode
       <<" inLabMom  "<<inLabMom<<" pTotLabMomentum "<<pTotLabMomentum<<G4endl;

    if(Amass>1)
    {    

      G4int  Step;
//  ..................................
    ElasticData * ElD1 = 0;
    //    ElasticData * ElD  = 0;

    G4String  hadrName = p->GetParticleName();

      size_t SizeData = SetOfElasticData.size();

      G4int NumberOfRecord = 0;

  if(verboselevel == 2)
  G4cout<<" A>1  sizeElD "<<SizeData<<G4endl;
//  .........................................
      if( SizeData!= 0)
      {
  if(verboselevel == 2)
  G4cout<<"  SampleT: SizeData!=0"<<G4endl;
        for(size_t kk = 0; kk<SizeData; kk++)
        {
  if(verboselevel == 2)
  G4cout<<"  SampleT: SetData kk "<<kk<<" HadrName "<<
                  SetOfElasticData[kk].hadrName
             <<"  "<< hadrName<<G4endl;

          if( SetOfElasticData[kk].AtomicWeight == Amass &&
                  SetOfElasticData[kk].hadrName == hadrName)
              {
                NumberOfRecord = kk;

  if(verboselevel == 2)
  G4cout<<"  SampleT: Hadron "<<SetOfElasticData[kk].hadrName<<G4endl;
              }
        }

  if(verboselevel == 1)
  G4cout<<"  SampleT: After the looking for: NumberOfRecord = "
        <<NumberOfRecord<<G4endl;

        if(NumberOfRecord > 0)
        {
  if(verboselevel == 1)
  G4cout<<"  SampleT: The ElD if found in Set: NumbRecord > 0 "
        <<NumberOfRecord<<G4cout;

          ElD1 = &SetOfElasticData[NumberOfRecord-0];
        }
        Step = 1;
      }       // if SizeData!=0
//  ...........................................
      else    //  SizeData==0 
      {
         ElasticData      ElD2(hName,Amass);
         SetOfElasticData.push_back(ElD2);
         SizeData       = SetOfElasticData.size();
         NumberOfRecord = SizeData-1;
         ElD1 = &SetOfElasticData[NumberOfRecord-0];
         Step = 0;

  if(verboselevel == 1)
  {
   G4cout<<" ### New SetELD ElD2 R1 "
       <<SetOfElasticData[0].R1<<"  "<<ElD2.R1<<G4endl;
   G4cout<<" ### New SetELD ElD2 R2 "
       <<SetOfElasticData[0].R2<<"  "<<ElD2.R2<<G4endl; 
   G4cout<<" ### Set R2 "<<SetOfElasticData[0].R2
       <<"   NumbRecords "<<SizeData<<G4endl;
   G4cout<<" ### Set R2 "<<ElD2.R2
       <<"   NumbRecords "<<SizeData<<G4endl;
   G4cout<<" ### Set dQ2 "<<SetOfElasticData[0].dQ2
       <<"  NumbRecords "<<SizeData<<G4endl;
   G4cout<<" ### Set HadrName "<<SetOfElasticData[0].hadrName
       <<"  TableQ2 "<<SetOfElasticData[0].TableQ2[3]<<G4endl;
  }
     }

  if(verboselevel == 1)
  G4cout<<"  SampleT: The number of members in DataSet  "
        <<SetOfElasticData.size()<<G4endl;
 
      R1    = ElD1->R1;
      R2    = ElD1->R2;
      Aeff  = ElD1->Aeff;
      Pnucl = ElD1->Pnucl;

  if(verboselevel == 1)
  G4cout<<" SampleT. R1 R2 Aeff Pnucl "<<R1<<"  "<<R2<<"  "
        <<Aeff<<"  "<<Pnucl<<" Momentum "<<pTotLabMomentum<<G4endl;

      Ran = G4UniformRand();

      Q2 = HadronNucleusQ2_2(p, Amass, pTotLabMomentum, Ran, Step, ElD1);

// /*const_cast<G4ParticleDefinition*>*/
   }   //  if Amass

    else  Q2 = HadronProtonQ2(p, pTotLabMomentum);

    return  Q2;
  }
/*
//  ++++++++++++++++++  ApplayYourself  ++++++++++++++++++++
    G4HadFinalState *
    G4ElasticHadrNucleusHE::ApplyYourself(
                          const  G4HadProjectile  &aHadron,
                                 G4Nucleus        &aNucleus)
 {
    G4int    nN, nZ;

    G4IonTable  * MyIonTable = const_cast<G4IonTable *>
             (G4ParticleTable::GetParticleTable()->GetIonTable());

 if(verboselevel == 1)
 {
 G4cout<<G4endl<<"----------- Begin Applay ------------"<<G4endl;
 G4cout<<" Total Energy "<<aHadron.GetTotalEnergy()<<G4endl;
 }
    G4Nucleus * aNucl      = &aNucleus;

    G4double    hadrMass;
    G4String    hadrName;

    theParticleChange.Clear();

//  ---------------  Hadron Definition  ---------------
    G4ParticleDefinition * hadrDef1 =
        const_cast<G4ParticleDefinition*>(aHadron.GetDefinition());

    G4ThreeVector  hadrMomentum = aHadron.Get4Momentum().vect();

    G4double   inLabMom  = aHadron.GetTotalMomentum()/1000; // GeV/c
    G4double   inEnHadr  = aHadron.GetTotalEnergy()/1000;   // GeV

 if(verboselevel == 1)
 G4cout<<" Hadr: Energy TotMom  "<<inLabMom<<"  "<<inEnHadr<<G4endl;
//  ---------------  Nucleus Definition  ---------------
    nN   = (int) aNucl->GetN();
    nZ   = (int) aNucl->GetZ();

    G4ParticleDefinition * secNuclDef;
    G4DynamicParticle    * secNuclDyn = new G4DynamicParticle();

    secNuclDef  =   MyIonTable->GetIon(nZ,  nN);
    secNuclDyn->SetDefinition(secNuclDef);
//  ----------------------------------------------------
// ----------------_ Randomization of Q2 as dSig/dt --------------

    G4double ranQ2;

    if(nN>1)  ranQ2 = SampleT(hadrDef1, inLabMom*1000,nN-nZ, nZ);

    else ranQ2 = HadronProtonQ2(aHadron.GetDefinition(), inLabMom*1000);

  if(verboselevel == 1)
  G4cout<<" Applay: Q2 = "<<ranQ2<<G4endl;
//  ----------------  Hadron kinematics ----------------
    G4double   MassHadr  = aHadron.GetDefinition()->GetPDGMass()/1000; // GeV
    G4double   MassNucl  = aNucleus.GetN()*0.93827;        // GeV
    G4double   sqrMass   = MassNucl*MassNucl+MassHadr*MassHadr;
// ---------------    For final state of hadron    ------------------
    G4double   invS      = sqrMass+2*inEnHadr*MassNucl;   // GeV^2
    G4double   invU      = 2*sqrMass-invS+ranQ2;          // GeV^2
    G4double   outEnHadr = (sqrMass-invU)/2/MassNucl;     // GeV
    G4double   outMomHdr = std::sqrt(outEnHadr*outEnHadr- // GeV
                                    MassHadr*MassHadr);
    G4double   cosHadr   = (-ranQ2-2*MassHadr*MassHadr+
                             2*inEnHadr*outEnHadr)/2
                             /inLabMom/outMomHdr;        
//----------------     For final state of nucleus   -----------------
    G4double   outEnNucl = (2*MassNucl*MassNucl+ranQ2)/2/MassNucl;
    G4double   outMomNcl = std::sqrt(outEnNucl*outEnNucl-  //  GeV
                                  MassNucl*MassNucl);
    G4double   cosNucl   = (invU-sqrMass+2*inEnHadr*outEnNucl)
                              /2/inLabMom/outMomNcl;   
    G4double   ranFi     = 6.2832*G4UniformRand();

  if(verboselevel == 1)
  G4cout<<" Apply:  Nucleus Momentum "<<outMomNcl<<" Cos "<<cosNucl<<G4endl;
// ------ The transformation from Hadron frame to absolute one -------
    G4double   RotMatrix[3][3];
    G4double   NewMomHdrSys[3];
    G4double   NewMomOldSys[3];
    G4double   NucMomHdrSys[3];
    G4double   NucMomOldSys[3];
//   ------------ The hadron angles in its own system ------------
    G4double   sinHadr   = std::sqrt(1-cosHadr*cosHadr);
    G4double   sinFiHadr = std::sin(ranFi);
    G4double   cosFiHadr = std::cos(ranFi);
//   ---- The unit vector of a hadron momentum in its own system ---
    NewMomHdrSys[0]    = sinHadr*cosFiHadr; 
    NewMomHdrSys[1]    = sinHadr*sinFiHadr;
    NewMomHdrSys[2]    = cosHadr;          
    G4double  sinNucl  = std::sqrt(1-cosNucl*cosNucl);
//   -----------  The unit vector of  a nucleus momentum  -----------
    NucMomHdrSys[0] = -sinNucl*cosFiHadr;
    NucMomHdrSys[1] = -sinNucl*sinFiHadr;
    NucMomHdrSys[2] =  cosNucl;
//  __________ The angles of unit vector in absolute system ----------
    G4double    rhoInMom = aHadron.Get4Momentum().vect().mag()/1000;
    G4double   cosTet  = aHadron.Get4Momentum().vect().z()/1000/rhoInMom;
    G4double   sinTet;
    G4double   sinFi;
    G4double   cosFi;
    if(cosTet>0.99999 || cosTet<-0.99999) 
    {
     sinTet = 0;
     sinFi  = 0;
     cosFi  = 1;
    }
    else
    {
     sinTet  = std::sqrt(1-cosTet*cosTet);
     sinFi   = aHadron.Get4Momentum().vect().y()/1000/sinTet/rhoInMom;
     cosFi   = aHadron.Get4Momentum().vect().x()/1000/sinTet/rhoInMom;
    }
//  -----------------  The rotation matrix  ---------------
     RotMatrix[0][0]   =  cosFi*cosTet;
     RotMatrix[0][1]   =  sinFi;
     RotMatrix[0][2]   = -cosFi*sinTet;
     RotMatrix[1][0]   = -sinFi*cosTet;
     RotMatrix[1][1]   =  cosFi;
     RotMatrix[1][2]   =  sinFi*sinTet;
     RotMatrix[2][0]   =  sinTet;
     RotMatrix[2][1]   =  0;
     RotMatrix[2][2]   =  cosTet;

  if(verboselevel==1)
  {
  G4cout<<" cosTet sinTet cosFi sinFi "<<cosTet<<"  "<<sinTet<<"  "
        <<cosFi<<"  "<<sinFi<<"  Rho "<<rhoInMom<<G4endl;
  G4cout<<" Rotation: "<<RotMatrix[0][0]<<"  "<<RotMatrix[0][1]<<"  "
        <<RotMatrix[0][2]<<G4endl;
  G4cout<<" Rotation: "<<RotMatrix[1][0]<<"  "<<RotMatrix[1][1]<<"  "
        <<RotMatrix[1][2]<<G4endl;
  G4cout<<" Rotation: "<<RotMatrix[2][0]<<"  "<<RotMatrix[2][1]<<"  "
        <<RotMatrix[2][2]<<G4endl;
  }
// --------------  Rotation of momenta  ----------------
     for(G4int ii=0; ii<3; ii++)
     {
       NewMomOldSys[ii] = 0;                        
       NucMomOldSys[ii] = 0;                        
       for(G4int ll=0; ll<3; ll++)
       {
//  --------------  Rotation of Hadron  ------------------
         NewMomOldSys[ii] = NewMomOldSys[ii]+
                  RotMatrix[ii][ll]*NewMomHdrSys[ll];
//  --------------  Rotation of Nucleus ------------------
         NucMomOldSys[ii] = NucMomOldSys[ii]+
                  RotMatrix[ii][ll]*NucMomHdrSys[ll];
       }
     }
//  --------------  The Nucleus Momentum  -------------
         G4ThreeVector  aNuclMom( NucMomOldSys[0]*outMomNcl*1000,
                                  NucMomOldSys[1]*outMomNcl*1000,
                                  NucMomOldSys[2]*outMomNcl*1000);

         secNuclDyn->SetMomentum(aNuclMom);

   if(verboselevel==1)
   G4cout<<" Hucleus: x y z "<<NucMomOldSys[0]*outMomNcl<<"  "
                        <<NucMomOldSys[1]*outMomNcl<<"  "
                        <<NucMomOldSys[2]*outMomNcl<<G4endl;
//  --------------  The Direction of Hadron Momentum  ---------------
         G4double pxnew = NewMomOldSys[0];  
         G4double pynew = NewMomOldSys[1];
         G4double pznew = NewMomOldSys[2];
//   ---------------------------------------------
     G4int NumbPart;
     NumbPart = theParticleChange.GetNumberOfSecondaries();

   if(verboselevel==1)
   G4cout<<"  Apply (end1): Q2 "<<ranQ2<<" NumbPart "<<NumbPart <<G4endl;

     theParticleChange.AddSecondary(secNuclDyn);

     theParticleChange.SetEnergyChange(outEnHadr*1000);
     theParticleChange.SetMomentumChange(pxnew, pynew, pznew);

     NumbPart = theParticleChange.GetNumberOfSecondaries();

   if(verboselevel==1)
   G4cout<<"  Apply (end2): Q2 "<<ranQ2<<" NumbPart "<<NumbPart <<G4endl;

//  ----------------------------------------------------
   if(verboselevel==1)
   G4cout<<G4endl<<"----------- End Applay ------------"<<G4endl;

  return &theParticleChange;
 }
*/
//  ########################################################
G4HadFinalState * G4ElasticHadrNucleusHE::ApplyYourself(
                          const  G4HadProjectile  &aHadron,
                                 G4Nucleus        &aNucleus)
{
  if(verboselevel > -1)
  {
   G4cout<<G4endl<<"----------- Begin Applay ------------"<<G4endl;
   G4cout<<" Total Energy "<<aHadron.GetTotalEnergy()<<G4endl;
  }

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
  else secNuclDef = G4ParticleTable::
                            GetParticleTable()->FindIon(nZ,nA,0,nZ);
//  ----------------------------------------------------
 if(verboselevel > -1)
 G4cout<<" Hadr: TotMom  "<<aHadron.GetTotalMomentum()
       <<" Energy  "<<aHadron.GetTotalEnergy()
       <<"  A  "<<A<<G4endl;
// ----------------_ Randomization of Q2 as dSig/dt --------------

    G4double ranQ2;

    if(A>1)  ranQ2 = SampleT(hadrDef, dParticle->GetTotalMomentum()/1000, 
                nA-nZ, nZ);
    else ranQ2 = HadronProtonQ2(aHadron.GetDefinition(), 
                  dParticle->GetTotalMomentum()/1000);

    ranQ2 = ranQ2*1000*1000;

    Q2res = ranQ2;

  if(verboselevel > -1)
  G4cout<<" Applay: Q2 = "<<ranQ2<<G4endl;
//  ----------------  Hadron kinematics ----------------
  G4double m1 = hadrDef->GetPDGMass();
  G4double m2 = secNuclDef->GetPDGMass();
  G4LorentzVector lv1 = dParticle->Get4Momentum();
  G4LorentzVector lv0(0.0,0.0,0.0,m2);
  G4LorentzVector lv1test = dParticle->Get4Momentum();
  G4LorentzVector lv0test(0.0,0.0,0.0,m2);

  if(verboselevel > -1)
  G4cout << "lv0= " << lv0 << "  lv1= " << lv1 << "   m1= "
         << m1 << " ranQ2= " << ranQ2<<G4endl;

  G4LorentzVector lv  = lv0 + lv1;
  G4ThreeVector bst = lv.boostVector();
  lv1.boost(-bst);
  lv0.boost(-bst);
  G4ThreeVector p1 = lv1.vect();
  G4double ptot = p1.mag();

  if (verboselevel > -1)
  G4cout<<" p: H, nucl "<<p1<<"  "<<lv0.vect()<<G4endl;  

  // Sampling in CM system
  G4double phi  = G4UniformRand()*twopi;
  G4double cost = 1. - 0.5*ranQ2/(ptot*ptot);

  if(std::abs(cost) > 1.0) cost = 2.0*G4UniformRand() - 1.0;
  G4double sint = std::sqrt((1.0-cost)*(1.0+cost));

  if (verboselevel > -1)
  G4cout << "cos(t)=" << cost << " std::sin(t)=" << sint << G4endl;

  G4ThreeVector v1(sint*std::cos(phi),sint*std::sin(phi),cost);
  p1 = p1.unit();
  v1.rotateUz(p1);
  v1 *= ptot;
  G4LorentzVector nlv1(v1.x(),v1.y(),v1.z(),
                            std::sqrt(ptot*ptot + m1*m1));
  G4LorentzVector nlv0 = lv0 + lv1 - nlv1;
  nlv0.boost(bst);
  nlv1.boost(bst);


  lv1 = nlv1-lv1test;
  lv0 = nlv0-lv0test;

  if(verboselevel > -1)
  G4cout<<" Q2 directly:  1 "<<lv1.m2()<<"  0 "<<lv0.m2()<<G4endl; 

  Q2res = -lv0.m2();

  G4double eFinal = nlv1.e() - m1;

  if (verboselevel > -1)
    G4cout << " P0= "<< nlv0 << "   P1= "
           << nlv1<<" m= " << m1 << " ekin0= " << eFinal
           << " ekin1= " << nlv0.e() - m2
           <<G4endl;
//  if(eFinal < 0.0) 
  if(verboselevel > -1)
  {
    G4cout << "G4HadronElastic WARNING ekin= " << eFinal
           << " after scattering of "
           << hadrDef->GetParticleName()
           << " p(GeV/c)= " << dParticle->GetTotalMomentum()/GeV
           << " on " << secNuclDef->GetParticleName()
           << G4endl;
///    eFinal = 0.0;
  }

  if(verboselevel > -1)
  {
  G4cout<<" nlv0 "<<nlv0<<" eFinal "<<eFinal<<"  m1 "<<m1<<G4endl;
  G4cout<<" nlv1 "<<nlv1<<" eFinal "<<eFinal<<"  m1 "<<m1<<G4endl;
  }

  theParticleChange.SetMomentumChange(nlv1.vect().unit());
  theParticleChange.SetEnergyChange(eFinal);

  G4DynamicParticle * aSec = new G4DynamicParticle(secNuclDef, nlv0);
  theParticleChange.AddSecondary(aSec);

  if(verboselevel > -1)
  {
  G4cout<<" nlv0 "<<nlv0<<" eFinal "<<eFinal<<"  m1 "<<m1<<G4endl;
  G4cout<<" nlv1 "<<nlv1<<" eFinal "<<eFinal<<"  m1 "<<m1<<G4endl;
  }

   if(verboselevel > -1)
   G4cout<<G4endl<<"----------- End Applay ------------"<<G4endl;

  return &theParticleChange;
}
//  ########################################################
  G4double G4ElasticHadrNucleusHE::
               HadronNucleusQ2_2(const G4ParticleDefinition * aHadron,
                                 G4int                  AWeight,
                                 G4double               LabMom,
                                 G4double               Rand,
                                 G4int                  Step,
                                 ElasticData          * pElD)
  {
    G4int        nucN, Step1=0, ii;

    G4int      kk=0, NumbOnE, iNumbQ2;
    G4double * dNumbQ2, * dNumbFQ2;
    G4double   Q2=0.0, Buf=0.0, Energy;

    Nstep = pElD->ONQ2;

    iContr = 2;

    G4String  hadrName = aHadron->GetParticleName();
    G4double  hadrMass = aHadron->GetPDGMass()/1000;

 if(verboselevel == 1)
 G4cout<<" Q2_2: hName "<<hadrName<<" Weight "<<AWeight<<" Mom "
      <<LabMom<<"  Mass "<<hadrMass<<" Rand "<<Rand<<G4endl; 

    Energy = std::sqrt(hadrMass*hadrMass+
                         LabMom*LabMom)-hadrMass;

    G4ParticleDefinition * pHadron = const_cast<G4ParticleDefinition *>
              (aHadron);
    G4DynamicParticle    * dHadron = new G4DynamicParticle::
                               G4DynamicParticle();

    G4ThreeVector HadrMomentum(0.0, 0.0, LabMom);
    
    dHadron->SetDefinition(pHadron);
    dHadron->SetMomentum(HadrMomentum);

    NumbOnE = pElD->GetNumberE(Energy);

 if(verboselevel == 1)
 G4cout<<"  Q2_2: Energy  NumbOnE "<<Energy<<"  "<<NumbOnE<<G4endl;

    nucN    = AWeight;
    dNumbQ2 = pElD->TableQ2;
    if(pElD->TableCrossSec[NumbOnE*pElD->ONQ2+2]<0.00001) Step=0;

 if(verboselevel == 1)
 G4cout<<" Q2_2. Before Step: Step = "<<Step<<"  NumbE "<<NumbOnE<<"  "
       <<*(dNumbQ2+3)<<"  "<<pElD->TableQ2[3]<<G4endl; 
//  .............................................
    if(Step==0)
    {
     GetHadronValues(aHadron, LabMom);
     Q2 = pElD->maxQ2;

 if(verboselevel == 1)
 G4cout<<" Q2_2. Step=0 (1):  R1 ElD->R1 "<<R1<<"  "<<pElD->R1
       <<"  Q2  "<<Q2<<G4endl;

     Weight = GetLightFq2(nucN, Q2, Step1);

     pElD->TableCrossSec[NumbOnE*pElD->ONQ2] = 0;

 if(verboselevel == 1)
 G4cout<<" Q2_2. Step=0 (1): maxQ2 "<<pElD->maxQ2<<"  dQ2 "<<pElD->dQ2
       <<" N "<<Nstep<<G4endl;
 
 // ------------------------------------------
        for(ii=1; ii<Nstep/*3*/; ii++)
        {
          Q2      = pElD->TableQ2[ii];

          Buf = GetLightFq2(nucN, Q2, Step1);
          pElD->TableCrossSec[NumbOnE*pElD->ONQ2+ii] = Buf/Weight;

 if(verboselevel == 1)
 G4cout<<"  Q2_2: Step 0: "<<ii<<"  Q2 "<<"   "<<Q2
       <<"  Buf  "<<Buf<<" FQ2 "<<Buf/Weight<<G4endl;
        }   // for ii

      RandMax  = Buf/Weight;
      CurrentN = 3;
    }  //  if Step
//  ......................................
 if(verboselevel == 1)
 G4cout<<"  Q2_2: After Step 0: RandMax Rand "<<RandMax<<"  "<<Rand
       <<"   CurrentN  "<<CurrentN<<G4endl;
 
    ii = 0;
    while(pElD->TableCrossSec[NumbOnE*pElD->ONQ2+ii]>0.0001) 
    {
       RandMax  = pElD->TableCrossSec[NumbOnE*pElD->ONQ2+ii];
       CurrentN = ii;
       ii++;
    }

    dNumbFQ2 = &pElD->TableCrossSec[NumbOnE*pElD->ONQ2];

 if(verboselevel == 1)
 G4cout<<"  Q2_2:  RandMax CurrentN  "<<RandMax<<"  "<<CurrentN<<" ii "<<ii<<G4endl;

    if(RandMax > Rand)
    {
     kk = 0;
     while(Rand>pElD->TableCrossSec[NumbOnE*pElD->ONQ2+kk])
     {
 if(verboselevel == 1)
 G4cout<<"  Q2_2: while Rand: Rand kk "<<Rand<<"  "<<kk<<" TableCrSec "
       <<pElD->TableCrossSec[NumbOnE*pElD->ONQ2+kk]<<G4endl;     
       kk++;

     }
     iNumbQ2 = kk;

 if(verboselevel == 1)
 G4cout<<"  Q2_2: Before GetQ2 kk iNumbQ2 "<<kk<<"  "<<iNumbQ2
       <<" dNumbFQ2 Rand "<<*dNumbFQ2<<"  "<<Rand<<G4endl;

     Q2 = GetQ2_2(iNumbQ2, dNumbQ2, dNumbFQ2, Rand);

 if(verboselevel == 1)
 G4cout<<"  Q2_2:  Max>Rand: Rand "<<Rand<<" Max "<<RandMax<<"  kk  "<<kk
       <<"  Q2  "<<Q2<<G4endl;        
     return Q2;
    }
    else
    {
      Step = 1;
      for(ii=CurrentN; ii<Nstep; ii++)
      { 

        Q2  = pElD->TableQ2[ii];  //(ii-0)*dQ2;
	Buf = GetLightFq2(nucN, Q2, Step);
        pElD->TableCrossSec[NumbOnE*pElD->ONQ2+ii] = Buf/Weight;

 if(verboselevel == 1)
 G4cout<<ii<<"  Q2_2:  Further:nucN  Q2 "<<nucN<<"   "<<Q2
       <<" Buf "<<Buf<<"  FQ2 "<<Buf/Weight
       <<"  Rand "<<Rand<<"  Max "<<RandMax<<G4endl;
        if(Buf/Weight>Rand) 
	{
          iNumbQ2 = ii;

 if(verboselevel == 1)
 G4cout<<"  Q2_2:  Befor GetQ2 ii iNumbQ2 "<<kk<<"  "<<iNumbQ2
       <<" dNumbFQ2 Rand "<<dNumbFQ2<<"  "<<Rand<<G4endl;

          Q2 = GetQ2_2(iNumbQ2, dNumbQ2, dNumbFQ2, Rand);
          CurrentN = ii;
          RandMax  = Buf/Weight;

 if(verboselevel == 1)
 G4cout<<"  Q2_2:  if Rand  ii  "<<ii<<" Buf "<<Buf/Weight<<" Rand "<<Rand
       <<"  Max "<<RandMax<<"  Q2  "<<Q2<<G4endl;

          break;
	}  //  if Rand

      }    //  for ii

    }      //  else

 return  Q2;
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

   if(abs(D0) < 0.00000001 || kk == 2)
   { 
    ranQ2 = (*(Q+kk-1)+(ranUni-*(F+kk-1))*
                 (*(Q+kk)-*(Q+kk-1))/(*(F+kk)-*(F+kk-1)));
    }
   else    
    {
      G4double DA = X1*F2+X3*F1+X2*F3-X3*F2-X1*F3-X2*F1;
      G4double DB = X2*F12+X1*F32+X3*F22-X2*F32-X3*F12-X1*F22;
      G4double DC = X3*F2*F12+X2*F1*F32+X1*F3*F22
                          -X1*F2*F32-X2*F3*F12-X3*F1*F22;
      ranQ2 = (DA*ranUni*ranUni+DB*ranUni+DC)/D0;
    }

  if(verboselevel == 1)
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
        GetLightFq2(G4int Nucleus, G4double Q2, G4int)
  {

// ----------------  The preparing of probability function  ------------

    G4double prec = Nucleus > 208  ?  1.0e-7 : 1.0e-6;

    G4double    Stot     = HadrTot*MbToGeV2;     //  Gev^-2
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

    G4double    Norm     = (R12*R1-Pnucl*R22*R2); //*HP->Aeff;
///      G4double    NormP    = R12*R1-PnuclP*R22*R2;
    G4double    R13      = R12*R1/R12B;
    G4double    R23      = Pnucl*R22*R2/R22B;
    G4double    Unucl    = Stot/2/pi/Norm*R13;
////      G4double    Unclprod = Stot/2/Pi1/NormP*R13Ap;
    G4double    FiH      = std::asin(HadrReIm/Rho2);
    G4double    NN2      = R23/R13;

    G4double    dddd;
 
 if(verboselevel == 3)
 {
  G4cout<<" HadrPar: tot Slope Im "<<Stot<<"  "<<Bhad
        <<"  "<<Asq<<G4endl;
  G4cout<<"  Light 0: R1 R2 FiH "<<R1<<"  "<<R2
        <<"  "<<FiH<<" Unucl "<<Unucl<<G4endl;
  G4cout<<"  Light 0: R13 R23 norm "<<R13<<"  "<<R23
        <<"  "<<Norm<<" Unucl "<<Unucl<<G4endl;
  G4cout<<"  Light 0: Asq "<<Asq<<G4endl;
 }
///      G4double    DDSec1p  = (DDSect2+
///                    DDSect3*std::log(1.06*2*Ehad/R1/4));

///      G4double    DDSec2p  = (DDSect2+
///                    DDSect3*std::log(1.06*2*Ehad/
//                     std::sqrt((R12+R22)/2)/4));

///      G4double    DDSec3p  = (DDSect2+
///                    DDSect3*std::log(1.06*2*Ehad/R2/4));

///      G4double    R12ApdR22Ap = 0.5*(R12Apd+R22Apd);

//                 iIntgr[0] = 0;

    G4double Prod0    = 0;
    G4double N1       = -1;
    G4double Tot0     = 0;
    G4double exp1;

    G4double Prod3 ;
    G4double exp2  ;
    G4double N4, N5, N2, Prod1, Prod2;
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
                            (1-std::exp(-Q2*dddd/*(1/exp1+1/exp2)*//4))/
		            dddd/*(1/exp1+1/exp2)*/*4*SetBinom[i1][m1];
               }                                   // m1
               Prod2 = Prod2 +Prod3*N5*SetBinom[i2][m2];

 if(verboselevel == 2)
 {
 G4cout<<"  Binom i2 m2 "<<m2<<"  "<<i2<<"  "<<SetBinom[i2][m2]<<G4endl;
 G4cout<<"m2 Prod2 "<<m2<<"  "<<Prod2<<"::";
 }

           }                                      // m2
         Prod1 = Prod1 + Prod2*N2*std::cos(FiH*(i1-i2));
//         Tot0  = Tot0  + Prod2*N2*std::sin(FiH*(i1-i2));

////         if (std::fabs(Prod2*N2/Prod1)<1e-6) break;
         if (std::fabs(Prod2*N2/Prod1)<prec) break;
        }                                         // i2
//        ImDistr = Tot0  + Tot0*N1;
         Prod0   = Prod0 + Prod1*N1;
////         if(std::fabs(N1*Prod1/Prod0) < 1e-6) break;
         if(std::fabs(N1*Prod1/Prod0) < prec) break;
      }                                           // i1
       Prod0        = Prod0*pi/MbToGeV2/4;  //  This is in mb

 if(verboselevel == 2)
 G4cout<<"  Q2  "<<Q2<<" Prod0 "<<Prod0<<G4endl;
   return Prod0;
  }
//  #########################################################
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
//  =====================================================
   void  G4ElasticHadrNucleusHE::
       GetKinematics(const G4ParticleDefinition * aHadron,
                           G4double MomentumH)
   {
     HadronName = aHadron->GetParticleName();

     GetHadronValues(aHadron, MomentumH);

     SigTot = HadrTot;
     ReOnIm = HadrReIm;

     if(Slope==0) Slope  = HadrSlope;

     IntConst = (1-Coeff1)/Slope;

     ProtonM   = 0.93827;

     HadronM   = aHadron->GetPDGMass()/1000.0;      // GeV
     HadronE   = std::sqrt(MomentumH*MomentumH+
                           HadronM*HadronM);
     HdrE      = HadronE;

     PM2    = ProtonM*ProtonM;
     HM2    = HadronM*HadronM;
     Sh     = 2.0*ProtonM*HadronE+PM2+HM2;    // GeV
     SqrtS  = std::sqrt(Sh);
     ConstU = 2*PM2+2*HM2-Sh;

     EcmH   = (Sh+HM2-PM2)/2/SqrtS;
     EcmP   = (Sh-HM2+PM2)/2/SqrtS;

     Kcm    = sqrt(EcmH*EcmH-HM2);
     MaxT   = 4*Kcm*Kcm;

     BoundaryP[0]=9.0; BoundaryTG[0]=5.0;BoundaryTL[0]=MaxT/2.0;
     BoundaryP[1]=20.0;BoundaryTG[1]=1.5;BoundaryTL[1]=MaxT;
     BoundaryP[2]=5.0; BoundaryTG[2]=1.0;BoundaryTL[2]=1.5;
     BoundaryP[3]=8.0; BoundaryTG[3]=3.0;BoundaryTL[3]=MaxT;
     BoundaryP[4]=7.0; BoundaryTG[4]=3.0;BoundaryTL[4]=MaxT;
     BoundaryP[5]=5.0; BoundaryTG[5]=2.0;BoundaryTL[5]=MaxT;
     BoundaryP[6]=5.0; BoundaryTG[6]=1.5;BoundaryTL[6]=3.0;

     HadrCodes[0] =  2212; HadrCodes[1] =  2112;
     HadrCodes[2] = -2212; HadrCodes[3] =  211;
     HadrCodes[4] = -211;  HadrCodes[5] =  321;
     HadrCodes[6] = -321;

     HadrCode = aHadron->GetPDGEncoding();

     G4int NumberH=0;

     while(HadrCode!=HadrCodes[NumberH]) NumberH++;

     if(MomentumH<BoundaryP[NumberH]) MaxTR = BoundaryTL[NumberH];
     else MaxTR = BoundaryTG[NumberH];

    if(verboselevel == 1)
  G4cout<<" GetKinem. Code MaxT "<<HadrCode<<"  "<<MaxT
        <<"  MaxTR "<<MaxTR<<G4endl;

     GetParametersHP(aHadron, MomentumH);
   }
//  +++++++++++++++++++++++++++++++++++++++
//  ++++++++++++++++++++++++++++++++++++++++
  void  G4ElasticHadrNucleusHE::
       GetParametersHP(const G4ParticleDefinition * aHadron,
                             G4double HadrP)
  {
    G4int  iStep;

    HadrCode = aHadron->GetPDGEncoding();

    switch(HadrCode)
    {
      case 2212 :
      case 3122 :
      case 3222 :
      case 3112 :
      case 3212 :
      case 3312 :
      case 3322 :
      case 3334 :
        {
          G4double EnP[6]={1.5,3.0,5.0,9.0,14.0,19.0};
          G4double C0P[6]={0.15,0.02,0.06,0.08,0.0003,0.0002};
          G4double C1P[6]={0.05,0.02,0.03,0.025,0.0,0.0};
          G4double B0P[6]={1.5,2.5,3.0,4.5,1.4,1.25};
          G4double B1P[6]={5.0,1.0,3.5,4.0,4.8,4.8};

          if(BoundaryP[0]<HadrP)
          {
            Coeff0=Coeff1=Coeff2=0.0;
            Slope0=Slope1=Slope2=1.0;
          }
          else
          {
            iStep = 0;  while(HadrP>EnP[iStep])  iStep++;

            Coeff0 = LineInterpol(EnP[iStep], EnP[iStep-1],
                                  C0P[iStep], C0P[iStep-1],
                                  HadrP);

            Coeff1 = LineInterpol(EnP[iStep], EnP[iStep-1],
                                  C1P[iStep], C1P[iStep-1],
                                  HadrP);

            Slope0 = LineInterpol(EnP[iStep], EnP[iStep-1],
                                  B0P[iStep], B0P[iStep-1],
                                  HadrP);

            Slope1 = LineInterpol(EnP[iStep], EnP[iStep-1],
                                  B1P[iStep], B1P[iStep-1],
                                  HadrP);
          }       //  else
         Coeff2 = 0;
         Slope2 = 5;
        } break;  //  case 2212
      case  2112 :
      case -2112 :
        {
          G4double EnN[5]={1.5,5.0,10.0,14.0,20.0};
          G4double C0N[5]={0.0,0.0,0.02,0.02,0.01};
          G4double C1N[5]={0.06,0.008,0.0015,0.001,0.0003};
          G4double B0N[5]={1.5,2.5,3.8,3.8,3.5};
          G4double B1N[5]={1.5,2.2,3.6,4.5,4.8};

          if(BoundaryP[1]<HadrP)
          {
            Coeff0=Coeff1=0.0;
            Slope0=Slope1=Slope2=1.0;
          }
          else
          {
            iStep = 0;  while(HadrP>EnN[iStep])   iStep++;

            Coeff0 = LineInterpol(EnN[iStep], EnN[iStep-1],
                                  C0N[iStep], C0N[iStep-1],
                                  HadrP);

            Coeff1 = LineInterpol(EnN[iStep], EnN[iStep-1],
                                  C1N[iStep], C1N[iStep-1],
                                  HadrP);

            Slope0 = LineInterpol(EnN[iStep], EnN[iStep-1],
                                  B0N[iStep], B0N[iStep-1],
                                  HadrP);

            Slope1 = LineInterpol(EnN[iStep], EnN[iStep-1],
                                  B1N[iStep], B1N[iStep-1],
                                  HadrP);
            Coeff2 = 0.8/HadrP/HadrP;
            Slope2 = 5;
          }       //  else
        } break;  //  case 2112

      case -2212 :
      case -3122 :
      case -3222 :
      case -3112 :
      case -3212 :
      case -3312 :
      case -3322 :
      case -3334 :
        {
          G4double EnP[2]={1.5,4.0};
          G4double C0P[2]={0.001,0.0005};
          G4double C1P[2]={0.003,0.001};
          G4double B0P[2]={2.5,4.5};
          G4double B1P[2]={1.0,4.0};
          if(BoundaryP[2]<HadrP)
          {
            Coeff0=Coeff1=Coeff2=0.0;
            Slope0=Slope1=Slope2=1.0;
          }

          else
          {
            iStep = 0;  while(HadrP>EnP[iStep])  iStep++;

            Coeff0 = LineInterpol(EnP[iStep], EnP[iStep-1],
                                  C0P[iStep], C0P[iStep-1],
                                  HadrP);

            Coeff1 = LineInterpol(EnP[iStep], EnP[iStep-1],
                                  C1P[iStep], C1P[iStep-1],
                                  HadrP);

            Slope0 = LineInterpol(EnP[iStep], EnP[iStep-1],
                                  B0P[iStep], B0P[iStep-1],
                                  HadrP);

            Slope1 = LineInterpol(EnP[iStep], EnP[iStep-1],
                                  B1P[iStep], B1P[iStep-1],
                                  HadrP);
          }       //  else

            Coeff2=0.0;
            Slope2=5.0;
        } break;  //  case -2212

      case 211 :
        {
          G4double EnPP[4]={1.0,2.0,3.0,4.0};
          G4double C0PP[4]={0.0,0.0,0.0,0.0};
          G4double C1PP[4]={0.15,0.08,0.02,0.01};
          G4double B0PP[4]={1.5,2.8,3.8,3.8};
          G4double B1PP[4]={0.8,1.6,3.6,4.6};

          if(BoundaryP[3]<HadrP)
          {
            Coeff0=Coeff1=Coeff2=0.0;
            Slope0=Slope1=Slope2=1.0;
          }
          else
          {
            iStep = 0;  while(HadrP>EnPP[iStep])   iStep++;

            Coeff0 = LineInterpol(EnPP[iStep], EnPP[iStep-1],
                                  C0PP[iStep], C0PP[iStep-1],
                                  HadrP);

            Coeff1 = LineInterpol(EnPP[iStep], EnPP[iStep-1],
                                  C1PP[iStep], C1PP[iStep-1],
                                  HadrP);

            Slope0 = LineInterpol(EnPP[iStep], EnPP[iStep-1],
                                  B0PP[iStep], B0PP[iStep-1],
                                  HadrP);

            Slope1 = LineInterpol(EnPP[iStep], EnPP[iStep-1],
                                  B1PP[iStep], B1PP[iStep-1],
                                  HadrP);
            Coeff2 = 0.02/HadrP;
            Slope2 = 5;
          }       //  else
        } break;  //  case 211

      case -211 :
        {
          G4double EnPP[4]={1.0,2.0,3.0,4.0};
          G4double C0PP[4]={0.0,0.0,0.0,0.0};
//          G4double C1PP[4]={0.15,0.08,0.02,0.01};
//          G4double B0PP[4]={1.5,2.8,3.8,3.8};
          G4double B1PP[4]={0.8,1.6,3.6,4.6};

          if(BoundaryP[4]<HadrP)
          {
            Coeff0=Coeff1=Coeff2=0.0;
            Slope0=Slope1=Slope2=1.0;
          }
          else
          {
            iStep = 0;  while(HadrP>EnPP[iStep])    iStep++;

            Coeff0 = LineInterpol(EnPP[iStep], EnPP[iStep-1],
                                  C0PP[iStep], C0PP[iStep-1],
                                  HadrP);
            Slope1 = LineInterpol(EnPP[iStep], EnPP[iStep-1],
                                  B1PP[iStep], B1PP[iStep-1],
                                  HadrP);
            Coeff2 = 0.02/HadrP;
            Slope2 = 5;
          }       //  else
        } break;  //  case -211

      case 321 :
        {
          G4double EnK[4]={1.4,2.33,3.0,5.0};
          G4double C0K[4]={0.0,0.0,0.0,0.0};
          G4double C1K[4]={0.01,0.007,0.005,0.003};
          G4double B0K[4]={1.5,2.0,3.8,3.8};
          G4double B1K[4]={1.6,1.6,1.6,1.6};

          if(BoundaryP[5]<HadrP)
          {
            Coeff0=Coeff1=Coeff2=0.0;
            Slope0=Slope1=Slope2=1.0;
          }
          else
          {
            iStep = 0;  while(HadrP>EnK[iStep])   iStep++;

            Coeff0 = LineInterpol(EnK[iStep], EnK[iStep-1],
                                  C0K[iStep], C0K[iStep-1],
                                  HadrP);

            Coeff1 = LineInterpol(EnK[iStep], EnK[iStep-1],
                                  C1K[iStep], C1K[iStep-1],
                                  HadrP);

            Slope0 = LineInterpol(EnK[iStep], EnK[iStep-1],
                                  B0K[iStep], B0K[iStep-1],
                                  HadrP);

            Slope1 = LineInterpol(EnK[iStep], EnK[iStep-1],
                                  B1K[iStep], B1K[iStep-1],
                                  HadrP);
            Coeff2 = 0.34/HadrP/HadrP/HadrP;
            Slope2 = 5;
          }       //  else
        } break;  //  case 321
      case -321 :
        {
          G4double EnKM[2]={1.4,4.0};
          G4double C0KM[2]={0.006,0.002};
          G4double C1KM[2]={0.00,0.00};
          G4double B0KM[2]={2.5,3.5};
          G4double B1KM[2]={1.6,1.6};

          if(BoundaryP[6]<HadrP)
          {
            Coeff0=Coeff1=Coeff2=0.0;
            Slope0=Slope1=Slope2=0.0;
          }
          else
          {
            iStep = 0;  while(HadrP>EnKM[iStep])   iStep++;

            Coeff0 = LineInterpol(EnKM[iStep], EnKM[iStep-1],
                                  C0KM[iStep], C0KM[iStep-1],
                                  HadrP);

            Coeff1 = LineInterpol(EnKM[iStep], EnKM[iStep-1],
                                  C1KM[iStep], C1KM[iStep-1],
                                  HadrP);

            Slope0 = LineInterpol(EnKM[iStep], EnKM[iStep-1],
                                  B0KM[iStep], B0KM[iStep-1],
                                  HadrP);

            Slope1 = LineInterpol(EnKM[iStep], EnKM[iStep-1],
                                  B1KM[iStep], B1KM[iStep-1],
                                  HadrP);
            Coeff2 = 0.01/HadrP/HadrP/HadrP;
            Slope2 = 5;
          }       //  else
        } break;  //  case -321
    }
  }
//  ++++++++++++++++++++++++++++++++++++++++
  G4double G4ElasticHadrNucleusHE::
               LineInterpol(G4double p1, G4double p2,
                            G4double c1, G4double c2,
                            G4double p)
  {
    G4double c;

    c = c1+(p-p1)/(p2-p1)*(c2-c1);
    return c;
  }
//  ============================================================
  G4float G4ElasticHadrNucleusHE::GetFt(G4double Q2)
  {
    G4float Fdistr, SqrQ2 = std::sqrt(Q2);

    Fdistr = (1-Coeff1-Coeff0/*-0.0*Coeff2*std::exp(ConstU)*/)
                   /Slope*(1-std::exp(-Slope*Q2))

                  +Coeff0*(1-std::exp(-Slope0*Q2))

                  +Coeff2/Slope2*std::exp(Slope2*ConstU)*
                      (std::exp(Slope2*Q2)-1)

                  +2*Coeff1/Slope1*(1/Slope1-(1/Slope1+SqrQ2)*
                  std::exp(-Slope1*SqrQ2))
                   ;

  if(verboselevel == 1)
  {
  G4cout<<" GetFt: Fdistr "<<Fdistr<<G4endl;
  G4cout<<" Coef 0 1 2 "<<Coeff0<<"  "<<Coeff1<<"  "
        <<Coeff2<<"  Q2 "<<Q2<<"  ConstU "<<ConstU<<G4endl;
  G4cout<<" Slope Slope0 1 2 "<<Slope<<"  "<<Slope0<<"   "
        <<Slope1<<"  "<<Slope2<<G4endl;
  }
    return Fdistr;
  }
//  +++++++++++++++++++++++++++++++++++++++
  G4float G4ElasticHadrNucleusHE::GetDistrFun(G4double Q2)
  {
    G4float Fq2;

    Fq2   = GetFt(Q2);

  if(verboselevel == 1)
  G4cout<<" GetDistrFun: Fq2  "<<Fq2<<G4endl;

    return Fq2/FmaxT;
  }
//  +++++++++++++++++++++++++++++++++++++++
  G4double G4ElasticHadrNucleusHE::
              GetQ2(G4double Ran)
  {
    G4double DDD0=MaxTR*0.5, DDD1=0.0, DDD2=MaxTR, delta;
    G4double Q2;

    FmaxT = GetFt(MaxTR);
    delta = GetDistrFun(DDD0)-Ran;

  if(verboselevel == 1)
  G4cout<<" GetQ2 Max TR FmaxT  "<<MaxTR<<"  "<<FmaxT
        <<"  delta  "<<delta<<G4endl;

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

    Q2 = DDD0;

    return Q2;
  }
//  ++++++++++++++++++++++++++++++++++++++++++
  G4double G4ElasticHadrNucleusHE::
       HadronProtonQ2(const G4ParticleDefinition * aHadron,
                            G4double inLabMom)
  {
    G4double Ran = G4UniformRand(), Q2;

    GetKinematics(aHadron, inLabMom);

    Q2 = GetQ2(Ran);

  if(verboselevel == 1)
  G4cout<<"  HRQ2HE Ran "<<Ran<<" q2 "<<Q2<<G4endl;

    return Q2;
  }

//  ===========================================
//  ##########################################################
////  void G4ElasticHadrNucleusHE::GetNucleusParameters(G4int Nucleus)
  void ElasticData::GetNucleusParameters(G4int Nucleus)
  {

    if(Nucleus == 208)
    {
//      R1 = 20.73; R2 = 15.74.
      R1       = 4.1408*std::pow(static_cast<double>(Nucleus),0.3018);
      R2       = 3.806*std::pow(Nucleus-10.068,0.2685);
      Pnucl    = 0.9;
      Aeff     = 1.1;
      R1       = 19.5;
      R1       = 20.5;    //  26.09.06
      R2       = 15.74;
      Pnucl    = 0.4;
      Aeff     = 0.7;
    }

    else if(Nucleus == 90)
    {
      R1    = 16.50;
      R2    = 11.62;
      Pnucl = 0.4;
      Aeff  = 0.9;
      R1    = 16.5;
      R2    = 11.62;
      Pnucl = 0.4;
      Aeff  = 0.7;
      R1    = R1*1.1;
    }

    else if(Nucleus == 58)
    {
      R1    = 15.0;
      R2    = 9.9;
      Pnucl = 0.45;
      Aeff  = 0.85;
      R1    = R1*1.05;
    }
      else if(Nucleus == 16)
    {
      R1    = 10.50;
      R2    = 5.5;
      Pnucl = 0.7;
      Aeff  = 0.98;
//      R1    = 11.3;
//      R2    = 2.5;
//      Pnucl = 0.75;
//      Aeff  = 0.9;
    }
    else
    {
      R1    = 4.45*std::pow(static_cast<double>(Nucleus-1),0.309);
    if(Nucleus == 28)
      R1    = 4.25*std::pow(static_cast<double>(Nucleus-1),0.309);
      R2    = 2.3*std::pow(static_cast<double>(Nucleus),0.36);
      Pnucl = 0.176+0.00167*Nucleus+
                      8.69E-6*Nucleus*Nucleus;
      Aeff  = 0.9;
      R1    = R1*0.90;
    }

    if(Nucleus == 9)
    {
      R1    = 9.0;
      R2    = 7.0;
      Pnucl = 0.190;
      Aeff  = 0.9;
    }
/*
    if(Nucleus == 12)
   {
      R1    = 9.336;
      R2    = 5.63;
      Pnucl = 0.197;
      Aeff  = 01.0;
   }
*/
/*
   if(Nucleus == 11)
   {
      R1    = 10.8;
      R2    = 7.5;
      Pnucl = 0.85;
      Aeff  = 1.2;
   }
*/
   if(Nucleus == 4)
   {
     R1    = 5.5;
     R1    = 6.0;   //  26.09.06
     R2    = 3.7;
     Pnucl = 0.4;
                  Aeff  = 0.87;
   }

//   G4cout<<" Nucl.Par. "<<Nucleus<<"  R1  "<<R1<<G4endl;
  }
//  ##########################################################
 void  G4ElasticHadrNucleusHE::
  GetHadronValues(const G4ParticleDefinition * aHadron,
                        G4double HadrMoment)
  {

//       G4ParticleDefinition * dHadron = aHadron->GetDefinition();
       G4double protM  = 0.93827;
       G4double protM2 = protM*protM;

       G4int iHadron(-1), iHadrCode;
       iHadrCode = aHadron->GetPDGEncoding();

       if(  iHadrCode == 2212 ||
            iHadrCode == 2112 ||
            iHadrCode == 3122 ||
            iHadrCode == 3222 ||
            iHadrCode == 3112 ||
            iHadrCode == 3212 ||
            iHadrCode == 3312 ||
            iHadrCode == 3322 ||
            iHadrCode == 3334 )   iHadron = 0;

       else if(
            iHadrCode == -2212 ||
            iHadrCode == -2112 ||
            iHadrCode == -3122 ||
            iHadrCode == -3222 ||
            iHadrCode == -3112 ||
            iHadrCode == -3212 ||
            iHadrCode == -3312 ||
            iHadrCode == -3322 ||
            iHadrCode == -3334 )   iHadron = 1;

       else if(  iHadrCode ==  211)     iHadron = 2;
       else if(  iHadrCode == -211)     iHadron = 3;
       else if(  iHadrCode ==  321)     iHadron = 4;
       else if(  iHadrCode == -321)     iHadron = 5;

       else   G4Exception(" There is not method for this hadron ");

       G4double mHadr      = aHadron->GetPDGMass()/1000.;         // In GeV
       G4double HadrEnergy = sqrt(mHadr*mHadr+HadrMoment*HadrMoment);

//       G4double HadrEnergy = aHadron->GetTotalEnergy()/1000.;  // In GeV
//       G4double HadrMoment = aHadron->GetTotalMomentum()/1000.;  // In GeV

       G4double sHadr      = 2*HadrEnergy*protM+protM2+mHadr*mHadr;
       G4double sqrS       = std::sqrt(sHadr);
       G4double Ecm        = (sHadr-mHadr*mHadr+protM2)/2/sqrS;
                MomentumCM = std::sqrt(Ecm*Ecm-protM2);

   if(HadrEnergy-mHadr<1.0)
    {
     G4cout<<" GetHadronValues: The energy = "<<HadrEnergy
           <<"  Hadr. momentum = "<<HadrMoment<<G4endl;
     G4Exception(" The hadron Energy is very low for this method!");
    }

        switch (iHadron)
        {

         case 0:                  //  proton

        G4double Delta;

          Delta=1;

              if(HadrEnergy<40)
                  Delta = 0.916+0.0021*HadrEnergy;
              HadrTot   = 5.2+5.2*std::log(HadrEnergy)
                          +51*std::pow(HadrEnergy,-0.35);            //  mb
              HadrSlope = 6.44+0.88*std::log(sHadr)-1;               //  GeV-2
              HadrReIm  = 0.13*std::log(sHadr/350)*std::pow(sHadr,-0.18);
              DDSect2   = 11;                                    //mb*GeV-2
              DDSect3   = 3;                                     //mb*GeV-2

//    if(HadrEnergy>1000) HadrReIm=0.15;

        if( iHadrCode == 3122 || iHadrCode == 3222 ||
            iHadrCode == 3112 || iHadrCode == 3212 )
            {
              HadrTot   *=0.80;
              HadrSlope *=0.85;
            }

        if( iHadrCode == 3312 || iHadrCode == 3322 )
            {
              HadrTot   *=0.70;
              HadrSlope *=0.75;
            }

         if( iHadrCode == 3334)
            {
              HadrTot   *=0.60;
              HadrSlope *=0.65;
            }

             break;

         case 1:              //   antiproton

              sqrS      = std::sqrt(sHadr);
              HadrTot   = 5.2+5.2*std::log(HadrEnergy)
                          +123.2*std::pow(HadrEnergy,-0.5);           //  mb
              HadrSlope = 8.32+0.57*std::log(sHadr); //GeV-2
           if(HadrEnergy<1000)
              HadrReIm  =0.06*(sqrS-2.236)*(sqrS-14.14)*std::pow(sHadr,-0.8);
           else
              HadrReIm  = 0.6*std::log(sHadr/350)*std::pow(sHadr,-0.25);

              DDSect2   = 11;                                     //mb*GeV-2
              DDSect3   = 3;                                      //mb*GeV-2

//    if(HadrEnergy>1000) HadrReIm=0.15;

        if( iHadrCode == -3122 || iHadrCode == -3222 ||
            iHadrCode == -3112 || iHadrCode == -3212 )
            {
              HadrTot   *=0.75;
              HadrSlope *=0.85;
            }

        if( iHadrCode == -3312 || iHadrCode == -3322 )
            {
              HadrTot   *=0.65;
              HadrSlope *=0.75;
            }

         if( iHadrCode == -3334)
            {
              HadrTot   *=0.55;
              HadrSlope *=0.65;
            }

           break;

         case 2:             //   pi plus

            if(HadrMoment>2.0)
              HadrTot    = 10.6+2.*log(HadrEnergy)+
                              25*std::pow(HadrEnergy,-0.43);           // mb
            else HadrTot = 40-50*(HadrMoment-1.5)*(HadrMoment-1.7);
              HadrSlope = 7.28+0.245*std::log(sHadr);                  //GeV-2
              HadrReIm  = 0.2*std::log(sHadr/100)*std::pow(sHadr,-0.15);
              DDSect2   = 4.6;                                    //mb*GeV-2
              DDSect3   = 1.33;                                   //mb*GeV-2
            break;

         case 3:             //   pi minus

              HadrTot   = 10.6+2*std::log(HadrEnergy)+
                          30*std::pow(HadrEnergy,-0.43);             // mb

            if(HadrMoment<1.399)
              HadrTot = HadrTot+21.0/0.4*(1.4-HadrMoment);

              HadrSlope = 7.28+0.245*std::log(sHadr);               // GeV-2
              HadrReIm  = 0.2*std::log(sHadr/100)*std::pow(sHadr,-0.15);
              DDSect2   = 4.6;                                 //mb*GeV-2
              DDSect3   = 1.33;                                //mb*GeV-2
            break;

         case 4:            //  K plus

              HadrTot   = 10.6+1.8*std::log(HadrEnergy)+
                               9.0*std::pow(HadrEnergy,-0.55);     // mb
         if(HadrEnergy>100) HadrSlope = 15.0;
         else
//              HadrSlope = 5.28+1.76*std::log(sHadr)-
              HadrSlope = 1.0+1.76*std::log(sHadr)-
                              2.84*std::pow(sHadr,-0.5);           // GeV-2
              HadrReIm  = 0.4*(sHadr-20)*(sHadr-150)*std::pow(sHadr+50,-2.1);
              DDSect2   = 3.5;                                //mb*GeV-2
              DDSect3   = 1.03;                               //mb*GeV-2
            break;

         case 5:              //   K minus

              HadrTot   = 10+1.8*std::log(HadrEnergy)
                               +25*std::pow(HadrEnergy,-0.5);   // mb
              HadrSlope = 6.98+0.127*std::log(sHadr);           // GeV-2
//         if(HadrEnergy<8) HadrReIm = 0.7;
//         else
              HadrReIm  = 0.4*(sHadr-20)*(sHadr-20)*std::pow(sHadr+50,-2.1);
              DDSect2   = 3.5;                             //mb*GeV-2
              DDSect3   = 1.03;                            //mb*GeV-2
            break;
      }

  if(verboselevel == 2)
  {
  G4cout<<"  Hadron. Code Mass Energy Momentum  "<<iHadrCode<<"  "
        <<mHadr<<"  "<<HadrEnergy<<"  "<<HadrMoment<<G4endl;
  G4cout<<"  TotCrSec Slope Re/Im "<<HadrTot<<"  "<<HadrSlope
        <<"  "<<HadrReIm<<G4endl;
  }
  }
//  =========================================================

/*  End of file  */
