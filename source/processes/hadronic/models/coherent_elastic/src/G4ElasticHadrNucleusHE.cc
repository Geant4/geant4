//  G4ElasticHadrNucleusHE.cc

#include  "G4ElasticHadrNucleusHE.hh"
#include  "Randomize.hh"
#include  "G4ParticleChange.hh"
#include  "G4Track.hh"
#include  "G4Nucleus.hh"
#include  "G4IonTable.hh"
#include  "globals.hh"
#include  "G4ParticleTable.hh"
#include  "G4NucleiProperties.hh"
#include  "g4std/strstream"
#include  "g4std/iostream"
#include  "g4std/fstream"

//  +++++++++++++++++  The constructor 1 ++++++++++++++++++
     G4ElasticHadrNucleusHE::
        G4ElasticHadrNucleusHE(const G4DynamicParticle * aHadron,
                                     G4Nucleus         * aNucleus) : 
               G4DiffElasticHadrNucleus(),  G4HadronicInteraction()
    {  
//  -------------  here all variables are in MeV  ---------------
        iKindWork      = 0;
        MyIonTable     = new G4IonTable();
        Factorials1[0] = 1;

      G4int ii;
      for( ii = 1; ii<110; ii++)
         {
          Factorials1[ii] = Factorial1(ii);
          }

       G4int  iNnucl = (int) aNucleus->GetN();

  if(iNnucl<11)
         {
  G4Exception(" This nucleus is very light for this model !!!");
         }

  if(iNnucl>208)
         {
  G4Exception(" This nucleus is very heavy for this model !!!");
         }

       G4double  dHdrEnergy = aHadron->GetTotalEnergy()/1000;

  if(dHdrEnergy < 2) 
  G4Exception(" The hadron energy is very low for this model !!!");

        CreationArray (aHadron, aNucleus);

      iContr = 0;
   }  //  Constructor

//  +++++++++++++++++ The constructor 2  ++++++++++++++++++
     G4ElasticHadrNucleusHE::
        G4ElasticHadrNucleusHE(const G4DynamicParticle * aHadron,
                                     G4Nucleus         * aNucleus,
                                     G4double             dEbeg,
                                     G4double             dEend,
                                     G4int                iNpoE,
                                     G4String          sNameFile) :
          G4DiffElasticHadrNucleus(),   G4HadronicInteraction()
  {  
//  -------------  here all variables are in MeV  ---------------

             G4int ii, kk;

             iKindWork   = 1;
             MyIonTable  = new G4IonTable();
             iPoE        = iNpoE;

       G4int  iNnucl = aNucleus->GetN();

  if(iNnucl<11)
         {
  G4Exception(" This nucleus is very light for this model !!!");
         }

  if(iNnucl>208)
         {
  G4Exception(" This nucleus is very heavy for this model !!!");
         }

       G4double  dHdrEnergy = aHadron->GetTotalEnergy()/1000;

  if(dHdrEnergy < 2) 
  G4Exception(" The hadron energy is very low for this model !!!");

        Factorials1[0] = 1;
          for( ii = 1; ii<110; ii++)
              Factorials1[ii] = Factorial1(ii);

             pTableCrSec = new G4double[ONQ2XE];     //  pointer
             pTableE     = new G4double[iNpoE];      //  pointer
             G4double dE = (dEend-dEbeg)/(iNpoE-1);  //   MeV

          G4ThreeVector  inVector(0.,0.,0.);
          G4DynamicParticle * innerHadron = new G4DynamicParticle;

          G4double    curE, curP;         
          G4double    mHadr = aHadron->GetMass();    //   MeV

          innerHadron->SetDefinition(aHadron->GetDefinition());

               for( ii = 0; ii<iNpoE; ii++)
                 {
                    curE         = dEbeg+ii*dE;                  //  MeV
                    curP         = sqrt(curE*curE-mHadr*mHadr);  //  MeV
                  * (pTableE+ii) = curE;

                    inVector.setZ(curP);                //  MeV
//                    innerHadron->SetTotalEnergy(curE);
                    innerHadron->SetMomentum(inVector);

//  ------ The creating of the distribution function array ------
                    CreationArray(innerHadron, aNucleus);

                      for( kk = 0; kk<ONQ2; kk++)
                        {
                         * (pTableCrSec+ii*ONQ2+kk) = iIntgr[kk];


//   G4cout<<" kk iInt Pointer "<<kk<<" "<<iIntgr[kk]<<
//        " "<< * (pTableCrSec+ii*ONQ2+kk)<<endl;
                        }  //   kk step Q2
                 }         //   ii   step Ei

//  ---------   The writing of array into sName File  ---------

       G4std::ofstream TestFile(sNameFile, G4std::ios::out);
             TestFile.setf(G4std::ios::scientific, 
                           G4std::ios::floatfield);
             TestFile.precision(9);

     TestFile << dEbeg<<" "<<dEend<<"  "<<iNpoE<<endl;

    for(ii=0; ii<ONQ2; ii++)
      {
        for(kk=0; kk<iPoE; kk++) 
          {
            TestFile << *(pTableCrSec+kk*ONQ2+ii)<<" ";

          }    //   kk
       TestFile << endl;

      }        //   ii
      TestFile.close();
      iContr = 0;
   }   //  Constructor 2

//  +++++++++++++++++ The constructor 3  ++++++++++++++++++
   G4ElasticHadrNucleusHE::
         G4ElasticHadrNucleusHE(const G4DynamicParticle * aHadron,
                                      G4Nucleus         * aNucleus,
                                      G4String           sNameFile) :
        G4DiffElasticHadrNucleus(),   G4HadronicInteraction()
   {
        MyIonTable  = new G4IonTable();
        G4int ii, kk, ik;
        G4double    aaa;

        iKindWork   = 2;

        G4double aNuc  = aNucleus->GetN();

  if(aNuc<10.5)
         {
  G4Exception(" This nucleus is very light for this model !!!");
         }

  if(aNuc>208.5)
         {
  G4Exception(" This nucleus is very heavy for this model !!!");
         }

       G4double  dHdrEnergy = aHadron->GetTotalEnergy()/1000;

  if(dHdrEnergy < 2) 
  G4Exception(" The hadron energy is very low for this model !!!");

        G4std::ifstream TestFile(sNameFile);
        TestFile.setf(ios::scientific);

        TestFile >> dEbeg1 >> dEend1 >> iPoE;

        G4double dE    = (dEend1-dEbeg1)/(iPoE-1);

        GetNucleusParameters(aNucleus);

        Nstep          = ONQ2;
        maxQ2          = GetQ2limit(R1);        //   MeV^2

        pTableCrSec    = new G4double[ONQ2XE];  //  pointer
        pTableE        = new G4double[5];       //  pointer

     for(kk=0; kk<iPoE; kk++) *(pTableE + kk) = dEbeg1+kk*dE; 

         for(ii=0; ii<ONQ2; ii++)
           {

             for(kk=0; kk<iPoE; kk++) 
                {
                 TestFile >>  *(pTableCrSec+kk*ONQ2+ii);
                }  //  kk
            }       //  ii

      TestFile.close();
      iContr = 0;
   }  //  Constructor   

//  +++++++++++++++++ The constructor 4  ++++++++++++++++++
     G4ElasticHadrNucleusHE::
        G4ElasticHadrNucleusHE(const G4DynamicParticle * aHadron,
                                     G4Nucleus         * aNucleus,
                                     G4double             dEbeg,
                                     G4double             dEend,
                                     G4int                iNpoE) :
        G4DiffElasticHadrNucleus(),   G4HadronicInteraction()
  {  
//  -------------  here all variables are in MeV  ---------------

             G4int ii, kk, ik;
            
       G4double dEbegIn     = 2000;
       G4double dEendIn     = 10000;    
                iKindWork   = 3;
                iPoE        = 5;
                Nstep       = ONQ2;

       G4String sNameHdr    = GetHadronName(aHadron);
       G4int    iNnucl      = (int) aNucleus->GetN();

  if(iNnucl<11)
         {
  G4Exception(" This nucleus is very light for this model !!!");
         }

  if(iNnucl>208)
         {
  G4Exception(" This nucleus is very heavy for this model !!!");
         }

       G4double  dHdrEnergy = aHadron->GetTotalEnergy()/1000;

  if(dHdrEnergy < 2) 
  G4Exception(" The hadron energy is very low for this model !!!");

                MyIonTable  = new G4IonTable();
 
             Factorials1[0] = 1;

          for( ii = 1; ii<110; ii++)
              Factorials1[ii] = Factorial1(ii);

          pTableCrSec      = new G4double[ONQ2XE*AreaNumb]; //pointer
          pTableE          = new G4double[iNpoE*AreaNumb];  //pointer

          G4DynamicParticle * innerHadron = new G4DynamicParticle;

          G4ThreeVector  inVector(0.,0.,0.);

          G4double    dE      = (dEendIn-dEbegIn)/(iPoE-1);//   MeV
          G4double    mHadr   = aHadron->GetMass();        //   MeV
          G4double    curE, curP;         

          innerHadron->SetDefinition(aHadron->GetDefinition());

          G4double    dPower; 

           for(G4int ik = 0; ik<AreaNumb; ik++) 
             {
               dPower = pow(10,ik);
               for( ii = 0; ii<iNpoE; ii++)
                 {  curE         = (dEbeg+ii*dE)*dPower;         //  MeV
                    curP         = sqrt(curE*curE-mHadr*mHadr);  //  MeV
                  * (pTableE+ii+ik*iNpoE) = curE;

                    inVector.setZ(curP);                //  MeV
                    innerHadron->SetMomentum(inVector);

//  ------ The creating of the distribution function array ------
                    CreationArray(innerHadron, aNucleus);

                      for( kk = 0; kk<ONQ2; kk++)
                        {
                         * (pTableCrSec+ii*ONQ2+kk+ik*ONQ2XE) =
                                                      iIntgr[kk];

                        }  //   kk   step Q2
                 }         //   ii   step Ei
             }             //   ik   step power

  G4cout<<" The array for elastic scattering of the "<<sNameHdr<<G4endl;
  G4cout<<" on the Nuclues "<<iNnucl<<" has been prepared !"
        <<G4endl<<G4endl;

//  ---------   The writing of array into sName File  ---------

       G4std::ostrstream sNameFile;
       sNameFile<<"Elastic//"<<sNameHdr<<"_"<<iNnucl<<".dat"<<ends;

       G4std::ofstream TestFile(sNameFile.str(), G4std::ios::out);
             TestFile.setf(G4std::ios::scientific, 
                           G4std::ios::floatfield);
             TestFile.precision(9);

       TestFile << dEbeg<<" "<<dEend<<"  "<<iNpoE<<G4endl;

  for(G4int ik=0; ik<AreaNumb; ik++)
  {
    for(ii=0; ii<ONQ2; ii++)
      {

        for(kk=0; kk<iPoE; kk++) 
          {
            TestFile << *(pTableCrSec+kk*ONQ2+ii+ik*ONQ2XE)<<" ";

          }    //   kk
       TestFile << G4endl;

      }        //   ii
   }           //   ik

      TestFile.close();

  G4cout<<" The array for elastic scattering of "<<sNameHdr<<G4endl;
  G4cout<<" on the Nuclues "<<iNnucl<<" has been written !"<<G4endl;
  G4cout<<" The Name of File is "<<sNameFile.str()<<G4endl<<G4endl;

      iContr = 0;
 }   //  Constructor 4  

//  +++++++++++++++++ The constructor 5  ++++++++++++++++++
   G4ElasticHadrNucleusHE::
         G4ElasticHadrNucleusHE(const G4DynamicParticle * aHadron,
                                      G4Nucleus         * aNucleus,
                                      G4int                 iNpoE) :
        G4DiffElasticHadrNucleus(),   G4HadronicInteraction()
   {
        G4int       ii, kk, ik;
        G4double    aaa;

        G4String sNameHdr  = GetHadronName(aHadron);
        G4int    iNnucl    = (int) aNucleus->GetN();

  if(iNnucl<11)
         {
  G4Exception(" This nucleus is very light for this model !!!");
         }

  if(iNnucl>208)
         {
  G4Exception(" This nucleus is very heavy for this model !!!");
         }

       G4double  dHdrEnergy = aHadron->GetTotalEnergy()/1000;

  if(dHdrEnergy < 2) 
  G4Exception(" The hadron energy is very low for this model !!!");

        G4std::ostrstream sNameFile;
        sNameFile<<"Elastic//"<<sNameHdr<<"_"<<iNnucl<<".dat"<<ends;

  G4cout <<" In constr. 5. Hadron - "<<sNameHdr
         <<".  Nucleus - "<<iNnucl<<G4endl;
  G4cout <<" The Name of File is "<<sNameFile.str()<<G4endl;

        MyIonTable     = new G4IonTable();

        GetNucleusParameters(aNucleus);

        iContr         = iNpoE;
        Nstep          = ONQ2;
        iKindWork      = 4;
        maxQ2          = GetQ2limit(R1);               //   MeV^2

        pTableCrSec    = new G4double[ONQ2XE*AreaNumb]; //  pointer
        pTableE        = new G4double[iNpoE*AreaNumb];  //  pointer

        G4double aNuc  = aNucleus->GetN();

        G4std::ifstream TestFile(sNameFile.str(), G4std::ios::in);
             TestFile.setf(G4std::ios::scientific);

        TestFile >> dEbeg1 >> dEend1 >> iPoE;

        G4double dE    = (dEend1-dEbeg1)/(iPoE-1);

          G4double  dPower;

   for(ik=0; ik<AreaNumb; ik++)
    {
           dPower = pow(10,ik);

     for(kk=0; kk<iPoE; kk++) 
        {               
            *(pTableE + kk + ik*iPoE) = (dEbeg1+kk*dE)*dPower; 

        }    //  kk

         for(ii=0; ii<ONQ2; ii++)
           {

             for(kk=0; kk<iPoE; kk++) 
                {
                 TestFile >>  *(pTableCrSec+kk*ONQ2+ii+ik*ONQ2XE);

                }  //  kk

           }       //  ii

   }               //  ik

      TestFile.close();

  G4cout<<" The array for elastic scattering of "<<sNameHdr<<G4endl;
  G4cout<<" on the Nuclues "<<iNnucl<<G4endl;
  G4cout<<" from the file "<<sNameFile.str()<<" has been readout !"
        <<G4endl<<G4endl;

   }  //  Constructor   

//  ++++++++++++++++++++   GetQ2limit  +++++++++++++++++++
     G4double G4ElasticHadrNucleusHE::
                   GetQ2limit(G4double R1)
   {
        G4double maxQ2 = 35*1000*1000/R1/R1;
                 dQ2   = maxQ2/(Nstep*3./2.-1.);

        iQ2[0] = 1;
        for(G4int ii=1; ii<Nstep; ii++)
            iQ2[ii] = ii<Nstep/2 ? 
              iQ2[ii-1]+dQ2 : iQ2[ii-1]+2*dQ2;

        return maxQ2;
   }

//  ++++++++++++++++++   GetHadronName  +++++++++++++++++++
     G4String G4ElasticHadrNucleusHE::
                   GetHadronName(const G4DynamicParticle * aHadron)
   {
       G4ParticleDefinition * dHadron = aHadron->GetDefinition();

       G4String sHadron;

         if(dHadron == G4Proton::Proton()||
            dHadron == G4Neutron::Neutron())         
                                       sHadron="Nucleon";

   else  if(dHadron == G4AntiProton::AntiProton()||
            dHadron == G4AntiNeutron::AntiNeutron())
                                       sHadron="AntiNucleon";

   else  if(dHadron == G4PionPlus::PionPlus())       
                                       sHadron="PionPlus";

   else  if(dHadron == G4PionMinus::PionMinus())     
                                       sHadron="PionMinus";

   else  if(dHadron == G4KaonPlus::KaonPlus())       
                                       sHadron="KaonPlus";

   else  if(dHadron == G4KaonMinus::KaonMinus())     
                                       sHadron="KaonMinus";

   else {
         G4cout << " There is not method for the hadron " <<
                  dHadron->GetParticleName() << endl;
        }
    return sHadron;
   }

//  ++++++++++++++++++  ApplayYourself  ++++++++++++++++++++
    G4VParticleChange *
    G4ElasticHadrNucleusHE::ApplyYourself(
                          const  G4Track    &aTrack,
                                 G4Nucleus  &aNucleus)
 {
              G4int    nN, nZ;

        const G4DynamicParticle * aParticle = aTrack.GetDynamicParticle();
              G4Nucleus *         aNucl     = &aNucleus;

        G4double   aNuclZ     = aNucl->GetZ();    
                  aNucleon    = aNucl->GetN();

                         nN   = aNucleon;
                         nZ   = aNuclZ;

        G4ParticleDefinition * secNuclDef;

        secNuclDef  =   MyIonTable->GetIon(  nZ,  nN,  0,  nZ);

        G4DynamicParticle *   secNuclDyn = new G4DynamicParticle();

        secNuclDyn->SetDefinition(secNuclDef);

        G4ThreeVector      NuclPos;

        G4Track      secNuclTrack = G4Track(secNuclDyn,
                                          0.0, NuclPos);

        G4double   ranQ2;       //      ranQ2   -  MeV^2

//  ____________ Randomization of Q2 as dSis/dt ______________

        if(iKindWork == 0)  ranQ2 = RandomElastic0(aParticle, aNucl);
        else                ranQ2 = RandomElastic1(aParticle, aNucl);

        G4double   inLabMom  = aParticle->GetTotalMomentum(); // MeV
        G4double   inEnHadr  = aParticle->GetTotalEnergy();   // MeV
        G4double   MassHadr  = aParticle->GetMass();          // MeV
        G4double   MassNucl  = aNucleus.GetN()*938;           // MeV
        G4double   sqrMass   = MassNucl*MassNucl+MassHadr*MassHadr;


// ---------------    For final state of hadron    ------------------

        G4double   invS      = sqrMass+2*inEnHadr*MassNucl;   // MeV^2
        G4double   invU      = 2*sqrMass-invS+ranQ2;          // MeV^2
        G4double   outEnHadr = (sqrMass-invU)/2/MassNucl;     // MeV
        G4double   outMomHdr = sqrt(outEnHadr*outEnHadr-      // MeV
                                    MassHadr*MassHadr);
        G4double   cosHadr   = (-ranQ2-2*MassHadr*MassHadr+
                               2*inEnHadr*outEnHadr)/2
                               /inLabMom/outMomHdr;        

//----------------     For final state of nucleus   -----------------

        G4double   outEnNucl = (2*MassNucl*MassNucl+ranQ2)/2/MassNucl;

        G4double   outMomNcl = sqrt(outEnNucl*outEnNucl-  //  MeV
                                    MassNucl*MassNucl);
        G4double   cosNucl   = (invU-sqrMass+2*inEnHadr*outEnNucl)
                                /2/inLabMom/outMomNcl;   

        G4double   ranFi     = 6.2832*G4UniformRand();

// ------ The transformation from Hadron frame to absolute one -------

             G4double   RotMatrix[3][3];
             G4double   NewMomHdrSys[3];
             G4double   NewMomOldSys[3];
             G4double   NucMomHdrSys[3];
             G4double   NucMomOldSys[3];

//   ------------ The hadron angles in its own system ------------
             G4double   sinHadr   = sqrt(1-cosHadr*cosHadr);
             G4double   sinFiHadr = sin(ranFi);
             G4double   cosFiHadr = cos(ranFi);

//   ---- The unit vector of a hadron momentum in its own system ---
                   NewMomHdrSys[0] = sinHadr*cosFiHadr; 
                   NewMomHdrSys[1] = sinHadr*sinFiHadr;
                   NewMomHdrSys[2] = cosHadr;          
        G4double   sinNucl         = sqrt(1-cosNucl*cosNucl);

//   -----------  The unit vector of  a nucleus momentum  -----------
                   NucMomHdrSys[0] = -sinNucl*cosFiHadr;
                   NucMomHdrSys[1] = -sinNucl*sinFiHadr;
                   NucMomHdrSys[2] =  cosNucl;

//  __________ The angles of unit vector in absolute system ----------
        G4double   cosTet  = aParticle->GetMomentumDirection().z();
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
               sinTet  = sqrt(1-cosTet*cosTet);
               sinFi   = aParticle->GetMomentumDirection().y()/
                         sinTet;
               cosFi   = aParticle->GetMomentumDirection().x()/
                         sinTet;
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

if(iContr == 137 )
  {
 G4cout<<" Angles CosTet SinTet CosFi SinFi "<<cosTet<<" "<<sinTet
       <<" "<<cosFi<<" "<<sinFi<<endl;
 G4cout<<" RotMatr "<< RotMatrix[0][0]<<" "
       << RotMatrix[0][1]<<" "<<RotMatrix[0][2]<<G4endl;
 G4cout<<"         "<< RotMatrix[1][0]<<" "<<RotMatrix[1][1]<<" "
       << RotMatrix[1][2]<<G4endl;
 G4cout<<"         "<< RotMatrix[2][0]<<" "<<RotMatrix[2][1]<<" "
       << RotMatrix[2][2]<<G4endl;
  }

                   for(G4int ii=0; ii<3; ii++)
                     {
                       NewMomOldSys[ii] = 0;                        
                       NucMomOldSys[ii] = 0;                        
                         for(G4int ll=0; ll<3; ll++)
                          {
                         NewMomOldSys[ii] = NewMomOldSys[ii]+
                               RotMatrix[ii][ll]*NewMomHdrSys[ll];
                         NucMomOldSys[ii] = NucMomOldSys[ii]+
                               RotMatrix[ii][ll]*NucMomHdrSys[ll];
                          }
                     }

              G4ThreeVector  aNuclMom( NucMomOldSys[0],     // in MeV
                                       NucMomOldSys[1],     //    MeV
                                       NucMomOldSys[2]);    //    MeV

                secNuclDyn->SetMomentum(aNuclMom);
                secNuclTrack.SetKineticEnergy(outEnNucl-MassNucl);

                secNuclTrack.SetMomentumDirection(aNuclMom);

                G4double pxnew = NewMomOldSys[0];    //  In  MeV
                G4double pynew = NewMomOldSys[1];    //      Mev
                G4double pznew = NewMomOldSys[2];    //      Mev

     if(iContr == 137)
{
G4cout<<G4endl<<" Q2 "<< ranQ2<<" Hdr Mom Tet "<<outMomHdr
      <<" "<< acos(cosHadr)*180/3.1416
      <<" Nucl Mom Tet "<<sqrt(pow(NucMomOldSys[0],2)+
                               pow(NucMomOldSys[1],2)+
                               pow(NucMomOldSys[2],2))
                              *outMomNcl<<"  "
      <<acos(cosNucl)*180/3.1416<<G4endl<<G4endl;
G4cout<<" Rotation: old UnitVec "
      <<NewMomHdrSys[0]<<" "<<NewMomHdrSys[1]<<" "
<<NewMomHdrSys[2]<<" Fi "<<ranFi<<G4endl;
G4cout<<"           new UnitVec "<<pxnew<<" "<<pynew<<" "
      <<pznew<<G4endl<<G4endl;
}

     if(iKindWork>0)
              theParticleChange.AddSecondary( &secNuclTrack);

       theParticleChange.SetEnergyChange(outEnHadr);
       theParticleChange.SetMomentumDirectionChange(pxnew, pynew, pznew);

       return &theParticleChange;
  }

// +++++++  The randomization of one dimensional array +++++++
    G4double G4ElasticHadrNucleusHE::RandomElastic0(
                           const G4DynamicParticle * aHadron,
                                 G4Nucleus         * aNucleus)
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
                  {
               ranQ2 = (iQ2[kk-1]+(ranUni-iIntgr[kk-1])*
                        (iQ2[kk]-iQ2[kk-1])
                        /(iIntgr[kk]-iIntgr[kk-1]));
                  }   

          else    {
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
                                 G4Nucleus         *  aNucleus)
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

       while(curE > *(pTableE+iiE) ) 
         {
           iiE++;

///   G4cout<<" Rand1: iiE "<< iiE <<" E "<<*(pTableE+iiE)<<endl;

         }

         G4double ranQ2;
         G4double ranUni = G4UniformRand();

            G4int kk;

       for(G4int mm=0; mm<3; mm++)
          {
                kk = 0;

            while( ranUni > 
                   *(pTableCrSec+iiE*ONQ2+mm*ONQ2+kk-ONQ2))  
             {

//  G4cout<<"find kk  ranUni Tab "<<kk<<" " <<ranUni<< "  "<<
//           *(pTableCrSec+iiE*ONQ2+mm*ONQ2+kk-ONQ2)<<endl;

               kk ++;
               k[mm]=kk;
             }     //   while   

///   G4cout<<" ranUni kk iIntgr[kk] "<<ranUni<<" "<<
///      kk<<" "<< *(pTableCrSec+iiE*ONQ2+mm*ONQ2+kk-ONQ2)<<endl;

//           if(kk == 1 ) kk = 2;
              X1 = *(pTableCrSec+iiE*ONQ2+mm*ONQ2+kk-ONQ2-1);
              X2 = *(pTableCrSec+iiE*ONQ2+mm*ONQ2+kk-ONQ2);
              X3 = *(pTableCrSec+iiE*ONQ2+mm*ONQ2+kk-ONQ2+1);
              Y1 = iQ2[kk-1];            //  MeV^2
              Y2 = iQ2[kk];
              Y3 = iQ2[kk+1];

       forQ2[mm] = InterPol(X1, X2, X3, Y1, Y2, Y3, ranUni);

/// G4cout<<endl;
/// G4cout<<" a   mm X1 X2 X3 "<<mm<<" "<<X1<<" "<<X2<<" "<<X3<<endl;
/// G4cout<<" a Y1 Y2 Y3  " <<Y1<< "  " <<Y2<< "  " <<Y3<< " forQ2 "
///                      <<forQ2[mm]<<endl;  
/// G4cout<<endl;

           }     //  mm

              X1 = *(pTableE+iiE-1);
              X2 = *(pTableE+iiE);
              X3 = *(pTableE+iiE+1);
              Y1 = forQ2[0];             //  MeV^2
              Y2 = forQ2[1];
              Y3 = forQ2[2];
           ranQ2 = InterPol(X1, X2, X3, Y1, Y2, Y3, curE);

/*
  if(ranQ2<0||ranQ2>250000)
   {
    G4int l=k[0];
  G4cout<<" mm=0 Q2 -1 0 +1 "<<iQ2[l-1]<<"  "<<iQ2[l]<<"  "
        <<iQ2[l+1]<<G4endl;
          l=k[1];
  G4cout<<" mm=1 Q2 -1 0 +1 "<<iQ2[l-1]<<"  "<<iQ2[l]<<"  "
        <<iQ2[l+1]<<G4endl;
          l=k[2];
  G4cout<<" mm=2 Q2 -1 0 +1 "<<iQ2[l-1]<<"  "<<iQ2[l]<<"  "
        <<iQ2[l+1]<<G4endl;
 
  G4cout<<" forQ2 : "<<forQ2[0]<<" "<<forQ2[1]<<" "
        <<forQ2[2]<<G4endl;
  
  G4cout<<" 1-ranUni : "<<1-ranUni<<G4endl;
  G4cout<<" Differ 0 : "<<1-
          *(pTableCrSec+iiE*ONQ2+0*ONQ2+k[0]-1-ONQ2)
        <<" "<<1-
          *(pTableCrSec+iiE*ONQ2+0*ONQ2+k[0]-ONQ2)
        <<" "<<1-
          *(pTableCrSec+iiE*ONQ2+0*ONQ2+k[0]+1-ONQ2)<<endl;
  G4cout<<" Differ 1 : "<<1-
          *(pTableCrSec+iiE*ONQ2+0*ONQ2+k[1]-1-ONQ2)
        <<" "<<1-
          *(pTableCrSec+iiE*ONQ2+0*ONQ2+k[1]-ONQ2)
        <<" "<<1-
          *(pTableCrSec+iiE*ONQ2+0*ONQ2+k[1]+1-ONQ2)<<endl;
  G4cout<<" Differ 2 : "<<1-
          *(pTableCrSec+iiE*ONQ2+0*ONQ2+k[2]-1-ONQ2)
        <<" "<<1-
          *(pTableCrSec+iiE*ONQ2+0*ONQ2+k[2]-ONQ2)
        <<" "<<1-
          *(pTableCrSec+iiE*ONQ2+0*ONQ2+k[2]+1-ONQ2)<<endl;

  G4cout<<"find kk  ranUni Tab mm=0 "<<kk<<" " <<ranUni<< "  "<<
           *(pTableCrSec+iiE*ONQ2+0*ONQ2+k[0]-1-ONQ2)<<endl;
  G4cout<<"find kk  ranUni Tab mm=0 "<<kk<<" " <<ranUni<< "  "<<
           *(pTableCrSec+iiE*ONQ2+0*ONQ2+k[0]-ONQ2)<<endl;
  G4cout<<"find kk  ranUni Tab mm=0 "<<kk<<" " <<ranUni<< "  "<<
           *(pTableCrSec+iiE*ONQ2+0*ONQ2+k[0]+1-ONQ2)<<endl;
  G4cout<<"find kk  ranUni Tab mm=1 "<<kk<<" " <<ranUni<< "  "<<
           *(pTableCrSec+iiE*ONQ2+1*ONQ2+k[1]-1-ONQ2)<<endl;
  G4cout<<"find kk  ranUni Tab mm=1 "<<kk<<" " <<ranUni<< "  "<<
           *(pTableCrSec+iiE*ONQ2+1*ONQ2+k[1]-ONQ2)<<endl;
  G4cout<<"find kk  ranUni Tab mm=1 "<<kk<<" " <<ranUni<< "  "<<
           *(pTableCrSec+iiE*ONQ2+1*ONQ2+k[1]+1-ONQ2)<<endl;
  G4cout<<"find kk  ranUni Tab mm=2 "<<kk<<" " <<ranUni<< "  "<<
           *(pTableCrSec+iiE*ONQ2+2*ONQ2+k[2]-1-ONQ2)<<endl;
  G4cout<<"find kk  ranUni Tab mm=2 "<<kk<<" " <<ranUni<< "  "<<
           *(pTableCrSec+iiE*ONQ2+2*ONQ2+k[2]-ONQ2)<<endl;
  G4cout<<"find kk  ranUni Tab mm=2 "<<kk<<" " <<ranUni<< "  "<<
           *(pTableCrSec+iiE*ONQ2+2*ONQ2+k[2]+1-ONQ2)<<endl;

 G4cout<<" kk "<<k[0]<<" "<<k[1]<<" "
       <<k[2]<<" ranUni "<<ranUni<<G4endl;
 G4cout<<" X1 X2 X3 "<<X1<<" "<<X2<<" "<<X3;
 G4cout<<" Y1 Y2 Y3 "<<Y1<<" "<<Y2<<" "<<Y3<<" ranQ2 "<<ranQ2<<endl;
 G4cout<<endl;

   }            //  if
*/

        return ranQ2;          //   MeV^2
 }

//  ++++++++++++ The filling of the one-dimensional array +++++++++++++
void  G4ElasticHadrNucleusHE::CreationArray(
              const G4DynamicParticle * aHadron,
                    G4Nucleus         * aNucleus)
    {
            GetNucleusParameters(aNucleus);

            Nstep  = ONQ2;
            maxQ2  = GetQ2limit(R1);             //   Mev^2

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

       G4double aNuc  = aNucleus->GetN();
       G4double intgrS, intgStep;

               iIntgr[0] = 0;

             for(G4int ii=0; ii<Nstep-1; ii++) 
              {

// ----------------  The integration  -------------

                G4double  curQ2  = iQ2[ii];             //  MeV^2
                G4double  ddQ2   =(iQ2[ii+1]-curQ2)/20; //  MeV^2
                G4double  curSum = 0; 
                G4double  curSec = 0;

                G4double aSimp = 1;
                   for(G4int ll=0; ll<20; ll++)                     
                     {                                              
                       curSec = aDiffElHadNcls.
                                 HadrNuclDifferCrSec(aHadron, 
                                                     aNucleus, 
                                                     curQ2); 

                       curQ2  = curQ2 + ddQ2;                        
                       curSum = curSum + curSec*aSimp;        

                       if(aSimp > 3 ) aSimp = 2;
                       else  aSimp = 4;
                       if(ll == 0) aSimp = 4;
                     }      //  ll                           

                    intgStep = curSum*ddQ2/3;      
//  --------------------------------------------------------------------

                iIntgr[ii+1] = iIntgr[ii]+intgStep; 
                intgrS       = intgrS + intgStep;
              }   //  ii

//  -------------------  The normalization  ----------------
            for(G4int ii=0; ii<Nstep; ii++)
             {
               iIntgr[ii] = iIntgr[ii]/iIntgr[Nstep-1];
              }   //  ii

    }             //  end  subroutine

//  +++++++++++++++++  The  method for light nucleus  +++++++++++++++++++
   void G4ElasticHadrNucleusHE::ArrayForLight( 
                          const G4DynamicParticle * aHadron,
                                G4Nucleus         * aNucleus)
     {
// -------------- !!! Here the main unit is GeV !!!  ---------------
       G4HadronValues::GetHadronValues(aHadron);
       G4double Ehad     = aHadron->GetTotalEnergy()/1000;  //  GeV

       G4double aNuc     = aNucleus->GetN();
       G4int    Nucleus  = aNuc;
       G4double intgrS;

             GetNucleusParameters(aNucleus);

       G4double PnuclP   = 0.001;

                maxQ2    = GetQ2limit(R1)/1000/1000;        //  GeV^2

// ----------------  The preparing of probability function  ------------

      G4double    MbToB    = 2.568;             //  from mb to GeV^-2
      G4double    Pi1      = 3.1416;
      G4double    Stot     = HadrTot*MbToB;     //  Gev^-2
      G4double    Bhad     = HadrSlope;         //  GeV^-2
      G4double    Asq      = 1+HadrReIm*HadrReIm;
      G4double    Rho2     = sqrt(Asq);

      G4double    R12      = R1*R1;
      G4double    R22      = R2*R2;
      G4double    R12B     = R12+2*Bhad;
      G4double    R22B     = R22+2*Bhad;
      G4double    R12Bp    = R12+20;
      G4double    R22Bp    = R22+20;
      G4double    R13Bp    = R12*R1/R12Bp;
      G4double    R23Bp    = R22*R2/R22Bp;
      G4double    R12Ap    = R12+20;
      G4double    R22Ap    = R22+20;
      G4double    R13Ap    = R12*R1/R12Ap;
      G4double    R23Ap    = R22*R2/R22Ap*PnuclP;
      G4double    R23dR13  = R23Ap/R13Ap;
      G4double    R12Apd   = 2/R12Ap;
      G4double    R22Apd   = 2/R22Ap;

      G4double    Norm     = R12*R1-Pnucl*R22*R2;
      G4double    NormP    = R12*R1-PnuclP*R22*R2;
      G4double    R13      = R12*R1/R12B;
      G4double    R23      = Pnucl*R22*R2/R22B;
      G4double    Unucl    = Stot/2/Pi1/Norm*R13;
      G4double    Unclprod = Stot/2/Pi1/NormP*R13Ap;
      G4double    FiH      = asin(HadrReIm/Rho2);
      G4double    NN2      = R23/R13;

      G4double    DDSec1p  = (DDSect2+
                    DDSect3*log(1.06*2*Ehad/R1/4));

      G4double    DDSec2p  = (DDSect2+
                    DDSect3*log(1.06*2*Ehad/sqrt((R12+R22)/2)/4));

      G4double    DDSec3p  = (DDSect2+
                    DDSect3*log(1.06*2*Ehad/R2/4));

      G4double    R12ApdR22Ap = 0.5*(R12Apd+R22Apd);

                 iIntgr[0] = 0;

             for(G4int ii=0; ii<Nstep-1; ii++) 
              {

      G4double    Q2       = iQ2[ii+1]/1000/1000;// internal in GeV^2

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
                G4double Prod2  = 0; //exp(-Q2/i2*R12B/4)/i2*R12B;
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
                                   (1-exp(-Q2*(1/exp1+1/exp2)/4))/
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
                 Prod1 = Prod1 + Prod2*N2*cos(FiH*(i1-i2));
                 Tot0  = Tot0  + Prod2*N2*sin(FiH*(i1-i2));
            if (abs(Prod2*N2/Prod1)<1e-6) break;

                   }                                         // i2
                  // ImDistr = Tot0  + Tot0*N1;
                     Prod0   = Prod0 + Prod1*N1;
            if(abs(N1*Prod1/Prod0) < 1e-6) break;

                 }                                           // i1
             Prod0        = Prod0*Pi1/2.568/4;  //  This is in mb
             iIntgr[ii+1] = Prod0;  //iIntgr[ii-1]+Prod0;

              }                                            //   ii (Q2)

          for(G4int ii = 0; ii<Nstep; ii++)
                   {  
                  iIntgr[ii]=iIntgr[ii]/iIntgr[Nstep-1];

                   }

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

///     if(Y1<0 || Y1>250000) 
///        G4cout<<" F12 F22 F32 D0 : "<<F12<<" "<<F22
///          <<" "<<F32<<" "<<D0<<endl;

          if(abs(D0) < 1e-8 || D0 == 0) 
                  {
                  ranQ2 = (Y2+(X-X2)*
                           (Y3-Y2) /(X3-X2));            //   MeV^2
                  }   //  ii

          else    {
             G4double DA = Y1*X2+Y3*X1+Y2*X3-Y3*X2-Y1*X3-Y2*X1;
             G4double DB = Y2*F12+Y1*F32+Y3*F22-Y2*F32-Y3*F12-Y1*F22;
             G4double DC = Y3*X2*F12+Y2*X1*F32+Y1*X3*F22
                           -Y1*X2*F32-Y2*X3*F12-Y3*X1*F22;
                   ranQ2 = (DA*X*X+DB*X+DC)/D0;           //   MeV^2
                  }
           return  ranQ2;
    }

/*  End of file  */
