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
#include  <rw/tpordvec.h>
#include  <rw/cstring.h>

    G4ParticleChange *
    G4ElasticHadrNucleusHE::ApplyYourself(
                          const  G4Track    &aTrack,
                                 G4Nucleus  &aNucleus)
 {
              G4int    nN, nZ;
        const G4DynamicParticle * aParticle = aTrack.GetDynamicParticle();
              G4Nucleus *         aNucl     = &aNucleus;
//              G4IonTable         aIonTable; // = new G4IonTable();

        G4double   aNuclZ     = aNucl->GetZ();    
                  aNucleon    = aNucl->GetN();
//  G4cout<< " 1 "<<G4endl;

               nN   = aNucleon;
               nZ   = aNuclZ;

        G4ParticleDefinition * secNuclDef; 

//        secNuclDef     =  G4ParticleTable::GetParticleTable()->
//           GetIonTable()->GetParticle(1);

        secNuclDef     =  G4ParticleTable::GetParticleTable()->
           GetIonTable()->GetIon(  6,  12,  0,  6);

//  G4cout<<" Name "<<  secNuclDef->GetParticleName()<<G4endl;

//  G4cout<< " 2 "<<G4endl;

//  G4cout<< " 3 "<<G4endl;

        G4DynamicParticle *   secNuclDyn = new G4DynamicParticle();
//  G4cout<< " 4 "<<G4endl;

            secNuclDyn->SetDefinition(secNuclDef);
//  G4cout<< " 5 "<<G4endl;

        G4ThreeVector          NuclPos;
//  G4cout<< " 6 "<<G4endl;

        G4Track                secNuclTrack = G4Track(secNuclDyn,
                                                      0.0, NuclPos);

        G4double   ranQ2     = RandomElastic(aParticle, aNucl);

        G4double   inLabMom  = aParticle->GetTotalMomentum()/1000;
        G4double   inEnHadr  = aParticle->GetTotalEnergy()/1000;
        G4double   MassHadr  = aParticle->GetMass()/1000;
        G4double   MassNucl  = aNucleus.GetN()*0.938;
        G4double   sqrMass   = MassNucl*MassNucl+MassHadr*MassHadr;

//  G4cout<<" Nucleus Mass "<<MassNucl<<" inEnHadr "<<
//         inEnHadr<<" inLabMom "<<inLabMom<<" MassH "<<
//             MassHadr<<G4endl;

// ---------------    For final state of hadron    ------------------

        G4double   invS      = sqrMass+2*inEnHadr*MassNucl;
        G4double   invU      = 2*sqrMass-invS+ranQ2;
        G4double   outEnHadr = (sqrMass-invU)/2/MassNucl;
        G4double   outMomHdr = sqrt(outEnHadr*outEnHadr-
                                    MassHadr*MassHadr);
        G4double   cosHadr   = (-ranQ2-2*MassHadr*MassHadr+
                          2*inEnHadr*outEnHadr)/2/inLabMom/outMomHdr;        

//  G4cout << " The angle of output Hadron : "<<
//                        acos(cosHadr)*180/3.1416<<
//  " Energy "<<outEnHadr<< " Mom "<<outMomHdr<<" Q2 "<< ranQ2<< G4endl;

//----------------     For final state of nucleus   -----------------

//  G4cout<<" sqrMass "<<sqrMass<<" ranQ2 "<<ranQ2<<" MassNucl "<<
//                MassNucl<<G4endl;
        G4double   outEnNucl = (2*MassNucl*MassNucl+ranQ2)/2/MassNucl;
        G4double   outMomNcl = sqrt(outEnNucl*outEnNucl-
                                    MassNucl*MassNucl);
        G4double   cosNucl   = (invU-sqrMass+2*inEnHadr*outEnNucl)
                                /2/inLabMom/outMomNcl;

//   G4cout << " The angle of output Nucleus : "<<
//                        acos(cosNucl)*180/3.1416<<
//    " Energy " << outEnNucl << " Mom "<<outMomNcl<< G4endl;

        G4double   ranFi     = 6.2832*G4UniformRand();

// ------ The transfer from Hadron frame to absolute one -------

             G4double   RotMatrix[3][3];
             G4double   NewMomHdrSys[3];
             G4double   NewMomOldSys[3];
             G4double   NucMomHdrSys[3];
             G4double   NucMomOldSys[3];

             G4double   sinHadr   = sqrt(1-cosHadr*cosHadr);
             G4double   sinFiHadr = sin(ranFi);
             G4double   cosFiHadr = cos(ranFi);

                        NewMomHdrSys[0] = sinHadr*cosFiHadr;                        
                        NewMomHdrSys[1] = sinHadr*sinFiHadr;
                        NewMomHdrSys[2] = cosHadr;

             G4double   sinNucl   = sqrt(1-cosNucl*cosNucl);

                        NucMomHdrSys[0] = -sinNucl*cosFiHadr;                        
                        NucMomHdrSys[1] = -sinNucl*sinFiHadr;
                        NucMomHdrSys[2] =  cosNucl;

        G4double   cosTet  = aParticle->GetMomentumDirection().z();

        G4double   sinTet;
        G4double   sinFi;
        G4double   cosFi;

                if(cosTet>0.9999 || cosTet<-0.9999) 
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

//   G4cout<<" Apply InLabMom "<< inLabMom <<G4endl;
//   G4cout<<" Apply Cos(Tet) initial part. "<< cosTet <<
//          " sinTet "<<sinTet<< G4endl;

                        RotMatrix[0][0]   =  cosFi*cosTet;
                        RotMatrix[0][1]   =  sinFi;
                        RotMatrix[0][2]   = -cosFi*sinTet;
                        RotMatrix[1][0]   = -sinFi*cosTet;
                        RotMatrix[1][1]   =  cosFi;
                        RotMatrix[1][2]   =  sinFi*sinTet;
                        RotMatrix[2][0]   =  sinTet;
                        RotMatrix[2][1]   =  0;
                        RotMatrix[2][2]   =  cosTet;

//   G4cout<<" Apply The rotation matrix : "<<G4endl;
//   G4cout<<RotMatrix[0][0]<< 
//                    " "<<RotMatrix[0][1]<<" "<< RotMatrix[0][2]<<G4endl;
//   G4cout<<RotMatrix[1][0]<< 
//                    " "<<RotMatrix[1][1]<<" "<< RotMatrix[1][2]<<G4endl;
//   G4cout<<RotMatrix[2][0]<< 
//                    " "<<RotMatrix[2][1]<<" "<< RotMatrix[2][2]<<G4endl;

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

                G4ThreeVector  aNuclMom( NucMomOldSys[0],
                                         NucMomOldSys[1],
                                         NucMomOldSys[2]);
//  G4cout<<" 8 "<<G4endl;
                secNuclDyn->SetMomentum(aNuclMom);
                secNuclTrack.SetKineticEnergy(outEnNucl-MassNucl);
//  G4cout<<" 9 "<<G4endl;

                secNuclTrack.SetMomentumDirection(aNuclMom);

//  G4cout<<" 10 "<<G4endl;

                G4double pxnew = NewMomOldSys[0];
                G4double pynew = NewMomOldSys[1];
                G4double pznew = NewMomOldSys[2];

//   G4cout<<" Apply Mom "<< NewMomOldSys[0]<<" "<< NewMomOldSys[1]<<
//                " "<< NewMomOldSys[2]<<" Angle "<<
//                        acos(cosHadr)*180/3.1416<<G4endl;

//   G4cout<<" "<<outMomHdr<<" "<<
//                acos(cosHadr)*180/3.1416<<" ";

//   G4cout<<" "<<sqrt(pow(NucMomOldSys[0],2)+
//            pow(NucMomOldSys[1],2)+
//            pow(NucMomOldSys[2],2))*outMomNcl<<" "<<
//            acos(cosNucl)*180/3.1416<<G4endl;
     
//     G4cout<<G4endl;

          theParticleChange.AddSecondary( &secNuclTrack);
//  G4cout<<" 11 "<<G4endl;

          theParticleChange.SetEnergyChange(outEnHadr*1000);
          theParticleChange.SetMomentumDirectionChange(pxnew, pynew, pznew);

          return &theParticleChange;
  }

// --------------------- The randomization of Q2 ---------------------

    G4double G4ElasticHadrNucleusHE::RandomElastic(
                           const G4DynamicParticle *  aHadron,
                                 G4Nucleus *          aNucleus)
         {

                      Nstep = 100;

             G4double R0    = sqrt(0.84*25.68)*pow(aNucleon, 0.3333);
             G4double maxQ2 = 3/R0*1000*1000; 

                      iQ2[0]    = 0;
                      iIntgr[0] = 0;
             G4double dQ2       = maxQ2/Nstep;
             G4double intgrS    = 0;
             G4double intgStep;
             G4int    intgN;

// ----------------  The preparing of probability function  ------------

//  G4cout << " In Random 3: R0 " << R0<<" dQ2 "<<dQ2<<
//                         " Nstep "<<Nstep<< G4endl;

        G4double Mom1 = aHadron->GetTotalMomentum();
        G4double Enr1 = aHadron->GetTotalEnergy();

//  G4cout << "In Random : Mom1 "<<Mom1<<" Enr1 "<<Enr1<<G4endl;
//  G4cout << G4endl;

             for(G4int ii=0; ii<Nstep-1; ii++) 
              {
                iQ2[ii+1]    = iQ2[ii]+dQ2;

// ---------------- The nonoptimal integration for gebuging -------------

                   G4double  ddQ2   =  dQ2/100;                      //
                   G4double  curQ2  = iQ2[ii];                       //
                   G4double  curSum = 0;                             //
                   G4double  curSec;                                 //
// G4cout << " befor Intgr "<<"ii "<<ii<<" "<<curSec<<G4endl;

                   for(G4int ll=0; ll<99; ll++)                      //
                     {                                               //
              curSec = aDiffElHadNcls.HadrNuclDifferCrSec(           //
                                          aHadron, aNucleus, curQ2); //
                       curQ2 = curQ2 + ddQ2;                         //
                      curSum = curSum + curSec*coefSimp[ll];         //
//G4cout << " In Intgr "<<ll<<" "<<" "<<curQ2<<" "<<curSec
//        <<" "<<coefSimp[ll]<<G4endl;
                     }      //  ll                                   //

//G4cout << " In out Intgr "<<ii<<" "<<curSum<<curQ2<<G4endl;
                    intgStep = curSum*ddQ2/3/1000/1000;              //
//  --------------------------------------------------------------------

//           intgStep = integral.Simpson(G4DiffElasticHadrNucleus,
//                     G4DiffElasticHadrNucleus::HadrNuclDifferCrSec,
//                      iQ2[ii], iQ2[ii+1], intgN);
                          
                iIntgr[ii+1] = iIntgr[ii]+intgStep; 
                intgrS       = intgrS + intgStep;
              }   //  ii

            for(ii=0; ii<=Nstep; ii++)
             {
               iIntgr[ii] = iIntgr[ii]/intgrS;
//   G4cout <<" Q2 "<< iQ2[ii]/1000/1000<< 
//           " Probability function " << iIntgr[ii]<<G4endl;
              }   //  ii

      
// --------------  The randomization of momentum transfered  -----------

         G4double ranQ2;

         for(G4int kk=1; kk<2; kk++)  //  The cicle for debug
         {
             G4double ranUni = G4UniformRand();

                           // Here must be the solusion of equation
                           // ranUni = iIntgr(Q2);

               G4int curStep = 0;
               for(ii=1; ii<Nstep; ii++)
                  {
//  G4cout << " Prob. func. " << iIntgr[ii]<< 
//                 " Uni "<<ranUni<<G4endl;
                   if(iIntgr[ii] > ranUni) break;
                     curStep += 1;
                  }   //  ii

               ranQ2 = (iQ2[curStep]+(ranUni-iIntgr[curStep])*
                           (iQ2[curStep+1]-iQ2[curStep])
           /(iIntgr[curStep+1]-iIntgr[curStep]))/1000/1000;
// ---------------------------------------------------------------------

//  G4cout<<" Random Uni "<<ranUni<<" Q2 "<<ranQ2<<" curStep "<<
//            curStep<<   " Elastic "<<intgrS<<G4endl;
//  G4cout<<" Diff : "<<ranUni-iIntgr[curStep]<<
//                 " "<< iIntgr[curStep+1]-iIntgr[curStep+1]<<G4endl;

//    G4cout<<ranQ2<<G4endl;
          }    //  kk

             return ranQ2;

    }     
