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
//               FTF test: p+C interactions; Inclusive NA61/SHINE
//
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//
//      File name:     Test30
//
//      Author:        V.Ivanchenko
//
//      Creation date: 12 March 2002
//
//      Modifications:
//      14.11.03 Renamed to cascade
//      09.05.06 Return back to test30
// -------------------------------------------------------------------
#include "globals.hh"
#include "G4Version.hh"
#include "G4ios.hh"
#include <fstream>
#include <iomanip>

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "FTFtest1.icc"
#include "G4ChipsComponentXS.hh"                  // Uzhi 29.01.13

#include "UZHI_diffraction.hh"

int main(int argc, char** argv)
{
  CLHEP::RanluxEngine defaultEngine( 1234567, 4 );
  CLHEP::HepRandom::setTheEngine( &defaultEngine );

  G4cout << "========================================================" << G4endl;
  G4cout << "======              FTF Test Start              ========" << G4endl;
  G4cout << "========================================================" << G4endl;

  // -------------------------------------------------------------------
  // Control on input

  if(argc < 2) {
    G4cout << "Input file is not specified! Exit" << G4endl;
    exit(1);
  }

  std::ifstream* fin = new std::ifstream();

  G4String fname = argv[1];
  fin->open(fname.c_str());
  if( !fin->is_open()) {
    G4cout << "Input file <" << fname << "> does not exist! Exit" << G4endl;
    exit(1);
  }

//-----------------------------------------------------------------------
  #include "FTFtest2.icc"   // Initialization
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
//G4double sigTot = 0; 
//G4double sigEl  = 0;
G4double sigIn  = 0;

//-------------------------- Global histograms  -------------------------
std::ofstream Pip("Pip.dat",std::ios::out);
std::ofstream Pim("Pim.dat",std::ios::out);
std::ofstream Pro("Pro.dat",std::ios::out);
std::ofstream Kp("Kp.dat",std::ios::out);

//G4double SqrtS;

G4double Pip1[10][10];  for(G4int jj=0; jj<10; jj++) {for(G4int ii=0; ii< 10; ii++){Pip1[jj][ii]=0.;}}
G4double Pip2[10][21];  for(G4int jj=0; jj<10; jj++) {for(G4int ii=0; ii< 21; ii++){Pip2[jj][ii]=0.;}} 
G4double Pip3[10][12];  for(G4int jj=0; jj<10; jj++) {for(G4int ii=0; ii< 12; ii++){Pip3[jj][ii]=0.;}}
G4double Pip4[10][12];  for(G4int jj=0; jj<10; jj++) {for(G4int ii=0; ii< 12; ii++){Pip4[jj][ii]=0.;}}

G4double Pim1[10][10];  for(G4int jj=0; jj<10; jj++) {for(G4int ii=0; ii< 10; ii++){Pim1[jj][ii]=0.;}}
G4double Pim2[10][21];  for(G4int jj=0; jj<10; jj++) {for(G4int ii=0; ii< 21; ii++){Pim2[jj][ii]=0.;}} 
G4double Pim3[10][12];  for(G4int jj=0; jj<10; jj++) {for(G4int ii=0; ii< 12; ii++){Pim3[jj][ii]=0.;}}
G4double Pim4[10][12];  for(G4int jj=0; jj<10; jj++) {for(G4int ii=0; ii< 12; ii++){Pim4[jj][ii]=0.;}}

G4double Pro1[10][10];  for(G4int jj=0; jj<10; jj++) {for(G4int ii=0; ii< 10; ii++){Pro1[jj][ii]=0.;}}
G4double Pro2[10][21];  for(G4int jj=0; jj<10; jj++) {for(G4int ii=0; ii< 21; ii++){Pro2[jj][ii]=0.;}} 
G4double Pro3[10][12];  for(G4int jj=0; jj<10; jj++) {for(G4int ii=0; ii< 12; ii++){Pro3[jj][ii]=0.;}}
G4double Pro4[10][19];  for(G4int jj=0; jj<10; jj++) {for(G4int ii=0; ii< 19; ii++){Pro4[jj][ii]=0.;}}

G4double Kp1[2][10];    for(G4int jj=0; jj< 2; jj++) {for(G4int ii=0; ii< 10; ii++){Kp1[jj][ii]=0.;}}
G4double Kp2[2][21];    for(G4int jj=0; jj< 2; jj++) {for(G4int ii=0; ii< 21; ii++){Kp2[jj][ii]=0.;}} 
G4double Kp3[2][12];    for(G4int jj=0; jj< 2; jj++) {for(G4int ii=0; ii< 12; ii++){Kp3[jj][ii]=0.;}}
G4double Kp4[2][12];    for(G4int jj=0; jj< 2; jj++) {for(G4int ii=0; ii< 12; ii++){Kp4[jj][ii]=0.;}}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // -------- Loop over run

  G4String line, line1;
  G4bool end = true;

  for(G4int run=0; run<100; run++) {
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//-------------------------- Current histograms -------------------------

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  do {
    #include "FTFtest3.icc"  // -------- Read input file
    #include "FTFtest4.icc"  // -------- Start run processing  
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

    G4cout << "cross(mb)in= " << cross_sec*1000./barn << G4endl
           << "cross(mb)el= " << cross_secel*1000./barn<<G4endl<<G4endl;

    cross_inel=cross_sec-cross_secel; // +++++++++++++++++++++++++

    cross_sec/=millibarn;   // Inel Cross section in mb
    cross_secel/=millibarn; // Elas Cross section in mb
    cross_inel/=millibarn;  // Inel Cross section in mb

    G4cout<<"Element A Z N: "<<A<<" "<<Z<<" "<<A-Z<<G4endl;
    G4cout<<"Proposed Xs (mb): Tot El In: "
          <<cross_sec<<" "<<cross_secel<<" "<<cross_inel<<G4endl;

//---------------------------------------------------------------------------
// Kossov cross sections      ---------------------------
    G4double chipsTot, chipsEl, chipsIn;

    static G4ChipsComponentXS* _instance = new G4ChipsComponentXS();
    G4ChipsComponentXS* CHIPSxsManager = _instance;

    G4bool CHIPapplic=true;                          //false;   Uzhi 29.01.13
    if(CHIPapplic)
    {
     chipsTot=CHIPSxsManager->GetTotalElementCrossSection(part,energy,Z,A-Z);
     chipsEl =CHIPSxsManager->GetElasticElementCrossSection(part,energy,Z,A-Z);
     chipsIn =CHIPSxsManager->GetInelasticElementCrossSection(part,energy,Z,A-Z);
     chipsTot/=millibarn; chipsEl/=millibarn; chipsIn/=millibarn;

     G4cout<<"CHIPS cross sections are used:----------------------"<<G4endl<<
             "Plab          Total        Elastic      Inelastic"   <<G4endl;
     G4cout<<" "<<Plab/GeV<<" "<< chipsTot<<" "<<chipsEl<<" "<<chipsIn <<G4endl<<G4endl;

     /*sigTot=chipsTot; sigEl=chipsEl;*/ sigIn=chipsIn;
    } else
    {
     //sigTot = cross_sec; 
     //sigEl  = cross_secel;
     sigIn  = cross_inel;

     G4cout<<"Proposed Xs (mb) are used: Tot El In: "
           <<cross_sec<<" "<<cross_secel<<" "<<cross_inel<<G4endl;
    }

//+++++++++++++++++++++++++++++++++ For each energy +++++++++++++++++++++
   //G4double E=energy+part->GetPDGMass();
   //SqrtS=std::sqrt(sqr(part->GetPDGMass())+sqr(938.)+2.*938.*E);       
   //SqrtS=SqrtS;
//   Ybeam=0.5*std::log((E+Plab)/(E-Plab));

//  G4int Ntotal=nevt;
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++       
//-------------------------------------------------------

    const G4DynamicParticle* sec = 0;
    G4ParticleDefinition* pd;
    G4ThreeVector  mom;
    G4LorentzVector labv, fm;
    G4double e, /* px, py, pz, pt, pt2,*/ theta;
    G4VParticleChange* aChange = 0;

//  G4double E=energy+part->GetPDGMass();                                  // Elab Proj
//  G4double SqrtS=std::sqrt(sqr(part->GetPDGMass())+sqr(938.)+2.*938.*E); // per  Proj+N
//  G4double Ycms=0.5*std::log((E+Plab)/(E-Plab));                         //      Proj+N

    // -------- Event loop
    G4cout<<"Events start "<<nevt<<G4endl;
    G4int Ninelast=nevt;                              // Uzhi
//=================================================================
    for (G4int iter=0; iter<nevt; ++iter) {
//=================================================================
      if(verbose > 0) G4cout<<"Start events loop***********************"<<G4endl;

      if(verbose>=1 || iter == modu*(iter/modu)) { 
        G4cout << "### " << iter << "-th event start " <<Plab/GeV<<G4endl;
      }

      if(saverand) {defaultEngine.saveStatus("initial.conf");}

      G4double e0 = energy;
      do {
        if(sigmae > 0.0) e0 = G4RandGauss::shoot(energy,sigmae);
      } while (e0 < 0.0);

      dParticle.SetKineticEnergy(e0);

      gTrack->SetStep(step);
      gTrack->SetKineticEnergy(e0);
      G4double amass = phys->GetNucleusMass();
      // note: check of 4-momentum balance for CHIPS is not guranteed due to
      // unknown isotope      
      aChange = proc->PostStepDoIt(*gTrack,*step); 

      G4double mass = part->GetPDGMass();

      if ( ionParticle ) 
      {
       e0/=ionA; G4double mass_N=938.*MeV;                              // Init 4-mom
       labv = G4LorentzVector(0.0, 0.0, std::sqrt(e0*(e0 + 2.*mass_N)), //   NN
	 		      e0 + mass_N +  mass_N); 
      } else
      {                                             
       labv = G4LorentzVector(0.0, 0.0, std::sqrt(e0*(e0 + 2.*mass)),   //   hA
  		              e0 + mass + amass); 
      }

      G4ThreeVector bst = labv.boostVector();          // To CMS NN in AA or hA
//------------
      G4LorentzVector labNN(0.0, 0.0, std::sqrt(e0*(e0 + 2.*mass)),e0 + mass + amass);
      G4ThreeVector boostNN = labNN.boostVector();
//------------

      // take into account local energy deposit
      G4double de = aChange->GetLocalEnergyDeposit();
      G4LorentzVector dee = G4LorentzVector(0.0, 0.0, 0.0, de); 
      labv -= dee;

      G4int n = aChange->GetNumberOfSecondaries();     // Multiplicity of prod. part.

      if(verbose>=1) G4cout<<" Uzhi ------------ N prod. part "<<n<<G4endl;
//++++++++++++++++ Variables for each event +++++++++++++++++++++++++++++
      if((verbose > 0) && (n < 2)) {G4cout<<"Multiplicity of produced < 2!!!"<<G4endl;}

//    G4int Ncharged=0;                          // Uzhi
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      for(G4int i=0; i<n; ++i)              // Loop over produced particles
      {
        sec = aChange->GetSecondary(i)->GetDynamicParticle();
        pd  = sec->GetDefinition();
        G4String pname=pd->GetParticleName();

        if(verbose>=1) G4cout<<" Part  "<<i<<" "<<pname
                             <<" "<<sec->Get4Momentum()/GeV
                             <<sec->Get4Momentum().mag()/GeV<<G4endl;

        fm  = sec->Get4Momentum();
//        fm.boost(-bst);
//        G4double rapidity=fm.rapidity();

        mom = fm.vect();
//      G4double mas = pd->GetPDGMass();
//	G4double p = mom.mag();

        labv -= fm;   // For checking energy-momentum conservation

	// electron can come only from internal conversion
	// its mass should be added to initial state
        if(pd == electron) { 

	  labv += G4LorentzVector(0.0,0.0,0.0,electron_mass_c2); 
	}

        //px = mom.x();
        //py = mom.y();
        //pz = mom.z(); pz=pz;
        G4double Pmod=mom.mag();
        G4double p=Pmod/GeV;;
        //pt = std::sqrt(px*px +py*py); pt2=sqr(pt/GeV); pt2=pt2;
        e  = fm.e(); // - m;
        theta = mom.theta();

//        G4double Xplus=(e+pz)/(E+Plab);
//        G4double CosTheta=std::cos(theta);
//        theta=theta*180./pi;
//        G4double costcm = std::cos(fm.theta());

G4int Itheta(-1); 
if(( 0.   < theta) && (theta <= 0.02)) {Itheta=0;}
if(( 0.02 < theta) && (theta <= 0.04)) {Itheta=1;}
if(( 0.04 < theta) && (theta <= 0.06)) {Itheta=2;}
if(( 0.06 < theta) && (theta <= 0.10)) {Itheta=3;}
if(( 0.10 < theta) && (theta <= 0.14)) {Itheta=4;}
if(( 0.14 < theta) && (theta <= 0.18)) {Itheta=5;}
if(( 0.18 < theta) && (theta <= 0.24)) {Itheta=6;}
if(( 0.24 < theta) && (theta <= 0.30)) {Itheta=7;}
if(( 0.30 < theta) && (theta <= 0.36)) {Itheta=8;}
if(( 0.36 < theta) && (theta <= 0.42)) {Itheta=9;}

G4int IthetaKp(-1);
if(( 0.02 < theta) && (theta <= 0.14)) {IthetaKp=0;}
if(( 0.14 < theta) && (theta <= 0.24)) {IthetaKp=1;}

//+++++++++++++++++ For each particle in the event ++++++++++++++++++++++
if(Itheta >= 0)
{
	 if  (pname == "proton")     
         {
               if(            (p <  1.  )){ G4int Ip=G4int((p- 0. )/0.1); Pro1[Itheta][Ip]++;}
          else if((1. < p) && (p <  5.2 )){ G4int Ip=G4int((p- 1. )/0.2); Pro2[Itheta][Ip]++;}
          else if((5.2< p) && (p < 10.  )){ G4int Ip=G4int((p- 5.2)/0.4); Pro3[Itheta][Ip]++;}
          else if((10.< p) && (p < 25.10)){ G4int Ip=G4int((p-10. )/0.8); Pro4[Itheta][Ip]++;};
         };

	 if ( pname == "pi+" ) 
         {
               if(            (p <  1.  )){ G4int Ip=G4int((p- 0. )/0.1); Pip1[Itheta][Ip]++;}
          else if((1. < p) && (p <  5.2 )){ G4int Ip=G4int((p- 1. )/0.2); Pip2[Itheta][Ip]++;}
          else if((5.2< p) && (p < 10.  )){ G4int Ip=G4int((p- 5.2)/0.4); Pip3[Itheta][Ip]++;}
          else if((10.< p) && (p < 19.59)){ G4int Ip=G4int((p-10. )/0.8); Pip4[Itheta][Ip]++;};
         };  

	 if ( pname == "pi-" ) 
         {
               if(            (p <  1.  )){ G4int Ip=G4int((p- 0. )/0.1); Pim1[Itheta][Ip]++;}
          else if((1. < p) && (p <  5.2 )){ G4int Ip=G4int((p- 1. )/0.2); Pim2[Itheta][Ip]++;}
          else if((5.2< p) && (p < 10.  )){ G4int Ip=G4int((p- 5.2)/0.4); Pim3[Itheta][Ip]++;}
          else if((10.< p) && (p < 19.59)){ G4int Ip=G4int((p-10. )/0.8); Pim4[Itheta][Ip]++;};
         };  

	 if ( (pname == "kaon+") && (IthetaKp >= 0) ) 
         {
               if(            (p <  1.  )){ G4int Ip=G4int((p- 0. )/0.1); Kp1[IthetaKp][Ip]++;}
          else if((1. < p) && (p <  5.2 )){ G4int Ip=G4int((p- 1. )/0.2); Kp2[IthetaKp][Ip]++;}
          else if((5.2< p) && (p < 10.  )){ G4int Ip=G4int((p- 5.2)/0.4); Kp3[IthetaKp][Ip]++;}
          else if((10.< p) && (p < 19.59)){ G4int Ip=G4int((p-10. )/0.8); Kp4[IthetaKp][Ip]++;};
         };

}  // end of if(Itheta >= 0)
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	de += e;

	//	delete sec;       	 
        delete aChange->GetSecondary(i);

      } //     end of the loop on particles
//+++++++++++++++++ Store after each event ++++++++++++++++++++++++++++++

      if(n == 0) {Ninelast--; G4cout<<"n=0 !!! "<<G4endl;}
      if(n == 1) {Ninelast--; G4cout<<"n=1 !!! "<<G4endl;}
      if(n == 2) {Ninelast--; G4cout<<"n=2 !!! "<<G4endl;}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      if(verbose > 0) 
        G4cout << "Energy/Momentum balance= " << labv << G4endl;

      aChange->Clear();

      if(verbose > 0) 
      { 
       G4cout << "End event =====================================" <<Plab<< G4endl; // Uzhi
       G4int Uzhi_i;                                                        // Uzhi
       G4cin >> Uzhi_i;                                                     // Uzhi
      }
//

    }   // End of the event loop ------------------------------------

    timer->Stop();
    G4cout << "  "  << *timer << G4endl;
    delete timer;

//++++++++++++++++++++++ After each energy run ++++++++++++++++++++++++++
G4cout<< "***********************************************************"<< G4endl;

G4cout<<"nevt Ninel "<<nevt<<" "<<Ninelast<<G4endl;
G4cout<<"Plab "<<Plab/GeV<<" SigIn "<<sigIn<<G4endl;
//-------------------------------------------------------------------


    if(verbose > 0) {
      G4cout << "###### End of run # " << run << "     ######" << G4endl;
    }
//++++++++++++++++++++++ After each energy run ++++++++++++++++++++++++++ Uzhi 

//sigTot=sigTot; sigEl=sigEl;

//----------------------------- Protons distributions------------------// Uzhi ++++
G4cout<< "******** Prot ******* at Plab "<<Plab/GeV<<" Xin " << sigIn<< G4endl;
Pro<<G4Version<<G4endl;
Pro<<"P S0_20 S20_40 S40_60 S60_100 S100_140 S140_180 S180_240 S240_300 S300_360 S360_400"<< G4endl;

G4double Pc=-0.05;
for(G4int ii=0; ii < 10; ii++)
{
 Pc+=0.1;
 Pro<<Pc<<" ";
 for(G4int jj=0;jj <10; jj++)
 {                                          
  Pro1[jj][ii]=1000.*Pro1[jj][ii]/Ninelast/0.1; Pro<<Pro1[jj][ii]<<" ";
 }
 Pro<<G4endl;
}

Pc=0.9;
for(G4int ii=0; ii <21; ii++)
{
 Pc+=0.2;
 Pro<<Pc<<" ";
 for(G4int jj=0;jj <10; jj++)
 {                                          
  Pro2[jj][ii]=1000.*Pro2[jj][ii]/Ninelast/0.2; Pro<<Pro2[jj][ii]<<" ";
 }
 Pro<<G4endl;
}                                                                           

Pc=5.0;
for(G4int ii=0; ii <12; ii++)
{
 Pc+=0.4;
 Pro<<Pc<<" ";
 for(G4int jj=0;jj <10; jj++)
 {                                          
  Pro3[jj][ii]=1000.*Pro3[jj][ii]/Ninelast/0.4; Pro<<Pro3[jj][ii]<<" ";
 }
 Pro<<G4endl;
}    

Pc=9.6;
for(G4int ii=0; ii <19; ii++)
{
 Pc+=0.8;
 Pro<<Pc<<" ";
 for(G4int jj=0;jj <10; jj++)
 {                                          
  Pro4[jj][ii]=1000.*Pro4[jj][ii]/Ninelast/0.8; Pro<<Pro4[jj][ii]<<" ";
 }
 Pro<<G4endl;
}    

//----------------------------- Pi+ distributions------------------// Uzhi ++++
G4cout<< "******** Pi+ ******* at Plab "<<Plab/GeV<<" Xin " << sigIn<< G4endl;
Pip<<G4Version<<G4endl;
Pip<<"P S0_20 S20_40 S40_60 S60_100 S100_140 S140_180 S180_240 S240_300 S300_360 S360_400"<< G4endl;

Pc=-0.05;
for(G4int ii=0; ii < 10; ii++)
{
 Pc+=0.1;
 Pip<<Pc<<" ";
 for(G4int jj=0;jj <10; jj++)
 {                                          
  Pip1[jj][ii]=sigIn*Pip1[jj][ii]/Ninelast/0.1; Pip<<Pip1[jj][ii]<<" ";
 }
 Pip<<G4endl;
}

Pc=0.9;
for(G4int ii=0; ii <21; ii++)
{
 Pc+=0.2;
 Pip<<Pc<<" ";
 for(G4int jj=0;jj <10; jj++)
 {                                          
  Pip2[jj][ii]=sigIn*Pip2[jj][ii]/Ninelast/0.2; Pip<<Pip2[jj][ii]<<" ";
 }
 Pip<<G4endl;
}                                                                           

Pc=5.0;
for(G4int ii=0; ii <12; ii++)
{
 Pc+=0.4;
 Pip<<Pc<<" ";
 for(G4int jj=0;jj <10; jj++)
 {                                          
  Pip3[jj][ii]=sigIn*Pip3[jj][ii]/Ninelast/0.4; Pip<<Pip3[jj][ii]<<" ";
 }
 Pip<<G4endl;
}    

Pc=9.6;
for(G4int ii=0; ii <12; ii++)
{
 Pc+=0.8;
 Pip<<Pc<<" ";
 for(G4int jj=0;jj <10; jj++)
 {                                          
  Pip4[jj][ii]=sigIn*Pip4[jj][ii]/Ninelast/0.8; Pip<<Pip4[jj][ii]<<" ";
 }
 Pip<<G4endl;
}    

//----------------------------- Pi- distributions------------------// Uzhi ++++
G4cout<< "******** Pi- ******* at Plab "<<Plab/GeV<<" Xin " << sigIn<< G4endl;
Pim<<G4Version<<G4endl;
Pim<<"P S0_20 S20_40 S40_60 S60_100 S100_140 S140_180 S180_240 S240_300 S300_360 S360_400"<< G4endl;

Pc=-0.05;
for(G4int ii=0; ii < 10; ii++)
{
 Pc+=0.1;
 Pim<<Pc<<" ";
 for(G4int jj=0;jj <10; jj++)
 {                                          
  Pim1[jj][ii]=sigIn*Pim1[jj][ii]/Ninelast/0.1; Pim<<Pim1[jj][ii]<<" ";
 }
 Pim<<G4endl;
}

Pc=0.9;
for(G4int ii=0; ii <21; ii++)
{
 Pc+=0.2;
 Pim<<Pc<<" ";
 for(G4int jj=0;jj <10; jj++)
 {                                          
  Pim2[jj][ii]=sigIn*Pim2[jj][ii]/Ninelast/0.2; Pim<<Pim2[jj][ii]<<" ";
 }
 Pim<<G4endl;
}                                                                           

Pc=5.0;
for(G4int ii=0; ii <12; ii++)
{
 Pc+=0.4;
 Pim<<Pc<<" ";
 for(G4int jj=0;jj <10; jj++)
 {                                          
  Pim3[jj][ii]=sigIn*Pim3[jj][ii]/Ninelast/0.4; Pim<<Pim3[jj][ii]<<" ";
 }
 Pim<<G4endl;
}    

Pc=9.6;
for(G4int ii=0; ii <12; ii++)
{
 Pc+=0.8;
 Pim<<Pc<<" ";
 for(G4int jj=0;jj <10; jj++)
 {                                          
  Pim4[jj][ii]=sigIn*Pim4[jj][ii]/Ninelast/0.8; Pim<<Pim4[jj][ii]<<" ";
 }
 Pim<<G4endl;
}    

//----------------------------- K+ distributions------------------// Uzhi ++++
G4cout<< "******** K+ ******* at Plab "<<Plab/GeV<<" Xin " << sigIn<< G4endl;
Kp<<G4Version<<G4endl;
Kp<<"P S20_140 S140_240"<< G4endl;

Pc=-0.05;
for(G4int ii=0; ii < 10; ii++)
{
 Pc+=0.1;
 Kp<<Pc<<" ";
 for(G4int jj=0;jj <2; jj++)
 {                                          
  Kp1[jj][ii]=sigIn*Kp1[jj][ii]/Ninelast/0.1; Kp<<Kp1[jj][ii]<<" ";
 }
 Kp<<G4endl;
}

Pc=0.9;
for(G4int ii=0; ii <21; ii++)
{
 Pc+=0.2;
 Kp<<Pc<<" ";
 for(G4int jj=0;jj <2; jj++)
 {                                          
  Kp2[jj][ii]=sigIn*Kp2[jj][ii]/Ninelast/0.2; Kp<<Kp2[jj][ii]<<" ";
 }
 Kp<<G4endl;
}                                                                           

Pc=5.0;
for(G4int ii=0; ii <12; ii++)
{
 Pc+=0.4;
 Kp<<Pc<<" ";
 for(G4int jj=0;jj <2; jj++)
 {                                          
  Kp3[jj][ii]=sigIn*Kp3[jj][ii]/Ninelast/0.4; Kp<<Kp3[jj][ii]<<" ";
 }
 Kp<<G4endl;
}    

Pc=9.6;
for(G4int ii=0; ii <12; ii++)
{
 Pc+=0.8;
 Kp<<Pc<<" ";
 for(G4int jj=0;jj <2; jj++)
 {                                          
  Kp4[jj][ii]=sigIn*Kp4[jj][ii]/Ninelast/0.8; Kp<<Kp4[jj][ii]<<" ";
 }
 Kp<<G4endl;
}    
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ Uzhi 
//    G4cerr << "###### End of run # " << run << "     ######" << G4endl;

  } while(end);
  }  // End of job ----------------------------------------------
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  delete pFrame;
//  delete lFrame;
//  delete sFrame;
  delete mate;
  delete fin;
  delete phys;
  partTable->DeleteAllParticles();

  G4cout << "###### End of test #####" << G4endl;
}
