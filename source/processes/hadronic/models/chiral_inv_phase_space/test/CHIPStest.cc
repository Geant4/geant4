///////////////////////////////////////////////////////////////////////////////
// File: CHIPStest.cc
// Author: M.Kossov
// Last modification: 08/99
// Description: Main function for the Test of CHIPS generator for Geant4
///////////////////////////////////////////////////////////////////////////////

//#define debug
#define pdebug

#include "G4UIterminal.hh"
#include "G4ios.hh"
#include "G4ProcessManager.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

//#ifdef XX
//  #include "YY.hh"
//#endif

#include "G4Quasmon.hh"

#include "g4std/iostream"
#include "g4std/fstream"
#include "g4std/iomanip"

//int main(int argc, char** argv) 
int main()
{
  //if (argc == 1) 
  //// Interactive mode
  //{
  //  G4UIsession* session = new G4UIterminal;
  //  session->SessionStart();
  //  delete session;
  //}
  ////Batch mode
  //else {
  //  G4String command = "/control/execute ";
  //  G4String fileName = argv[1];
  //  UI->ApplyCommand(command+fileName);    
  //}

  //-------- set output format-------
   G4cout.setf( G4std::ios::scientific, G4std::ios::floatfield );

  //-------- write results onto a file --------
   G4std::ofstream outFile( "chipstest.out", G4std::ios::out);
   outFile.setf( G4std::ios::scientific, G4std::ios::floatfield );

  //--------- Example how to write in the outFile -----------
  //outFile << " " << G4endl;

  //@@... Temporary: initialize the G4ParticleTable (particle definition procedure)
  G4ParticleTable::GetParticleTable()->SetVerboseLevel(2);

  G4ParticleDefinition gam("gam",0.,0.,0.,0,0,0,0,0,0,"gauge",0,0,22,true,0.,NULL);
  G4ParticleDefinition eta("eta",547.3,.00118,0.,0,0,0,0,0,0,"meson",0,0,221,true,0.,NULL);
  G4ParticleDefinition piZ("piZ",134.9764,0.,0.,0,0,0,2,0,0,"meson",0,0,111,true,0.,NULL);
  G4ParticleDefinition piP("piP",139.56995,0.,1.,0,0,0,2,2,0,"meson",0,0,211,true,0.,NULL);
  G4ParticleDefinition piM("piM",139.56995,0.,-1.,0,0,0,2,-2,0,"meson",0,0,-211,true,0.,NULL);
  G4ParticleDefinition kaZ("kaZ",497.672,0.,0.,0,0,0,1,-1,0,"meson",0,0,311,true,0.,NULL);
  G4ParticleDefinition kaM("kaM",493.677,0.,1.,0,0,0,1,1,0,"meson",0,0,321,true,0.,NULL);
  G4ParticleDefinition akZ("akZ",497.672,0.,0.,0,0,0,1,-1,0,"meson",0,0,-311,true,0.,NULL);
  G4ParticleDefinition akP("akP",493.677,0.,-1.,0,0,0,1,1,0,"meson",0,0,-321,true,0.,NULL);
  G4ParticleDefinition etP("etP",957.78,.203,0.,0,0,0,0,0,0,"meson",0,0,331,false,0.,NULL);
  G4ParticleDefinition omega("omega",781.94,8.41,0.,2,0,0,0,0,0,"meson",0,0,223,false,0.,NULL);
  G4ParticleDefinition rhoZ("rhoZ",770.,150.7,0.,2,0,0,2,0,0,"meson",0,0,113,false,0.,NULL);
  G4ParticleDefinition rhoP("rhoP",770.,150.7,1.,2,0,0,2,1,0,"meson",0,0,213,false,0.,NULL);
  G4ParticleDefinition rhoM("rhoM",770.,150.7,-1.,2,0,0,2,-1,0,"meson",0,0,-213,false,0.,NULL);
  G4ParticleDefinition kSZ("kSZ",896.1,50.5,0.,2,0,0,1,-1,0,"meson",0,0,313,false,0.,NULL);
  G4ParticleDefinition kSC("kSC",891.66,50.8,1.,2,0,0,1,1,0,"meson",0,0,323,false,0.,NULL);
  G4ParticleDefinition akSZ("akSZ",896.1,50.5,0.,2,0,0,1,-1,0,"meson",0,0,-313,false,0.,NULL);
  G4ParticleDefinition akSC("akSC",891.66,50.8,-1.,2,0,0,1,1,0,"meson",0,0,-323,false,0.,NULL);
  G4ParticleDefinition phi("phi",1019.413,4.43,0.,2,0,0,0,0,0,"meson",0,0,333,false,0.,NULL);
  G4ParticleDefinition apr("apr",938.27231,0.,-1.,1,0,0,1,1,0,"baryon",0,-1,-2212,true,0.,NULL);
  //                       Name,M,Width,Charge,2S,P,C,2IS,2ISZ,G,"Type",LepN,BarN,code,stab,T,*)

  G4ParticleDefinition* proton=G4ParticleTable::GetParticleTable()->FindParticle(2212);
  G4ParticleDefinition* antiproton=G4ParticleTable::GetParticleTable()->FindParticle(-2212);
  G4LorentzVector prot4mom(0.,0.,0.,proton->GetPDGMass());
  G4LorentzVector apro4mom(0.,0.,0.,antiproton->GetPDGMass());  
  G4int rPDG=113;
  G4double maxM=1500.;
  for (G4int jr=0; jr<10000; jr++)
  {
    G4QHadron curRes(rPDG, maxM);
    G4LorentzVector tmp4Mom = curRes.Get4Momentum();
    outFile << tmp4Mom.m() << G4endl;
  }
  for (G4int ir=0; ir<0; ir++)
  {
    G4LorentzVector totSum = prot4mom + apro4mom;
#ifdef debug
    cout << "Main: p4mom =" << prot4mom << ", ap4mom =" << apro4mom << G4endl;
#endif
    G4Quasmon* pan = new G4Quasmon(2212, -2212, prot4mom, apro4mom);
#ifdef debug
    cout << "===>>> Now call HadronizeQuasmon hadronization function" << G4endl;
#endif
    G4QHadronVector output = pan->HadronizeQuasmon();
    G4int totNOfHadrons = output.entries();
#ifdef pdebug
    cout<<"DONE^^^^^^^^^^^^:ir="<<ir<<": A#of generated hadrons ="<<totNOfHadrons<<G4endl;
#endif
    for (G4int index=0; index<totNOfHadrons; index++)
    {
      G4double m = output[index]->GetMass();
      G4LorentzVector lorV = output[index]->Get4Momentum();
      totSum    -= lorV;
#ifdef pdebug
      cout<<"Hadron #"<<index<<", m="<<m <<", LV="<<lorV<<", M="<<totSum.m()<<G4endl;
#endif
    }
#ifdef pdebug
    cout << "CHECK: Lorentz Vector sum:" << totSum << G4endl;
#endif
    delete pan;
  }

  return EXIT_SUCCESS;
}



