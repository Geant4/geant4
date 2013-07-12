#include "globals.hh"
#include "G4ios.hh"
#include <fstream>
#include <iomanip>

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Proton.hh"
#include "G4PionMinus.hh"
#include "G4PionPlus.hh"
#include "G4KaonMinus.hh"
#include "G4KaonPlus.hh"
#include "G4AntiProton.hh"

#include "G4ChipsComponentXS.hh"                  // Uzhi 29.01.13

int main()
{
    G4cout<<"--------- Test CHIPS Xs for pp, pi-p, pi+p, k-p, k+p and Pbar P -------"<<G4endl;
    G4ParticleTable* partTable = G4ParticleTable::GetParticleTable();
    partTable->SetReadiness();

//  Kossov cross sections ---------------------------
    static G4ChipsComponentXS* _instance = new G4ChipsComponentXS();
    G4ChipsComponentXS* CHIPSxsManager = _instance;

    G4String  namePart;
    G4ParticleDefinition* part(0);
    G4double chipsTot, chipsEl, chipsIn, Plab, energy;

//------------------  Define target    -------------------------------
    G4int A = 1;       // Proton
    G4int Z = 1;       // Proton

//------------------ Define projectile -------------------------------
    namePart = "proton";
    G4cout<<"                      "<<namePart<<G4endl;

    std::ofstream ppCX("ppCX.dat",std::ios::out);
    part = G4Proton::Proton(); 

    ppCX<<"Plab Tkin Xtotal Xelastic Xinelastic "<<namePart<<G4endl;

    Plab=90.;
    for(G4int i=1; i<335; i++)
    {
     if(i < 292 )     Plab+= 10.;
     else if(i < 309) Plab+=1000.;
     else if(i < 317) Plab+=10000.;
     else if(i < 326) Plab+=100000.;
     else            Plab+=1000000.;

     energy=std::sqrt(sqr(Plab)+sqr(part->GetPDGMass())) - part->GetPDGMass();

     chipsTot=CHIPSxsManager->GetTotalElementCrossSection(part,energy,Z,A-Z);
     chipsEl =CHIPSxsManager->GetElasticElementCrossSection(part,energy,Z,A-Z);
     chipsIn =CHIPSxsManager->GetInelasticElementCrossSection(part,energy,Z,A-Z);

     chipsTot /=millibarn; chipsEl /=millibarn; chipsIn /=millibarn;

     ppCX<<" "<<Plab/GeV<<" "<<energy/GeV<<" "
         <<chipsTot<<" "<<chipsEl<<" "<<chipsIn<<G4endl;
    }

//------------------ Define projectile -------------------------------
    namePart = "pi-";
    G4cout<<"                      "<<namePart<<G4endl;

    std::ofstream pimpCX("pimpCX.dat",std::ios::out);
    part =G4PionMinus::PionMinus(); 

    pimpCX<<"Plab Tkin Xtotal Xelastic Xinelastic "<<namePart<<G4endl;

    Plab=90.;
    for(G4int i=1; i<335; i++)
    {
     if(i < 292 )     Plab+= 10.;
     else if(i < 309) Plab+=1000.;
     else if(i < 317) Plab+=10000.;
     else if(i < 326) Plab+=100000.;
     else            Plab+=1000000.;

     energy=std::sqrt(sqr(Plab)+sqr(part->GetPDGMass())) - part->GetPDGMass();

     chipsTot=CHIPSxsManager->GetTotalElementCrossSection(part,energy,Z,A-Z);
     chipsEl =CHIPSxsManager->GetElasticElementCrossSection(part,energy,Z,A-Z);
     chipsIn =CHIPSxsManager->GetInelasticElementCrossSection(part,energy,Z,A-Z);

     chipsTot /=millibarn; chipsEl /=millibarn; chipsIn /=millibarn;

     pimpCX<<" "<<Plab/GeV<<" "<<energy/GeV<<" "
           <<chipsTot<<" "<<chipsEl<<" "<<chipsIn<<G4endl;
    }

//------------------ Define projectile -------------------------------
    namePart = "pi+";
    G4cout<<"                      "<<namePart<<G4endl;

    std::ofstream pippCX("pippCX.dat",std::ios::out);
    part =G4PionPlus::PionPlus(); 

    pippCX<<"Plab Tkin Xtotal Xelastic Xinelastic "<<namePart<<G4endl;

    Plab=90.;
    for(G4int i=1; i<335; i++)
    {
     if(i < 292 )     Plab+= 10.;
     else if(i < 309) Plab+=1000.;
     else if(i < 317) Plab+=10000.;
     else if(i < 326) Plab+=100000.;
     else            Plab+=1000000.;

     energy=std::sqrt(sqr(Plab)+sqr(part->GetPDGMass())) - part->GetPDGMass();

     chipsTot=CHIPSxsManager->GetTotalElementCrossSection(part,energy,Z,A-Z);
     chipsEl =CHIPSxsManager->GetElasticElementCrossSection(part,energy,Z,A-Z);
     chipsIn =CHIPSxsManager->GetInelasticElementCrossSection(part,energy,Z,A-Z);

     chipsTot /=millibarn; chipsEl /=millibarn; chipsIn /=millibarn;

     pippCX<<" "<<Plab/GeV<<" "<<energy/GeV<<" "
           <<chipsTot<<" "<<chipsEl<<" "<<chipsIn<<G4endl;
    }

//------------------ Define projectile -------------------------------
    namePart = "kaon-";
    G4cout<<"                      "<<namePart<<G4endl;

    std::ofstream kmpCX("kmpCX.dat",std::ios::out);
    part =G4KaonMinus::KaonMinus(); 

    kmpCX<<"Plab Tkin Xtotal Xelastic Xinelastic "<<namePart<<G4endl;

    Plab=90.;
    for(G4int i=1; i<335; i++)
    {
     if(i < 292 )     Plab+= 10.;
     else if(i < 309) Plab+=1000.;
     else if(i < 317) Plab+=10000.;
     else if(i < 326) Plab+=100000.;
     else            Plab+=1000000.;

     energy=std::sqrt(sqr(Plab)+sqr(part->GetPDGMass())) - part->GetPDGMass();

     chipsTot=CHIPSxsManager->GetTotalElementCrossSection(part,energy,Z,A-Z);
     chipsEl =CHIPSxsManager->GetElasticElementCrossSection(part,energy,Z,A-Z);
     chipsIn =CHIPSxsManager->GetInelasticElementCrossSection(part,energy,Z,A-Z);

     chipsTot /=millibarn; chipsEl /=millibarn; chipsIn /=millibarn;

     kmpCX<<" "<<Plab/GeV<<" "<<energy/GeV<<" "
          <<chipsTot<<" "<<chipsEl<<" "<<chipsIn<<G4endl;
    }

//------------------ Define projectile -------------------------------
    namePart = "kaon+";
    G4cout<<"                      "<<namePart<<G4endl;

    std::ofstream kppCX("kppCX.dat",std::ios::out);
    part =G4KaonPlus::KaonPlus(); 

    kppCX<<"Plab Tkin Xtotal Xelastic Xinelastic "<<namePart<<G4endl;

    Plab=90.;
    for(G4int i=1; i<335; i++)
    {
     if(i < 292 )     Plab+= 10.;
     else if(i < 309) Plab+=1000.;
     else if(i < 317) Plab+=10000.;
     else if(i < 326) Plab+=100000.;
     else            Plab+=1000000.;

     energy=std::sqrt(sqr(Plab)+sqr(part->GetPDGMass())) - part->GetPDGMass();

     chipsTot=CHIPSxsManager->GetTotalElementCrossSection(part,energy,Z,A-Z);
     chipsEl =CHIPSxsManager->GetElasticElementCrossSection(part,energy,Z,A-Z);
     chipsIn =CHIPSxsManager->GetInelasticElementCrossSection(part,energy,Z,A-Z);

     chipsTot /=millibarn; chipsEl /=millibarn; chipsIn /=millibarn;

     kppCX<<" "<<Plab/GeV<<" "<<energy/GeV<<" "
          <<chipsTot<<" "<<chipsEl<<" "<<chipsIn<<G4endl;
    }

//------------------ Define projectile -------------------------------
    namePart = "anti_proton";
    G4cout<<"                      "<<namePart<<G4endl;

    std::ofstream antippCX("antippCX.dat",std::ios::out);
    part =G4AntiProton::AntiProton(); 

    antippCX<<"Plab Tkin Xtotal Xelastic Xinelastic "<<namePart<<G4endl;

    Plab=90.;
    for(G4int i=1; i<335; i++)
    {
     if(i < 292 )     Plab+= 10.;
     else if(i < 309) Plab+=1000.;
     else if(i < 317) Plab+=10000.;
     else if(i < 326) Plab+=100000.;
     else            Plab+=1000000.;

     energy=std::sqrt(sqr(Plab)+sqr(part->GetPDGMass())) - part->GetPDGMass();

     chipsTot=CHIPSxsManager->GetTotalElementCrossSection(part,energy,Z,A-Z);
     chipsEl =CHIPSxsManager->GetElasticElementCrossSection(part,energy,Z,A-Z);
     chipsIn =CHIPSxsManager->GetInelasticElementCrossSection(part,energy,Z,A-Z);

     chipsTot /=millibarn; chipsEl /=millibarn; chipsIn /=millibarn;

     antippCX<<" "<<Plab/GeV<<" "<<energy/GeV<<" "
              <<chipsTot<<" "<<chipsEl<<" "<<chipsIn<<G4endl;
    }

    partTable->DeleteAllParticles();
    G4cout << "###### End of test #####" << G4endl;
}
