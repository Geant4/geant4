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
#include "globals.hh"
#include "G4Version.hh"
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
#include "G4HadronNucleonXsc.hh"                  // Uzhi MK Cross-Sections
#include "G4VComponentCrossSection.hh"            // Uzhi
#include "G4ComponentAntiNuclNuclearXS.hh"        // Uzhi

int main()
{
    G4cout<<"--------- Test Xs for Pbar+A interactions -------------- ";
    G4ParticleTable* partTable = G4ParticleTable::GetParticleTable();
    partTable->SetReadiness();

//  Our AntiNucleus-Nucleus cross sections ----------
    G4VComponentCrossSection* cs =new G4ComponentAntiNuclNuclearXS();

//  Kossov cross sections ---------------------------
    static G4ChipsComponentXS* _instance = new G4ChipsComponentXS();
    G4ChipsComponentXS* CHIPSxsManager = _instance;

    G4String  namePart;
    G4ParticleDefinition* part(0);

//------------------ Define projectile -------------------------------
    namePart = "anti_proton";
    G4cout<<namePart<<G4endl;

    part =G4AntiProton::AntiProton(); 
    cs->BuildPhysicsTable(*part);
// -------------------------------------------------------------------

    G4double cross_sec, cross_secel, cross_inel;
    G4double chipsTot, chipsEl, chipsIn, Plab, energy;

    G4int A(0), Z(0);

//------------------  Define target    -------------------------------
    G4cout<<"                      Pbar+H                            "<<G4endl;
    A = 1;       // Proton
    Z = 1;       // Proton

    std::ofstream XpbarH("XpbarH.dat",std::ios::out);

    XpbarH<<G4Version<<G4endl; /////////////////////////////////////////////////////////////
    XpbarH<<"Plab Tkin Xtotal Xelastic Xinelastic ChipsTot ChipsEl ChipsInel"<<G4endl;

    Plab=90.;
    for(G4int i=1; i<335; i++)
    {
     if(i < 292 )     Plab+= 10.;
     else if(i < 309) Plab+=1000.;
     else if(i < 317) Plab+=10000.;
     else if(i < 326) Plab+=100000.;
     else            Plab+=1000000.;

     energy=std::sqrt(sqr(Plab)+sqr(part->GetPDGMass())) - part->GetPDGMass();

     cross_sec = cs->GetTotalElementCrossSection(part, energy, Z, (G4double) A);
     cross_inel=cs->GetInelasticElementCrossSection(part, energy, Z, (G4double)A);
     cross_secel=cs->GetElasticElementCrossSection(part, energy, Z, (G4double)A);

     chipsTot=CHIPSxsManager->GetTotalElementCrossSection(part,energy,Z,A-Z);
     chipsEl =CHIPSxsManager->GetElasticElementCrossSection(part,energy,Z,A-Z);
     chipsIn =CHIPSxsManager->GetInelasticElementCrossSection(part,energy,Z,A-Z);

     chipsTot /=millibarn; chipsEl /=millibarn; chipsIn /=millibarn;

     XpbarH<<" "<<Plab/GeV<<" "<<energy/GeV
           <<" "<<cross_sec/millibarn<<" "<<cross_secel/millibarn<<" "<<cross_inel/millibarn
           <<" "<<chipsTot           <<" "<<chipsEl              <<" "<<chipsIn<<G4endl;
    }

//------------------  Define target    -------------------------------
    G4cout<<"                      Pbar+H2                           "<<G4endl;
    A = 2;       // H-2
    Z = 1;       // H-2

    std::ofstream XpbarH2("XpbarH2.dat",std::ios::out);

    XpbarH2<<G4Version<<G4endl; /////////////////////////////////////////////////////////////
    XpbarH2<<"Plab Tkin Xtotal Xelastic Xinelastic"<<G4endl;

    Plab=90.;
    for(G4int i=1; i<335; i++)
    {
     if(i < 292 )     Plab+= 10.;
     else if(i < 309) Plab+=1000.;
     else if(i < 317) Plab+=10000.;
     else if(i < 326) Plab+=100000.;
     else            Plab+=1000000.;

     energy=std::sqrt(sqr(Plab)+sqr(part->GetPDGMass())) - part->GetPDGMass();

     cross_sec = cs->GetTotalElementCrossSection(part, energy, Z, (G4double) A);
     cross_inel=cs->GetInelasticElementCrossSection(part, energy, Z, (G4double)A);
     cross_secel=cs->GetElasticElementCrossSection(part, energy, Z, (G4double)A);

//     chipsTot=CHIPSxsManager->GetTotalElementCrossSection(part,energy,Z,A-Z);
//     chipsEl =CHIPSxsManager->GetElasticElementCrossSection(part,energy,Z,A-Z);
//     chipsIn =CHIPSxsManager->GetInelasticElementCrossSection(part,energy,Z,A-Z);

//     chipsTot /=millibarn; chipsEl /=millibarn; chipsIn /=millibarn;

     XpbarH2<<" "<<Plab/GeV<<" "<<energy/GeV
           <<" "<<cross_sec/millibarn<<" "<<cross_secel/millibarn<<" "<<cross_inel/millibarn<<G4endl;
//         <<" "<<chipsTot           <<" "<<chipsEl              <<" "<<chipsIn<<G4endl;
    }

//------------------  Define target    -------------------------------
    G4cout<<"                      Pbar+He                           "<<G4endl;
    A = 4;       // He
    Z = 2;       // He

    std::ofstream XpbarHe("XpbarHe4.dat",std::ios::out);

    XpbarHe<<G4Version<<G4endl; /////////////////////////////////////////////////////////////
    XpbarHe<<"Plab Tkin Xtotal Xelastic Xinelastic"<<G4endl;

    Plab=90.;
    for(G4int i=1; i<335; i++)
    {
     if(i < 292 )     Plab+= 10.;
     else if(i < 309) Plab+=1000.;
     else if(i < 317) Plab+=10000.;
     else if(i < 326) Plab+=100000.;
     else            Plab+=1000000.;

     energy=std::sqrt(sqr(Plab)+sqr(part->GetPDGMass())) - part->GetPDGMass();

     cross_sec = cs->GetTotalElementCrossSection(part, energy, Z, (G4double) A);
     cross_inel=cs->GetInelasticElementCrossSection(part, energy, Z, (G4double)A);
     cross_secel=cs->GetElasticElementCrossSection(part, energy, Z, (G4double)A);

//     chipsTot=CHIPSxsManager->GetTotalElementCrossSection(part,energy,Z,A-Z);
//     chipsEl =CHIPSxsManager->GetElasticElementCrossSection(part,energy,Z,A-Z);
//     chipsIn =CHIPSxsManager->GetInelasticElementCrossSection(part,energy,Z,A-Z);

//     chipsTot /=millibarn; chipsEl /=millibarn; chipsIn /=millibarn;

     XpbarHe<<" "<<Plab/GeV<<" "<<energy/GeV
           <<" "<<cross_sec/millibarn<<" "<<cross_secel/millibarn<<" "<<cross_inel/millibarn<<G4endl;
//         <<" "<<chipsTot           <<" "<<chipsEl              <<" "<<chipsIn<<G4endl;
    }

//------------------  Define target    -------------------------------
    G4cout<<"                      Pbar+Li-6                         "<<G4endl;
    A = 6;       // Li-6
    Z = 3;       // Li-6

    std::ofstream XpbarLi6("XpbarLi6.dat",std::ios::out);

    XpbarLi6<<G4Version<<G4endl; /////////////////////////////////////////////////////////////
    XpbarLi6<<"Plab Tkin Xtotal Xelastic Xinelastic"<<G4endl;

    Plab=90.;
    for(G4int i=1; i<335; i++)
    {
     if(i < 292 )     Plab+= 10.;
     else if(i < 309) Plab+=1000.;
     else if(i < 317) Plab+=10000.;
     else if(i < 326) Plab+=100000.;
     else            Plab+=1000000.;

     energy=std::sqrt(sqr(Plab)+sqr(part->GetPDGMass())) - part->GetPDGMass();

     cross_sec = cs->GetTotalElementCrossSection(part, energy, Z, (G4double) A);
     cross_inel=cs->GetInelasticElementCrossSection(part, energy, Z, (G4double)A);
     cross_secel=cs->GetElasticElementCrossSection(part, energy, Z, (G4double)A);

//     chipsTot=CHIPSxsManager->GetTotalElementCrossSection(part,energy,Z,A-Z);
//     chipsEl =CHIPSxsManager->GetElasticElementCrossSection(part,energy,Z,A-Z);
//     chipsIn =CHIPSxsManager->GetInelasticElementCrossSection(part,energy,Z,A-Z);

//     chipsTot /=millibarn; chipsEl /=millibarn; chipsIn /=millibarn;

     XpbarLi6<<" "<<Plab/GeV<<" "<<energy/GeV
           <<" "<<cross_sec/millibarn<<" "<<cross_secel/millibarn<<" "<<cross_inel/millibarn<<G4endl;
//         <<" "<<chipsTot           <<" "<<chipsEl              <<" "<<chipsIn<<G4endl;
    }

//------------------  Define target    -------------------------------
    G4cout<<"                      Pbar+Be-9                         "<<G4endl;
    A = 9;       // Be-9
    Z = 4;       // Be-9

    std::ofstream XpbarBe9("XpbarBe9.dat",std::ios::out);

    XpbarBe9<<G4Version<<G4endl; /////////////////////////////////////////////////////////////
    XpbarBe9<<"Plab Tkin Xtotal Xelastic Xinelastic"<<G4endl;

    Plab=90.;
    for(G4int i=1; i<335; i++)
    {
     if(i < 292 )     Plab+= 10.;
     else if(i < 309) Plab+=1000.;
     else if(i < 317) Plab+=10000.;
     else if(i < 326) Plab+=100000.;
     else            Plab+=1000000.;

     energy=std::sqrt(sqr(Plab)+sqr(part->GetPDGMass())) - part->GetPDGMass();

     cross_sec = cs->GetTotalElementCrossSection(part, energy, Z, (G4double) A);
     cross_inel=cs->GetInelasticElementCrossSection(part, energy, Z, (G4double)A);
     cross_secel=cs->GetElasticElementCrossSection(part, energy, Z, (G4double)A);

//     chipsTot=CHIPSxsManager->GetTotalElementCrossSection(part,energy,Z,A-Z);
//     chipsEl =CHIPSxsManager->GetElasticElementCrossSection(part,energy,Z,A-Z);
//     chipsIn =CHIPSxsManager->GetInelasticElementCrossSection(part,energy,Z,A-Z);

//     chipsTot /=millibarn; chipsEl /=millibarn; chipsIn /=millibarn;

     XpbarBe9<<" "<<Plab/GeV<<" "<<energy/GeV
           <<" "<<cross_sec/millibarn<<" "<<cross_secel/millibarn<<" "<<cross_inel/millibarn<<G4endl;
//         <<" "<<chipsTot           <<" "<<chipsEl              <<" "<<chipsIn<<G4endl;
    }


//------------------  Define target    -------------------------------
    G4cout<<"                      Pbar+C-12                         "<<G4endl;
    A = 12;       // C-12
    Z =  6;       // C-12

    std::ofstream XpbarC12("XpbarC12.dat",std::ios::out);

    XpbarC12<<G4Version<<G4endl; /////////////////////////////////////////////////////////////
    XpbarC12<<"Plab Tkin Xtotal Xelastic Xinelastic"<<G4endl;

    Plab=90.;
    for(G4int i=1; i<335; i++)
    {
     if(i < 292 )     Plab+= 10.;
     else if(i < 309) Plab+=1000.;
     else if(i < 317) Plab+=10000.;
     else if(i < 326) Plab+=100000.;
     else            Plab+=1000000.;

     energy=std::sqrt(sqr(Plab)+sqr(part->GetPDGMass())) - part->GetPDGMass();

     cross_sec = cs->GetTotalElementCrossSection(part, energy, Z, (G4double) A);
     cross_inel=cs->GetInelasticElementCrossSection(part, energy, Z, (G4double)A);
     cross_secel=cs->GetElasticElementCrossSection(part, energy, Z, (G4double)A);

//     chipsTot=CHIPSxsManager->GetTotalElementCrossSection(part,energy,Z,A-Z);
//     chipsEl =CHIPSxsManager->GetElasticElementCrossSection(part,energy,Z,A-Z);
//     chipsIn =CHIPSxsManager->GetInelasticElementCrossSection(part,energy,Z,A-Z);

//     chipsTot /=millibarn; chipsEl /=millibarn; chipsIn /=millibarn;

     XpbarC12<<" "<<Plab/GeV<<" "<<energy/GeV
           <<" "<<cross_sec/millibarn<<" "<<cross_secel/millibarn<<" "<<cross_inel/millibarn<<G4endl;
//         <<" "<<chipsTot           <<" "<<chipsEl              <<" "<<chipsIn<<G4endl;
    }


//------------------  Define target    -------------------------------
    G4cout<<"                      Pbar+Al-27                        "<<G4endl;
    A = 27;       // Al-27
    Z = 13;       // Al-27

    std::ofstream XpbarAl27("XpbarAl27.dat",std::ios::out);

    XpbarAl27<<G4Version<<G4endl; /////////////////////////////////////////////////////////////
    XpbarAl27<<"Plab Tkin Xtotal Xelastic Xinelastic"<<G4endl;

    Plab=90.;
    for(G4int i=1; i<335; i++)
    {
     if(i < 292 )     Plab+= 10.;
     else if(i < 309) Plab+=1000.;
     else if(i < 317) Plab+=10000.;
     else if(i < 326) Plab+=100000.;
     else            Plab+=1000000.;

     energy=std::sqrt(sqr(Plab)+sqr(part->GetPDGMass())) - part->GetPDGMass();

     cross_sec = cs->GetTotalElementCrossSection(part, energy, Z, (G4double) A);
     cross_inel=cs->GetInelasticElementCrossSection(part, energy, Z, (G4double)A);
     cross_secel=cs->GetElasticElementCrossSection(part, energy, Z, (G4double)A);

//     chipsTot=CHIPSxsManager->GetTotalElementCrossSection(part,energy,Z,A-Z);
//     chipsEl =CHIPSxsManager->GetElasticElementCrossSection(part,energy,Z,A-Z);
//     chipsIn =CHIPSxsManager->GetInelasticElementCrossSection(part,energy,Z,A-Z);

//     chipsTot /=millibarn; chipsEl /=millibarn; chipsIn /=millibarn;

     XpbarAl27<<" "<<Plab/GeV<<" "<<energy/GeV
           <<" "<<cross_sec/millibarn<<" "<<cross_secel/millibarn<<" "<<cross_inel/millibarn<<G4endl;
//         <<" "<<chipsTot           <<" "<<chipsEl              <<" "<<chipsIn<<G4endl;
    }


//------------------  Define target    -------------------------------
    G4cout<<"                      Pbar+Cu-64                        "<<G4endl;
    A = 64;       // Cu-64
    Z = 29;       // Cu-64

    std::ofstream XpbarCu64("XpbarCu64.dat",std::ios::out);

    XpbarCu64<<G4Version<<G4endl; /////////////////////////////////////////////////////////////
    XpbarCu64<<"Plab Tkin Xtotal Xelastic Xinelastic"<<G4endl;

    Plab=90.;
    for(G4int i=1; i<335; i++)
    {
     if(i < 292 )     Plab+= 10.;
     else if(i < 309) Plab+=1000.;
     else if(i < 317) Plab+=10000.;
     else if(i < 326) Plab+=100000.;
     else            Plab+=1000000.;

     energy=std::sqrt(sqr(Plab)+sqr(part->GetPDGMass())) - part->GetPDGMass();

     cross_sec = cs->GetTotalElementCrossSection(part, energy, Z, (G4double) A);
     cross_inel=cs->GetInelasticElementCrossSection(part, energy, Z, (G4double)A);
     cross_secel=cs->GetElasticElementCrossSection(part, energy, Z, (G4double)A);

//     chipsTot=CHIPSxsManager->GetTotalElementCrossSection(part,energy,Z,A-Z);
//     chipsEl =CHIPSxsManager->GetElasticElementCrossSection(part,energy,Z,A-Z);
//     chipsIn =CHIPSxsManager->GetInelasticElementCrossSection(part,energy,Z,A-Z);

//     chipsTot /=millibarn; chipsEl /=millibarn; chipsIn /=millibarn;

     XpbarCu64<<" "<<Plab/GeV<<" "<<energy/GeV
           <<" "<<cross_sec/millibarn<<" "<<cross_secel/millibarn<<" "<<cross_inel/millibarn<<G4endl;
//         <<" "<<chipsTot           <<" "<<chipsEl              <<" "<<chipsIn<<G4endl;
    }


//------------------  Define target    -------------------------------
    G4cout<<"                      Pbar+Sn-119                       "<<G4endl;
    A = 119;       // Sn-119
    Z =  50;       // Sn-119

    std::ofstream XpbarSn("XpbarSn119.dat",std::ios::out);

    XpbarSn<<G4Version<<G4endl; /////////////////////////////////////////////////////////////
    XpbarSn<<"Plab Tkin Xtotal Xelastic Xinelastic"<<G4endl;

    Plab=90.;
    for(G4int i=1; i<335; i++)
    {
     if(i < 292 )     Plab+= 10.;
     else if(i < 309) Plab+=1000.;
     else if(i < 317) Plab+=10000.;
     else if(i < 326) Plab+=100000.;
     else            Plab+=1000000.;

     energy=std::sqrt(sqr(Plab)+sqr(part->GetPDGMass())) - part->GetPDGMass();

     cross_sec = cs->GetTotalElementCrossSection(part, energy, Z, (G4double) A);
     cross_inel=cs->GetInelasticElementCrossSection(part, energy, Z, (G4double)A);
     cross_secel=cs->GetElasticElementCrossSection(part, energy, Z, (G4double)A);

//     chipsTot=CHIPSxsManager->GetTotalElementCrossSection(part,energy,Z,A-Z);
//     chipsEl =CHIPSxsManager->GetElasticElementCrossSection(part,energy,Z,A-Z);
//     chipsIn =CHIPSxsManager->GetInelasticElementCrossSection(part,energy,Z,A-Z);

//     chipsTot /=millibarn; chipsEl /=millibarn; chipsIn /=millibarn;

     XpbarSn<<" "<<Plab/GeV<<" "<<energy/GeV
           <<" "<<cross_sec/millibarn<<" "<<cross_secel/millibarn<<" "<<cross_inel/millibarn<<G4endl;
//         <<" "<<chipsTot           <<" "<<chipsEl              <<" "<<chipsIn<<G4endl;
    }

//------------------  Define target    -------------------------------
    G4cout<<"                      Pbar+Pb-208                       "<<G4endl;
    A = 208;       // Pb-208
    Z =  82;       // Pb-208

    std::ofstream XpbarPb("XpbarPb208.dat",std::ios::out);

    XpbarPb<<G4Version<<G4endl; /////////////////////////////////////////////////////////////
    XpbarPb<<"Plab Tkin Xtotal Xelastic Xinelastic"<<G4endl;

    Plab=90.;
    for(G4int i=1; i<335; i++)
    {
     if(i < 292 )     Plab+= 10.;
     else if(i < 309) Plab+=1000.;
     else if(i < 317) Plab+=10000.;
     else if(i < 326) Plab+=100000.;
     else            Plab+=1000000.;

     energy=std::sqrt(sqr(Plab)+sqr(part->GetPDGMass())) - part->GetPDGMass();

     cross_sec = cs->GetTotalElementCrossSection(part, energy, Z, (G4double) A);
     cross_inel=cs->GetInelasticElementCrossSection(part, energy, Z, (G4double)A);
     cross_secel=cs->GetElasticElementCrossSection(part, energy, Z, (G4double)A);

//     chipsTot=CHIPSxsManager->GetTotalElementCrossSection(part,energy,Z,A-Z);
//     chipsEl =CHIPSxsManager->GetElasticElementCrossSection(part,energy,Z,A-Z);
//     chipsIn =CHIPSxsManager->GetInelasticElementCrossSection(part,energy,Z,A-Z);

//     chipsTot /=millibarn; chipsEl /=millibarn; chipsIn /=millibarn;

     XpbarPb<<" "<<Plab/GeV<<" "<<energy/GeV
           <<" "<<cross_sec/millibarn<<" "<<cross_secel/millibarn<<" "<<cross_inel/millibarn<<G4endl;
//         <<" "<<chipsTot           <<" "<<chipsEl              <<" "<<chipsIn<<G4endl;
    }

//------------------  Define target    -------------------------------
    G4cout<<"                      Pbar+U-238                        "<<G4endl;
    A = 238;       // U-238
    Z =  92;       // U-238

    std::ofstream XpbarU("XpbarU238.dat",std::ios::out);

    XpbarU<<G4Version<<G4endl; /////////////////////////////////////////////////////////////
    XpbarU<<"Plab Tkin Xtotal Xelastic Xinelastic"<<G4endl;

    Plab=90.;
    for(G4int i=1; i<335; i++)
    {
     if(i < 292 )     Plab+= 10.;
     else if(i < 309) Plab+=1000.;
     else if(i < 317) Plab+=10000.;
     else if(i < 326) Plab+=100000.;
     else            Plab+=1000000.;

     energy=std::sqrt(sqr(Plab)+sqr(part->GetPDGMass())) - part->GetPDGMass();

     cross_sec = cs->GetTotalElementCrossSection(part, energy, Z, (G4double) A);
     cross_inel=cs->GetInelasticElementCrossSection(part, energy, Z, (G4double)A);
     cross_secel=cs->GetElasticElementCrossSection(part, energy, Z, (G4double)A);

//     chipsTot=CHIPSxsManager->GetTotalElementCrossSection(part,energy,Z,A-Z);
//     chipsEl =CHIPSxsManager->GetElasticElementCrossSection(part,energy,Z,A-Z);
//     chipsIn =CHIPSxsManager->GetInelasticElementCrossSection(part,energy,Z,A-Z);

//     chipsTot /=millibarn; chipsEl /=millibarn; chipsIn /=millibarn;

     XpbarU<<" "<<Plab/GeV<<" "<<energy/GeV
           <<" "<<cross_sec/millibarn<<" "<<cross_secel/millibarn<<" "<<cross_inel/millibarn<<G4endl;
//         <<" "<<chipsTot           <<" "<<chipsEl              <<" "<<chipsIn<<G4endl;
    }


    partTable->DeleteAllParticles();
    G4cout << "###### End of test #####" << G4endl;
}
