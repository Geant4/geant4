
#include <iostream>
#include <fstream>

#ifdef USEROOT
#include "TFile.h"
#include "TTree.h"
#endif // USEROOT

using namespace std;

int main(int argc, char *argv[])
{
#ifdef USEROOT
  // Verbosity for debugging:
  Int_t verboseLevel = 0;

  Int_t event;
  Int_t bulletType;
  Double_t bulletE;

  Int_t massini, mzini; 
  Double_t exini;
  Int_t mulncasc, mulnevap,mulntot;
  Double_t bimpact;
  Int_t jremn,kfis;
  Double_t estfis;
  Int_t izfis,iafis,ntrack,itypcasc,avv, zvv;
  Double_t  enerj,plab,tetlab,philab;

  Double_t pxlab,pylab,pzlab, momAnglelab;

  Int_t last_event = -1;

  const Int_t maxpart = 255;

  // Variables for Ntuple (TTree)
  Int_t Event;
  Int_t BulletType;
  Double_t BulletE;

  Int_t Massini, Mzini; 
  Double_t Exini;
  Int_t Mulncasc, Mulnevap,Mulntot;

  Double_t Bimpact;
  Int_t Jremn, Kfis;
  Double_t Estfis;
  Int_t Izfis, Iafis;
  Int_t Ntrack;
  Int_t Itypcasc[maxpart], Avv[maxpart], Zvv[maxpart];
  Double_t Enerj[maxpart], Plab[maxpart], Tetlab[maxpart], Philab[maxpart];

  ifstream in;

  if(argc != 3) {
    cout <<"Usage: inclout2root outFileName.out rootFileName.root" << endl;
    return -1;
  }

  // Input filename:
  in.open(argv[1]);
  cout <<"Reading data from file: " << argv[1] << endl;

  // Output filename:
  TFile *datafile = new TFile(argv[2], "RECREATE");

  cout <<"Creating ROOT file: " << argv[2] << endl;

  // Create a new tree:
  TTree *h101 = new TTree("h101", "Output from INCL4 thin target simulation");

  // With these branches (variables)
  h101->Branch("Event", &Event, "Event/I");
  h101->Branch("BulletType", &BulletType, "BulletType/I");
  h101->Branch("BulletE", &BulletE, "BulletE/D");

  h101->Branch("Massini", &Massini, "Massini/I");
  h101->Branch("Mzini", &Mzini, "Mzini/I"); 
  h101->Branch("Exini", &Exini, "Exini/D");
  h101->Branch("Mulncasc", &Mulncasc, "Mulncasc/I");
  h101->Branch("Mulnevap", &Mulnevap, "Mulnevap/I");
  h101->Branch("Mulntot", &Mulntot, "Mulntot/I");

  h101->Branch("Bimpact", &Bimpact, "Bimpact/D");
  h101->Branch("Jremn", &Jremn, "Jremn/I");
  h101->Branch("Kfis", &Kfis, "Kfis/I");
  h101->Branch("Estfis", &Estfis, "Estfis/D");
  h101->Branch("Izfis", &Izfis, "Izfis/I");
  h101->Branch("Iafis", &Iafis, "Iafis/I");

  h101->Branch("Ntrack", &Ntrack, "Ntrack/I");
  h101->Branch("Itypcasc", Itypcasc, "Itypcasc[Ntrack]/I");
  h101->Branch("Avv", Avv, "Avv[Ntrack]/I"); 
  h101->Branch("Zvv", Zvv, "Zvv[Ntrack]/I");
  h101->Branch("Enerj", Enerj, "Enerj[Ntrack]/D");
  h101->Branch("Plab", Plab, "Plab[Ntrack]/D");
  h101->Branch("Tetlab", Tetlab, "Tetlab[Ntrack]/D");
  h101->Branch("Philab", Philab, "Philab[Ntrack]/D");

  Int_t eventNumber;
  Int_t particle;
  Bool_t beginOfConversion = true;

  particle = 0;

  while(1) {
    particle = 0;

    in >> event >> bulletType >> bulletE >> massini >> mzini >> exini >> mulncasc >> mulnevap >> mulntot >> bimpact >> jremn >> kfis >> estfis >> izfis >> iafis >> ntrack >> itypcasc >> avv >> zvv >> enerj >> plab >> tetlab >> philab;

    if(!in.good()) break;

    Event = event;
    BulletType = bulletType;
    BulletE = bulletE;
    Massini = massini;
    Mzini = mzini; 
    Exini = exini;
    Mulncasc = mulncasc;
    Mulnevap = mulnevap;
    Mulntot = mulntot;
    Bimpact = bimpact;
    Jremn = jremn;
    Kfis = kfis;
    Estfis = estfis;
    Izfis = izfis; 
    Iafis = iafis;
    Ntrack = ntrack;
    Itypcasc[particle] = itypcasc;
    Avv[particle] = avv;
    Zvv[particle] = zvv;
    Enerj[particle] = enerj;
    Plab[particle] = plab;
    Tetlab[particle] = tetlab;
    Philab[particle] = philab;

    for(particle = 1; particle < Ntrack; particle++) {
      in >> event >> bulletType >> bulletE >> massini >> mzini >> exini >> mulncasc >> mulnevap >> mulntot >> bimpact >> jremn >> kfis >> estfis >> izfis >> iafis >> ntrack >> itypcasc >> avv >> zvv >> enerj >> plab >> tetlab >> philab;
      Event = event;
      BulletType = bulletType;
      BulletE = bulletE;
      Massini = massini;
      Mzini = mzini; 
      Exini = exini;
      Mulncasc = mulncasc;
      Mulnevap = mulnevap;
      Mulntot = mulntot;
      Bimpact = bimpact;
      Jremn = jremn;
      Kfis = kfis;
      Estfis = estfis;
      Izfis = izfis; 
      Iafis = iafis;
      Ntrack = ntrack;
      Itypcasc[particle] = itypcasc;
      Avv[particle] = avv;
      Zvv[particle] = zvv;
      Enerj[particle] = enerj;
      Plab[particle] = plab;
      Tetlab[particle] = tetlab;
      Philab[particle] = philab;
    }

    h101->Fill();
  }

  datafile->Write();

  in.close();

  cout <<"Done." << endl;

#endif //USEROOT

#ifndef USEROOT // If ROOT headers and libraries not available:
  cout <<"This probram needs ROOT (http://root.cern.ch) libraries to work." << endl;
  cout <<"If you wish to use it please install ROOT and enable ROOT integration in ";
  cout << "the GNUmakefile." << endl;
#endif // USEROOT

  return 0;
}
