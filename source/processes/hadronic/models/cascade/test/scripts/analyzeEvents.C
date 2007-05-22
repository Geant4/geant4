{
   gROOT->Reset();
  gROOT->SetStyle("clearRetro");

#include "Riostream.h";

  enum particleType { nuclei = 0, proton = 1, neutron = 2, pionPlus = 3, pionMinus = 5, pionZero = 7, foton = 10 };

  float MeV=0.001;
  ifstream in;
  in.open("../data/data.out");

  int runId, eventId, particleId, fragmentA, fragmentZ;
  float kineticEnergy, momX, momY, momZ, exitationEnergy;
  int nlines = 0;
  TFile *f = new TFile("analyzeEvents.root","RECREATE");

  TH1F *hpE = new TH1F("hpE","p kinetic energy",100,0,1);
  TH1F *hnE = new TH1F("hnE","n kinetic energy",100,0,1);
  TH1F *hgE = new TH1F("hgE","gamma energy",100,0,1);

  TH1F *hppE = new TH1F("hppE","pi+ energy",100,0,1);

  
  TNtuple *ntuple = new TNtuple("ntuple","foorified data from cascade.cc  ouput","runId:eventId:particleId:kineticEnergy:momX:momY:momZ:fragmentA:fragmentZ:exitationEnergy");
  // read first comment line and second line with parameters

  char commentLine[80], executable[15];
  //   fgets(&line,80,in);

  int runId, nCollisions, bulletType, targetA, targetZ;  
  float bulletMomZ;
  in >> runId >> nCollisions >> bulletType >> bulletMomZ >> targetA >> targetZ;

  cout << "runId        : " << runId       << endl; 
  cout << "# collisions : " << nCollisions << endl;  
  cout << "bullet type  : " << bulletType  << endl; 
  cout << "bullet Zmom  : " << bulletMomZ  << " MeV"<< endl; 
  cout << "target A     : " << targetA     << endl; 
  cout << "target Z     : " << targetZ     <<endl;

  double eMax = bulletMomZ;

  while (1) {
    in >> runId >> eventId >> particleId >> kineticEnergy >> momX >> momY >> momZ >> fragmentA >> fragmentZ >> exitationEnergy;
    if (!in.good()) break;
    //    if (nlines < 5) printf("x=%8f, y=%8f, z=%8f\n",runId, eventId, particleId);
    if (particleId==proton) hpE->Fill(kineticEnergy/eMax);
    if (particleId==neutron) hnE->Fill(kineticEnergy/eMax);
    if (particleId==foton) hgE->Fill(kineticEnergy/eMax);
    if (particleId==pionPlus) hppE->Fill(kineticEnergy/eMax);

    ntuple->Fill(runId, eventId, particleId, kineticEnergy/eMax, momX/eMax, momY/eMax, momZ/eMax, fragmentA, fragmentZ, exitationEnergy/eMax);
    nlines++;
  };
  printf(" found %d lines \n",nlines);

  
  in.close();
  f->Write();

}
