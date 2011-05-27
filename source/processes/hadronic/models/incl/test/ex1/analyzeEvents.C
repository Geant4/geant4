#include "Riostream.h";

void analyzeEvents() {

   gROOT->Reset();
  gROOT->SetStyle("clearRetro");



  enum particleType { nuclei = 0, proton = 1, neutron = 2, pionPlus = 3, pionMinus = 5, pionZero = 7, foton = 10 };

  float MeV=0.001;
  ifstream in;
  in.open("data/data.out");

  int runId, eventId, particleId, modelId,fragmentA, fragmentZ;
  float kineticEnergy, momX, momY, momZ, exitationEnergy;
  int nlines = 0;
  TFile *f = new TFile("data/analyzeEvents.root","RECREATE");

  TH1F *hpE = new TH1F("hpE","p kinetic energy",100,0,1);
  TH1F *hnE = new TH1F("hnE","n kinetic energy",100,0,1);
  TH1F *hgE = new TH1F("hgE","gamma energy",100,0,1);

  TH1F *hppE = new TH1F("hppE","pi+ energy",100,0,1);

  
  TNtuple *ntuple = new TNtuple("ntuple","foorified data from cascade.cc  ouput","runId:eventId:particleId:modelId:kineticEnergy:momX:momY:momZ:momAngle:fragmentA:fragmentZ:exitationEnergy:coulombOK");
  // read first comment line and second line with parameters

  char commentLine[80], executable[15];
  //   fgets(&line,80,in);

  int runId, nCollisions, bulletType, targetA, targetZ;  
  int coulombOK;
  float mMax;
  in >> runId >> nCollisions >> bulletType >>  mMax >> targetA >> targetZ;

  cout << "runId        : " << runId       << endl; 
  cout << "# collisions : " << nCollisions << endl;  
  cout << "bullet type  : " << bulletType  << endl; 
  cout << "bullet Zmom  : " << mMax  << " MeV"<< endl; 
  cout << "target A     : " << targetA     << endl; 
  cout << "target Z     : " << targetZ     <<endl;

  double eMax = mMax;

  while (1) {
    in >> runId >> eventId >> particleId >> modelId >> kineticEnergy >> momX >> momY >> momZ >> fragmentA >> fragmentZ >> exitationEnergy >> coulombOK;
    if (!in.good()) break;
    //    if (nlines < 5) printf("x=%8f, y=%8f, z=%8f\n",runId, eventId, particleId);
    if (particleId==proton) hpE->Fill(kineticEnergy/eMax);
    if (particleId==neutron) hnE->Fill(kineticEnergy/eMax);
    if (particleId==foton) hgE->Fill(kineticEnergy/eMax);
    if (particleId==pionPlus) hppE->Fill(kineticEnergy/eMax);
    double t=sqrt(mMax*938.27/1000);
    double GeV=1000;
    double momAngle = acos(momZ/(sqrt(momX*momX + momY*momY + momZ*momZ)))*(180.0/3.14);
    ntuple->Fill(runId, eventId, particleId, modelId, kineticEnergy, momX, momY, momZ, momAngle, fragmentA, fragmentZ, exitationEnergy, coulombOK);

    nlines++;
  };
  printf(" found %d lines \n",nlines);

  
  in.close();
  f->Write();
  if (gROOT->IsBatch()) return;
}
