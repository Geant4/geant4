
#include "Riostream.h"
void rcsv() {

  ifstream in;
  in.open("out.csv");

  unsigned int index;
  double rgauss,rbw;

  unsigned int nlines = 0;
  TH1D* hrg = new TH1D("hrg","rgauss distribution",100,-4,4);

  char sep;
  while(true) {
    in >> index >> sep >> rgauss >> sep >> rbw;
    if(!in.good()) break;
    if(nlines<5) printf("index=%d, rgauss=%g, rbw=%g\n",index,rgauss,rbw);
    hrg->Fill(rgauss);
    nlines++;
  }

  in.close();

  /////////////////////////////////////////////////
  /// plotting : //////////////////////////////////
  /////////////////////////////////////////////////
  TCanvas* plotter = new TCanvas("canvas","",10,10,800,600);
  plotter->Divide(1,1);  
  plotter->cd(1);
  hrg->Draw();
  plotter->Update();

    
}
