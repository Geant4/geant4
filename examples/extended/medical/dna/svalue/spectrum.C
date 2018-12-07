// *********************************************************************
// To execute this macro under ROOT after your simulation ended, 
//   1 - launch ROOT (usually type 'root' at your machine's prompt)
//   2 - type '.X spectrum.C' at the ROOT session prompt
//   3 - OR type directly 'root spectrum.C'
// this will generate an output energy spectrum spectrum.txt, 
// according to selected statistics and distribution
// (energy unit assumed by svalue example: eV)
// *********************************************************************

{
gROOT->Reset();

// 1) SELECT number of values to shoot
int nbValues = 100;
//

double value;
TNtuple * ntuple = new TNtuple ("","","value");

FILE * fp = fopen("spectrum.txt","w");

TRandom r;

for (int i=0;i<nbValues;i++) 
{
// 2) SELECT distribution
//    here Gaussian with mean, sigma
//    see https://root.cern.ch/doc/master/classTRandom.html
  value = r.Gaus(1.,0.1);
//
  ntuple->Fill(value);
  fprintf(fp,"%f \n",value);
}
fclose(fp);

ntuple->Draw("value");
}
