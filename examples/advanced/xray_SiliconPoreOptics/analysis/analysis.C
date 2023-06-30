#include "TROOT.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TBranch.h"
#include <TFile.h>
#include <TH1D.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <math.h>
#include "TLegend.h"
#include "TApplication.h"
#include "TString.h"
#include <algorithm>
#include "TPie.h"

using namespace std;
TApplication myapp("app",NULL,NULL);

const double pi = 3.14159;
double to_deg(double angle)
{
  return angle * 180 / pi;
}

//////////////////////////////////////////////////////

int main(int argc, char *argv[]){
  // input = ASCII data filename
  string ifilename = "";

  if(argc == 1)
  {
    cout << "No input file selected. Please use: ./read_tree_spectrum [input_file.root]" << endl;
    return 0;
  }
  else
  {
    ifilename = argv[1];
  }
  TString ifilenameforroot = ifilename;

  // create the file, the Tree and a few branches
  TFile *f = new TFile(ifilenameforroot, "READ");
  TTree *mytree = (TTree*)f->Get("XraySPO");  // name of the ntuple in the rootfile

  // create variables to store column values
  int eventID = 0;
  int trackID = 0;
  double x = 0;
  double y = 0;
  double z = 0;
  double theta = 0;
  double phi = 0;
  int parentID = 0;
  char vol_name[500];
  char proc_name[500];
  int num_reflections = 0;

  mytree->SetBranchAddress("eventID", &eventID);
  mytree->SetBranchAddress("vol_name", &vol_name);
  mytree->SetBranchAddress("trackID", &trackID);
  mytree->SetBranchAddress("x", &x);
  mytree->SetBranchAddress("y", &y);
  mytree->SetBranchAddress("z", &z);
  mytree->SetBranchAddress("theta", &theta);
  mytree->SetBranchAddress("phi", &phi);
  mytree->SetBranchAddress("proc_name", &proc_name);
  mytree->SetBranchAddress("parentID", &parentID);
  mytree->SetBranchAddress("num_reflections", &num_reflections);

  int bin_theta = 50;
  int bin_phi = 100;

  int min_val_theta = 0;
  int max_val_theta = 5;
  int min_val_phi = -100;
  int max_val_phi = 100;

  double bin_width_theta = (max_val_theta-min_val_theta)/(double)bin_theta;
  double bin_width_phi = (max_val_phi-min_val_phi)/(double)bin_phi;

  TH1D *angle_theta = new TH1D ("theta", "Efficiency - theta [deg]", bin_theta, min_val_theta, max_val_theta);
  TH1D *angle_phi = new TH1D ("phi", "Efficiency - phi [deg]", bin_phi, min_val_phi, max_val_phi);

  TH1I *reflections = new TH1I ("reflections", "Number of reflections", 8, -1.5, 6.5);

  // Read all entries and fill the histograms
  Long64_t nentries = mytree->GetEntries();
  cout << "Reading " << nentries << " from the TTree." << endl;

  int trackID_old = 0;

  bool passed_entrance = false;
  bool passed_exit = false;
  int eventID_old = 0;

  int check = -999;

  int count_enter = 0;
  int count_exit = 0;

  double theta_exit = 0 ;
  double phi_exit = 0;

  // std::string new_particle_event = "";
  for (Long64_t i=0; i<nentries; i++)
  {
    mytree->GetEntry(i);
    if (eventID != eventID_old || eventID == 0)
    {
      passed_entrance = false;
      passed_exit = false;

      if ((string)vol_name == "pDummyEntrance")
      {
        passed_entrance = true;
        passed_exit = false;
        count_enter += 1;
      }
    }
    else
    {
      if ((string)vol_name == "pDummyExit" && passed_entrance == true)
      {
        // Plot the angle at the
        theta_exit = theta;
        phi_exit = phi;
        passed_exit = true;
      }

      if (passed_exit)
      {
        if ((string)vol_name == "pDummySphere")
        {
          passed_entrance = false;
          passed_exit = false;

          if (num_reflections < 1) num_reflections = -1;
          reflections->Fill(num_reflections);

          // Check if the eventID has not already been counted before
          if (eventID == check) cout << "Evento uguale! " << endl;

          theta_exit = pi-theta_exit;
          if (phi_exit <= 0) phi_exit = -pi-phi_exit;
          else phi_exit = pi-phi_exit;

          if (theta_exit >= min_val_theta && theta_exit <= max_val_theta)
          {
            if (phi_exit >= min_val_phi && phi_exit <= max_val_phi)
            {
              count_exit+=1;
              angle_theta->Fill(to_deg(theta_exit));
              angle_phi->Fill(to_deg(phi_exit));
            }
          }
          check = eventID;
        }
      }
    }
    eventID_old = eventID;
  }

  // pie chart
  Float_t val0 = reflections->GetBinContent(1);  // number of entries with 0 reflections
  Float_t val2 = reflections->GetBinContent(3);  // number of entries with 1 reflection
  Float_t val3 = reflections->GetBinContent(4);
  Float_t val4 = reflections->GetBinContent(5);
  Float_t val5 = reflections->GetBinContent(6) + reflections->GetBinContent(7); // value for 4+ reflections
  cout << "Values are: " << val0 << " " << val2 << " " << val3 << " " << val4 << " " << val5 << " " << endl;


  Float_t vals[] = {val0,val2,val3,val4,val5};
  Int_t colors[] = {0,1,2,3,4};
  Int_t nvals = sizeof(vals)/sizeof(vals[0]);

  TCanvas *cpie = new TCanvas("cpie", "cpie", 0, 0, 900, 900);
  TPie *pie1 = new TPie("pie1", "100 keV, 1deg., SS model",nvals,vals,colors);
  pie1->SetLabelsOffset(.01);
  pie1->SetRadius(.2);

  pie1->SetEntryLabel(0, "0");
  pie1->SetEntryLabel(1, "1");
  pie1->SetEntryLabel(2, "2");
  pie1->SetEntryLabel(3, "3");
  pie1->SetEntryLabel(4, ">3");
  pie1->SetLabelFormat("%txt (%perc)");
  pie1->Draw();

  // Normalization
  cout << "Count particle: <entered: " << count_enter << "> <exited: " << count_exit << ">" << endl;
  angle_theta->Scale(1.0/((double)count_enter));
  angle_phi->Scale(1.0/((double)count_enter));
  angle_theta->Scale(1.0/(bin_width_theta));
  angle_phi->Scale(1.0/(bin_width_phi));

  TCanvas *c1 = new TCanvas("c1", "c1", 0, 50, 900, 900);
  int c1_id = c1->GetCanvasID();
  cout << "The first canvas ID is: " << c1_id << endl;
  c1->SetLogy();
  angle_theta->SetFillColor(kBlack);
  angle_theta->SetMarkerColor(kBlack);
  angle_theta->Draw("");
  angle_theta->GetXaxis()->SetTitle("Theta [deg]");
  angle_theta->GetYaxis()->SetTitle("Normalized efficiency [deg^-1]");
  angle_theta->GetXaxis()->SetTitleSize(0.03);
  angle_theta->GetYaxis()->SetTitleSize(0.03);
  cout << "integral of angle_theta from histo:" << angle_theta->Integral() << " and from fraction: " << count_exit/count_enter << endl;

  TCanvas *c2 = new TCanvas("c2", "c2", 0, 100, 900, 900);
  int c2_id = c2->GetCanvasID();
  cout << "The second canvas ID is: " << c2_id << endl;
  c2->SetLogy();
  angle_phi->SetFillColor(kBlack);
  angle_phi->SetMarkerColor(kBlack);
  angle_phi->Draw("");
  angle_phi->GetXaxis()->SetTitle("Phi [deg]");
  angle_phi->GetYaxis()->SetTitle("Normalized efficiency [deg^-1]");
  angle_phi->GetXaxis()->SetTitleSize(0.03);
  angle_phi->GetYaxis()->SetTitleSize(0.03);
  cout << "integral of angle_phi is:" << angle_theta->Integral() << " and from fraction: " << count_exit/count_enter << endl;

  myapp.Run(true);

  f->Close();
  return 0;
}
