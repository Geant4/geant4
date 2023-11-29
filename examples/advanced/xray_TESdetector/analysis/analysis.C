#include <TROOT.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TString.h>
#include <TBranch.h>
#include <TStyle.h>
#include <TPie.h>
#include <TFile.h>
#include <TLegend.h>
#include <TH1D.h>
#include <TH2D.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <algorithm>

using namespace std;
TApplication myapp("app",NULL,NULL);

//////////////////////////////////////////////////////

int main(int argc, char *argv[]){
  string ifilename = "";
  int no_events = 10000000; // number of primaries

  if(argc == 1)
  {
      cout << "No input file selected. Please use: ./read_tree_spectrum [input_file.root]" << endl;
      return 0;
  }
  else
  {
      ifilename = argv[1];
      if (argc == 3) no_events = atoi(argv[2]);
  }

  TString ifilenameforroot = ifilename;
  TFile *f = new TFile(ifilenameforroot, "READ");

  int no_pixel = 317;

  // Histogram definition
  //  Detector deposits
  TH1D *pixel_dep = new TH1D ("pixel", "Avg. energy deposit per pixel [MeV]", no_pixel, 0, no_pixel);  //(nbinx, xdown, xup, nbiny, ylow, yup) for drawing the hexagon
  TH1D *part_per_pixel = new TH1D ("Part per pixel", "No. particles per pixel", no_pixel, 0, no_pixel);

  // Spectra for energies (initial, incident and deposited)
  int bin_number = 500;
  int max_val = 10000;
  int min_val = 0;
  double bin_width = (max_val-min_val)/(double)bin_number;
  TH1D *initial_energies = new TH1D ("energies4", "GCR Protons impacting spectra on the detector", bin_number, min_val, max_val);
  TH1D *incident_energies = new TH1D ("energies2", " ", bin_number, min_val, max_val);
  TH1D *dep_energies = new TH1D ("energies3", " ", bin_number, min_val, max_val);

  // Spectra for energies (initial, incident and deposited) - for primaries (protons)
  TH1D *initial_energies_P = new TH1D ("energies4_P", " ", bin_number, min_val, max_val);
  TH1D *incident_energies_P = new TH1D ("energies2_P", " ", bin_number, min_val, max_val);
  TH1D *dep_energies_P = new TH1D ("energies3_P", "Primary GCR Protons impacting spectra on the detector", bin_number, min_val, max_val);

  // Incident particles
  TH1D *inc_particles = new TH1D ("energies6", "GCR Protons particle fluxes on the detector composition", 5, 0, 5);

  // Variables to store the tuples' values
  int eventID = 0;
  char vol_name[500];
  int trackID = 0;
  double x = 0;
  double y = 0;
  double z = 0;
  double theta = 0;
  double phi = 0;
  int parentID = 0;
  int pixel_number = 0;
  double step_energy_dep = 0;
  int step_number = 0;
  double init_kinetic_energy = 0;
  double kinetic_energy = 0;
  char particle_name[500];
  char pre_step_name[500];
  char post_step_name[500];

  // Define tuple elements
  TTree *mytree = (TTree*)f->Get("TES_Tuple");
  Long64_t nentries = mytree->GetEntries();
  mytree->SetBranchAddress("eventID", &eventID);
  mytree->SetBranchAddress("vol_name", &vol_name);
  mytree->SetBranchAddress("trackID", &trackID);
  mytree->SetBranchAddress("x", &x);
  mytree->SetBranchAddress("y", &y);
  mytree->SetBranchAddress("z", &z);
  mytree->SetBranchAddress("theta", &theta);
  mytree->SetBranchAddress("phi", &phi);
  mytree->SetBranchAddress("parentID", &parentID);
  mytree->SetBranchAddress("pixel_number", &pixel_number);
  mytree->SetBranchAddress("step_energy_dep", &step_energy_dep);
  mytree->SetBranchAddress("step_number", &step_number);
  mytree->SetBranchAddress("init_kinetic_energy", &init_kinetic_energy);
  mytree->SetBranchAddress("kinetic_energy", &kinetic_energy);
  mytree->SetBranchAddress("particle_name", &particle_name);
  mytree->SetBranchAddress("pre_step_name", &pre_step_name);
  mytree->SetBranchAddress("post_step_name", &post_step_name);

  cout << "Reading " << nentries << " from the TTree." << endl;
  string vol_name_old = "";

  int eventID_old = 0;
  int trackID_old = 0;
  double total_step_dep = 0;
  double total_step_dep_P = 0;
  double kine = 0;        // kinetic energy
  double kine_P = 0;      // kinetic energy for protons
  double in_kine = 0;     // initial kinetic energy
  double in_kine_P = 0;   // initial kinetic energy for protons
  bool entered = false;
  bool entered_P = false;
  char arr[50];
  string str = "";
  int pixel_number_old = 0;

  for (Long64_t i=0; i<nentries; i++)
  {
    mytree->GetEntry(i);

    if (entered || entered_P)
    {
      if (eventID != eventID_old || (eventID == eventID_old && trackID != trackID_old))
      {
        if (entered && total_step_dep > 0)
        {
          initial_energies->Fill(in_kine);
          incident_energies->Fill(kine);
          dep_energies->Fill(total_step_dep);
          inc_particles->Fill(arr, 1);
          total_step_dep = 0;
        }
        if (entered_P && total_step_dep_P > 0)
        {
          initial_energies_P->Fill(in_kine_P);
          incident_energies_P->Fill(kine_P);
          dep_energies_P->Fill(total_step_dep_P);
          total_step_dep_P = 0;
        }
        entered = false;
        entered_P = false;
      }
    }
    str = (string)particle_name;

    if((string)vol_name == "Bipxl")
    {
      if (vol_name_old != "Bipxl")
      {
        if (!entered)
        {
          entered = true;
          kine = kinetic_energy;
          in_kine = init_kinetic_energy;

          if (str != "proton" && str != "gamma" && str != "e+" && str != "e-")
          {
            str = "Other particles";
          }
          strcpy(arr, str.c_str());

          // Same but onyl for primaries
          if (!entered_P)
          {
            if (str == "proton" && parentID == 0)
            {
              entered_P = true;
              kine_P = kinetic_energy;
              in_kine_P = init_kinetic_energy;
            }
          }
        }
      }

      // Sum the total energy deposit from all particles
      total_step_dep += step_energy_dep;

      // Sum the total energy deposit from the primaries
      if (str == "proton" && parentID == 0) total_step_dep_P += step_energy_dep;
      // hist_det_count->Fill(x, y, step_energy_dep);

      if (pixel_number != pixel_number_old)
      {
        part_per_pixel->Fill(-pixel_number);
      }

      pixel_dep->Fill(-pixel_number, step_energy_dep);
    }
    eventID_old = eventID;
    trackID_old = trackID;
    vol_name_old = (string)vol_name;
    pixel_number_old = pixel_number;
  }

  cout << "Starts plotting..." << endl;

  // Draw histograms
  TCanvas *c1 = new TCanvas ("c1", "Energy dep", 0, 0, 1000, 900);
  double radius = 27.0;
  double pi = 3.1415;
  double T = no_events/(0.407*4*pi*pi*radius*radius);     //equivalent time
  double S = 2.1; //cm2
  dep_energies->Scale(1.0/(T*S*bin_width));               //scale based on the equivalent time and detector surface
  initial_energies->Scale(1.0/(T*S*bin_width));
  incident_energies->Scale(1.0/(T*S*bin_width));
  gStyle->SetOptStat(0);
  c1->SetLogx();
  c1->SetLogy();
  dep_energies->SetLineColor(kRed);
  initial_energies->Draw("hist");
  initial_energies->SetLineColor(kBlack);
  incident_energies->Draw("histsame");
  dep_energies->Draw("histsame");
  incident_energies->SetLineColor(kBlue);
  initial_energies->GetYaxis()->SetRangeUser(1.0e-6, 1.0);
  initial_energies->GetXaxis()->SetTitle("Energy [MeV]");
  initial_energies->GetYaxis()->SetTitle("Counts/cm2/s/MeV");
  auto legend1 = new TLegend(0.75,0.8,0.9,0.9);
  legend1->AddEntry(dep_energies,"Edep");
  legend1->AddEntry(incident_energies, "Einc");
  legend1->AddEntry(initial_energies, "Ei");
  legend1->Draw("same");

  TCanvas *c2 = new TCanvas ("c2", "Proton energy dep", 50, 0, 1000, 900);
  dep_energies_P->Scale(1.0/(T*S*bin_width));
  initial_energies_P->Scale(1.0/(T*S*bin_width));
  incident_energies_P->Scale(1.0/(T*S*bin_width));
  gStyle->SetOptStat(0);
  c2->SetLogx();
  c2->SetLogy();
  dep_energies_P->Draw("hist");
  dep_energies_P->SetLineColor(kRed);
  initial_energies_P->Draw("histsame");
  initial_energies_P->SetLineColor(kBlack);
  incident_energies_P->Draw("histsame");
  incident_energies_P->SetLineColor(kBlue);
  dep_energies_P->GetYaxis()->SetRangeUser(1.0e-6, 1.0);
  dep_energies_P->GetXaxis()->SetTitle("Energy [MeV]");
  dep_energies_P->GetYaxis()->SetTitle("Counts/cm2/s/MeV");
  auto legend2 = new TLegend(0.75,0.8,0.9,0.9);
  legend2->AddEntry(dep_energies,"Edep");
  legend2->AddEntry(incident_energies, "Einc");
  legend2->AddEntry(initial_energies, "Ei");
  legend2->Draw("same");

  TCanvas *cpie = new TCanvas("Particless distribution", "Particles distribution", 100, 0, 1000, 1000);
  string label1 = inc_particles->GetXaxis()->GetBinLabel(1);
  string label2 = inc_particles->GetXaxis()->GetBinLabel(2);
  string label3 = inc_particles->GetXaxis()->GetBinLabel(3);
  string label4 = inc_particles->GetXaxis()->GetBinLabel(4);
  string label5 = inc_particles->GetXaxis()->GetBinLabel(5);
  Float_t val1 = inc_particles->GetBinContent(1);
  Float_t val2 = inc_particles->GetBinContent(2);
  Float_t val3 = inc_particles->GetBinContent(3);
  Float_t val4 = inc_particles->GetBinContent(4);
  Float_t val5 = inc_particles->GetBinContent(5);
  Float_t vals[] = {val1,val2,val3,val4,val5};
  Int_t colors[] = {1,2,3,4,5};
  Int_t nvals = sizeof(vals)/sizeof(vals[0]);
  TPie *pie1 = new TPie("pie1", "Particles distribution",nvals,vals,colors);
  pie1->SetLabelsOffset(.01);
  pie1->SetRadius(.2);
  pie1->SetEntryLabel(0, label1.c_str());
  pie1->SetEntryLabel(1, label2.c_str());
  pie1->SetEntryLabel(2, label3.c_str());
  pie1->SetEntryLabel(3, label4.c_str());
  pie1->SetEntryLabel(4, label5.c_str());
  pie1->SetLabelFormat("%txt (%perc)");
  pie1->Draw();

  TCanvas *c4_num = new TCanvas ("c4_num", "Particle count per pixel", 150, 0, 1000, 900);
  part_per_pixel->GetXaxis()->SetTitle("Pixel number");
  part_per_pixel->GetYaxis()->SetTitle("Counts");
  part_per_pixel->SetFillColor(kBlack);
  part_per_pixel->SetMarkerColor(kBlack);
  part_per_pixel->Draw();

  // Draw the histogram of the single pixel deposition
  TCanvas *c5 = new TCanvas ("c5", "Avg. energy deposit per pixel", 250, 0, 1000, 900);
  pixel_dep->Divide(part_per_pixel);
  pixel_dep->GetXaxis()->SetTitle("Pixel number");
  pixel_dep->GetYaxis()->SetTitle("Avg. energy deposit [MeV]");
  pixel_dep->SetFillColor(kBlack);
  pixel_dep->SetMarkerColor(kBlack);
  pixel_dep->Draw("");

  myapp.Run(true);
  f->Close();
  return 0;
}