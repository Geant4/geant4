//*********************************************************************************
// Modified by Sara Zein to calculate the damage probability per plasmid
//_________________________________________________________________________________
//*********************************************************************************

// This macro requires the molecular-dna.root file generated from molecularDNA example
// To run this file just insert this command to the terminal:
// root .X plasmid.C

//***************************************//
// Please define the parameters below    //
// ifile, r3, Nbp (as shown in terminal) //
//***************************************//

{
  //*******************************************************************************//
  // If you need to add multiple root outputs, by multithreading, use this command:
  system("hadd -O -f molecular-dna.root molecular-dna_t*.root");

  // Define these parameters of the simulation
  char ifile[256] = "molecular-dna.root";  // input filepath
  const Int_t numberOfPlasmids = 10144;
  Double_t r3 = 4.42e-6 * 4.42e-6 * 4.84e-6;  // r * r * r  (world cube side m)
  Double_t Nbp = 4.367 * 0.001 * numberOfPlasmids;  // Mbp // Length of the DNA chain in Mbp
  // multiplying by 10142 since we are testing 10142 k plasmids
  Double_t mass = 997 * r3;  // density * r3

  //*******************************************************************************//

  typedef std::pair<int64_t, int64_t> ipair;
  bool greaterPair(const ipair& l, const ipair& r);
  bool smallerPair(const ipair& l, const ipair& r);

  void BinLogX(TH1 * h);

  gROOT->Reset();
  gStyle->SetPalette(1);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(00000);

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  // Initialize output histograms
  TCanvas* cdamage = new TCanvas("cdamage", "DNA Damage Distribution", 900, 120, 600, 400);
  cdamage->SetLogy();
  TH1F* h1damage = new TH1F("h1damage", "h1damage", 28, 0, 14);
  TCanvas* c4damage = new TCanvas("c4damage", "Direct damage per plasmid", 900, 120, 600, 400);
  // c4damage->SetLogy();
  TH1F* h4damage = new TH1F("h4damage", "h4damage", 28, 0, 14);
  TCanvas* cSSB = new TCanvas("cSSB", "SSB Distribution", 900, 120, 600, 400);
  cSSB->SetLogy();
  TH1F* h1SSB = new TH1F("h1SSB", "h1SSB", 28, 0, 14);
  TCanvas* cDSB = new TCanvas("cDSB", "DSB Distribution", 900, 120, 600, 400);
  cDSB->SetLogy();
  TH1F* h1DSB = new TH1F("h1DSB", "h1DSB", 40, 0, 10);
  TCanvas* ccount = new TCanvas("cCount", "Damage per plasmid Distribution", 900, 120, 600, 400);
  TH1F* h1count = new TH1F("h1count", "h1count", 11000, 0, 11000);
  TGraph2D* gr1 = new TGraph2D();

  // Open root file
  TFile* f = TFile::Open(ifile);

  // Initialize Variables
  Int_t EB, ES, OHB, OHS, HB, HS, FL;
  Int_t total_EB, total_ES, total_OHB, total_OHS, total_HB, total_HS, total_FL;
  Float_t total_EB2, total_ES2, total_OHB2, total_OHS2, total_HB2, total_HS2, total_FL2;
  Float_t SD_EB, SD_ES, SD_OHB, SD_OHS, SD_HB, SD_HS;
  Float_t SD_SSB, SD_SSBp, SD_SSB2p, SD_sSSB, SD_SSBd, SD_SSBi, SD_SSBm;
  Float_t SD_DSB, SD_DSBp, SD_DSBpp, SD_sDSB, SD_DSBd, SD_DSBi, SD_DSBm, SD_DSBh;

  Int_t SSB, SSBp, SSB2p;
  Int_t total_SSB, total_SSBp, total_SSB2p;
  Float_t total_SSB2, total_SSBp2, total_SSB2p2;
  Int_t DSB, DSBp, DSBpp;
  Int_t total_DSB, total_DSBp, total_DSBpp;
  Float_t total_DSB2, total_DSBp2, total_DSBpp2;

  Int_t SSBd, SSBi, SSBm;
  Int_t total_sSSB, total_SSBd, total_SSBi, total_SSBm;
  Float_t total_sSSB2, total_SSBd2, total_SSBi2, total_SSBm2;
  Int_t DSBd, DSBi, DSBm, DSBh;
  Int_t total_sDSB, total_DSBd, total_DSBi, total_DSBm, total_DSBh;
  Float_t total_sDSB2, total_DSBd2, total_DSBi2, total_DSBm2, total_DSBh2;

  Double_t dose = 0;
  Double_t SD_dose = 0;

  Double_t EB_yield = 0;
  Double_t ES_yield = 0;
  Double_t OHB_yield = 0;
  Double_t OHS_yield = 0;
  Double_t HB_yield = 0;
  Double_t HS_yield = 0;
  Double_t SD_EB_yield = 0;
  Double_t SD_ES_yield = 0;
  Double_t SD_OHB_yield = 0;
  Double_t SD_OHS_yield = 0;
  Double_t SD_HB_yield = 0;
  Double_t SD_HS_yield = 0;

  Double_t SSB_yield = 0;
  Double_t SSBp_yield = 0;
  Double_t SSB2p_yield = 0;
  Double_t SD_SSB_yield = 0;
  Double_t SD_SSBp_yield = 0;
  Double_t SD_SSB2p_yield = 0;
  Double_t DSB_yield = 0;
  Double_t DSBp_yield = 0;
  Double_t DSBpp_yield = 0;
  Double_t SD_DSB_yield = 0;
  Double_t SD_DSBp_yield = 0;
  Double_t SD_DSBpp_yield = 0;

  Double_t sSSB_yield = 0;
  Double_t SSBi_yield = 0;
  Double_t SSBd_yield = 0;
  Double_t SSBm_yield = 0;
  Double_t SD_sSSB_yield = 0;
  Double_t SD_SSBi_yield = 0;
  Double_t SD_SSBd_yield = 0;
  Double_t SD_SSBm_yield = 0;
  Double_t sDSB_yield = 0;
  Double_t DSBi_yield = 0;
  Double_t DSBd_yield = 0;
  Double_t DSBm_yield = 0;
  Double_t DSBh_yield = 0;
  Double_t SD_sDSB_yield = 0;
  Double_t SD_DSBi_yield = 0;
  Double_t SD_DSBd_yield = 0;
  Double_t SD_DSBm_yield = 0;
  Double_t SD_DSBh_yield = 0;

  total_EB = 0;
  total_ES = 0;
  total_OHB = 0;
  total_OHS = 0;
  total_HB = 0;
  total_HS = 0;

  total_SSB = 0;
  total_SSBp = 0;
  total_SSB2p = 0;
  total_SSB2 = 0;
  total_SSBp2 = 0;
  total_SSB2p2 = 0;
  total_DSB = 0;
  total_DSBp = 0;
  total_DSBpp = 0;
  total_DSB2 = 0;
  total_DSBp2 = 0;
  total_DSBpp2 = 0;

  total_sSSB = 0;
  total_SSBd = 0;
  total_SSBi = 0;
  total_SSBm = 0;
  total_sSSB2 = 0;
  total_SSBd2 = 0;
  total_SSBi2 = 0;
  total_SSBm2 = 0;
  total_sDSB = 0;
  total_DSBd = 0;
  total_DSBi = 0;
  total_DSBm = 0;
  total_DSBh = 0;
  total_sDSB2 = 0;
  total_DSBd2 = 0;
  total_DSBi2 = 0;
  total_DSBm2 = 0;
  total_DSBh2 = 0;

  Double_t eVtoJ = 1.60218e-19;
  Double_t EnergyDeposited_eV = 0;
  Double_t acc_edep = 0;
  Double_t acc_edep2 = 0;

  Double_t Energy;
  Double_t BPID;
  Int_t Strand;
  Int_t StrandDamage;
  Int_t Event1, directB;
  Char_t Primary;
  char* primaryName = new char[32];
  char* type = new char[256];
  char* sourceClassification = new char[256];

  Double_t xx, yy, zz;

  // Read trees and leaves from root file, and give values to variables
  TTree* tree = (TTree*)f->Get("tuples/primary_source");
  Float_t number = (Float_t)tree->GetEntries();

  vector<pair<int, int64_t>> DSBBPID;

  // For reading species production
  tree = (TTree*)f->Get("tuples/damage");
  tree->SetBranchAddress("Primary", &Primary);
  tree->SetBranchAddress("Energy", &Energy);
  tree->SetBranchAddress("EaqBaseHits", &EB);
  tree->SetBranchAddress("EaqStrandHits", &ES);
  tree->SetBranchAddress("OHBaseHits", &OHB);
  tree->SetBranchAddress("OHStrandHits", &OHS);
  tree->SetBranchAddress("HBaseHits", &HB);
  tree->SetBranchAddress("HStrandHits", &HS);
  tree->SetBranchAddress("TypeClassification", type);
  tree->SetBranchAddress("BasePair", &BPID);
  tree->SetBranchAddress("Event", &Event1);
  tree->SetBranchAddress("Strand", &Strand);
  tree->SetBranchAddress("StrandDamage", &StrandDamage);
  tree->SetBranchAddress("Position_x_um", &xx);
  tree->SetBranchAddress("Position_y_um", &yy);
  tree->SetBranchAddress("Position_z_um", &zz);
  tree->SetBranchAddress("DirectBreaks", &directB);

  int primea = 0;
  Long64_t nentries = tree->GetEntries();
  int damagecount = 0;
  int hitcount = 0;
  // Double_t plasmidsD [25992];
  Double_t plasmidsD[2599200];
  for (int i = 0; i < nentries; i++) {
    tree->GetEntry(i);
    total_EB += EB;
    total_EB2 += pow(EB, 2);
    total_ES += ES;
    total_ES2 += pow(ES, 2);
    total_OHB += OHB;
    total_OHB2 += pow(OHB, 2);
    total_OHS += OHS;
    total_OHS2 += pow(OHS, 2);
    total_HB += HB;
    total_HB2 += pow(HB, 2);
    total_HS += HS;
    total_HS2 += pow(HS, 2);
    if ((string)type == "DSB" || (string)type == "DSB+" || (string)type == "DSB++") {
      DSBBPID.push_back(make_pair(i, (int64_t)BPID));
    }
    if (StrandDamage != 0) {
      {
        damagecount++;
        h1count->Fill(Strand);
      }
      primea++;
    }
  }
  int unbrokenP = 0;
  int brokenP = 0;
  for (int i = 0; i < 11000; i++) {
    plasmidsD[i] = h1count->GetBinContent(i);
    if (plasmidsD[i] == 0 && i < numberOfPlasmids) unbrokenP++;
  }

  Double_t totalDs = 0.0;
  for (int i = 0; i < 11000; i++) {
    if (plasmidsD[i] != 0) {
      h4damage->Fill(plasmidsD[i]);
      totalDs += plasmidsD[i];
      brokenP++;
    }
  }

  Double_t X;
  for (int i = 0; i < 28; i++) {
    X = h4damage->GetBinContent(i);
    h4damage->SetBinContent(i, X * 100 / totalDs);
  }

  // Sort DSBs from the one with lower ID value to the one with higher ID value
  // Then find the number of fragments that have been produced
  sort(DSBBPID.begin(), DSBBPID.end(), smallerPair);
  for (int ie = 0; ie < DSBBPID.size() - 1; ie++) {
    int64_t dsbfragment = DSBBPID[ie + 1].second - DSBBPID[ie].second;
  }

  // Calculate the standard deviation of species
  SD_EB = sqrt(((total_EB2 / number) - pow(total_EB / number, 2)) / (number - 1));
  SD_ES = sqrt(((total_ES2 / number) - pow(total_ES / number, 2)) / (number - 1));
  SD_OHB = sqrt(((total_OHB2 / number) - pow(total_OHB / number, 2)) / (number - 1));
  SD_OHS = sqrt(((total_OHS2 / number) - pow(total_OHS / number, 2)) / (number - 1));
  SD_HB = sqrt(((total_HB2 / number) - pow(total_HB / number, 2)) / (number - 1));
  SD_HS = sqrt(((total_HS2 / number) - pow(total_HS / number, 2)) / (number - 1));

  // Read damage classification SSB, SSB+, 2SSB, DSB, DSB+, DSB++
  // As they have been defined in: Nikjoo, H., O’Neill, O., Goodhead, T., & Terrissol, M. 1997,
  // Computational modelling of low-energy electron-induced DNA damage by early physical
  // and chemical events, International Journal of Radiation Biology, 71, 467.
  tree = (TTree*)f->Get("tuples/classification");
  tree->SetBranchAddress("Primary", &Primary);
  tree->SetBranchAddress("Energy", &Energy);
  tree->SetBranchAddress("SSB", &SSB);
  tree->SetBranchAddress("SSBp", &SSBp);
  tree->SetBranchAddress("2SSB", &SSB2p);
  tree->SetBranchAddress("DSB", &DSB);
  tree->SetBranchAddress("DSBp", &DSBp);
  tree->SetBranchAddress("DSBpp", &DSBpp);

  Long64_t nentriesC = tree->GetEntries();
  for (int i = 0; i < nentriesC; i++) {
    tree->GetEntry(i);

    total_SSBp += SSBp;
    total_SSBp2 += pow(SSBp, 2);
    total_SSB2p += SSB2p;
    total_SSB2p2 += pow(SSB2p, 2);
    total_SSB += SSB;
    total_SSB2 += pow(SSB, 2);

    total_DSBp += DSBp;
    total_DSBp2 += pow(DSBp, 2);
    total_DSBpp += DSBpp;
    total_DSBpp2 += pow(DSBpp, 2);
    total_DSB += DSB;
    total_DSB2 += pow(DSB, 2);
  }

  // Calculate the standard deviation
  SD_SSB = sqrt(((total_SSB2 / number) - pow(total_SSB / number, 2)) / (number - 1));
  SD_SSBp = sqrt(((total_SSBp2 / number) - pow(total_SSBp / number, 2)) / (number - 1));
  SD_SSB2p = sqrt(((total_SSB2p2 / number) - pow(total_SSB2p / number, 2)) / (number - 1));

  SD_DSB = sqrt(((total_DSB2 / number) - pow(total_DSB / number, 2)) / (number - 1));
  SD_DSBp = sqrt(((total_DSBp2 / number) - pow(total_DSBp / number, 2)) / (number - 1));
  SD_DSBpp = sqrt(((total_DSBpp2 / number) - pow(total_DSBpp / number, 2)) / (number - 1));

  // Read damage classification SSBd, SSBi, SSBm, DSBd, DSBi, DSBm, DSBh
  // As they have been defined in: Nikjoo, H., O’Neill, O., Goodhead, T., & Terrissol, M. 1997,
  // Computational modelling of low-energy electron-induced DNA damage by early physical
  // and chemical events, International Journal of Radiation Biology, 71, 467.
  tree = (TTree*)f->Get("tuples/source");
  tree->SetBranchAddress("Primary", primaryName);
  tree->SetBranchAddress("Energy", &Energy);
  tree->SetBranchAddress("SSBd", &SSBd);
  tree->SetBranchAddress("SSBi", &SSBi);
  tree->SetBranchAddress("SSBm", &SSBm);
  tree->SetBranchAddress("DSBd", &DSBd);
  tree->SetBranchAddress("DSBi", &DSBi);
  tree->SetBranchAddress("DSBm", &DSBm);
  tree->SetBranchAddress("DSBh", &DSBh);

  int iprime = 0;

  Long64_t nentriesS = tree->GetEntries();
  for (int i = 0; i < nentriesS; i++) {
    tree->GetEntry(i);

    total_SSBd += SSBd;
    total_SSBd2 += pow((SSBd), 2);
    total_SSBi += SSBi;
    total_SSBi2 += pow((SSBi), 2);
    total_SSBm += SSBm;
    total_SSBm2 += pow((SSBm), 2);
    total_sSSB += SSBd + SSBi + SSBm;
    total_sSSB2 += pow((SSBd + SSBi + SSBm), 2);

    total_DSBd += DSBd;
    total_DSBd2 += pow(DSBd, 2);
    total_DSBi += DSBi;
    total_DSBi2 += pow(DSBi, 2);
    total_DSBm += DSBm;
    total_DSBm2 += pow(DSBm, 2);
    total_DSBh += DSBh;
    total_DSBh2 += pow(DSBh, 2);
    total_sDSB += DSBd + DSBi + DSBm + DSBh;
    total_sDSB2 += pow((DSBd + DSBi + DSBm + DSBh), 2);

    if (SSBd != 0) h1SSB->Fill(SSBd);
    if (DSBd != 0) h1DSB->Fill(DSBd);
    if (SSBd != 0) h1damage->Fill(SSBd);
    if (DSBd != 0) h1damage->Fill(DSBd);
    if (SSB != 0 || DSBd != 0) {
      iprime++;
    }
  }

  Double_t Y;
  for (int i = 0; i < 28; i++) {
    Y = h1damage->GetBinContent(i);
    h1damage->SetBinContent(i, Y * 100 / totalDs);
  }

  // Calculate the standard deviation
  SD_sSSB = sqrt(((total_sSSB2 / number) - pow(total_sSSB / number, 2)) / (number - 1));
  SD_SSBd = sqrt(((total_SSBd2 / number) - pow(total_SSBd / number, 2)) / (number - 1));
  SD_SSBi = sqrt(((total_SSBi2 / number) - pow(total_SSBi / number, 2)) / (number - 1));
  SD_SSBm = sqrt(((total_SSBm2 / number) - pow(total_SSBm / number, 2)) / (number - 1));

  SD_sDSB = sqrt(((total_sDSB2 / number) - pow(total_sDSB / number, 2)) / (number - 1));
  SD_DSBd = sqrt(((total_DSBd2 / number) - pow(total_DSBd / number, 2)) / (number - 1));
  SD_DSBi = sqrt(((total_DSBi2 / number) - pow(total_DSBi / number, 2)) / (number - 1));
  SD_DSBm = sqrt(((total_DSBm2 / number) - pow(total_DSBm / number, 2)) / (number - 1));
  SD_DSBh = sqrt(((total_DSBh2 / number) - pow(total_DSBh / number, 2)) / (number - 1));

  // Measure the Deposited Energy in the whole volume that includes DNA chain

  tree = (TTree*)f->Get("tuples/chromosome_hits");
  tree->SetBranchAddress("e_chromosome_kev", &EnergyDeposited_eV);
  nentries = tree->GetEntries();
  for (int i = 0; i < nentries; i++) {
    tree->GetEntry(i);
    acc_edep += EnergyDeposited_eV * 1e3;
    acc_edep2 += EnergyDeposited_eV * EnergyDeposited_eV * 1e6;
  }
  tree->SetBranchAddress("e_dna_kev", &EnergyDeposited_eV);
  nentries = tree->GetEntries();
  for (int i = 0; i < nentries; i++) {
    tree->GetEntry(i);
    acc_edep += (EnergyDeposited_eV * 1e3);
    acc_edep2 += (EnergyDeposited_eV * EnergyDeposited_eV * 1e6);
  }

  // Close the root file to free space
  f->Close();
  // Calculate the absorbed dose
  dose = acc_edep * eVtoJ / mass;

  double norm = 1;
  // Calculate the yields, together with their standard deviation
  EB_yield = (Double_t)total_EB / dose / Nbp;
  ES_yield = (Double_t)total_ES / dose / Nbp;
  OHB_yield = (Double_t)total_OHB / dose / Nbp;
  OHS_yield = (Double_t)total_OHS / dose / Nbp;
  HB_yield = (Double_t)total_HB / dose / Nbp;
  HS_yield = (Double_t)total_HS / dose / Nbp;

  SD_EB_yield = SD_EB / dose / Nbp;
  SD_ES_yield = SD_ES / dose / Nbp;
  SD_OHB_yield = SD_OHB / dose / Nbp;
  SD_OHS_yield = SD_OHS / dose / Nbp;
  SD_HB_yield = SD_HB / dose / Nbp;
  SD_HS_yield = SD_HS / dose / Nbp;

  SSB_yield = (Double_t)norm * total_SSB / dose / Nbp;
  SSBp_yield = (Double_t)norm * total_SSBp / dose / Nbp;
  SSB2p_yield = (Double_t)norm * total_SSB2p / dose / Nbp;

  DSB_yield = (Double_t)norm * total_DSB / dose / Nbp;
  DSBp_yield = (Double_t)norm * total_DSBp / dose / Nbp;
  DSBpp_yield = (Double_t)norm * total_DSBpp / dose / Nbp;

  SD_SSB_yield = norm * SD_SSB / dose / Nbp;
  SD_SSBp_yield = norm * SD_SSBp / dose / Nbp;
  SD_SSB2p_yield = norm * SD_SSB2p / dose / Nbp;

  SD_DSB_yield = norm * SD_DSB / dose / Nbp;
  SD_DSBp_yield = norm * SD_DSBp / dose / Nbp;
  SD_DSBpp_yield = norm * SD_DSBpp / dose / Nbp;

  sSSB_yield = (Double_t)norm * total_sSSB / dose / Nbp;
  SSBi_yield = (Double_t)norm * total_SSBi / dose / Nbp;
  SSBd_yield = (Double_t)norm * total_SSBd / dose / Nbp;
  SSBm_yield = (Double_t)norm * total_SSBm / dose / Nbp;

  sDSB_yield = (Double_t)norm * total_sDSB / dose / Nbp;
  DSBi_yield = (Double_t)norm * total_DSBi / dose / Nbp;
  DSBd_yield = (Double_t)norm * total_DSBd / dose / Nbp;
  DSBm_yield = (Double_t)norm * total_DSBm / dose / Nbp;
  DSBh_yield = (Double_t)norm * total_DSBh / dose / Nbp;

  SD_sSSB_yield = norm * SD_sSSB / dose / Nbp;
  SD_SSBi_yield = norm * SD_SSBi / dose / Nbp;
  SD_SSBd_yield = norm * SD_SSBd / dose / Nbp;
  SD_SSBm_yield = norm * SD_SSBm / dose / Nbp;

  SD_sDSB_yield = norm * SD_sDSB / dose / Nbp;
  SD_DSBi_yield = norm * SD_DSBi / dose / Nbp;
  SD_DSBd_yield = norm * SD_DSBd / dose / Nbp;
  SD_DSBm_yield = norm * SD_DSBm / dose / Nbp;
  SD_DSBh_yield = norm * SD_DSBh / dose / Nbp;

  // Print output in terminal

  float total_SSB_totalYield = SSB_yield + SSBp_yield + SSB2p_yield;
  float total_DSB_totalYield = DSB_yield + DSBp_yield + DSBpp_yield;

  cout << "\n"
       << " Output file: " << ifile << '\n'
       << "\nDose Absorbed (Gy): " << dose << '\n'
       << "Particle : " << primaryName << '\t' << "Energy (MeV) : " << Energy << '\t'
       << "Number of Primaries : " << number << '\n'

       << "  Output Damage : " << '\n'

       << '\t' << "Number of plasmids:         " << numberOfPlasmids << "   \t" << '\n'
       << '\t' << "Unbroken plasmids:          " << unbrokenP << "   \t" << '\n'
       << '\t' << "Broken plasmids:            " << brokenP << "   \t" << '\n'
       << '\t' << "Rate of damages per plasmid " << totalDs / numberOfPlasmids << "   \t" << '\n'

       << '\t' << "Rate of damages per plasmid per dose " << '\t'
       << totalDs / numberOfPlasmids / dose << " plasmid-1Gy-1" << '\n'
       << '\t' << "Rate of damages per dose per Mbp " << totalDs / dose / Nbp << " Gy-1Mbp-1"
       << '\n'

       << '\t' << "     Species Hits " << '\n'
       << '\t' << "EaqBaseHits    " << EB_yield * dose * Nbp << "   \t" << '\n'
       << '\t' << "EaqStrandHits  " << ES_yield * dose * Nbp << "   \t" << '\n'
       << '\t' << "OHBaseHits     " << OHB_yield * dose * Nbp << "   \t" << '\n'
       << '\t' << "OHStrandHits   " << OHS_yield * dose * Nbp << "   \t" << '\n'
       << '\t' << "HBaseHits      " << HB_yield * dose * Nbp << "   \t" << '\n'
       << '\t' << "HStrandHits    " << HS_yield * dose * Nbp << "   \t" << '\n'
       << '\n'
       << '\t' << "     Damage number" << '\n'
       << '\t' << "SSB            " << SSB_yield * dose * Nbp << "   \t" << '\n'
       << '\t' << "SSB+           " << SSBp_yield * dose * Nbp << "   \t" << '\n'
       << '\t' << "2SSB           " << SSB2p_yield * dose * Nbp << "   \t" << '\n'
       << '\t' << "SSB total      " << total_SSB_totalYield * dose * Nbp << '\n'
       << '\t' << "DSB            " << DSB_yield * dose * Nbp << "   \t" << '\n'
       << '\t' << "DSB+           " << DSBp_yield * dose * Nbp << "   \t" << '\n'
       << '\t' << "DSB++          " << DSBpp_yield * dose * Nbp << "   \t" << '\n'
       << '\t' << "DSB total      " << total_DSB_totalYield * dose * Nbp << '\n'
       << '\n'
       << '\t' << "     Breaks number " << '\n'
       << '\t' << "SSB direct     " << SSBd_yield * dose * Nbp << "   \t" << '\n'
       << '\t' << "SSB indirect   " << SSBi_yield * dose * Nbp << "   \t" << '\n'
       << '\t' << "SSB mixed      " << SSBm_yield * dose * Nbp << "   \t" << '\n'
       << '\t' << "SSB total      " << sSSB_yield * dose * Nbp << "   \t" << '\n'
       << '\t' << "DSB direct     " << DSBd_yield * dose * Nbp << "   \t" << '\n'
       << '\t' << "DSB indirect   " << DSBi_yield * dose * Nbp << "   \t" << '\n'
       << '\t' << "DSB mixed      " << DSBm_yield * dose * Nbp << "   \t" << '\n'
       << '\t' << "DSB hybrid     " << DSBh_yield * dose * Nbp << "   \t" << '\n'
       << '\t' << "DSB total      " << sDSB_yield * dose * Nbp << "   \t" << '\n'
       << '\n'
       << '\t' << "SSB/DSB        " << sSSB_yield / sDSB_yield << '\n'
       << '\n';

  // Plot Histograms

  cSSB->GetCanvas()->cd();
  h1SSB->SetStats(false);
  h1SSB->SetMarkerSize(0.1);
  h1SSB->SetMarkerColor(kBlue);
  h1SSB->SetLineColor(kBlue);
  h1SSB->SetTitle("");
  h1SSB->SetYTitle(" ");
  h1SSB->SetXTitle("Number of direct SSB per event");
  h1SSB->SetFillColor(kBlue);
  h1SSB->Draw();

  cDSB->GetCanvas()->cd();
  h1DSB->SetStats(false);
  h1DSB->SetMarkerSize(0.1);
  h1DSB->SetMarkerColor(kGreen + 2);
  h1DSB->SetLineColor(kGreen + 2);
  h1DSB->SetTitle("");
  h1DSB->SetYTitle(" ");
  h1DSB->SetXTitle("Number of direct DSB per event");
  h1DSB->SetFillColor(kGreen + 2);
  h1DSB->Draw();

  cdamage->GetCanvas()->cd();
  h1damage->SetStats(false);
  h1damage->SetMarkerSize(0.1);
  h1damage->SetMarkerColor(kMagenta);
  h1damage->SetLineColor(kMagenta);
  h1damage->SetTitle("");
  h1damage->SetYTitle("Percentage % ");
  h1damage->SetXTitle("Number of damages per event");
  h1damage->SetFillColor(kMagenta);
  h1damage->Draw();

  c4damage->GetCanvas()->cd();
  h4damage->SetStats(false);
  h4damage->SetMarkerSize(0.1);
  h4damage->SetMarkerColor(kRed - 4);
  h4damage->SetLineColor(kRed - 4);
  h4damage->SetTitle("");
  h4damage->SetYTitle("Percentage % ");
  h4damage->SetXTitle("Number of damages per plasmid");
  h4damage->SetFillColor(kRed - 4);
  h4damage->Draw();

  ccount->GetCanvas()->cd();
  h1count->SetStats(false);
  h1count->SetMarkerSize(0.1);
  h1count->SetMarkerColor(kCyan);
  h1count->SetLineColor(kCyan);
  h1count->SetTitle("");
  h1count->SetYTitle("Number of damages");
  h1count->SetXTitle("Plasmid ID");
  h1count->Draw();
}

// Some important bools that are needed to run the root macro file
bool greaterPair(const ipair& l, const ipair& r)
{
  return l.second > r.second;
}
bool smallerPair(const ipair& l, const ipair& r)
{
  return l.second < r.second;
}

void BinLogX(TH1* h)
{
  TAxis* axis = h->GetXaxis();
  int bins = axis->GetNbins();
  Axis_t from = axis->GetXmin();
  Axis_t to = axis->GetXmax();
  Axis_t width = (to - from) / bins;
  Axis_t* new_bins = new Axis_t[bins + 1];
  for (int i = 0; i <= bins; i++) {
    new_bins[i] = TMath::Power(10, from + i * width);
  }
  axis->Set(bins, new_bins);
  delete[] new_bins;
}
