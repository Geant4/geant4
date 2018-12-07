void ReadASCII(TString source_file_ascii, TString output_file_root)
{

// This macro converts the results of the simulation stored in ASCII files to ROOT files. 
  ifstream in;
  in.open(source_file_ascii);

  Float_t x,y,z, edep;
  Int_t nlines = 0;
  TFile *f = new TFile(output_file_root,"RECREATE");
  TH2F *h20 = new TH2F("h20","h20",801,-100.125,100.125, 801, -100.125, 100.125);

   while (1) {
      in >> x >> y >> z >> edep;
      if (!in.good()) break;
    //  if (edep !=0.) printf("x=%8f, y=%8f, edep=%8f\n",x,y,edep);

      h20->Fill(x,y,edep);
      nlines++;
   }
   printf(" found %d points\n",nlines);

   in.close();

   f->Write();
}
