struct SpeciesInfo
{
   SpeciesInfo()
   {
   }
   SpeciesInfo(const SpeciesInfo& right) :
      fG(right.fG),
      fGerr(right.fGerr),
      fLET(right.fLET),
      fLETerr(right.fLETerr),
      fName(right.fName)
   {}
   SpeciesInfo& operator=(const SpeciesInfo& right)
   {
      if(this == &right) return *this;
      fG = right.fG;
      fGerr = right.fGerr;
      fLET = right.fLET;
      fLETerr = right.fLETerr;
      fName = right.fName;
      return *this;
   }

   std::vector<Double_t> fG;
   std::vector<Double_t> fGerr;
   std::vector<Double_t> fLET;
   std::vector<Double_t> fLETerr;
   string fName;
};

const char* filetypes[] = {
   "PostScript", "*.ps",
   "Encapsulated PostScript", "*.eps",
   "PDF files", "*.pdf",
   "Gif files", "*.gif",
   "PNG files", "*.png",
   "All files", "*",
   0, 0
};

TGTab *gTab = nullptr;

void Save()
{
  TGFileInfo fi;
  fi.fFileTypes = filetypes;

  new TGFileDialog(gClient->GetRoot(),gClient->GetRoot(),kFDSave,&fi);
  gROOT->GetListOfCanvases()->At(gTab->GetCurrent())->SaveAs(fi.fFilename);
}

void plotG_LET()
{

   std::map<Int_t, SpeciesInfo> speciesInfo;

   gROOT->SetStyle("Plain");
   gStyle->SetPalette(1);
   gStyle->SetCanvasBorderMode(0);
   gStyle->SetFrameBorderMode(0);
   gStyle->SetPadTickX(1);
   gStyle->SetPadTickY(1);

   TGMainFrame *main = new TGMainFrame(gClient->GetRoot(), 200, 200);
   gTab = new TGTab(main, 200, 200);

   Int_t ncols, tag, runID;
   Double_t LET, LET_sigma, Gvalue, Gvalue_sigma;
   string name;
   string dummy;
   string string_key, string_value;

   ifstream file;
   file.open("Species.txt",std::ios::in);

   runID = 0;

   while(1) {
      // Read LET values
      file >> dummy >> LET >> dummy >> LET_sigma;
      if (file.eof()) break;

      std::getline(file,dummy);

      if (!std::getline(file,string_key)) break;
      std::istringstream key(string_key); // Read keys (name of molecule)
      if (!std::getline(file,string_value)) break;
      std::istringstream value(string_value); // Read G value
      while(1) {
         key >> name >> tag;
         value >> Gvalue >> Gvalue_sigma;
         if (!key || !value) break;
         speciesInfo[tag].fName = name;
         speciesInfo[tag].fG.resize(runID+1);
         speciesInfo[tag].fGerr.resize(runID+1);
         speciesInfo[tag].fLET.resize(runID+1);
         speciesInfo[tag].fLETerr.resize(runID+1);

         speciesInfo[tag].fG[runID] = Gvalue;
         speciesInfo[tag].fGerr[runID] = Gvalue_sigma;
         speciesInfo[tag].fLET[runID] = LET;
         speciesInfo[tag].fLETerr[runID] = LET_sigma;
      }
      runID++;
   }
   file.close();

   for (auto it_map : speciesInfo) {

      auto map = it_map.second;
      TGraphErrors* gSpecies = new TGraphErrors(map.fG.size(),
                                                map.fLET.data(),
                                                map.fG.data(),
                                                map.fLETerr.data(),
                                                map.fGerr.data());

      Int_t color = (2+it_map.first)%TColor::GetNumberOfColors();
      if (color == 5 || color == 10 || color == 0) ++color;


      TGCompositeFrame *tf = gTab->AddTab(map.fName.c_str());
      TGCompositeFrame *frame = new TGCompositeFrame(tf, 60, 60,
                                                     kHorizontalFrame);

      tf->AddFrame(frame, new TGLayoutHints(kLHintsExpandX|kLHintsExpandY,
                   10,10,10,2));

      TRootEmbeddedCanvas *c1 = new TRootEmbeddedCanvas(map.fName.c_str(),
                                                        frame, 700, 500);
      frame->AddFrame(c1, new TGLayoutHints(kLHintsExpandX|kLHintsExpandY,
                      10,10,10,2));
      c1->GetCanvas()->SetLogx();

      TGHorizontalFrame* hframe = new TGHorizontalFrame(tf, 200, 40);

      TGTextButton* save = new TGTextButton(hframe, "&Save as ...",
                                            "Save()");
      hframe->AddFrame(save, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

      TGTextButton *exit = new TGTextButton(hframe, "&Exit ",
                                            "gApplication->Terminate()");
      hframe->AddFrame(exit, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

      tf->AddFrame(hframe, new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2));

      gSpecies->SetTitle(map.fName.c_str());
      gSpecies->SetMarkerStyle(20+it_map.first);
      gSpecies->SetMarkerColor(color);
      gSpecies->GetXaxis()->SetTitle("LET (keV/um)");
      gSpecies->GetXaxis()->SetTitleOffset(1.1);
      gSpecies->GetYaxis()->SetTitle("G value (molecules/100 eV)");
      gSpecies->GetYaxis()->SetTitleOffset(1.2);
      gSpecies->Draw("AP");
   }

   main->AddFrame(gTab, new TGLayoutHints(kLHintsBottom | kLHintsExpandX |
                                          kLHintsExpandY, 2, 2, 5, 1));

   main->MapSubwindows();
   main->Resize();   // resize to default size
   main->MapWindow();
}

