struct SpeciesInfoAOS
{
   SpeciesInfoAOS()
   {
      fNEvent = 0;
      fNumber = 0;
      fG = 0.;
      fG2 = 0.;
   }

   SpeciesInfoAOS(const SpeciesInfoAOS& right) // Species A(B);
   {
      fNEvent = right.fNEvent;
      fNumber = right.fNumber;
      fG = right.fG;
      fG2 = right.fG2;
      fName = right.fName;
   }

   SpeciesInfoAOS& operator=(const SpeciesInfoAOS& right) // A = B
   {
      if(&right == this) return *this;
      fNEvent = right.fNEvent;
      fNumber = right.fNumber;
      fG = right.fG;
      fG2 = right.fG2;
      fName = right.fName;
      return *this;
   }

   Int_t fNEvent;
   Int_t fNumber;
   Double_t fG;
   Double_t fG2;
   string fName;
};

//------------------------------------------------------------------------

struct SpeciesInfoSOA
{
   SpeciesInfoSOA()
   {
      fRelatErr = 0;
   }

   SpeciesInfoSOA(const SpeciesInfoSOA& right) :
      fG(right.fG),
      fGerr(right.fGerr),
      fTime(right.fTime),
      fRelatErr(right.fRelatErr),
      fName(right.fName)
   {}

   SpeciesInfoSOA& operator=(const SpeciesInfoSOA& right)
   {
      if(this == &right) return *this;
      fG = right.fG;
      fGerr = right.fGerr;
      fTime = right.fTime;
      fRelatErr = right.fRelatErr;
      fName = right.fName;
      return *this;
   }

   std::vector<Double_t> fG;
   std::vector<Double_t> fGerr;
   std::vector<Double_t> fTime;
   Double_t fRelatErr;
   Double_t fMin, fMax;
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


void plotG()
{
   gROOT->SetStyle("Plain");
   gStyle->SetPalette(1);
   gStyle->SetCanvasBorderMode(0);
   gStyle->SetFrameBorderMode(0);
   gStyle->SetPadTickX(1);
   gStyle->SetPadTickY(1);

   TGMainFrame *main = new TGMainFrame(gClient->GetRoot(), 200, 200);
   gTab = new TGTab(main, 200, 200);
   Double_t timeA, sumG, sumG2;
   Int_t speciesID, number, nEvent;
   char speciesName[500];

   TFile* file = new TFile;
   file = TFile::Open("scorer.root");

   TTree* tree = (TTree*)file->Get("species");
   tree->SetBranchAddress("speciesID", &speciesID);
   tree->SetBranchAddress("number", &number);
   tree->SetBranchAddress("nEvent", &nEvent);
   tree->SetBranchAddress("speciesName", &speciesName);
   tree->SetBranchAddress("time", &timeA);
   tree->SetBranchAddress("sumG", &sumG);
   tree->SetBranchAddress("sumG2", &sumG2);

   Long64_t nentries = tree->GetEntries();

   if(nentries == 0) {
      cout << "No entries found in the tree species contained in the file "
      << file->GetPath() << endl;
      exit(1);
   }

   std::map<int, std::map<double, SpeciesInfoAOS>> speciesTimeInfo;

   for (Int_t j=0; j < nentries; j++) {
     tree->GetEntry(j);

     SpeciesInfoAOS& infoAOS = speciesTimeInfo[speciesID][timeA];

     infoAOS.fNumber += number;
     infoAOS.fG += sumG;
     infoAOS.fG2 += sumG2;
     infoAOS.fNEvent += nEvent;
     infoAOS.fName = speciesName;
   }

   std::map<Int_t, SpeciesInfoSOA> speciesInfo;

   auto it_SOA = speciesTimeInfo.begin();
   auto end_SOA = speciesTimeInfo.end();

   for (; it_SOA!=end_SOA;++it_SOA) {
      const Int_t _speciesID = it_SOA->first;
      SpeciesInfoSOA& info = speciesInfo[_speciesID];

      auto it2 = it_SOA->second.begin();
      auto end2 = it_SOA->second.end();

      info.fName = it2->second.fName;
      const size_t size2 = it_SOA->second.size();
      info.fG.resize(size2);
      info.fGerr.resize(size2);
      info.fTime.resize(size2);

      Int_t color = (2+_speciesID)%TColor::GetNumberOfColors();
      if(color == 5 || color == 10 || color == 0) ++color;

      for (int i2 = 0 ;it2!=end2;++it2, ++i2) {
         SpeciesInfoAOS& infoAOS = it2->second;

         Double_t _SumG2 = infoAOS.fG2;
         Double_t _MeanG = infoAOS.fG/infoAOS.fNEvent;
         Double_t _Gerr = (infoAOS.fNEvent > 1) ? sqrt((_SumG2/infoAOS.fNEvent - pow(_MeanG,2))
                                                    /(infoAOS.fNEvent-1) ) : 0.;

         info.fG[i2] = _MeanG;
         info.fGerr[i2] = _Gerr;
         info.fTime[i2] = it2->first;

        info.fRelatErr += _Gerr/(_MeanG + 1e-30);
        if(info.fG[i2] != 0)
        {
          std::cout<<info.fTime[i2]<<" "<<info.fG[i2]<<" "<<info.fGerr[i2]<<" "<<info.fName<<std::endl;
        }
  }

      TGraphErrors* gSpecies = new TGraphErrors(info.fG.size(),
                                                info.fTime.data(),
                                                info.fG.data(),
                                                0,
                                                info.fGerr.data());

      TGCompositeFrame *tf = gTab->AddTab(info.fName.c_str());
      TGCompositeFrame *frame = new TGCompositeFrame(tf, 60, 60,
                                                     kHorizontalFrame);

      tf->AddFrame(frame, new TGLayoutHints(kLHintsExpandX|kLHintsExpandY,
                                            10,10,10,2));

      TRootEmbeddedCanvas *c1 = new TRootEmbeddedCanvas(info.fName.c_str(),
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

      gSpecies->SetTitle(info.fName.c_str());
      gSpecies->SetMarkerStyle(20+_speciesID);
      gSpecies->SetMarkerColor(color);
      gSpecies->SetLineColor(color);
      gSpecies->GetXaxis()->SetTitle("Time (ns)");
      gSpecies->GetXaxis()->SetTitleOffset(1.1);
      gSpecies->GetYaxis()->SetTitle("G value (molecules/100 eV)");
      gSpecies->GetYaxis()->SetTitleOffset(1.2);
      gSpecies->Draw("ALP");
   }

   main->AddFrame(gTab, new TGLayoutHints(kLHintsBottom | kLHintsExpandX |
                                          kLHintsExpandY, 2, 2, 5, 1));

   main->MapSubwindows();
   main->Resize();   // resize to default size
   main->MapWindow();

}
