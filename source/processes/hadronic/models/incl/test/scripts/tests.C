{
   gROOT->Reset();

   //Add the tutorials directory to the macro path
   //This is necessary in case this macro is executed from another user directory
   TString dir = gSystem->UnixPathName(TCint::GetCurrentMacroName());
   dir.ReplaceAll("demos.C","");
   dir.ReplaceAll("/./","");
   const char *current = gROOT->GetMacroPath();
   gROOT->SetMacroPath(Form("%s:%s",current,dir.Data()));
   
   TControlBar *bar = new TControlBar("vertical", "Demos",10,10);
   bar->AddButton("test1: dd.C - Double-differential plots",".x scripts/dd.C", "Double-Differential (p(1.2GeV) + Pb)");
   bar->AddButton("test2: fragments.C - Nucleus fragment plots",".x scripts/fragments.C", "Nucleus fragments (p(1.0GeV) + Pb)");
   bar->AddButton("browser", "new TBrowser;", "Start the ROOT Browser");
   bar->AddButton("Exit", ".q", "Exit");
   bar->SetButtonWidth(90);
   bar->Show();
   gROOT->SaveContext();
}
