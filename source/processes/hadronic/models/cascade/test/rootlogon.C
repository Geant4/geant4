{

//
// This macro generates a Controlbar menu: To see the output, click begin_html <a href="gif/demos.gif" >here</a> end_html
// To execute an item, click with the left mouse button.
// To see the HELP of a button, click on the right mouse button.

   printf("\n Geant4 cascade benchmarking\n\n");
   printf("\nType \".x [scriptname].C\" to run the spesific script.\n");

   gROOT->Reset();
   gStyle->SetScreenFactor(1); //if you have a large screen, select 1,2 or 1.4

   bar = new TControlBar("vertical", "INC",10,10);

   bar->AddButton("ROOT Browser",     "new TBrowser;",  "Start the ROOT Browser");
//   bar->AddButton("Benchmarks",".x benchmarks.C", "Cascade benchmarking.");
   bar->AddButton("Cascade interface",   ".x cascade.C", "General test for cascade codes.");

   bar->Show();
   gROOT->SaveContext();
}

