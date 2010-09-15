{
  gROOT->ProcessLine(".x ./scripts/rootlogon.C");
  gROOT->SetStyle("clearRetro");

  gROOT->LoadMacro("./analysis/libAnalysis.C++");

  HistoFactory *hf = new HistoFactory();
  CalculationAnalysis *fortranAnalysis = new CalculationAnalysis(hf, "./tmp/run23ref.root",
								 "h101", 3780.6172,
								 100000);

  CalculationAnalysis *cppAnalysis = new CalculationAnalysis(hf, "./tmp/run23.root", //"./tmp/run49_XFOISA8.root",
							     "h101", 3780.62, //5664.46,
							     100000);

  
  ViewManager *vm = new ViewManager(3,2);
  Double_t Emax = 1000.0;
  Double_t Emin = 1.0;
  Bool_t logE = true;

  vm->plot(cppAnalysis->fillNeutronDoubleDiffXS(1.0, 1.0, logE, Emin, Emax));
  vm->plot(fortranAnalysis->fillNeutronDoubleDiffXS(1.0, 1.0, logE, Emin, Emax), true);

  vm->plot(cppAnalysis->fillNeutronDoubleDiffXS(15.0, 2.0, logE, Emin, Emax));
  vm->plot(fortranAnalysis->fillNeutronDoubleDiffXS(15.0, 2.0, logE, Emin, Emax), true);

  vm->plot(cppAnalysis->fillNeutronDoubleDiffXS(20.0, 2.0, logE, Emin, Emax));
  vm->plot(fortranAnalysis->fillNeutronDoubleDiffXS(20.0, 2.0, logE, Emin, Emax), true);

  vm->plot(cppAnalysis->fillNeutronDoubleDiffXS(30.0, 2.0, logE, Emin, Emax));
  vm->plot(fortranAnalysis->fillNeutronDoubleDiffXS(30.0, 2.0, logE, Emin, Emax), true);

  vm->plot(cppAnalysis->fillNeutronDoubleDiffXS(50.0, 2.0, logE, Emin, Emax));
  vm->plot(fortranAnalysis->fillNeutronDoubleDiffXS(50.0, 2.0, logE, Emin, Emax), true);

  vm->plot(cppAnalysis->fillNeutronDoubleDiffXS(80.0, 2.0, logE, Emin, Emax));
  vm->plot(fortranAnalysis->fillNeutronDoubleDiffXS(80.0, 2.0, logE, Emin, Emax), true);

  ViewManager *vm2 = new ViewManager(3,2);

  vm2->plot(cppAnalysis->fillNeutronDoubleDiffXS(90.0, 1.0, logE, Emin, Emax));
  vm2->plot(fortranAnalysis->fillNeutronDoubleDiffXS(90.0, 1.0, logE, Emin, Emax), true);

  vm2->plot(cppAnalysis->fillNeutronDoubleDiffXS(100.0, 2.0, logE, Emin, Emax));
  vm2->plot(fortranAnalysis->fillNeutronDoubleDiffXS(100.0, 2.0, logE, Emin, Emax), true);

  vm2->plot(cppAnalysis->fillNeutronDoubleDiffXS(110.0, 2.0, logE, Emin, Emax));
  vm2->plot(fortranAnalysis->fillNeutronDoubleDiffXS(110.0, 2.0, logE, Emin, Emax), true);

  vm2->plot(cppAnalysis->fillNeutronDoubleDiffXS(140.0, 2.0, logE, Emin, Emax));
  vm2->plot(fortranAnalysis->fillNeutronDoubleDiffXS(140.0, 2.0, logE, Emin, Emax), true);

  vm2->plot(cppAnalysis->fillNeutronDoubleDiffXS(150.0, 2.0, logE, Emin, Emax));
  vm2->plot(fortranAnalysis->fillNeutronDoubleDiffXS(150.0, 2.0, logE, Emin, Emax), true);

  vm2->plot(cppAnalysis->fillNeutronDoubleDiffXS(160.0, 2.0, logE, Emin, Emax));
  vm2->plot(fortranAnalysis->fillNeutronDoubleDiffXS(160.0, 2.0, logE, Emin, Emax), true);

  ViewManager *vm3 = new ViewManager(3,2);

  vm3->plot(cppAnalysis->fillHisto("Masp", 20, 0.5, 20.5));
  vm3->plot(fortranAnalysis->fillHisto("Masp", 20, 0.5, 20.5), true);

  vm3->plot(cppAnalysis->fillHisto("Mzsp", 20, 0.5, 20.5));
  vm3->plot(fortranAnalysis->fillHisto("Mzsp", 20, 0.5, 20.5), true);

  vm3->plot(cppAnalysis->fillHisto("Exsp", 100, 0.0, 1000));
  vm3->plot(fortranAnalysis->fillHisto("Exsp", 100, 0.0, 1000), true);

  vm3->plot(cppAnalysis->fillHisto("Massini", 210, 0.5, 210.5));
  vm3->plot(fortranAnalysis->fillHisto("Massini", 210, 0.5, 210.5), true);

  vm3->plot(cppAnalysis->fillHisto("Mzini", 210, 0.5, 210.5));
  vm3->plot(fortranAnalysis->fillHisto("Mzini", 210, 0.5, 210.5), true);

  vm3->plot(cppAnalysis->fillHisto("Exini", 100, 0.0, 5000));
  vm3->plot(fortranAnalysis->fillHisto("Exini", 100, 0.0, 5000), true);

  ViewManager *vm4 = new ViewManager(3,2);

  vm4->plot(cppAnalysis->fillHisto("Pxsp", 100, -500.0, 500.0));
  vm4->plot(fortranAnalysis->fillHisto("Pxsp", 100, -500.0, 500.0), true);

  vm4->plot(cppAnalysis->fillHisto("Pysp", 100, -500.0, 500.0));
  vm4->plot(fortranAnalysis->fillHisto("Pysp", 100, -500.0, 500.0), true);

  vm4->plot(cppAnalysis->fillHisto("Pzsp", 100, -500.0, 5000.0));
  vm4->plot(fortranAnalysis->fillHisto("Pzsp", 100, -500.0, 5000.0), true);

  vm4->plot(cppAnalysis->fillHisto("Pxrem", 100, -500.0, 500.0));
  vm4->plot(fortranAnalysis->fillHisto("Pxrem", 100, -500.0, 500.0), true);

  vm4->plot(cppAnalysis->fillHisto("Pyrem", 100, -500.0, 500.0));
  vm4->plot(fortranAnalysis->fillHisto("Pyrem", 100, -500.0, 500.0), true);

  vm4->plot(cppAnalysis->fillHisto("Pzrem", 100, -500.0, 5000.0));
  vm4->plot(fortranAnalysis->fillHisto("Pzrem", 100, -500.0, 5000.0), true);

  ViewManager *vm5 = new ViewManager(3,2);
  vm5->plot(cppAnalysis->fillHisto("Bimpact", 100, 0.0, 20.0));
  vm5->plot(fortranAnalysis->fillHisto("Bimpact", 100, 0.0, 20.0), true);

  vm5->plot(cppAnalysis->fillHisto("Enerj", 100, 1.0, 4000.0, true, "Avv==1 && Zvv == 0 && Ityp==1"));
  vm5->plot(fortranAnalysis->fillHisto("Enerj", 100, 1.0, 4000.0, true, "Avv==1 && Zvv==0 && Ityp==1"), true);

  vm5->plot(cppAnalysis->fillHisto("Tetlab", 180, 0.0, 180.0, false, "Avv==1 && Zvv == 0 && Ityp==1"));
  vm5->plot(fortranAnalysis->fillHisto("Tetlab", 180, 0.0, 180.0, false, "Avv==1 && Zvv==0 && Ityp==1"), true);

  vm5->plot(cppAnalysis->fillHisto("Philab", 100, -180.0, 180.0, false, "Avv==1 && Zvv == 0 && Ityp==1"));
  vm5->plot(fortranAnalysis->fillHisto("Philab", 100, -180.0, 180.0, false, "Avv==1 && Zvv==0 && Ityp==1"), true);

  ViewManager *vm6 = new ViewManager(3,2);
  Double_t Emax = 1000.0;
  Double_t Emin = 1.0;
  Bool_t logE = true;

  vm6->plot(cppAnalysis->fillCascadeNeutronDoubleDiffXS(1.0, 1.0, logE, Emin, Emax));
  vm6->plot(fortranAnalysis->fillCascadeNeutronDoubleDiffXS(1.0, 1.0, logE, Emin, Emax), true);

  vm6->plot(cppAnalysis->fillCascadeNeutronDoubleDiffXS(15.0, 2.0, logE, Emin, Emax));
  vm6->plot(fortranAnalysis->fillCascadeNeutronDoubleDiffXS(15.0, 2.0, logE, Emin, Emax), true);

  vm6->plot(cppAnalysis->fillCascadeNeutronDoubleDiffXS(20.0, 2.0, logE, Emin, Emax));
  vm6->plot(fortranAnalysis->fillCascadeNeutronDoubleDiffXS(20.0, 2.0, logE, Emin, Emax), true);

  vm6->plot(cppAnalysis->fillCascadeNeutronDoubleDiffXS(30.0, 2.0, logE, Emin, Emax));
  vm6->plot(fortranAnalysis->fillCascadeNeutronDoubleDiffXS(30.0, 2.0, logE, Emin, Emax), true);

  vm6->plot(cppAnalysis->fillCascadeNeutronDoubleDiffXS(50.0, 2.0, logE, Emin, Emax));
  vm6->plot(fortranAnalysis->fillCascadeNeutronDoubleDiffXS(50.0, 2.0, logE, Emin, Emax), true);

  vm6->plot(cppAnalysis->fillCascadeNeutronDoubleDiffXS(80.0, 2.0, logE, Emin, Emax));
  vm6->plot(fortranAnalysis->fillCascadeNeutronDoubleDiffXS(80.0, 2.0, logE, Emin, Emax), true);

  ViewManager *vm7 = new ViewManager(3,2);
  Double_t Emax = 1000.0;
  Double_t Emin = 1.0;
  Bool_t logE = true;

  vm7->plot(cppAnalysis->fillCascadeNeutronDoubleDiffXS(90.0, 2.0, logE, Emin, Emax));
  vm7->plot(fortranAnalysis->fillCascadeNeutronDoubleDiffXS(90.0, 2.0, logE, Emin, Emax), true);

  vm7->plot(cppAnalysis->fillCascadeNeutronDoubleDiffXS(100.0, 2.0, logE, Emin, Emax));
  vm7->plot(fortranAnalysis->fillCascadeNeutronDoubleDiffXS(100.0, 2.0, logE, Emin, Emax), true);

  vm7->plot(cppAnalysis->fillCascadeNeutronDoubleDiffXS(110.0, 2.0, logE, Emin, Emax));
  vm7->plot(fortranAnalysis->fillCascadeNeutronDoubleDiffXS(110.0, 2.0, logE, Emin, Emax), true);

  vm7->plot(cppAnalysis->fillCascadeNeutronDoubleDiffXS(140.0, 2.0, logE, Emin, Emax));
  vm7->plot(fortranAnalysis->fillCascadeNeutronDoubleDiffXS(140.0, 2.0, logE, Emin, Emax), true);

  vm7->plot(cppAnalysis->fillCascadeNeutronDoubleDiffXS(150.0, 2.0, logE, Emin, Emax));
  vm7->plot(fortranAnalysis->fillCascadeNeutronDoubleDiffXS(150.0, 2.0, logE, Emin, Emax), true);

  vm7->plot(cppAnalysis->fillCascadeNeutronDoubleDiffXS(160.0, 2.0, logE, Emin, Emax));
  vm7->plot(fortranAnalysis->fillCascadeNeutronDoubleDiffXS(160.0, 2.0, logE, Emin, Emax), true);

  ViewManager *vm8 = new ViewManager(3,2);
  Double_t Emax = 1000.0;
  Double_t Emin = 1.0;
  Bool_t logE = true;

  vm8->plot(cppAnalysis->fillEvaporationNeutronDoubleDiffXS(1.0, 1.0, logE, Emin, Emax));
  vm8->plot(fortranAnalysis->fillEvaporationNeutronDoubleDiffXS(1.0, 1.0, logE, Emin, Emax), true);

  vm8->plot(cppAnalysis->fillEvaporationNeutronDoubleDiffXS(15.0, 2.0, logE, Emin, Emax));
  vm8->plot(fortranAnalysis->fillEvaporationNeutronDoubleDiffXS(15.0, 2.0, logE, Emin, Emax), true);

  vm8->plot(cppAnalysis->fillEvaporationNeutronDoubleDiffXS(20.0, 2.0, logE, Emin, Emax));
  vm8->plot(fortranAnalysis->fillEvaporationNeutronDoubleDiffXS(20.0, 2.0, logE, Emin, Emax), true);

  vm8->plot(cppAnalysis->fillEvaporationNeutronDoubleDiffXS(30.0, 2.0, logE, Emin, Emax));
  vm8->plot(fortranAnalysis->fillEvaporationNeutronDoubleDiffXS(30.0, 2.0, logE, Emin, Emax), true);

  vm8->plot(cppAnalysis->fillEvaporationNeutronDoubleDiffXS(50.0, 2.0, logE, Emin, Emax));
  vm8->plot(fortranAnalysis->fillEvaporationNeutronDoubleDiffXS(50.0, 2.0, logE, Emin, Emax), true);

  vm8->plot(cppAnalysis->fillEvaporationNeutronDoubleDiffXS(80.0, 2.0, logE, Emin, Emax));
  vm8->plot(fortranAnalysis->fillEvaporationNeutronDoubleDiffXS(80.0, 2.0, logE, Emin, Emax), true);

  ViewManager *vm9 = new ViewManager(3,2);
  Double_t Emax = 1000.0;
  Double_t Emin = 1.0;
  Bool_t logE = true;

  vm9->plot(cppAnalysis->fillEvaporationNeutronDoubleDiffXS(90.0, 2.0, logE, Emin, Emax));
  vm9->plot(fortranAnalysis->fillEvaporationNeutronDoubleDiffXS(90.0, 2.0, logE, Emin, Emax), true);

  vm9->plot(cppAnalysis->fillEvaporationNeutronDoubleDiffXS(100.0, 2.0, logE, Emin, Emax));
  vm9->plot(fortranAnalysis->fillEvaporationNeutronDoubleDiffXS(100.0, 2.0, logE, Emin, Emax), true);

  vm9->plot(cppAnalysis->fillEvaporationNeutronDoubleDiffXS(110.0, 2.0, logE, Emin, Emax));
  vm9->plot(fortranAnalysis->fillEvaporationNeutronDoubleDiffXS(110.0, 2.0, logE, Emin, Emax), true);

  vm9->plot(cppAnalysis->fillEvaporationNeutronDoubleDiffXS(140.0, 2.0, logE, Emin, Emax));
  vm9->plot(fortranAnalysis->fillEvaporationNeutronDoubleDiffXS(140.0, 2.0, logE, Emin, Emax), true);

  vm9->plot(cppAnalysis->fillEvaporationNeutronDoubleDiffXS(150.0, 2.0, logE, Emin, Emax));
  vm9->plot(fortranAnalysis->fillEvaporationNeutronDoubleDiffXS(150.0, 2.0, logE, Emin, Emax), true);

  vm9->plot(cppAnalysis->fillEvaporationNeutronDoubleDiffXS(160.0, 2.0, logE, Emin, Emax));
  vm9->plot(fortranAnalysis->fillEvaporationNeutronDoubleDiffXS(160.0, 2.0, logE, Emin, Emax), true);

  ViewManager *vm10 = new ViewManager(3,2);
  Double_t Emax = 1000.0;
  Double_t Emin = 1.0;
  Bool_t logE = true;

  vm10->plot(cppAnalysis->fillCascadeSpectatorNeutronDoubleDiffXS(1.0, 1.0, logE, Emin, Emax));
  vm10->plot(fortranAnalysis->fillCascadeSpectatorNeutronDoubleDiffXS(1.0, 1.0, logE, Emin, Emax), true);

  vm10->plot(cppAnalysis->fillCascadeSpectatorNeutronDoubleDiffXS(15.0, 2.0, logE, Emin, Emax));
  vm10->plot(fortranAnalysis->fillCascadeSpectatorNeutronDoubleDiffXS(15.0, 2.0, logE, Emin, Emax), true);

  vm10->plot(cppAnalysis->fillCascadeSpectatorNeutronDoubleDiffXS(20.0, 2.0, logE, Emin, Emax));
  vm10->plot(fortranAnalysis->fillCascadeSpectatorNeutronDoubleDiffXS(20.0, 2.0, logE, Emin, Emax), true);

  vm10->plot(cppAnalysis->fillCascadeSpectatorNeutronDoubleDiffXS(30.0, 2.0, logE, Emin, Emax));
  vm10->plot(fortranAnalysis->fillCascadeSpectatorNeutronDoubleDiffXS(30.0, 2.0, logE, Emin, Emax), true);

  vm10->plot(cppAnalysis->fillCascadeSpectatorNeutronDoubleDiffXS(50.0, 2.0, logE, Emin, Emax));
  vm10->plot(fortranAnalysis->fillCascadeSpectatorNeutronDoubleDiffXS(50.0, 2.0, logE, Emin, Emax), true);

  vm10->plot(cppAnalysis->fillCascadeSpectatorNeutronDoubleDiffXS(80.0, 2.0, logE, Emin, Emax));
  vm10->plot(fortranAnalysis->fillCascadeSpectatorNeutronDoubleDiffXS(80.0, 2.0, logE, Emin, Emax), true);

  ViewManager *vm11 = new ViewManager(3,2);
  Double_t Emax = 1000.0;
  Double_t Emin = 1.0;
  Bool_t logE = true;

  vm11->plot(cppAnalysis->fillCascadeSpectatorNeutronDoubleDiffXS(90.0, 2.0, logE, Emin, Emax));
  vm11->plot(fortranAnalysis->fillCascadeSpectatorNeutronDoubleDiffXS(90.0, 2.0, logE, Emin, Emax), true);

  vm11->plot(cppAnalysis->fillCascadeSpectatorNeutronDoubleDiffXS(100.0, 2.0, logE, Emin, Emax));
  vm11->plot(fortranAnalysis->fillCascadeSpectatorNeutronDoubleDiffXS(100.0, 2.0, logE, Emin, Emax), true);

  vm11->plot(cppAnalysis->fillCascadeSpectatorNeutronDoubleDiffXS(110.0, 2.0, logE, Emin, Emax));
  vm11->plot(fortranAnalysis->fillCascadeSpectatorNeutronDoubleDiffXS(110.0, 2.0, logE, Emin, Emax), true);

  vm11->plot(cppAnalysis->fillCascadeSpectatorNeutronDoubleDiffXS(140.0, 2.0, logE, Emin, Emax));
  vm11->plot(fortranAnalysis->fillCascadeSpectatorNeutronDoubleDiffXS(140.0, 2.0, logE, Emin, Emax), true);

  vm11->plot(cppAnalysis->fillCascadeSpectatorNeutronDoubleDiffXS(150.0, 2.0, logE, Emin, Emax));
  vm11->plot(fortranAnalysis->fillCascadeSpectatorNeutronDoubleDiffXS(150.0, 2.0, logE, Emin, Emax), true);

  vm11->plot(cppAnalysis->fillCascadeSpectatorNeutronDoubleDiffXS(160.0, 2.0, logE, Emin, Emax));
  vm11->plot(fortranAnalysis->fillCascadeSpectatorNeutronDoubleDiffXS(160.0, 2.0, logE, Emin, Emax), true);

  ViewManager *vm12 = new ViewManager(3,2);
  Double_t Emax = 1000.0;
  Double_t Emin = 1.0;
  Bool_t logE = true;

  vm12->plot(cppAnalysis->fillSpectatorDeExcitationNeutronDoubleDiffXS(1.0, 1.0, logE, Emin, Emax));
  vm12->plot(fortranAnalysis->fillSpectatorDeExcitationNeutronDoubleDiffXS(1.0, 1.0, logE, Emin, Emax), true);

  vm12->plot(cppAnalysis->fillSpectatorDeExcitationNeutronDoubleDiffXS(15.0, 2.0, logE, Emin, Emax));
  vm12->plot(fortranAnalysis->fillSpectatorDeExcitationNeutronDoubleDiffXS(15.0, 2.0, logE, Emin, Emax), true);

  vm12->plot(cppAnalysis->fillSpectatorDeExcitationNeutronDoubleDiffXS(20.0, 2.0, logE, Emin, Emax));
  vm12->plot(fortranAnalysis->fillSpectatorDeExcitationNeutronDoubleDiffXS(20.0, 2.0, logE, Emin, Emax), true);

  vm12->plot(cppAnalysis->fillSpectatorDeExcitationNeutronDoubleDiffXS(30.0, 2.0, logE, Emin, Emax));
  vm12->plot(fortranAnalysis->fillSpectatorDeExcitationNeutronDoubleDiffXS(30.0, 2.0, logE, Emin, Emax), true);

  vm12->plot(cppAnalysis->fillSpectatorDeExcitationNeutronDoubleDiffXS(50.0, 2.0, logE, Emin, Emax));
  vm12->plot(fortranAnalysis->fillSpectatorDeExcitationNeutronDoubleDiffXS(50.0, 2.0, logE, Emin, Emax), true);

  vm12->plot(cppAnalysis->fillSpectatorDeExcitationNeutronDoubleDiffXS(80.0, 2.0, logE, Emin, Emax));
  vm12->plot(fortranAnalysis->fillSpectatorDeExcitationNeutronDoubleDiffXS(80.0, 2.0, logE, Emin, Emax), true);

  ViewManager *vm13 = new ViewManager(3,2);
  Double_t Emax = 1000.0;
  Double_t Emin = 1.0;
  Bool_t logE = true;

  vm13->plot(cppAnalysis->fillSpectatorDeExcitationNeutronDoubleDiffXS(90.0, 2.0, logE, Emin, Emax));
  vm13->plot(fortranAnalysis->fillSpectatorDeExcitationNeutronDoubleDiffXS(90.0, 2.0, logE, Emin, Emax), true);

  vm13->plot(cppAnalysis->fillSpectatorDeExcitationNeutronDoubleDiffXS(100.0, 2.0, logE, Emin, Emax));
  vm13->plot(fortranAnalysis->fillSpectatorDeExcitationNeutronDoubleDiffXS(100.0, 2.0, logE, Emin, Emax), true);

  vm13->plot(cppAnalysis->fillSpectatorDeExcitationNeutronDoubleDiffXS(110.0, 2.0, logE, Emin, Emax));
  vm13->plot(fortranAnalysis->fillSpectatorDeExcitationNeutronDoubleDiffXS(110.0, 2.0, logE, Emin, Emax), true);

  vm13->plot(cppAnalysis->fillSpectatorDeExcitationNeutronDoubleDiffXS(140.0, 2.0, logE, Emin, Emax));
  vm13->plot(fortranAnalysis->fillSpectatorDeExcitationNeutronDoubleDiffXS(140.0, 2.0, logE, Emin, Emax), true);

  vm13->plot(cppAnalysis->fillSpectatorDeExcitationNeutronDoubleDiffXS(150.0, 2.0, logE, Emin, Emax));
  vm13->plot(fortranAnalysis->fillSpectatorDeExcitationNeutronDoubleDiffXS(150.0, 2.0, logE, Emin, Emax), true);

  vm13->plot(cppAnalysis->fillSpectatorDeExcitationNeutronDoubleDiffXS(160.0, 2.0, logE, Emin, Emax));
  vm13->plot(fortranAnalysis->fillSpectatorDeExcitationNeutronDoubleDiffXS(160.0, 2.0, logE, Emin, Emax), true);
}
