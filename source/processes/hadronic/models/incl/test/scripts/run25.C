{
  gROOT->ProcessLine(".x ./scripts/rootlogon.C");
  gROOT->SetStyle("clearRetro");

  gROOT->LoadMacro("./analysis/libAnalysis.C++");

  HistoFactory *hf = new HistoFactory();
  CalculationAnalysis *fortranAnalysis = new CalculationAnalysis(hf, "./tmp/run25ref.root",
							     "h101", 2107.26,
							     100000);

  CalculationAnalysis *cppAnalysis = new CalculationAnalysis(hf, "./tmp/run25.root",
							     "h101", 2107.26,
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
}
