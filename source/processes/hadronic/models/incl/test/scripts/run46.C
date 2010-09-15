{
  gROOT->LoadMacro("./analysis/libAnalysis.C++");
  gROOT->ProcessLine(".x ./analysis/rootlogon.C");

  HistoFactory *hf = new HistoFactory();
  CalculationAnalysis *fortranAnalysis = new CalculationAnalysis(hf, "./tmp/run46ref.root", //"../tmp/c135_c_i43_dres_fb_INIT.root",
							     "h101", 2056.1111,
							     500000);

  CalculationAnalysis *cppAnalysis = new CalculationAnalysis(hf, "./tmp/run46.root",
							     "h101", 2056.1111,
							     500000);

  
  ViewManager *vm = new ViewManager(3,2);
  Double_t Emax = 1000.0;
  Double_t Emin = 1.0;
  Bool_t logE = true;
  vm->plot(cppAnalysis->fillNeutronDoubleDiffXS(1.0, 1.0, logE, Emin, Emax));
  vm->plot(fortranAnalysis->fillNeutronDoubleDiffXS(1.0, 1.0, logE, Emin, Emax), true);
  ThreeColumnReader *r0deg =
    new ThreeColumnReader("./data/carbon/c_c_135MeV/C_C_0deg");
  vm->plot(r0deg->getGraph(), true);

  vm->plot(cppAnalysis->fillNeutronDoubleDiffXS(15.0, 2.0, logE, Emin, Emax));
  vm->plot(fortranAnalysis->fillNeutronDoubleDiffXS(15.0, 2.0, logE, Emin, Emax), true);
  ThreeColumnReader *r15deg =
    new ThreeColumnReader("./data/carbon/c_c_135MeV/C_C_15deg");
  vm->plot(r15deg->getGraph(), true);

  vm->plot(cppAnalysis->fillNeutronDoubleDiffXS(30.0, 2.0, logE, Emin, Emax));
  vm->plot(fortranAnalysis->fillNeutronDoubleDiffXS(30.0, 2.0, logE, Emin, Emax), true);
  ThreeColumnReader *r30deg =
    new ThreeColumnReader("./data/carbon/c_c_135MeV/C_C_30deg");
  vm->plot(r30deg->getGraph(), true);

  vm->plot(cppAnalysis->fillNeutronDoubleDiffXS(50.0, 2.0, logE, Emin, Emax));
  vm->plot(fortranAnalysis->fillNeutronDoubleDiffXS(50.0, 2.0, logE, Emin, Emax), true);
  ThreeColumnReader *r50deg =
    new ThreeColumnReader("./data/carbon/c_c_135MeV/C_C_50deg");
  vm->plot(r50deg->getGraph(), true);

  vm->plot(cppAnalysis->fillNeutronDoubleDiffXS(80.0, 2.0, logE, Emin, Emax));
  vm->plot(fortranAnalysis->fillNeutronDoubleDiffXS(80.0, 2.0, logE, Emin, Emax), true);
  ThreeColumnReader *r80deg =
    new ThreeColumnReader("./data/carbon/c_c_135MeV/C_C_80deg");
  vm->plot(r80deg->getGraph(), true);

  vm->plot(cppAnalysis->fillNeutronDoubleDiffXS(110.0, 2.0, logE, Emin, Emax));
  vm->plot(fortranAnalysis->fillNeutronDoubleDiffXS(110.0, 2.0, logE, Emin, Emax), true);
  ThreeColumnReader *r110deg =
    new ThreeColumnReader("./data/carbon/c_c_135MeV/C_C_110deg");
  vm->plot(r110deg->getGraph(), true);
}
