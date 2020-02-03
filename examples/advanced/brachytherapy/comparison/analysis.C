{
// For Flexi source
// eventually convert from ASCII output files to ROOT files
/*
gROOT->ProcessLine(".L convert.C");
ReadASCII("EnergyDeposition_Flexi_livermore.out", "brachytherapy_Flexi_livermore.root");
ReadASCII("EnergyDeposition_Flexi_penelope.out", "brachytherapy_Flexi_penelope.root");
ReadASCII("EnergyDeposition_Flexi_opt0.out", "brachytherapy_Flexi_opt0.root");
ReadASCII("EnergyDeposition_Flexi_opt3.out", "brachytherapy_Flexi_opt3.root");
ReadASCII("EnergyDeposition_Flexi_opt4.out", "brachytherapy_Flexi_opt4.root");
*/

// g(r) is calculated from the root files 
gROOT->ProcessLine(".L TG43_relative_dose.C");
Read("Flexi", "livermore");
Read("Flexi", "penelope");
Read("Flexi", "opt0");
Read("Flexi", "opt3");
Read("Flexi", "opt4");
//gROOT->ProcessLine(".x compare.C");
// or...
// plot all results with the alternative EM physics constructors
gROOT->ProcessLine(".x compare_all.C");
// use compare.C when only one physics approach is used

// For Oncura source
/*gROOT->ProcessLine(".L TG43_relative_dose.C");
Read("Oncura", "livermore");
Read("Oncura", "penelope");
Read("Oncura", "opt0");
Read("Oncura", "opt3");
Read("Oncura", "opt4");
gROOT->ProcessLine(".x compare_6711.C");
//or ...
// gROOT->ProcessLine(".x compare_6711_all.C");
*/
}

