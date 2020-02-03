{
gROOT->ProcessLine(".L TG43_relative_dose.C");
Read("Oncura", "livermore");
Read("Oncura", "penelope");
Read("Oncura", "opt0");
Read("Oncura", "opt4");
gROOT->ProcessLine(".x compare_6711.C");
}

