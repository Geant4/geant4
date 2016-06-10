{
gROOT->ProcessLine(".x $G4INSTALL/examples/extended/hadronic/Hadr00/scripts/Style.C");

for(ipart = 0; ipart<4; ipart++) {
for(ixs=0; ixs<7; ixs++) {
  bool yes = true;
  if(ipart != 1 && (ixs == 4 || ixs == 5)) yes = false;
  if(yes)gROOT->ProcessLine(".x $G4INSTALL/examples/extended/hadronic/Hadr00/scripts/PlotSingle.C");
}
}
}
