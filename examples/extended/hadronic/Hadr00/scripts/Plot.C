{
gROOT->ProcessLine(".x $G4INSTALL/examples/extended/hadronic/Hadr00/scripts/Style.C");

for(ipart = 0; ipart<4; ipart++) {
for(ixs=0; ixs<4; ixs++) {
gROOT->ProcessLine(".x $G4INSTALL/examples/extended/hadronic/Hadr00/scripts/PlotSingle.C");
}
}
}
