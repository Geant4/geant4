{
gROOT->ProcessLine(".x $G4INSTALL/tests/test46/scripts/Style.C");
TCanvas c1("c1"," ",0.5, 5, 800, 600);

iac[1] = 0;
iac[2] = 0;

//for(iener = 0; iener<3; iener++) {
iener = 2;
c1.Divide(2,2);
for(iplot=0; iplot<nplot; iplot++) {
gROOT->ProcessLine(".x $G4INSTALL/tests/test46/scripts/PlotSingle.C");
}
c1.Print("a"+fil[iener] + ".gif");
//}
}
