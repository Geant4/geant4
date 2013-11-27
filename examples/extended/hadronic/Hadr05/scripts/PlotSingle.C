{ 

c1 = new TCanvas("c1"," ",0.5, 5, 800, 600);

gtit = hed[ixs] + " Cross Section for " + part[ipart];
hh[ixs] = gPad->DrawFrame(x1[ixs], y1[ixs], x2[ixs], y2[ixs], gtit);
hh[ixs]->GetYaxis()->SetTitle(axtit[1]);
gPad->SetLogy();
if(ixs==0 || ixs==2) hh[ixs]->GetXaxis()->SetTitle(axtit[2]);
else  hh[ixs]->GetXaxis()->SetTitle(axtit[0]);

leg[ixs] = new TLegend(0.9, 0.6, 1.1, 0.9);
 
for(itarg=0; itarg<ntarg; itarg++) {
gROOT->ProcessLine(".x $G4INSTALL/examples/extended/hadronic/Hadr00/scripts/AddMC.C");
cout << "Target# " << itarg << "  max# " << ntarg << endl;
}
cout << "Loop is completed ixs= " << ixs << endl;
leg[ixs]->Draw("SAME");
cout << "Legend is done " << endl;
c1->Update();
c1->Print("a"+filp[ipart] + fil0[ixs] + ".gif");
}
