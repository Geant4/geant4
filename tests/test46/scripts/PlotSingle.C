{ 

cout << "PlotSingle: " << iplot << endl;
c1.cd(iplot+1);

hh[iplot] = gPad->DrawFrame(0.,0.,x1[iplot]*ener[iener],y1[iplot],tit[iplot]);
hh[iplot]->GetYaxis()->SetTitle(axtit[1]);
//gPad->SetLogy();
hh[iplot]->GetXaxis()->SetTitle(axtit[0]);
//else  hh[ixs]->GetXaxis()->SetTitle(axtit[0]);

leg[iplot] = new TLegend(0.6, 0.6, 0.9, 0.9);
leg[iplot]->SetHeader(part[iener]);
 
for(idir=0; idir<ndir; idir++) {
if(iac[idir]>0)gROOT->ProcessLine(".x $G4INSTALL/tests/test46/scripts/AddMC.C");
}
//if(iplot==3) 
leg[iplot]->Draw("SAME");
//c1.Update();
cout << "PlotSingle done " << endl;
}
