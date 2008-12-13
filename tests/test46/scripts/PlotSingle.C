{ 

c1.cd(iplot+1);

hh[iplot] = gPad->DrawFrame(0.,0.,x1[iplot]*ener[iener],y1[iplot],tit[iplot]);
hh[iplot]->GetYaxis()->SetTitle(axtit[1]);
//gPad->SetLogy();
hh[iplot]->GetXaxis()->SetTitle(axtit[0]);
//else  hh[ixs]->GetXaxis()->SetTitle(axtit[0]);

leg[iplot] = new TLegend(0.65, 0.7, 0.9, 0.9);
 
for(idir=0; idir<ndir; idir++) {
gROOT->ProcessLine(".x $G4INSTALL/tests/test46/scripts/AddMC.C");
}
leg[iplot]->Draw("SAME");
cout << "Legend is done " << endl;
//c1->Update();
}
