{
  // run first analyzeEvents.C
  TCanvas *c =new TCanvas("c","c",0,0,600,400);
  c->Divide(3,2);

  c->cd(1);ntuple.Draw("momX");
  c->cd(2); ntuple.Draw("momX:momZ","","surf");
  c->cd(3); ntuple.Draw("momX:momY:momZ","particleId==10","P");
  c->cd(4); ntuple.Draw("momX:momY:momZ","particleId==1");

 c->cd(5);ntuple.Draw("momZ >>hmomz(100,0,1)","particleId==1");
  ntuple.Draw("momZ >>+hmomz","particleId==2");
  c->cd(6);ntuple.Draw("fragmentA:fragmentZ >>haz(20,0,20,30,0,30)");

  TCanvas *e =new TCanvas("e","e",0,0,600,400);
  e->Divide(3,2);
  e->cd(1); ntuple.Draw("particleId","",""); ntuple.Draw("particleId","","E1 same");
  e->cd(2); ntuple.Draw("momX:particleId","","box"); 
  e->cd(3); ntuple.Draw("momZ:particleId","","box"); //or cont1
  e->cd(4); ntuple.Draw("kineticEnergy:particleId","","box");   
  e->cd(5); ntuple.Draw("sqrt(momZ*momZ)/sqrt(momX*momX+momY*momY):particleId","","box");   

   TH1F *h1 = new TH1F("h1","h1",100,0,1.1);

  e->cd(6); gPad->SetLogy(); 
  // h1->Scale(0.1);
ntuple.Draw("kineticEnergy >>h1","particleId==2"); 
ntuple.Draw("kineticEnergy","particleId==2","E1 same");// neutron spectrum  
// h1->Scale(0.1);
ntuple->SetLineStyle(3);
  ntuple.Draw("kineticEnergy","particleId==1","same"); // neutron spectrum  

  TCanvas *d =new TCanvas("d","d",0,0,600,400);
  d->Divide(3,2);
ntuple->SetLineStyle(2);
ntuple->SetLineWidth(2);
 d->cd(1); gPad->SetLogy(); ntuple.Draw("kineticEnergy"); 
  d->cd(2); ntuple.Draw("kineticEnergy:momX");
  d->cd(3); ntuple.Draw("kineticEnergy:momZ");

  ntuple->Draw(">>myList", "particleId==1");
  TEventList *list = (TEventList*)gDirectory->Get("myList");
  ntuple->SetEventList(list);

ntuple->SetLineStyle(0);
ntuple->SetLineWidth(2);
  d->cd(4); ntuple->Draw("kineticEnergy "); // now draws only protons
  d->cd(5); ntuple->Draw("kineticEnergy:momX ");
  d->cd(6); ntuple->Draw("kineticEnergy:momZ ");

  TCanvas *g =new TCanvas("g","g",0,0,600,400);
  g->cd(1); ntuple.Draw("sqrt(momX*momX+momY*momY+momZ*momZ)");   

}
