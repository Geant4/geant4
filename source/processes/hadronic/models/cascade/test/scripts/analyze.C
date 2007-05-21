{
  // run first analyzeEvents.C
  TCanvas *c =new TCanvas("c","c",0,0,600,400);
  c->Divide(3,2);

  c->cd(1);ntuple.Draw("momX");
  c->cd(2); ntuple.Draw("momX:momZ");
  c->cd(3); ntuple.Draw("momX:momY:momZ","particleId==10");
  c->cd(4); ntuple.Draw("momX:momY:momZ","particleId==1");

  c->cd(5);ntuple.Draw("momZ >>hmomz(100,0,0.01)","particleId==1");
  ntuple.Draw("momZ >>+hmomz","particleId==2");
  c->cd(6);ntuple.Draw("fragmentA:fragmentZ >>haz(20,0,20,30,0,30)");

  TCanvas *d =new TCanvas("d","d",0,0,600,400);
  d->Divide(3,2);

  d->cd(1); ntuple.Draw("kineticEnergy");
  d->cd(2); ntuple.Draw("kineticEnergy:momX");
  d->cd(3); ntuple.Draw("kineticEnergy:momZ");

  ntuple->Draw(">>myList", "particleId==1");
  TEventList *list = (TEventList*)gDirectory->Get("myList");
  ntuple->SetEventList(list);

  d->cd(4); ntuple->Draw("kineticEnergy "); // now draws only protons
  d->cd(5); ntuple->Draw("kineticEnergy:momX ");
  d->cd(6); ntuple->Draw("kineticEnergy:momZ ");
}
