{
  // run first analyzeEvents.C
  TCanvas *c =new TCanvas("c","c",0,0,600,400);
  c->Divide(3,2);

  c->cd(1);ntuple.Draw("momX");
  c->cd(2); ntuple.Draw("momX:momZ","","surf");
  c->cd(3); ntuple.Draw("momX:momY:momZ","particleId==10");
  c->cd(4); ntuple.Draw("momX:momY:momZ","particleId==1");

  c->cd(5);ntuple.Draw("momZ >>hmomz(100,0,0.01)","particleId==1");
  ntuple.Draw("momZ >>+hmomz","particleId==2");
  c->cd(6);ntuple.Draw("fragmentA:fragmentZ >>haz(20,0,20,30,0,30)");

  TCanvas *e =new TCanvas("e","e",0,0,600,400);
  e->Divide(3,2);
  e->cd(1); ntuple.Draw("particleId");
  e->cd(2); ntuple.Draw("momZ:particleId","","box"); //or cont1
  e->cd(3); ntuple.Draw("kineticEnergy:particleId","","box");   

  e->cd(4); ntuple.Draw("sqrt(momZ*momZ)/sqrt(momX*momX+momY*momY):particleId","","box");   


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
