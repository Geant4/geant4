{
  // run first analyzeEvents.C
  //  Interactivly: ntuple->Scan("modelId:particleId")

  ntuple->Print();
  ntuple->Show(1); // print first

  TCanvas *c =new TCanvas("c","c",0,0,600,400);
  c->Divide(3,2);

  c->cd(1);ntuple.Draw("momX","particleId==2");
  TLatex l;
  l.SetTextAlign(23);
  l.SetTextSize(0.05);
  l.DrawLatex(1,1,"e^{+}e^{-}#rightarrowZ^{0}#rightarrowI#bar{I}, q#bar{q}");

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


  TCanvas *g =new TCanvas("g","g",0,0,600,400);
  g->Divide(3,2);
  ntuple->SetLineStyle(1);
  g->cd(1); ntuple.Draw("sqrt(momX*momX+momY*momY+momZ*momZ)");   
  g->cd(2); ntuple.Draw("modelId");   
  g->cd(3); ntuple.Draw("modelId:particleId", "", "box");   

  TH1F *h2 = new TH1F("h2","h2",100,0,1.1);

  g->cd(4); gPad->SetLogy(); 
  // h1->Scale(0.1);
  ntuple.Draw("kineticEnergy >>h2","particleId==1"); 
  ntuple->SetLineStyle(3);
  ntuple.Draw("kineticEnergy","particleId==1 && modelId<6","same"); // proton spectrum  

  TH1F *h3 = new TH1F("h3","h3",100,0,1.1);
  g->cd(5); gPad->SetLogy(); 
  // h1->Scale(0.1);
  ntuple.Draw("kineticEnergy >>h3","particleId==2"); 
  ntuple->SetLineStyle(3);
  ntuple.Draw("kineticEnergy","particleId==2 && modelId<6","same"); // neutron spectrum  

  ntuple->SetLineStyle(2);
  ntuple.Draw("kineticEnergy","particleId==2 && modelId>=6","same"); // evaporation
  l.SetTextAlign(23); // center
  l.SetTextSize(0.05);
  l.DrawLatex(0.5,20000, "Secondary neutrons from");
  l.DrawLatex(0.5,10000, "Bi (p, X n) 90 MeV");

  l.SetTextAlign(1); // left corner
  l.SetTextSize(0.04);
  l.DrawLatex(0.8, 80, "Total");
  l.DrawLatex(0.3, 5 , "Evaporation");
  l.DrawLatex(0.05, 20, "INC with exitons");

  TLine line(0.8,80,0.6,20);
  line.SetLineWidth(1);
  line.Draw();

  h3.GetYaxis()->SetLabelOffset(0.00);

  h3.GetYaxis()->SetTitle("d#sigma/dE (mb/MeV)");
  h3.GetXaxis()->SetTitle("E_{kin}/E_{90 MeV}");

  g_5->Print("nFromSubModels.eps");
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
  l.DrawLatex(0,1,"e^{+}e^{-}#rightarrowZ^{0}#rightarrowI#bar{I}, q#bar{q}");
  d->cd(5); ntuple->Draw("kineticEnergy:momX ");
  l.DrawLatex(1,3,"e^{+}e^{-}#rightarrowZ^{0}#rightarrowI#bar{I}, q#bar{q}");
  d->cd(6); ntuple->Draw("kineticEnergy:momZ ");
  l.DrawLatex(3,1,"e^{+}e^{-}#rightarrowZ^{0}#rightarrowI#bar{I}, q#bar{q}");
}
