void run45 () {
  TFile *ff = new TFile("tmp/run45ref.root");
  TFile *cf = new TFile("tmp/run45.root");

  TCut noZeroMass("Masp > 0");

  TTree *c = (TTree *) cf->Get("h101");
  c->SetLineColor(kRed);
  TTree *ref = (TTree *) ff->Get("h101");

  TCanvas *c1 = new TCanvas();
  c1->SetTitle("Projectile spectator");
  c1->Divide(2,2);
  c1->cd(1);
  ref->Draw("Masp", noZeroMass);
  c->Draw("Masp", noZeroMass, "same");

  c1->cd(2);
  ref->Draw("Mzsp", noZeroMass);
  c->Draw("Mzsp", noZeroMass, "same");

  c1->cd(3);
  ref->Draw("Exsp", noZeroMass);
  c->Draw("Exsp", noZeroMass, "same");

  c1->cd(4);
  ref->Draw("Bimpact");
  c->Draw("Bimpact", "", "same");

  c1->Print("run45.ps(");

  TCanvas *c3 = new TCanvas();
  c3->SetTitle("Projectile spectator kinematics");
  c3->Divide(2, 2);
  
  c3->cd(1);
  ref->Draw("Tsp", noZeroMass);
  c->Draw("Tsp", noZeroMass, "same");

  c3->cd(2);
  ref->Draw("Pxsp", noZeroMass);
  c->Draw("Pxsp", noZeroMass, "same");

  c3->cd(3);
  ref->Draw("Pysp", noZeroMass);
  c->Draw("Pysp", noZeroMass, "same");

  c3->cd(4);
  ref->Draw("Pzsp", noZeroMass);
  c->Draw("Pzsp", noZeroMass, "same");

  c3->Print("run45.ps");

  TCanvas *c2 = new TCanvas();
  c2->SetTitle("Cascade remnant");
  c2->Divide(2,2);
  c2->cd(1);
  ref->Draw("Massini");
  c->Draw("Massini", "", "same");

  c2->cd(2);
  ref->Draw("Mzini");
  c->Draw("Mzini", "", "same");

  c2->cd(3);
  ref->Draw("Exini");
  c->Draw("Exini", "", "same");

  c2->cd(4);
  ref->Draw("Jremn");
  c->Draw("Jremn", "", "same");

  c2->Print("run45.ps");

  TCanvas *c4 = new TCanvas();
  c4->SetTitle("Debugging plots");
  c4->Divide(2,2);
  TCut zeroes("Pxsp == 0 && Pysp == 0 && Pzsp == 0");
  TCut zeroesForMassZero("Masp == 0");
  c4->cd(1);
  c->Draw("Tsp", zeroes, "");
  ref->Draw("Tsp", zeroes, "same");
  c4->cd(2);
  ref->Draw("Tsp", zeroesForMassZero, "");
  c->Draw("Tsp", zeroesForMassZero, "same");

  c4->Print("run45.ps");

  TCanvas *c5 = new TCanvas();
  c5->SetTitle("Debugging plots");
  c5->Divide(2,2);
  TCut zeroes("Pxsp == 0 && Pysp == 0 && Pzsp == 0");
  TCut zeroesForMassZero("Masp == 0");
  c5->cd(1);
  c->Draw("Ntrack", "", "");
  ref->Draw("Ntrack", "", "same");
  c5->cd(2);
  ref->Draw("Zvv", "Avv == 1", "");
  c->Draw("Zvv", "Avv == 1", "same");
  
  c5->Print("run45.ps");

  TCanvas *c6 = new TCanvas();
  c6->SetTitle("Ityp test (neutrons)");
  c6->Divide(2,2);

  c6->cd(1);
  TCut spectatorNeutronsFromCascade("Avv == 1 && Zvv == 0 && Ityp==-1");
  c->Draw("Enerj", spectatorNeutronsFromCascade);
  ref->Draw("Enerj", spectatorNeutronsFromCascade, "same");

  c6->cd(3);
  c->Draw("Tetlab", spectatorNeutronsFromCascade);
  ref->Draw("Tetlab", spectatorNeutronsFromCascade, "same");

  c6->cd(2);
  TCut nonSpectatorNeutronsFromCascade("Avv == 1 && Zvv == 0 && Ityp==1");
  c->Draw("Enerj", nonSpectatorNeutronsFromCascade);
  ref->Draw("Enerj", nonSpectatorNeutronsFromCascade, "same");

  c6->cd(4);
  c->Draw("Tetlab", nonSpectatorNeutronsFromCascade);
  ref->Draw("Tetlab", nonSpectatorNeutronsFromCascade, "same");
  
  c6->Print("run45.ps");

  TCut neutronsFromFermiBreakUpOfSpectator("Avv == 1 && Zvv == 0 && Ityp==-2");
  TCut neutronsFromFermiBreakUpOfRemnant("Avv == 1 && Zvv == 0 && Ityp==0");
  c6->cd(1);
  ref->Draw("Enerj", neutronsFromFermiBreakUpOfSpectator);
  c->Draw("Enerj", neutronsFromFermiBreakUpOfSpectator, "same");
  
  c6->cd(3);
  ref->Draw("Tetlab", neutronsFromFermiBreakUpOfSpectator);
  c->Draw("Tetlab", neutronsFromFermiBreakUpOfSpectator, "same");

  c6->cd(2);
  ref->Draw("Enerj", neutronsFromFermiBreakUpOfRemnant);
  c->Draw("Enerj", neutronsFromFermiBreakUpOfRemnant, "same");

  c6->cd(4);
  c->Draw("Tetlab", neutronsFromFermiBreakUpOfRemnant);
  ref->Draw("Tetlab", neutronsFromFermiBreakUpOfRemnant, "same");

  c6->Modified();
  c6->Update();

  c6->Print("run45.ps");

  TCut protonsFromFermiBreakUpOfSpectator("Avv == 1 && Zvv == 1 && Ityp==-2");
  TCut protonsFromFermiBreakUpOfRemnant("Avv == 1 && Zvv == 1 && Ityp==0");
  c6->cd(1);
  ref->Draw("Enerj", protonsFromFermiBreakUpOfSpectator);
  c->Draw("Enerj", protonsFromFermiBreakUpOfSpectator, "same");
  
  c6->cd(3);
  ref->Draw("Tetlab", protonsFromFermiBreakUpOfSpectator);
  c->Draw("Tetlab", protonsFromFermiBreakUpOfSpectator, "same");

  c6->cd(2);
  ref->Draw("Enerj", protonsFromFermiBreakUpOfRemnant);
  c->Draw("Enerj", protonsFromFermiBreakUpOfRemnant, "same");

  c6->cd(4);
  c->Draw("Tetlab", protonsFromFermiBreakUpOfRemnant);
  ref->Draw("Tetlab", protonsFromFermiBreakUpOfRemnant, "same");

  c6->Modified();
  c6->Update();

  c6->Print("run45.ps");

  TCanvas *c7 = new TCanvas();
  c7->SetTitle("Ityp test (protons)");
  c7->Divide(2,2);

  c7->cd(1);
  TCut spectatorProtonsFromCascade("Avv == 1 && Zvv == 1 && Ityp==-1");
  c->Draw("Enerj", spectatorProtonsFromCascade);
  ref->Draw("Enerj", spectatorProtonsFromCascade, "same");

  c7->cd(3);
  c->Draw("Tetlab", spectatorProtonsFromCascade);
  ref->Draw("Tetlab", spectatorProtonsFromCascade, "same");

  c7->cd(2);
  TCut nonSpectatorProtonsFromCascade("Avv == 1 && Zvv == 1 && Ityp==1");
  c->Draw("Enerj", nonSpectatorProtonsFromCascade);
  ref->Draw("Enerj", nonSpectatorProtonsFromCascade, "same");

  c7->cd(4);
  c->Draw("Tetlab", nonSpectatorProtonsFromCascade);
  ref->Draw("Tetlab", nonSpectatorProtonsFromCascade, "same");

  c7->Print("run45.ps");

  TCanvas *c666 = new TCanvas();
  c666->SetTitle("This page has been intentionally left empty");

  c666->Print("run45.ps)");
}
