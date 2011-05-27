int sx=500, sy=300; // canvas size
int gx=550, gy=350; // grid for canvas locationSo
enum color {OLD=1, DEV=8, IND=9, };  // black for stable old, gree for development
void coulomb()
{
  // 0 default
  // 1 bullet
  // 2 target
  // 3 G4ElementaryParticleCollider
  // 4 G4IntraNucleiCascader
  // 5 G4NonEquilibriumEvaporator
  // 6 G4EquilibriumEvaporator
  // 7 G4Fissioner
  // 8 G4BigBanger
  // enum particleType { nuclei = 0, proton = 1, neutron = 2, pionPlus = 3, pionMinus = 5, pionZero = 7, foton = 10 };
  // runId, eventId, particleId, modelId, kineticEnergy/eMax, momX/eMax, momY/eMax, momZ/eMax, fragmentA, fragmentZ, exitationEnergy/eMax);

  // run first analyzeEvents.C
  //  Interactivly: ntuple->Scan("modelId:particleId")
  gROOT->SetStyle("clearRetro");

  gSystem->Exec("root -q analyzeEvents.C");           // run batch and return to continue 
  TFile *f = new TFile("../data/analyzeEvents.root"); // get standard file
  TNtuple *ntuple = (TNtuple*)f->Get("ntuple");       // get standard tuple
  ntuple->Print();
  ntuple->Show(1); // show first event
 
  // subroutines
  compare();
  standard();
  fast();
 
}

void compare(){  // compare interface vs cascade
  //______________________________________________________________________
  TCanvas *c6 =new TCanvas("c6","c6", 2*gx, 3*gy, sx, sy);
  c6->Divide(3,4);
  ntuple->SetLineWidth(1);

  c6->cd(1); gPad->SetLogy(); 
  ntuple->SetLineWidth(2); ntuple->SetLineColor(OLD); // old tested  physics with  black (1 and thik line 2)
     ntuple->Draw("particleId", "runId==1"); 
  ntuple->SetLineWidth(1); ntuple->SetLineColor(DEV); 
  ntuple->Draw("particleId", "runId==2", "same"); 

  c6->cd(2); gPad->SetLogy(); 
  ntuple->SetLineWidth(2); ntuple->SetLineColor(OLD); 
  ntuple->Draw("kineticEnergy","particleId==1 && runId==1"); 

  ntuple->SetLineWidth(1); ntuple->SetLineColor(DEV); 
  ntuple->Draw("kineticEnergy","particleId==1 && runId==2", "same"); 
  TLatex l61;   l61.SetTextSize(0.1); l61.DrawLatex(1,1,"Proton");

  c6->cd(3); gPad->SetLogy(); 
  ntuple->SetLineWidth(2); ntuple->SetLineColor(OLD); 
  ntuple->Draw("kineticEnergy","particleId==2 && runId==1"); 
  ntuple->SetLineWidth(1); ntuple->SetLineColor(DEV); 
  ntuple->Draw("kineticEnergy","particleId==2 && runId==2", "same"); 
  TLatex l61;   l61.SetTextSize(0.1); l61.DrawLatex(1,1,"Neutron");


  c6->cd(4); gPad->SetLogy(); 
    ntuple->SetLineWidth(2); ntuple->SetLineColor(OLD); 
  ntuple->Draw("kineticEnergy","particleId==10 && runId==1"); 
  
  ntuple->SetLineWidth(1); ntuple->SetLineColor(DEV);
  ntuple->Draw("kineticEnergy","particleId==10 && runId==2", "same"); 
  TLatex l64;   l64.SetTextSize(0.1); l64.DrawLatex(1,1,"Gamma");

  c6->cd(5); gPad->SetLogy(); 
  ntuple->SetLineWidth(2); ntuple->SetLineColor(OLD);
  ntuple->Draw("kineticEnergy","particleId==3 && runId==1"); 
    ntuple->SetLineWidth(1); ntuple->SetLineColor(DEV); 
   ntuple->Draw("kineticEnergy","particleId==3 && runId==2","same"); 
  TLatex l65;   l65.SetTextSize(0.1); l65.DrawLatex(1,1,"#pi^{+}");

  c6->cd(6); gPad->SetLogy(); 
  ntuple->SetLineWidth(2); ntuple->SetLineColor(OLD);
  ntuple->Draw("kineticEnergy","particleId==5 && runId==1"); 
    ntuple->SetLineWidth(1); ntuple->SetLineColor(DEV); 
   ntuple->Draw("kineticEnergy","particleId==5 && runId==2","same"); 
  TLatex l66;   l66.SetTextSize(0.1); l66.DrawLatex(1,1,"#pi^{-}");

  c6->cd(7); gPad->SetLogy(); 
  ntuple->SetLineWidth(2); ntuple->SetLineColor(OLD);
  ntuple->Draw("kineticEnergy","particleId==7 && runId==1"); 
    ntuple->SetLineWidth(1); ntuple->SetLineColor(DEV); 
   ntuple->Draw("kineticEnergy","particleId==7 && runId==2","same"); 
  TLatex l67;   l67.SetTextSize(0.1); l67.DrawLatex(1,1,"#pi^{0}");

    //  TH1F *h6 = new TH1F("h6","h6",100,0,100);  
  //  ntuple->Draw("kineticEnergy >>h6","particleId==1 && runId==1"); 

    c6->cd(8); gPad->SetLogy(); 
  ntuple->SetLineWidth(2); ntuple->SetLineColor(OLD);
  ntuple->Draw("kineticEnergy","particleId==2 && runId==1"); 
  TLatex l68;
  l68.DrawLatex(1,1,"Neutron");


  c6->cd(9);
  ntuple->SetLineWidth(2); ntuple->SetMarkerColor(OLD);
  ntuple->Draw("momX*momX +momY*momY:kineticEnergy","particleId==1 && runId==1", "box"); 
    ntuple->SetLineWidth(1); ntuple->SetLineColor(DEV);
   ntuple->Draw("momX*momX +momY*momY:kineticEnergy","particleId==1 && runId==2"," box same"); 
  TLatex l69;   l69.SetTextSize(0.1); l69.DrawLatex(1,1,"Ekin vs. pT");


  c6->cd(10);
ntuple->SetLineWidth(2); ntuple->SetLineColor(OLD);
  ntuple->SetLineWidth(2); ntuple->SetMarkerColor(OLD);
   ntuple->Draw("momZ*momZ:kineticEnergy","particleId==1 && runId==1","box"); 
  ntuple->SetLineWidth(1); ntuple->SetMarkerColor(DEV);
  ntuple->SetLineWidth(1); ntuple->SetLineColor(DEV);
    ntuple->Draw("momZ*momZ:kineticEnergy","particleId==1 && runId==2","box same"); 
  TLatex l610;   l610.SetTextSize(0.1); l610.DrawLatex(1,1,"Ekin vs. pII");

  c6->cd(11);
ntuple->SetLineWidth(2); ntuple->SetLineColor(OLD);
  ntuple->SetLineWidth(2); ntuple->SetMarkerColor(OLD);
   ntuple->Draw("momZ*momZ:momX*momX +momY*momY","particleId==1 && runId==1","box"); 
  ntuple->SetLineWidth(1); ntuple->SetMarkerColor(DEV);
  ntuple->SetLineWidth(1); ntuple->SetLineColor(DEV);
    ntuple->Draw("momZ*momZ:momX*momX +momY*momY","particleId==1 && runId==2","box same"); 
  TLatex l611;   l611.SetTextSize(0.1); l611.DrawLatex(1,1,"pT vs. pII");


  ntuple->SetLineWidth(1); ntuple->SetLineColor(1); // Reset to default
};


void fast(){  // analyze kinetic energy spectrum of protons for different p-A confiqurations, using target c9 
  //______________________________________________________________________
  TCanvas *c0 =new TCanvas("c0","c0", 0*gx, 3*gy, sx, sy);
  c0->Divide(3,4);
  ntuple->SetLineWidth(1);
  c0->cd(1); gPad->SetLogy(); 
  ntuple->Draw("kineticEnergy","particleId==1 && runId==5"); 

  ntuple->SetLineWidth(2);
  ntuple->Draw("kineticEnergy","particleId==1 && runId==1","same"); 
  ntuple->Draw("kineticEnergy","particleId==1 && runId==1 && coulombOK==1","same"); 
  ntuple->SetLineWidth(1);
  ntuple->Draw("kineticEnergy","particleId==1 && runId==2","same"); 
  ntuple->Draw("kineticEnergy","particleId==1 && runId==3","same"); 
  ntuple->Draw("kineticEnergy","particleId==1 && runId==4","same"); 
  ntuple->Draw("kineticEnergy","particleId==1 && runId==6","same"); 
  ntuple->Draw("kineticEnergy","particleId==1 && runId==7","same"); 
  ntuple->Draw("kineticEnergy","particleId==1 && runId==8","same"); 
  ntuple->Draw("kineticEnergy","particleId==1 && runId==9","same"); 
  ntuple->Draw("kineticEnergy","particleId==1 && runId==10","same"); 
  ntuple->Draw("kineticEnergy","particleId==1 && runId==11","same"); 

  c0->cd(2); gPad->SetLogy(); 
  ntuple->Draw("kineticEnergy","particleId==1 && runId==2","same"); 
  ntuple->Draw("kineticEnergy","particleId==1 && runId==2 && coulombOK==1","same"); 
  c0->cd(3); gPad->SetLogy(); 
  ntuple->Draw("kineticEnergy","particleId==1 && runId==3","same"); 
  ntuple->Draw("kineticEnergy","particleId==1 && runId==3 && coulombOK==1","same"); 
  c0->cd(4); gPad->SetLogy(); 
  ntuple->Draw("kineticEnergy","particleId==1 && runId==4","same"); 
  ntuple->Draw("kineticEnergy","particleId==1 && runId==4 && coulombOK==1","same"); 
  c0->cd(5); gPad->SetLogy(); 
  ntuple->Draw("kineticEnergy","particleId==1 && runId==5","same"); 
  ntuple->Draw("kineticEnergy","particleId==1 && runId==5 && coulombOK==1","same"); 
  c0->cd(6); gPad->SetLogy(); 
  ntuple->Draw("kineticEnergy","particleId==1 && runId==6","same"); 
  ntuple->Draw("kineticEnergy","particleId==1 && runId==6 && coulombOK==1","same"); 
  c0->cd(7); gPad->SetLogy(); 
  ntuple->Draw("kineticEnergy","particleId==1 && runId==7","same"); 
  ntuple->Draw("kineticEnergy","particleId==1 && runId==7 && coulombOK==1","same"); 
  c0->cd(8); gPad->SetLogy(); 
  ntuple->Draw("kineticEnergy","particleId==1 && runId==8","same"); 
  ntuple->Draw("kineticEnergy","particleId==1 && runId==8 && coulombOK==1","same"); 
  c0->cd(9); gPad->SetLogy(); 
  ntuple->Draw("kineticEnergy","particleId==1 && runId==9","same"); 
  ntuple->Draw("kineticEnergy","particleId==1 && runId==9 && coulombOK==1","same"); 
  c0->cd(10); gPad->SetLogy(); 
  ntuple->Draw("kineticEnergy","particleId==1 && runId==10","same"); 
  ntuple->Draw("kineticEnergy","particleId==1 && runId==10 && coulombOK==1","same"); 
  c0->cd(11); gPad->SetLogy(); 
  ntuple->Draw("kineticEnergy","particleId==1 && runId==11","same"); 
  ntuple->Draw("kineticEnergy","particleId==1 && runId==11 && coulombOK==1","same"); 

  c0->cd(12); gPad->SetLogy(); 
  ntuple->Draw("kineticEnergy","kineticEnergy<0.05 && particleId==1 && runId==11"); 
  ntuple->Draw("kineticEnergy","kineticEnergy<0.05 && particleId==1 && runId==11 && coulombOK==1","same"); 
};


void standard() {
  //______________________________________________________________________
  TCanvas *c1 =new TCanvas("c1","c1", 0*gx,0*gy, sx, sy);
  c1->Divide(3,2);
 ntuple->SetLineWidth(1);
  c1->cd(1); gPad->SetLogy(); ntuple->Draw("modelId","","");

  TLatex l;
  l.SetTextAlign(23);
  l.SetTextSize(0.03);
  l.DrawLatex(1,1,"e^{+}e^{-}#rightarrowZ^{0}#rightarrowI#bar{I}, q#bar{q}");

  c1->cd(2);ntuple->Draw("kineticEnergy","particleId==1");
  //  c1->cd(2); ntuple->Draw("kineticEnergy","particleId==1 && modelId==5");


  c1->cd(3); gPad->SetLogy(); ntuple->Draw("fragmentA");

  //  c1->cd(4);ntuple->Draw("modelId:fragmentA >>haz(20,0,20,30,0,30)");
  c1->cd(4);ntuple->Draw("modelId:fragmentA");

  //  c1->cd(5);ntuple->Draw("momZ >>hmomz(100,0,1)","particleId==1");
  //ntuple->Draw("momZ >>+hmomz","particleId==2");

  //  c1->cd(5);ntuple->Draw("fragmentA:kineticEnergy >>haz(20,0,20,30,0,30)");
  c1->cd(5);ntuple->Draw("kineticEnergy:fragmentA","","box");

  c1->cd(6);ntuple->Draw("fragmentZ:fragmentA","","box");
  //  c1->cd(6);ntuple->Draw("fragmentA:fragmentZ >>haz(20,0,20,30,0,30)");


  //__________________________________________________________________________
  TCanvas *c2 =new TCanvas("c2","c2", 1*gx, 0*gy, sx, sy);
  c2->Divide(3,2);
  c2->cd(1); ntuple->Draw("particleId","",""); ntuple->Draw("particleId","","E1 same");

  c2->cd(2); ntuple->Draw("modelId:particleId", "", "box");   

  c2->cd(3); ntuple->Draw("kineticEnergy:particleId","","box");   

  c2->cd(4); ntuple->Draw("momX:particleId","","box"); 

  c2->cd(5); ntuple->Draw("momY:particleId","","box"); 

  c2->cd(6); ntuple->Draw("momZ:particleId","","box"); 

  //______________________________________________________________________________
  TCanvas *c3 =new TCanvas("c3","c3",0*gx, 1*gy, sx, sy);
  c3->Divide(3,2);
  ntuple->SetLineStyle(1);

  c3->cd(1); gPad->SetLogy(); ntuple->Draw("kineticEnergy >>h0","particleId==1","");   
ntuple->Draw("kineticEnergy","particleId==2","same");

  c3->cd(2);gPad->SetLogy(); ntuple->Draw("kineticEnergy","particleId==2"," ");

  c3->cd(3); gPad->SetLogy(); 
  // h1->Scale(0.1);
  ntuple->Draw("kineticEnergy >>h1","particleId==2"); 
  ntuple->Draw("kineticEnergy","particleId==2","E1 same");
  // h1->Scale(0.1);
  ntuple->SetLineStyle(3);
  ntuple->Draw("kineticEnergy","particleId==1","same"); 


  TH1F *h2 = new TH1F("h2","h2",100,0,1.1);

  c3->cd(4); gPad->SetLogy(); 
  // h1->Scale(0.1);
  ntuple->SetLineStyle(1);
  ntuple->Draw("kineticEnergy >>h2","particleId==1"); 
  ntuple->SetLineStyle(3);
  ntuple->Draw("kineticEnergy","particleId==1 && modelId<6","same"); 

  c3->cd(5); gPad->SetLogy(); 
  // h1->Scale(0.1);
  ntuple->SetLineStyle(1);
  ntuple->Draw("kineticEnergy","particleId==2"); 
  ntuple->SetLineStyle(3);
  ntuple->Draw("kineticEnergy","particleId==2 && coulombOK==1","same");
  TLatex l35;  l35.SetTextSize(0.1); l35.DrawLatex(0,1,"Proton");

  c3->cd(6);  gPad->SetLogx();  gPad->SetLogy(); 
  // h1->Scale(0.1);
  ntuple->SetLineStyle(1);
  ntuple->Draw("kineticEnergy","particleId==2");
  ntuple->SetLineStyle(3);
  ntuple->Draw("kineticEnergy","particleId==2 && coulombOK==1", "same");
  TLatex l36; l36.SetTextSize(0.1); l36.DrawLatex(0,1,"Proton");

  //___________________________________________________________________________________
  TCanvas *c4 =new TCanvas("c4","c4",1*gx,1*gy, sx, sy);
  c4->Divide(3,2);
  TH1F *h3 = new TH1F("h3","h3",100,0,1.1);

  // enum particleType { nuclei = 0, proton = 1, neutron = 2, pionPlus = 3, pionMinus = 5, pionZero = 7, foton = 10 };
  c4->cd(1); gPad->SetLogy(); 

  ntuple->SetLineStyle(1);
  ntuple->Draw("kineticEnergy >>h3","particleId==1"); 
  ntuple->SetLineStyle(2);
  ntuple->SetLineColor(6);
  ntuple->Draw("kineticEnergy","particleId==1 && modelId>5","same"); 
  ntuple->SetLineStyle(3);
  //  ntuple->Draw("kineticEnergy","particleId==1 && modelId>5","same"); 

  c4->cd(2); gPad->SetLogy(); 
  TH1F *h4 = new TH1F("h4","h4",100,0,1.1);
  ntuple->SetLineStyle(1);
  ntuple->SetLineColor(3);
  ntuple->Draw("kineticEnergy >>h4","particleId==2"); 
  ntuple->SetLineStyle(2);
  ntuple->SetLineColor(6);
   ntuple->Draw("kineticEnergy","particleId==2 && modelId>5","same"); 
  ntuple->SetLineStyle(2);
  ntuple->SetLineColor(2);
  //  ntuple->Draw("kineticEnergy","particleId==2 && coulombOK==1","same"); 

  c4->cd(3); gPad->SetLogy(); 
  //  TH1F *h4 = new TH1F("h4","h4",100,0,1.1);
  //ntuple->Draw("kineticEnergy >>h4","particleId==0"); 
ntuple->Draw("kineticEnergy","particleId==0"); 
  ntuple->SetLineStyle(3);
  //ntuple->Draw("kineticEnergy","particleId==0 && coulombOK==1","same"); 

  h3->GetYaxis()->SetLabelOffset(0.00);
  h3->GetYaxis()->SetTitle("d#sigma/dT (mb/MeV)");
    h3->GetXaxis()->SetTitle("T / T Max MeV");
    //h3.GetXaxis()->SetTitle("T / 35 MeV");
 
  c4->cd(4); gPad->SetLogy(); 
  ntuple->Draw("kineticEnergy:fragmentA","particleId==0","box"); 

  c4->cd(5); gPad->SetLogy(); 
  ntuple->Draw("kineticEnergy","fragmentA==2 && particleId==0"); 
  ntuple->SetLineColor(7);
  ntuple->Draw("kineticEnergy","fragmentA==4 && particleId==0","same"); 
  ntuple->SetLineColor(8);
  ntuple->Draw("kineticEnergy","particleId==0 && fragmentA>20 && fragmentA<70","same"); 


  //    c4_1->Print("bi90bert.png");
    // cc_1->Print("sn35bert.png");
    //____________________________________________________________________
  TCanvas *c5 =new TCanvas("c5","c5", 1*gx, 3*gy, sx , sy);
  c5->Divide(3,2);
  ntuple->SetLineStyle(1);
  ntuple->SetLineWidth(2);
  ntuple->SetMarkerStyle(2);
  c5->cd(1); ntuple->Draw("momX:momY:momZ","particleId==1");

  c5->cd(2); ntuple->Draw("kineticEnergy:momX");

  c5->cd(3); ntuple->Draw("kineticEnergy:momZ");

  ntuple->Draw(">>myList", "particleId==1");
  TEventList *list = (TEventList*)gDirectory->Get("myList");
  ntuple->SetEventList(list);
  ntuple->SetLineStyle(0);
  ntuple->SetLineWidth(2);

  c5->cd(4); ntuple->Draw("kineticEnergy:momX:momZ"); // now draws only protons


  //  ntuple->Draw("kineticEnergy ","runId==25","same"); // now draws only protons
  l.DrawLatex(0,1,"e^{+}e^{-}#rightarrowZ^{0}#rightarrowI#bar{I}, q#bar{q}");

  c5->cd(5); ntuple->Draw("kineticEnergy:momX ");
  l.DrawLatex(1,3,"e^{+}e^{-}#rightarrowZ^{0}#rightarrowI#bar{I}, q#bar{q}");

  c5->cd(6); ntuple->Draw("kineticEnergy:momZ ");
  l.DrawLatex(3,1,"e^{+}e^{-}#rightarrowZ^{0}#rightarrowI#bar{I}, q#bar{q}");

}
