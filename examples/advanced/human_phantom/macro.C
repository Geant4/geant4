{
TFile f("g4humanphantom.root");

TDirectory* dir = (TDirectory*)f.Get("ntuple");
TTree* ntuple = (TTree*)dir->Get("ntuple1");
ntuple -> Print();   

Double_t organ_id;
Double_t edep;
   
ntuple->GetBranch("ID")->SetAddress(&organ_id);   
ntuple->GetBranch("Edep")->SetAddress(&edep);   

// Print the content of the ntuple   
Int_t nevent = ntuple->GetEntries();

   for ( Int_t i=0; i<nevent; i++ ) {
     ntuple->GetEvent(i);
     cout << "organ id, edep (MeV): " 
          << organ_id << ", " << edep << endl;
   }
}
