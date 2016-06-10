{
gROOT -> Reset();
TFile f("human_phantom.root");

TDirectory* dir = (TDirectory*)f.Get("human_phantom_ntuple");
TTree* ntuple = (TTree*)dir->Get("1");
ntuple -> Print();   
 
// Print the content of the ntuple  
Int_t nevent = Int_t(ntuple->GetEntries());

Double_t xx;
Double_t edep;
ntuple->GetBranch("organID")->SetAddress(&xx);      
ntuple->GetBranch("edep")->SetAddress(&edep);   
 
for ( Int_t i=0; i<nevent; i++ ) {
     ntuple->GetEvent(i);
     cout << "organ ID, edep (MeV): " 
          << xx << ", " << edep << endl;
   }


}
