{
gROOT -> Reset();
TFile f("brachytherapy.root");

ntuple -> Print();   
 
Int_t index;
Double_t xx;
Double_t yy;
Double_t zz;
Double_t edep;
ntuple->GetBranch("xx")->SetAddress(&xx);   
ntuple->GetBranch("yy")->SetAddress(&yy);   
ntuple->GetBranch("zz")->SetAddress(&zz);   
ntuple->GetBranch("edep")->SetAddress(&edep);   
 
// Print the content of the ntuple  
/*Int_t nevent = Int_t(ntuple->GetEntries());

for ( Int_t i=0; i<nevent; i++ ) {
     ntuple->GetEvent(i);
     cout << "xx, yy, zz, edep: " 
          << xx << ", " << yy << ", " << zz << ", " << edep << endl;
   }
*/

// The phantom is 30 cm wide along x, y, z
// the voxel size is 1 mm. The number of voxels is 300 along x, y, z

// Plot the energy deposition in the phantom in 3D
TCanvas* c1 = new TCanvas("c1", " ");

TH3F* edepDDistribution3D = new TH3F("h30", "3Dedepxyz", 
				     300, -150, 150, // binning, xmin, xmax, along x direction
				     300, -150, 150, // binning, xmin, xmax, along y direction
				     300, -150, 150);// binning, xmin, xmax, along z direction
						     
gStyle->SetPalette(1); 

ntuple.Draw("xx:yy:zz:edep>>h30", "", "colz");
}
