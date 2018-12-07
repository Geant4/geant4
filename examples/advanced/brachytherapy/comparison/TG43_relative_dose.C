void Read(TString source, TString physics_list){
// Create output file geant4_dose.txt with the dose rate distribution, calculated
// with the simulation results containted in brachytherapy.root

gROOT -> Reset();
TString fileName="brachytherapy_"+source+"_"+physics_list+".root";
std::cout<< "Reading " << fileName << std::endl;
//const char * c = fileName.c_str();
TFile f(fileName);
					     
Double_t Seed_length = 0.35; //seed length in cm

Double_t EnergyMap[401]; //2D map of total energy in "radial distance (mm)" and "angle (5 degrees)"
Int_t Voxels[401]; //the number of voxels used to provide dose to each element of the energy map
Double_t normDose[401]; //Energy map divided by voxels used to make cell, normalised to energy deposition at 1cm, 90 degrees
Double_t GeomFunction[401]; //Geometry Function, normalised to the geometry function at the reference point
Double_t GeometryFunctionZero;  //Geometry function at reference point, 1cm and 90 degrees
Double_t beta;  //beta angle for Geometry Function calculation
Double_t R;     //radial distance in cm
Double_t K;     //polar angle in radians
Double_t Radial[401]; //radial dose function
Double_t radius; //radius (mm)
Int_t radInt; //nearest integer of radius (mm)
Int_t numberOfBins=801;

for (int i=0; i <401; i++)
 {
 EnergyMap[i]=0.;
 Voxels[i]=0.;
}

//Build Energy Deposition Map
for (int k=0; k< numberOfBins; k++)
 {
   for (int m=0; m< numberOfBins; m++) 
 {
   Double_t xx_histo = h20->GetXaxis()->GetBinCenter(k);
   Double_t yy_histo = h20->GetYaxis()->GetBinCenter(m);
   Double_t edep_histo=h20->GetBinContent(k, m);
   radius = sqrt(xx_histo*xx_histo+yy_histo*yy_histo);
 //  if ((edep_histo!=0) && radius < 12. && radius > 9) std::cout << "histo: " << xx_histo << ", " << yy_histo 
   //                                                             << ", radius: " << radius <<", edep: "<< edep_histo << std::endl;

    if (radius != 0){
		      radInt = TMath::Nint(4*radius);
		      if ((radInt>0)&&(radInt<=400))
			{
			 EnergyMap[radInt]+= edep_histo;
			 Voxels[radInt]+= 1;
                      //   if (radius < 12. && radius > 9 && edep_histo!=0)std::cout<< "Radius: " << radius << ", radInt:"<<radInt << ", EnergyMap: "<< EnergyMap[radInt]<< ", voxels: " << Voxels[radInt]<< std::endl;
                         
				}
			}

}}

//Create Normalised Dose Map
std::cout << "The energy deposition at the reference point is " << EnergyMap[40] << std::endl;
Double_t tempNormValue = EnergyMap[40]/Voxels[40]; 
//value at 1cm, 90 degrees, the normalisation point
std::cout << "Dose rate ditribution (distances in cm)" << std::endl;

ofstream myfile;
TString outputFileName ="geant4_dose_"+ source+"_"+physics_list+".txt";
//const char * cOutputFileName  = fileName.c_str();
myfile.open(outputFileName);
std::cout << "file " << outputFileName << " is created "<<std::endl; 

for (int i=0; i<=400; i++)
{
 R = double(i)/40; //distance in CM!!!
 if (Voxels[i]>0) normDose[i] = EnergyMap[i]/Voxels[i]/tempNormValue;
    else normDose[i] = 0;

 
            
 if (R>  0.05)
    {
   // cout << R << "     " << normDose[i] << endl;  
    myfile << R <<  "     " << normDose[i] << "\n";                     
    }
}

myfile.close();
}

