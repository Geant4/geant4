#include "../cascade/src/G4NucleusModel.cc"
#include "../cascade/src/G4VRegionModel.cc"
#include "../cascade/src/G4RegionModel.cc"
#include "../cascade/src/G4BertiniData.cc"
#include <fstream>


//test for G4NucleusModel::GetBindingEnergy(). A is kept constant and Z is chan//ged. A is read from file "bindingEnergy.in".

int main(){

  int Z, A;
  ifstream in;
  G4int verboseLevel = 0;
  G4double bindingEnergy;

  if(verboseLevel > 0) cout << "testing G4NucleusModel::GetAtomicMass()" <<endl;
  
  in.open("bindingEnergy.in", ios::in);
  
  G4NucleusModel testNucleus;
  testNucleus.CreateModel(10,5);

  in >> A;
  in.close();

  if(verboseLevel > 0) cout << "a test nucleus created " <<endl;
  
  
  for(Z=1; Z <= A; Z++){
     if(A < 273){
       testNucleus.SetParameters(A, Z);
       bindingEnergy = testNucleus.GetBindingEnergy(); 
       if(bindingEnergy > 0){
	 if(verboseLevel > 0) cout << "A: " << A << " Z: " << Z << " ";
	 if(verboseLevel == 0) cout << Z <<" ";
	 
	 cout << bindingEnergy << endl;
       }
     }
    else break;
  }

  return 0;
}
