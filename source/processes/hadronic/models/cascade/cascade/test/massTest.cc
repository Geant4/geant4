#include "../cascade/src/G4NucleusModel.cc"
#include "../cascade/src/G4VRegionModel.cc"
#include "../cascade/src/G4RegionModel.cc"
#include "../cascade/src/G4BertiniData.cc"
#include <fstream>


//Test for G4NucleusModel::GetAtomicMass().
//Test keeps the ratio N/Z constant and varies Z & N.

int main(){

  int i, j;
  ifstream in;
  G4double neutronProtonRatio;
  G4int verboseLevel = 0;
  
  if(verboseLevel > 0) cout << "testing G4NucleusModel::GetAtomicMass()" <<endl;
  in.open("massTest.in", ios::in);

  
  G4NucleusModel testNucleus;
  testNucleus.CreateModel(10,5);
  in >> neutronProtonRatio;

  if(verboseLevel > 0) cout << "a test nucleus created " <<endl;
  

  for(i=1; i < 100; i++){
    j = static_cast<int>(neutronProtonRatio * i);
    if((i+j) < 273){
    
      if(verboseLevel > 0)
	cout << "A: " << j + i << " Z: " << i << " ";
      if(verboseLevel == 0) cout << i <<" ";

      testNucleus.SetParameters((i+j), i);
    
      cout << testNucleus.GetAtomicMass() << endl;
    }
    else
      break;
  }
  return 0;
}






