#include "../cascade/src/G4NucleusModel.cc"
#include "../cascade/src/G4VRegionModel.cc"
#include "../cascade/src/G4RegionModel.cc"
#include "../cascade/src/G4BertiniData.cc"
#include <fstream>

//Test for G4NucleusModel::GetAtomicMass().
// A = N + Z is a constant and is read from file. 

int main(){

  int z, A;
  ifstream in;
  G4int verboseLevel = 0;
  
  if(verboseLevel > 0) cout << "testing G4NucleusModel::GetAtomicMass()" <<endl;
  
  in.open("massTest3.in", ios::in);

  
  G4NucleusModel testNucleus;
  testNucleus.CreateModel(10,5);

  in >> A;
  //  cout << Z <<" " << Nlow <<" " <<Nhigh << endl;

  if(verboseLevel > 0) cout << "a test nucleus created " <<endl;
  
  
  for(z=1; z <= A; z++){
     if(A < 273){
    
      if(verboseLevel > 0)
	cout << "A: " << A << " Z: " << z << " ";
      if(verboseLevel == 0) cout << z <<" ";

      testNucleus.SetParameters(A, z);
    
      cout << testNucleus.GetAtomicMass() << endl;
     }
    else
      break;
  }
  return 0;
}






