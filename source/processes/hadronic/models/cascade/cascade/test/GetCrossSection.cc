#include <fstream>
#include "../src/G4BertiniData.cc" 

// test method for G4double G4BertiniData::GetCrossSection()
//interaction = (p,p) single pion production

int main() {

  G4BertiniData data;
  G4int i, stopIndex;
  G4double *energy, *result, initialEnergy, finalEnergy, step;
  ifstream in;
    

  in.open("test/GetCrossSection.in", ios::in);
  //should contain 3 parameter values 
  
  in >> initialEnergy >> step >> finalEnergy;

  if(initialEnergy > finalEnergy || initialEnergy < 0 || step < 0
     || finalEnergy < 0){
    cout << initialEnergy <<" "<< step <<" "<< finalEnergy << endl;
    cout << "Not valid parameter values!" << endl;
    return 1;
  }

  
  stopIndex = static_cast<G4int>((finalEnergy - initialEnergy)/step) + 1;

  energy = new G4double[stopIndex];
  result = new G4double[stopIndex];

 for(i=0; i<stopIndex; i++) energy[i] = initialEnergy + step*i;

 for(i=0; i<stopIndex; i++){



       result[i]=data.GetCrossSection(G4ProtonProtonSingleProd, energy[i]);
 }
   
 for(i=0; i<stopIndex; i++){
   cout << energy[i] <<"\t"<< result[i] << endl;
 
 }

 in.close();
 
   
  return 0;
}


