#include "../src/G4BertiniData.cc"
//test for (p,p) cross sections in class G4BertiniData
// for energy region 0-20 MeV

int main(){

  G4BertiniData data;
  
  for(int j=1; j < 200; j++) 
    cout << 1+0.1*j <<"\t" << data.GetCrossSection(G4ProtonProtonTotal, 1 + 0.1*j) << endl;


  return 0;

} 
