#include "../cascade/src/G4BertiniData.cc"

//test method for G4HETCData::Instance()

main(){
  G4BertiniData *dataBase;
    
  G4BertiniData *data1 = dataBase->Instance();  //client
  data1->SetVerboseLevel(1);
  data1->ProtonProtonTotalXSec[0] = 1;
  cout <<"Value in first instance: "<< data1->ProtonProtonTotalXSec[0] << endl;
  
  G4BertiniData *data2 = dataBase->Instance(); //uses the first and only inctance of G4BertiniData!
  cout <<"Value in second instance: "<< data2->ProtonProtonTotalXSec[0] << endl;
}






