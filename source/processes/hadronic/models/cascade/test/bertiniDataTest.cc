#include "G4ios.hh"
#include "G4BertiniData.hh"

int main()
{
   G4cout << G4endl << "testing G4BertiniData" << G4endl;
   G4BertiniData *db = new G4BertiniData();                 // sever

   G4BertiniData *data1 = db->Instance();                   // client 
   G4BertiniData *data2 = db->Instance();                   // old instance used

   delete db;
   delete data1;
   delete data2;
}    
    

