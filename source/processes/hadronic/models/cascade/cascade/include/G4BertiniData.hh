#ifndef G4BertiniData_h
#define G4BertiniData_h 1

#include "globals.hh"

class G4BertiniData //: public G4VIntraNuclearTransportModel 
{
  public:
      G4BertiniData(){ 
	G4cout << " G4BertiniData constructor" << G4endl;
}
      ~G4BertiniData(){
	G4cout << " G4BertiniData destructor" << G4endl;
}
 
   static G4BertiniData* Instance()
   {
      if (!theInstance) theInstance = new G4BertiniData();
      else G4cout << "singleton instance already created" << G4endl;
      return theInstance;
   }

  private:
    static G4BertiniData* theInstance;
};

#endif










