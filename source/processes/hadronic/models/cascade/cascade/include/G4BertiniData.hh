//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
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










