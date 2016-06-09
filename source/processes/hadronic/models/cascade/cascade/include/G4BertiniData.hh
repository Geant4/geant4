//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
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










