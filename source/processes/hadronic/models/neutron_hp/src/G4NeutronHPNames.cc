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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPNames.hh"
#include "G4SandiaTable.hh"

  const G4String G4NeutronHPNames::theString[99] = {"Hydrogen", "Helium",
 "Lithium", "Berylium", "Boron", "Carbon", "Nitrogen", "Oxygen", "Fluorine",
 "Neon", "Sodium", "Magnesium", "Aluminum", "Silicon", "Phosphorous", 
 "Sulfur", "Chlorine", "Argon", "Potassium", "Calcium", "Scandium",
 "Titanium", "Vanadium", "Chromium", "Manganese", "Iron", "Cobalt", "Nickel",
 "Copper", "Zinc", "Gallium", "Germanium", "Arsenic", "Selenium", "Bromine",
 "Krypton", "Rubidium", "Strontium", "Yttrium", "Zirconium", "Niobium",
 "Molybdenum", "Technetium", "Ruthenium", "Rhodium", "Palladium", "Silver",
 "Cadmium", "Indium", "Tin", "Antimony", "Tellurium", "Iodine", "Xenon",
 "Cesium", "Barium", "Lanthanum", "Cerium", "Praseodymium", "Neodymium",
 "Promethium", "Samarium", "Europium", "Gadolinium", "Terbium", "Dysprosium",
 "Holmium", "Erbium", "Thulium", "Ytterbium", "Lutetium", "Hafnium",
 "Tantalum", "Tungsten", "Rhenium", "Osmium", "Iridium", "Platinium", "Gold",
 "Mercury", "Thallium", "Lead", "Bismuth", "Polonium", "Astatine", "Radon", 
 "Francium", "Radium", "Actinium ", "Thorium", "Protactinium", "Uranium", 
 "Neptunium", "Plutonium", "Americium", "Curium", "Berkelium", "Californium",
 "Einsteinium"};


  G4NeutronHPDataUsed G4NeutronHPNames::GetName(G4int A, G4int Z, G4String base, G4String rest, G4bool & aFlag)
  {
    G4NeutronHPDataUsed result;
    aFlag = true;
//    G4cout << "Names::GetName entered"<<G4endl;
    G4int myA = A;
    G4int myZ = Z;
    G4String * theName = NULL;
    G4String theFileName("");
    G4int offA = 0, offZ = 0, inc = 1;
    
    G4std::ifstream * check = new G4std::ifstream(".dummy");
    G4bool first = true;
//    G4cout << "entered GetName!!!"<<G4endl;
     do   
     {
      aFlag = true;
      G4String * biff = new G4String(); // delete here as theName
      *biff = base+"/"+"CrossSection/"+itoa(myZ)+"_"+itoa(myA)+"_"+theString[myZ-1];
      
      if(theName!=NULL) delete theName;
      theName = biff;
      result.SetName(*theName);
      result.SetA(myA);
      result.SetZ(myZ);
//  G4cout <<"HPWD 1 "<<*theName<<G4endl;
#ifdef G4USE_STD_NAMESPACE
      check = new G4std::ifstream(*theName);
#else
      check = new G4std::ifstream(*theName,G4std::ios::in|G4std::ios::nocreate);
#endif
      if(!(*check)) 
      {
	check->close();
	delete check;
        aFlag = false;
        if(first)
        {
          aFlag = true;
          first = false;
          biff = new G4String(); // delete here as theName
          *biff = base+"/"+"CrossSection/"+itoa(myZ)+"_"+"nat"+"_"+theString[myZ-1];
          if(theName!=NULL) delete theName;
          theName = biff;
//      G4cout <<"HPWD 2 "<<*theName<<G4endl;
          result.SetName(*theName);
          G4double natA = myZ/G4SandiaTable::GetZtoA(myZ);
          result.SetA(natA);
          result.SetZ(myZ);
#ifdef G4USE_STD_NAMESPACE
      check = new G4std::ifstream(*theName);
#else
      check = new G4std::ifstream(*theName,G4std::ios::in|G4std::ios::nocreate);
#endif
          if (!(*check)) 
          {
	    check->close();
	    delete check;
            aFlag = false;
          }
          else
          {
            biff = new G4String(); // delete here as theName
            if(theName!=NULL) delete theName;
            *biff = base+"/"+rest+itoa(myZ)+"_"+"nat"+"_"+theString[myZ-1];  
            theName = biff;
//      G4cout <<"HPWD 3 "<<*theName<<G4endl;
            result.SetName(*theName);
            G4double natA = myZ/G4SandiaTable::GetZtoA(myZ);
            result.SetA(natA);
            result.SetZ(myZ);
          }
        }
      }
      else
      {
        biff = new G4String(); // delete here as theName
        *biff = base+"/"+rest+itoa(myZ)+"_"+itoa(myA)+"_"+theString[myZ-1];  
        if(theName!=NULL) delete theName;
        theName = biff;
//      G4cout <<"HPWD 4 "<<*theName<<G4endl;
        result.SetName(*theName);
        result.SetA(myA);
        result.SetZ(myZ);
      }
      if (abs(myZ-Z)>theMaxOffSet||myZ==0||myA==0)
        if(inc>0)
        {
          inc*= -1;
          myZ = Z;
          myA = A;
        }else{
          G4cout <<"G4NeutronHPNames: Sorry, this material does not come near to any data."<<G4endl;
          G4cout <<"G4NeutronHPNames: Please make sure NeutronHPCrossSections points to the" << G4endl;
          G4cout <<"                  directory, the neutron scattering data are located in." << G4endl;
          G4cout << "G4NeutronHPNames: The material was: A="<<A<<", Z="<<Z<<G4endl;
          G4Exception("In case the data sets are at present not available in the neutron data library, please contact Hans-Peter.Wellisch@cern.ch");
          delete theName;
          theFileName = "";
          return result;
        }
      if (abs(myA-A)>theMaxOffSet)
      {
        first = true;
        myA = A;
        myZ+=inc;
      }else{
        myA+=inc;
      }
    }
    while(!(*check));
//    G4cout << "Names::GetName: last theName proposal = "<< *theName <<" "<<A<<" "<<Z<<G4endl;
//    G4cout << "File-name: "<<*theName<<G4endl;
    if(getenv("NeutronHPNamesLogging")) G4cout << "Names::GetName: last theName proposal = "<< *theName <<" "<<A<<" "<<Z<<G4endl;
    delete theName;
    return result;
  }
