// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPNames.hh"
#include "G4SandiaTable.hh"

  G4NeutronHPNames::G4NeutronHPNames(){theMaxOffSet = 5;}
  G4NeutronHPNames::G4NeutronHPNames(G4int maxOffSet){theMaxOffSet = maxOffSet;}
  G4NeutronHPNames::~G4NeutronHPNames(){}

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
//    G4cout << "Names::GetName entered"<<endl;
    G4int myA = A;
    G4int myZ = Z;
    G4String * theName = NULL;
    G4String theFileName("");
    G4int offA = 0, offZ = 0, inc = 1;
    
    ifstream check;
    G4bool first = true;
//    G4cout << "entered GetName!!!"<<endl;
     do   
     {
      aFlag = true;
      char the1[100] = {""};
      ostrstream ost1(the1, 100, ios::out);
      ost1 <<base<<"/"<<"CrossSection/"<<myZ<<"_"<<myA<<"_"<<theString[myZ-1];
      G4String * biff = new G4String(the1); // delete here as theName
      if(theName!=NULL) delete theName;
      theName = biff;
      result.SetName(*theName);
      result.SetA(myA);
      result.SetZ(myZ);
      check.open(*theName);
      if(!(check)) 
      {
        aFlag = false;
        if(first)
        {
          aFlag = true;
          first = false;
          char the1[100] = {""};
          ostrstream ost1(the1, 100, ios::out);
          ost1 <<base<<"/"<<"CrossSection/"<<myZ<<"_"<<"nat"<<"_"<<theString[myZ-1];
          biff = new G4String(the1); // delete here as theName
          if(theName!=NULL) delete theName;
          theName = biff;
          result.SetName(*theName);
          G4double natA = myZ/G4SandiaTable::GetZtoA(myZ);
          result.SetA(natA);
          result.SetZ(myZ);
          check.open(*theName);
          if (!check) 
          {
            aFlag = false;
          }
          else
          {
            char the1[100] = {""};
            ostrstream ost1(the1, 100, ios::out);
            ost1 <<base<<"/"<<rest<<myZ<<"_"<<"nat"<<"_"<<theString[myZ-1];  
            biff = new G4String(the1); // delete here as theName
            if(theName!=NULL) delete theName;
            theName = biff;
            result.SetName(*theName);
            G4double natA = myZ/G4SandiaTable::GetZtoA(myZ);
            result.SetA(natA);
            result.SetZ(myZ);
          }
        }
      }
      else
      {
        char the1[100] = {""};
        ostrstream ost1(the1, 100, ios::out);
        ost1 <<base<<"/"<<rest<<myZ<<"_"<<myA<<"_"<<theString[myZ-1];  
        biff = new G4String(the1); // delete here as theName
        if(theName!=NULL) delete theName;
        theName = biff;
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
          G4cout <<"G4NeutronHPNames: Sorry, this material does not come near to any data."<<endl;
          G4cout <<"G4NeutronHPNames: Please make sure NeutronHPCrossSections points to the" << endl;
          G4cout <<"                  directory, the neutron scattering data are located in." << endl;
          G4cout << "G4NeutronHPNames: The material was: A="<<A<<", Z="<<Z<<endl;
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
    while(!(check));
//    G4cout << "Names::GetName: last theName proposal = "<< *theName <<" "<<A<<" "<<Z<<endl;
//    G4cout << "File-name: "<<*theName<<endl;
    delete theName;
    return result;
  }
