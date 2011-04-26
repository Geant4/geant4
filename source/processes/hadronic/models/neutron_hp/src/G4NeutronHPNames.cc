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
// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// 24-Jan-07 Enable to use exact data only and add warnig when substitute file is used T. Koi
// 30-Jan-07 Modified method of searching substitute isotope data by T. Koi
// 07-06-12 fix memory leaking by T. Koi
// 07-06-25 Change data selection logic when G4NEUTRONHP_SKIP_MISSING_ISOTOPES is turn on
//          Natural Abundance data are allowed. by T. Koi
// 07-07-06 Allow _nat_ final state even for isotoped cross sections by T. Koi
// 08-09-01 Add protection that deuteron data do not selected for hydrogen and so on by T. Koi
//
#include "G4NeutronHPNames.hh"
#include "G4SandiaTable.hh"
#include "G4HadronicException.hh"
#include <fstream>

  const G4String G4NeutronHPNames::theString[100] = {"Hydrogen", "Helium",
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
 "Francium", "Radium", "Actinium", "Thorium", "Protactinium", "Uranium", 
 "Neptunium", "Plutonium", "Americium", "Curium", "Berkelium", "Californium",
 "Einsteinium","Fermium"};


  G4String G4NeutronHPNames::GetName(G4int i) { return theString[i]; }

G4NeutronHPDataUsed G4NeutronHPNames::GetName(G4int A, G4int Z, G4String base, G4String rest, G4bool & aFlag)
{

   G4NeutronHPDataUsed result;
   aFlag = true;
if(getenv("NeutronHPNames")) G4cout << "Names::GetName entered for Z = " << Z << ", A = " << A <<G4endl;

    G4int myA = A;
    G4int myZ = Z;

    if(Z>92.5&&!getenv("AllowForHeavyElements") ) 
    {
      G4cerr << "Please contact Hans-Peter.Wellisch@cern.ch"<<G4endl;
      throw G4HadronicException(__FILE__, __LINE__, "G4NeutronHPNames::GetName - data with Z>92 are not provided");
    }

    G4String * theName = 0;
    G4String theFileName("");

//    G4int inc = 1;

    G4int flip_Z = 1;
    G4int delta_Z = 0;

    G4int flip_A = 1;
    G4int delta_A = 0;
    
    std::ifstream * check = new std::ifstream(".dummy");
    G4bool first = true;
if(getenv("NeutronHPNames"))  G4cout << "entered GetName!!!"<<G4endl;
    do   
    {
       aFlag = true;
       G4String * biff = new G4String(); // delete here as theName
       *biff = base+"/"+"CrossSection/"+itoa(myZ)+"_"+itoa(myA)+"_"+theString[myZ-1];
      
       if(theName!=0) delete theName;
       theName = biff;
       result.SetName(*theName);
       result.SetA(myA);
       result.SetZ(myZ);
if(getenv("NeutronHPNames")) G4cout <<"HPWD 1 "<<*theName<<G4endl;

     // T.K. debug for memory leak
     if ( check != 0 )
     {
        check->close();
        delete check;
     } 
       check = new std::ifstream(*theName);
       if ( !(*check) ) 
       {
	  check->close();
	  delete check;
          check = 0;
          aFlag = false;
          if ( first )
          {
             aFlag = true;
             first = false;
             biff = new G4String(); // delete here as theName
             *biff = base+"/"+"CrossSection/"+itoa(myZ)+"_"+"nat"+"_"+theString[myZ-1];
             if(theName!=0) delete theName;
             theName = biff;
if(getenv("NeutronHPNames"))    G4cout <<"HPWD 2 "<<*theName<<G4endl;
             result.SetName(*theName);
             G4double natA = myZ/G4SandiaTable::GetZtoA(myZ);
             result.SetA(natA);
             result.SetZ(myZ);
             check = new std::ifstream(*theName);
             if ( !(*check) ) 
             {
                check->close();
	        delete check;
                check = 0;
                aFlag = false;
             }
             else
             {
                biff = new G4String(); // delete here as theName
                if(theName!=0) delete theName;
                *biff = base+"/"+rest+itoa(myZ)+"_"+"nat"+"_"+theString[myZ-1];  
                theName = biff;
if(getenv("NeutronHPNames"))    G4cout <<"HPWD 3 "<<*theName<<G4endl;
                result.SetName(*theName);
                G4double natA = myZ/G4SandiaTable::GetZtoA(myZ);
                result.SetA(natA);
                result.SetZ(myZ);
                result.SetNaturalAbundanceFlag();
             }
          }
       }
       else
       {
// 070706 T. Koi Modified 
/*
          biff = new G4String(); // delete here as theName
          *biff = base+"/"+rest+itoa(myZ)+"_"+itoa(myA)+"_"+theString[myZ-1];  
          if(theName!=0) delete theName;
          theName = biff;
if(getenv("NeutronHPNames"))    G4cout <<"HPWD 4 "<<*theName<<G4endl;
          result.SetName(*theName);
          result.SetA(myA);
          result.SetZ(myZ);
*/

          G4double tmpA = myA;
          std::ifstream* file = NULL;
          G4String fileName;

          if ( rest == "/CrossSection/" )
          {

             fileName = base+"/"+rest+itoa(myZ)+"_"+itoa(myA)+"_"+theString[myZ-1];
if(getenv("NeutronHPNames"))    G4cout <<"HPWD 4a "<<*theName<<G4endl;

          }
          else
          {

// For FS
             fileName = base+"/"+rest+itoa(myZ)+"_"+itoa(myA)+"_"+theString[myZ-1];
             file = new std::ifstream(fileName);

             if ( *file )
             {

// isotope FS
if(getenv("NeutronHPNames"))    G4cout <<"HPWD 4b1 "<<*theName<<G4endl;
             }
             else
             {

// _nat_ FS
                fileName  = base+"/"+rest+itoa(myZ)+"_"+"nat"+"_"+theString[myZ-1];

                delete file;
                file = new std::ifstream(fileName);
                if ( *file )
                {

// FS neither isotope nor _nat_
if(getenv("NeutronHPNames"))    G4cout <<"HPWD 4b2a "<<*theName<<G4endl;
                   G4double natA = myZ/G4SandiaTable::GetZtoA(myZ);
                   tmpA = natA;
                }
                else
                {
if(getenv("NeutronHPNames"))    G4cout <<"HPWD 4b2c "<<*theName<<G4endl;
                }
             }

             delete file;

          }

          result.SetName(fileName);
          result.SetA(tmpA);
          result.SetZ(myZ);

       }

       do 
       {
//        if (std::abs(myZ-Z)>theMaxOffSet||myZ==0||myA==0)
          if ( delta_Z > theMaxOffSet )
          { 
             //if ( inc > 0 )
             //{
             //   inc*= -1;
             //   myZ = Z;
             //   myA = A;
             //}
             //else
             //{
                G4cout <<"G4NeutronHPNames: Sorry, this material does not come near to any data."<<G4endl;
                G4cout <<"G4NeutronHPNames: Please make sure G4NEUTRONHPDATA points to the" << G4endl;
                G4cout <<"                  directory, the neutron scattering data are located in." << G4endl;
                G4cout << "G4NeutronHPNames: The material was: A="<<A<<", Z="<<Z<<G4endl;
                throw G4HadronicException(__FILE__, __LINE__, "In case the data sets are at present not available in the neutron data library, please contact Hans-Peter.Wellisch@cern.ch");
                delete theName;
                theFileName = "";
                return result;
             //}
          }

          //if ( std::abs( myA - A ) > theMaxOffSet )
          if ( delta_A > 2*theMaxOffSet )
          {
             delta_A = 0;
             flip_A = 1;

             first = true;

             if ( flip_Z > 0 ) 
             {
                delta_Z +=1; 
             }
             myZ = Z + flip_Z * delta_Z;
             flip_Z *= -1;
             
             myA = A;
             if ( myZ > 99 ) 
             {
                myZ = 99;
             }
             if ( myZ < 1 ) 
             {
                myZ = 1;
             }
              
//             myZ += inc;
          }
          else
          {
             if ( flip_A > 0 )
             {
                delta_A += 1;
             }
             myA = A + flip_A * delta_A; 
             flip_A *= -1;

             if ( myA < 1 ) 
             {
                myA = 1;
             }
              
//             myA += inc;
          }

       }
       while( myZ == 0 || myA == 0 );  // No meaning 

    }
    while((!check) || (!(*check)));

    if(getenv("NeutronHPNamesLogging") || getenv("NeutronHPNames")) 
    {
      G4cout << "Names::GetName: last theName proposal = "<< G4endl;
      G4cout << *theName <<" "<<A<<" "<<Z<<" "<<result.GetName()<<G4endl;
    }

// administration and anouncement for lacking of exact data in NDL 
    if ( Z != result.GetZ() || A != result.GetA() )
    {
       if ( rest == "/CrossSection/" )
       {
          G4String reac = base;
          G4String dir = getenv("G4NEUTRONHPDATA"); 
          reac.erase ( 0 , dir.length() );
          if ( getenv ( "G4NEUTRONHP_SKIP_MISSING_ISOTOPES" ) && !( Z == result.GetZ() && result.IsThisNaturalAbundance() ) )
          {
             G4cout << "NeutronHP: " << reac << " file for Z = " << Z << ", A = " << A << " is not found and CrossSection set to 0." << G4endl;
             G4String new_name = base+"/"+rest+"0_0_Zero";  
             result.SetName( new_name );
          }
          else
          { 
             //080901 Add protection that deuteron data do not selected for hydrogen and so on by T. Koi
             if ( (reac.find("Inelastic") != reac.size() && 
                   ((Z == 1 && A == 1) || (Z == 2 && A == 4) ) ) 
                 ||   
                  (reac.find("Capture") != reac.size() && (Z == 2 && A == 4) ) )
             {
                G4String new_name = base+"/"+rest+"0_0_Zero";
                result.SetName( new_name );
             }
             else
             {
                G4cout << "NeutronHP: " << reac << " file for Z = " << Z << ", A = " << A << " is not found and NeutronHP will use " << result.GetName() << G4endl;
             }
          }
       }
    }

    delete theName;
    if(aFlag)
    {
      check->close();
      delete check;
      check = 0;
    }
    return result;
  }
