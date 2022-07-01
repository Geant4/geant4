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
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
// June-2019 - E. Mendoza --> Modification to allow using an incomplete data library if the G4NEUTRONHP_SKIP_MISSING_ISOTOPES environmental flag is defined. The missing XS are set to 0.
// Oct-2019 - E. Mendoza --> remove restriction of using isotopes with Z>92

#include "G4ParticleHPNames.hh"
#include "G4ParticleHPManager.hh"
#include "G4SandiaTable.hh"
#include "G4HadronicException.hh"
#include "G4HadronicParameters.hh"
#include <fstream>

  const G4String G4ParticleHPNames::theString[100] = {"Hydrogen", "Helium",
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


G4String G4ParticleHPNames::GetName(G4int i) { return theString[i]; }

G4ParticleHPDataUsed G4ParticleHPNames::GetName(G4int A, G4int Z, G4int M, G4String base, G4String rest, G4bool & aFlag)
{
  #ifdef G4VERBOSE
  G4int verboseLevel = G4ParticleHPManager::GetInstance()->GetVerboseLevel();
  #endif

   //G4cout << Z << " " << A << " " << M << " " << base << " " << rest << G4endl;

  //Excited isomer indicator
  std::stringstream ss;
  G4String sM;
  if (M > 0) {
    ss << "m";
    ss << M;
    ss >> sM;
    ss.clear();
  }

  G4ParticleHPDataUsed result;
  aFlag = true;
  
  #ifdef G4VERBOSE
  if (std::getenv("NeutronHPNames") && G4HadronicParameters::Instance()->GetVerboseLevel() > 0)
    G4cout << "Names::GetName entered for Z = " << Z << ", A = " << A <<G4endl;
  #endif
  
  G4int myA = A;
  G4int myZ = Z;

  G4String * theName = 0;
  G4String theFileName("");

  //G4int inc = 1;

  G4int flip_Z = 1;
  G4int delta_Z = 0;

  G4int flip_A = 1;
  G4int delta_A = 0;
    
  //std::ifstream * check = new std::ifstream(".dummy");
  std::istringstream* check = NULL;
  G4bool first = true;

  #ifdef G4VERBOSE
  if (std::getenv("NeutronHPNames") && G4HadronicParameters::Instance()->GetVerboseLevel() > 0)
    G4cout << "entered GetName!!!"<<G4endl;
  #endif

  do   
    {
       aFlag = true;
       G4String * biff = new G4String(); // delete here as theName
       *biff = base+"/CrossSection/"+itoa(myZ)+"_"+itoa(myA)+sM+"_"+theString[myZ-1];
      
       if(theName!=0) delete theName;
       theName = biff;
       result.SetName(*theName);
       result.SetA(myA);
       result.SetZ(myZ);
       result.SetM(M);
       //if(std::getenv("NeutronHPNames")) G4cout <<"HPWD 1 "<<*theName<<G4endl;

       // T.K. debug for memory leak
       if ( check != NULL ) {
         //check->close();
         delete check;
       } 

       //check = new std::ifstream(*theName);
       check = new std::istringstream(std::ios::in);
       G4ParticleHPManager::GetInstance()->GetDataStream2(*theName,*check);
       if ( !(*check) ) 
       {
	  //check->close();
	  delete check;
          check = 0;
          aFlag = false;
          if ( first )
          {
             aFlag = true;
             first = false;
             biff = new G4String(); // delete here as theName
             *biff = base+"/CrossSection/"+itoa(myZ)+"_"+"nat"+"_"+theString[myZ-1];
             delete theName;
             theName = biff;
             //if(std::getenv("NeutronHPNames")) G4cout <<"HPWD 2 "<<*theName<<G4endl;
             result.SetName(*theName);
             G4double natA = myZ/G4SandiaTable::GetZtoA(myZ);
             result.SetA(natA);
             result.SetZ(myZ);
             result.SetM(M);
             //check = new std::ifstream(*theName);
             check = new std::istringstream(std::ios::in);
             G4ParticleHPManager::GetInstance()->GetDataStream2(*theName,*check);
             if ( !(*check) ) 
             {
                //check->close();
	        delete check;
                check = 0;
                aFlag = false;
             }
             else
             {
                biff = new G4String(); // delete here as theName
                *biff = base+"/"+rest+"/"+itoa(myZ)+"_"+"nat"+"_"+theString[myZ-1];  
                if ( rest=="/CrossSection" ) *biff = base+rest+"/"+itoa(myZ)+"_"+"nat"+"_"+theString[myZ-1];  
                delete theName;
                theName = biff;
                //if(std::getenv("NeutronHPNames")) G4cout <<"HPWD 3 "<<*theName<<G4endl;
                result.SetName(*theName);
                natA = myZ/G4SandiaTable::GetZtoA(myZ);
                result.SetA(natA);
                result.SetZ(myZ);
                result.SetM(M);
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
          if(std::getenv("NeutronHPNames")) G4cout <<"HPWD 4 "<<*theName<<G4endl;
          result.SetName(*theName);
          result.SetA(myA);
          result.SetZ(myZ);
          */

          G4double tmpA = myA;
          //std::ifstream* file = NULL;
          std::istringstream* file = NULL;
          G4String fileName;

          if ( rest == "/CrossSection" )
          {

            //fileName = base+"/"+rest+"/"+itoa(myZ)+"_"+itoa(myA)+sM+"_"+theString[myZ-1];
            fileName = base+rest+"/"+itoa(myZ)+"_"+itoa(myA)+sM+"_"+theString[myZ-1];
            //if(std::getenv("NeutronHPNames")) G4cout <<"HPWD 4a "<<*theName<<G4endl;

          } else {

            // For FS
            fileName = base+"/"+rest+"/"+itoa(myZ)+"_"+itoa(myA)+sM+"_"+theString[myZ-1];
            file = new std::istringstream(std::ios::in);
            G4ParticleHPManager::GetInstance()->GetDataStream2(fileName,*file);

            if (*file) {

              // isotope FS
              //if(std::getenv("NeutronHPNames")) G4cout <<"HPWD 4b1 "<<*theName<<G4endl;
            } else {

              // _nat_ FS
              fileName = base+"/"+rest+"/"+itoa(myZ)+"_"+"nat"+"_"+theString[myZ-1];

              delete file;
              //file = new std::ifstream(fileName);
              file = new std::istringstream(std::ios::in);
              G4ParticleHPManager::GetInstance()->GetDataStream2(fileName,*file);
              if (*file) {

                // FS neither isotope nor _nat_
                //if(std::getenv("NeutronHPNames")) G4cout <<"HPWD 4b2a "<<*theName<<G4endl;
                G4double natA = myZ/G4SandiaTable::GetZtoA(myZ);
                tmpA = natA;
              } else {
                //if(std::getenv("NeutronHPNames")) G4cout <<"HPWD 4b2c "<<*theName<<G4endl;
                fileName="INVALID";
              }
            }

            delete file;
          }

          result.SetName(fileName);
          result.SetA(tmpA);
          result.SetZ(myZ);
          result.SetM(M);
       }

       do 
       {
         if (delta_Z > theMaxOffSet) {
           if (!G4ParticleHPManager::GetInstance()->GetSkipMissingIsotopes()) {
             #ifdef G4VERBOSE
             if ( G4HadronicParameters::Instance()->GetVerboseLevel() > 0 ) {
	       G4cout << "G4ParticleHPNames: There are no data available for some isotopes in this material " << G4endl;
	       G4cout << "G4ParticleHPNames: nor are there data for nearby isotopes." << G4endl;
	       G4cout << "G4ParticleHPNames: Please make sure G4NEUTRONHPDATA points to the directory " << G4endl;
	       G4cout << "G4ParticleHPNames: in which the neutron scattering data are located." << G4endl;
	       G4cout << "G4ParticleHPNames: The material was A = " << A << ", Z = " << Z << G4endl;
	     }
             #endif
	     throw G4HadronicException(__FILE__, __LINE__, "In case the data sets are at present not available in the neutron data library, please contact Hadron Group Coordinator");
           } else {
             check = new std::istringstream(std::ios::in);
             break;
           }
         }

         //if ( std::abs( myA - A ) > theMaxOffSet )
         if (delta_A > 2*theMaxOffSet) {
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
             if ( myZ > 100 ) 
             {
                myZ = 100;
             }
             if ( myZ < 1 ) 
             {
                myZ = 1;
             }
              
             //myZ += inc;
         } else {
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
              
             //myA += inc;
         }

       }
       while( myZ == 0 || myA == 0 );  // No meaning // Loop checking, 11.05.2015, T. Koi

    }
    while((!check) || (!(*check))); // Loop checking, 11.05.2015, T. Koi

    #ifdef G4VERBOSE  
    if( (std::getenv("NeutronHPNamesLogging") || std::getenv("NeutronHPNames")) &&
        G4HadronicParameters::Instance()->GetVerboseLevel() > 0 )
    {
      G4cout << "Names::GetName: last theName proposal = "<< G4endl;
      G4cout << *theName <<" "<<A<<" "<<Z<<" "<<result.GetName()<<G4endl;
    }
    #endif
    
    // administration and anouncement for lacking of exact data in NDL 
    if ( Z != result.GetZ() || A != result.GetA() )
    {
       if ( rest == "/CrossSection" )
       {
          G4String reac = base;
          G4String dir = G4FindDataDir("G4NEUTRONHPDATA");
          reac.erase ( 0 , dir.length() );
          if ( G4ParticleHPManager::GetInstance()->GetSkipMissingIsotopes() && !( Z == result.GetZ() && result.IsThisNaturalAbundance() ) )
          {
             #ifdef G4VERBOSE
             if ( verboseLevel > 0 && G4HadronicParameters::Instance()->GetVerboseLevel() > 0 ) {
                G4cout << "NeutronHP: " << reac << " file for Z = " << Z << ", A = " << A << " is not found and CrossSection set to 0." << G4endl;
             }
             #endif
             G4String new_name = base+"/"+rest+"/"+"0_0_Zero";  
             result.SetName( new_name );
          }
          else
          { 
             //080901 Add protection that deuteron data do not selected for hydrogen and so on by T. Koi
             //160216 Increase protencted isotopes for fixing problem on charged particle HP
             if ( ( reac.find("Inelastic") != reac.size() && ( (Z == 1 && A == 1) || (Z == 1 && A == 2) || (Z == 1 && A == 3) || (Z == 2 && A == 3) || (Z == 2 && A == 4) ) ) 
               || ( reac.find("Capture")   != reac.size() && ( (Z == 1 && A == 3) || (Z == 2 && A == 4) ) )  
               || ( reac.find("Fission")   != reac.size() && ( (Z == 88 && A == 224) || (Z == 88 && A == 225) || (Z == 89 && A == 225) || (Z == 88 && A == 226) ) ) ) 
                   
             {
                G4String new_name = base+"/"+rest+"/"+"0_0_Zero";
                result.SetName( new_name );
             }
             else
             {
	        #ifdef G4VERBOSE
                if ( verboseLevel > 0 && G4HadronicParameters::Instance()->GetVerboseLevel() > 0 ) {
                   G4cout << "NeutronHP: " << reac << " file for Z = " << Z << ", A = " << A << " is not found and NeutronHP will use " << result.GetName() << G4endl;
                }
		#endif
             }
          }
       }
    }

    delete theName;
    if(aFlag)
    {
      //check->close();
      delete check;
      check = NULL;
    }
    return result;
}
