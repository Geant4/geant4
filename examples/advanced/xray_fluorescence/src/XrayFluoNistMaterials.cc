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
//
//
// Author: Alfonso Mantero (Alfonso.Mantero@ge.infn.it)
//
// History:
// -----------
//  20 Aug 2001  Alfonso Mantero   Created
//
// -------------------------------------------------------------------

#include "XrayFluoNistMaterials.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

XrayFluoNistMaterials::XrayFluoNistMaterials()
{ CreateMaterials();}

XrayFluoNistMaterials::~XrayFluoNistMaterials()
{
  delete    dolorite;        
  delete    HPGe;            
  delete    SiLi;            
  delete    mars1;           
  delete    anorthosite;     
  delete    basalt;          
  delete    gabbro;
  delete    gabbroWD;
  delete    gabbroRF;
  delete    Air;             
  delete    Sci;             
  delete    Vacuum;          
  delete    madaBasalt;      
  delete    icelandicBasalt; 
  delete    icelandicWD; 
  delete    icelandicRF; 
  delete    GaAs;            
  delete    galactic;        
  delete    copper; 
  //  delete    nickel;
  delete    hawaiianRF;
  delete    hawaiianWD;         



}
XrayFluoNistMaterials* XrayFluoNistMaterials::instance = 0;

XrayFluoNistMaterials* XrayFluoNistMaterials::GetInstance()
{
  if (instance == 0)
    {
      instance = new XrayFluoNistMaterials;
     
    }
  return instance;
}

G4Material* XrayFluoNistMaterials::GetMaterial(G4String material)  
{

  //instancing G4NistManager
  nistMan = G4NistManager::Instance();
  nistMan->SetVerbose(0);

  //If not available at NIST, look for the local version of it
  G4Material* mat =  nistMan->FindOrBuildMaterial(material);
  if (!mat) {
    mat = G4Material::GetMaterial(material);
  }
  if (!mat) {G4cout << material << "Not Found, Please Retry"<< G4endl;}
  return mat;
}


void XrayFluoNistMaterials::CreateMaterials()
{

  G4double density;            
  std::vector<G4int>  natoms;
  std::vector<G4double> fractionMass;
  std::vector<G4String> elements;


  //instancing G4NistManager
  nistMan = G4NistManager::Instance();
  nistMan->SetVerbose(1);

  // Materials Definitions


  ///////////////////////
  // Madagascar Basalt //
  ///////////////////////


  // Define Madagascar Basalt main components  0054.PP.0044 sample
  density = 3*g/cm3;
  elements.push_back("Si");  fractionMass.push_back(0.1992);  // 0.007093 mol/g(mat)    
  elements.push_back("Ti");  fractionMass.push_back(0.02027); // 4.235e-4               
  elements.push_back("Al");  fractionMass.push_back(0.04758); // 0.001763               
  elements.push_back("Fe");  fractionMass.push_back(0.1303);  // 0.002333               
  elements.push_back("Mn");  fractionMass.push_back(0.001549);// 2.820e-5               
  elements.push_back("Mg");  fractionMass.push_back(0.08141); // 0.003350               
  elements.push_back("Ca");  fractionMass.push_back(0.06468); // 0.001614               
  elements.push_back("Na");  fractionMass.push_back(0.01692); // 7.360e-4               
  elements.push_back("K");   fractionMass.push_back(0.008576);// 2.193e-4               
  elements.push_back("P");   fractionMass.push_back(0.001977);// 6.383e-5               
  elements.push_back("O");   fractionMass.push_back(0.427538);// 0.02672                

   // sum is 0.04434383 total number of moles of atoms in one gram of material
  // 2.248766e8 g per 10.000.000 moles.

  G4Material* madaBasaltMain= nistMan->ConstructNewMaterial("MadaBasaltMain",elements, fractionMass, density);
  elements.clear();
  fractionMass.clear();

// Define Madagascar Basalt traces components  0054.PP.0044 sample
  density = 3*g/cm3;

  
  elements.push_back("Ti");  natoms.push_back(33);  
  elements.push_back("Ba");  natoms.push_back(4131);
  elements.push_back("Ce");  natoms.push_back(694); 
  elements.push_back("Co");  natoms.push_back(965); 
  elements.push_back("Cr");  natoms.push_back(5584);
  elements.push_back("La");  natoms.push_back(269); 
  elements.push_back("Nb");  natoms.push_back(259); 
  elements.push_back("Nd");  natoms.push_back(410); 
  elements.push_back("Ni");  natoms.push_back(389); 
  elements.push_back("Rb");  natoms.push_back(227); 
  elements.push_back("Sc");  natoms.push_back(212); 
  elements.push_back("Sr");  natoms.push_back(8686);
  elements.push_back("V");   natoms.push_back(4203);
  elements.push_back("Y");   natoms.push_back(272); 
  elements.push_back("Zn");  natoms.push_back(1440);
  elements.push_back("Th");  natoms.push_back(19);  
  elements.push_back("Sm");  natoms.push_back(93);  
  elements.push_back("Eu");  natoms.push_back(32);  
  elements.push_back("Gd");  natoms.push_back(89);  
  elements.push_back("Tb");  natoms.push_back(13);  
  elements.push_back("Yb");  natoms.push_back(15);  
  elements.push_back("Lu");  natoms.push_back(2);   
  elements.push_back("Ta");  natoms.push_back(15);  
  elements.push_back("Hf");  natoms.push_back(62);  

  //tot 28114/10e7  weight: 2335253.28 g per 10e6 moles

  G4Material* madaBasaltTraces= nistMan->ConstructNewMaterial("MadaBasaltTraces", elements, natoms, density);
  elements.clear();
  natoms.clear();

  // Define Madacagascar Basalt complete material  0054.PP.0044 sample
  density = 3*g/cm3;
  
  madaBasalt= new G4Material("MadaBasalt", density, 2);
  madaBasalt->AddMaterial(madaBasaltMain,    0.9897);
  madaBasalt->AddMaterial(madaBasaltTraces,  0.0103);



  ///////////////////////////////////////////
  // Iceland    Basalt 0029.PP.0035 sample //
  ///////////////////////////////////////////

  elements.push_back("Si");  fractionMass.push_back(0.2313); 
  elements.push_back("Ti");  fractionMass.push_back(0.0127); 
  elements.push_back("Al");  fractionMass.push_back(0.0702); 
  elements.push_back("Fe");  fractionMass.push_back(0.1134); 
  elements.push_back("Mn");  fractionMass.push_back(0.0019); 
  elements.push_back("Mg");  fractionMass.push_back(0.0349); 
  elements.push_back("Ca");  fractionMass.push_back(0.0756); 
  elements.push_back("Na");  fractionMass.push_back(0.0892); 
  elements.push_back("K");   fractionMass.push_back(0.0032); 
  elements.push_back("P");   fractionMass.push_back(0.00096);
  elements.push_back("S");   fractionMass.push_back(0.0004); 
  elements.push_back("O");   fractionMass.push_back(0.36624);

   // Define Icelandic Basalt main components  0029.PP.0035 sample
  density = 3*g/cm3;
  G4Material* icelandicBasaltMain= nistMan->ConstructNewMaterial("IceBasaltMain",elements, fractionMass, density);
  elements.clear();
  fractionMass.clear();

  // Define Icelandic Basalt traces components  0029.PP.0035 sample
  density = 3*g/cm3;
 
  elements.push_back("Ba");  natoms.push_back(756);  
  elements.push_back("Ce");  natoms.push_back(328);  
  elements.push_back("Co");  natoms.push_back(643);  
  elements.push_back("Cr");  natoms.push_back(1000); 
  elements.push_back("Cu");  natoms.push_back(1396); 
  elements.push_back("Ga");  natoms.push_back(190);  
  elements.push_back("La");  natoms.push_back(103);  
  elements.push_back("Mo");  natoms.push_back(9);    
  elements.push_back("Nb");  natoms.push_back(114);  
  elements.push_back("Nd");  natoms.push_back(104);  
  elements.push_back("Ni");  natoms.push_back(544);  
  elements.push_back("Rb");  natoms.push_back(78);   
  elements.push_back("S");   natoms.push_back(5550); 
  elements.push_back("Sc");  natoms.push_back(531);  
  elements.push_back("Sr");  natoms.push_back(1353); 
  elements.push_back("U");   natoms.push_back(22);   
  elements.push_back("V");   natoms.push_back(4533); 
  elements.push_back("Y");   natoms.push_back(408);  
  elements.push_back("Zn");  natoms.push_back(1259); 
  elements.push_back("Zr");  natoms.push_back(1274); 

  G4Material* icelandicBasaltTraces= nistMan->ConstructNewMaterial("IceBasaltTraces", elements, natoms, density);

  elements.clear();
  natoms.clear();

  // Define Icelandic Basalt complete material  0029.PP.0035 sample
  density = 3*g/cm3;
  icelandicBasalt= new G4Material("IceBasalt", density, 2);
  icelandicBasalt->AddMaterial(icelandicBasaltMain,    0.9978);
  icelandicBasalt->AddMaterial(icelandicBasaltTraces,  0.0022);


  ///////////////////////
  //    Dolorite       //
  ///////////////////////

  // Define dolorite main components 0055.PP.0038 sample

  density = 3*g/cm3;

  elements.push_back("Fe");    fractionMass.push_back(0.1750);
  elements.push_back("Ti");    fractionMass.push_back(0.0082);
  elements.push_back("Ca");    fractionMass.push_back(0.0753);
  elements.push_back("Si");    fractionMass.push_back(0.2188);
  elements.push_back("Al");    fractionMass.push_back(0.0676);
  elements.push_back("Mg");    fractionMass.push_back(0.0008);
  elements.push_back("O");     fractionMass.push_back(0.4377);
  elements.push_back("Mn");    fractionMass.push_back(0.0015);
  elements.push_back("Na");    fractionMass.push_back(0.0134);
  elements.push_back("K");     fractionMass.push_back(0.0011);
  elements.push_back("P");     fractionMass.push_back(0.0006);


  G4Material* doloriteMain = nistMan->ConstructNewMaterial("Dolorite", elements, fractionMass, density);

  elements.clear();
  fractionMass.clear();

  // define traces in dolorite 0055.PP.0038 sample

  density = 3*g/cm3;

  elements.push_back("Nb");    natoms.push_back(5);   
  elements.push_back("Zr");    natoms.push_back(91);  
  elements.push_back("Y");     natoms.push_back(29);  
  elements.push_back("Sr");    natoms.push_back(140); 
  elements.push_back("Rb");    natoms.push_back(3);   
  elements.push_back("Ga");    natoms.push_back(20);  
  elements.push_back("Zn");    natoms.push_back(99);  
  elements.push_back("Ni");    natoms.push_back(77);  
  elements.push_back("Sc");    natoms.push_back(32);  
  elements.push_back("V");     natoms.push_back(314); 
  elements.push_back("Cr");    natoms.push_back(130); 
  elements.push_back("Co");    natoms.push_back(56);  
  elements.push_back("Cu");    natoms.push_back(119); 
  elements.push_back("Ba");    natoms.push_back(38);  
  elements.push_back("Ce");    natoms.push_back(15);  
  elements.push_back("Nd");    natoms.push_back(9);   

  G4Material* tracesOfDolorite= nistMan->ConstructNewMaterial("TracesOfDolorite", elements, natoms, density);

  elements.clear();
  natoms.clear();

  // define dolorite (full) --  0055.PP.0038 sample

  density = 3*g/cm3;
  dolorite = new G4Material("Dolorite", density, 2);
  dolorite->AddMaterial(tracesOfDolorite, 0.0027842352);
  dolorite->AddMaterial(doloriteMain, 0.9972157648);

  ///////////////////////
  //       Mars1       //
  ///////////////////////


 // define mars1 --  01.PP.0030 sample

  density = 3*g/cm3;

  elements.push_back("Fe");    fractionMass.push_back(0.100916);  
  elements.push_back("Ti");    fractionMass.push_back(0.0186804); 
  elements.push_back("Ca");    fractionMass.push_back(0.0404091); 
  elements.push_back("Si");    fractionMass.push_back(0.196378);  
  elements.push_back("Al");    fractionMass.push_back(0.103282);  
  elements.push_back("Mg");    fractionMass.push_back(0.0241622); 
  elements.push_back("Mn");    fractionMass.push_back(0.00184331);
  elements.push_back("Na");    fractionMass.push_back(0.0177908); 
  elements.push_back("K");     fractionMass.push_back(0.00574498);
  elements.push_back("P");     fractionMass.push_back(0.00280169);
  elements.push_back("O");     fractionMass.push_back(0.48799152);


  G4Material* mars1Main = nistMan->ConstructNewMaterial("Mars1 Main components", elements, fractionMass, density);

  elements.clear();
  fractionMass.clear();

  elements.push_back("Nb");    natoms.push_back(55);   
  elements.push_back("Zr");    natoms.push_back(433);  
  elements.push_back("Y");     natoms.push_back(58);   
  elements.push_back("Sr");    natoms.push_back(968);  
  elements.push_back("Rb");    natoms.push_back(16);   
  elements.push_back("Ga");    natoms.push_back(24);   
  elements.push_back("Zn");    natoms.push_back(109);  
  elements.push_back("Ni");    natoms.push_back(70);   
  elements.push_back("Sc");    natoms.push_back(21);   
  elements.push_back("V");     natoms.push_back(134);  
  elements.push_back("Cr");    natoms.push_back(141);  
  elements.push_back("Co");    natoms.push_back(30);   
  elements.push_back("Cu");    natoms.push_back(19);   
  elements.push_back("Ba");    natoms.push_back(580);  
  elements.push_back("Pb");    natoms.push_back(4);    
  elements.push_back("S");     natoms.push_back(444);  
  elements.push_back("U");     natoms.push_back(2);   

  density = 3*g/cm3;
  G4Material* tracesOfMars1 = nistMan->ConstructNewMaterial("TracesOfMars1", elements, natoms, density);

  elements.clear();
  natoms.clear();

  density = 3*g/cm3;
  mars1 = new G4Material("Mars1", density, 2);
  mars1->AddMaterial(tracesOfMars1, 0.0044963163);
  mars1->AddMaterial(mars1Main, 0.9955036837);

  /////////////////////////////////
  //  Hawaiian -- WD coposition  //
  /////////////////////////////////

  density = 3*g/cm3;

  elements.push_back("Fe");    fractionMass.push_back(1.1819860E-01);  
  elements.push_back("Ti");    fractionMass.push_back(2.2781000E-02); 
  elements.push_back("Ca");    fractionMass.push_back(4.5026100E-02); 
  elements.push_back("Si");    fractionMass.push_back(2.0518860E-01);  
  elements.push_back("Al");    fractionMass.push_back(1.3285430E-01);  
  elements.push_back("Mg");    fractionMass.push_back(2.4120000E-03); 
  elements.push_back("Na");    fractionMass.push_back(2.2257000E-02); 
  elements.push_back("K");     fractionMass.push_back(4.9812000E-03);
  elements.push_back("O");     fractionMass.push_back(4.4630120E-01);

  hawaiianWD = nistMan->ConstructNewMaterial("HawaiianWD", elements, fractionMass, density);

  elements.clear();
  fractionMass.clear();

  //////////////////////////////////
  //  Hawaiian -- RF composition  //
  //////////////////////////////////

  density = 3*g/cm3;


  elements.push_back("Fe");    fractionMass.push_back(1.1120460E-01);  
  elements.push_back("Ti");    fractionMass.push_back(2.1582000E-02); 
  elements.push_back("Ca");    fractionMass.push_back(4.3596700E-02); 
  elements.push_back("Si");    fractionMass.push_back(2.1313440E-01);  
  elements.push_back("Al");    fractionMass.push_back(1.0374280E-01);  
  elements.push_back("Mg");    fractionMass.push_back(1.9296000E-02); 
  elements.push_back("Na");    fractionMass.push_back(2.8192200E-02); 
  elements.push_back("K");     fractionMass.push_back(5.8114000E-03);
  elements.push_back("P");     fractionMass.push_back(4.8004000E-03);
  elements.push_back("Mn");    fractionMass.push_back(2.3235000E-03);
  elements.push_back("O");     fractionMass.push_back(4.4531600E-01);

  hawaiianRF = nistMan->ConstructNewMaterial("HawaiianRF", elements, fractionMass, density);

  elements.clear();
  fractionMass.clear();

  //////////////////////////////////
  //  Icelandic -- WD composition  //
  //////////////////////////////////

  density = 3*g/cm3;


  elements.push_back("Si");    fractionMass.push_back(2.2949340E-01);    
  elements.push_back("Ti");    fractionMass.push_back(1.1990000E-02); 
  elements.push_back("Al");    fractionMass.push_back(7.0396900E-02);  
  elements.push_back("Fe");    fractionMass.push_back(1.1330280E-01);  
  elements.push_back("Mg");    fractionMass.push_back(3.4974000E-02); 
  elements.push_back("Ca");    fractionMass.push_back(7.5758200E-02); 
  elements.push_back("Na");    fractionMass.push_back(1.8547500E-02);  
  elements.push_back("K");     fractionMass.push_back(3.3208000E-03);
  elements.push_back("O");     fractionMass.push_back(4.4121640E-01);
  
  icelandicWD = nistMan->ConstructNewMaterial("IcelandicWD", elements, fractionMass, density);
  
  elements.clear();
  fractionMass.clear();


  //////////////////////////////////
  //  Icelandic -- RF composition  //
  //////////////////////////////////

  density = 3*g/cm3;


  elements.push_back("Si");    fractionMass.push_back(2.4304800E-01);    
  elements.push_back("Ti");    fractionMass.push_back(1.3788500E-02); 
  elements.push_back("Al");    fractionMass.push_back(6.5103900E-02);  
  elements.push_back("Fe");    fractionMass.push_back(1.1819860E-01);  
  elements.push_back("Mn");    fractionMass.push_back(2.3235000E-03);
  elements.push_back("Mg");    fractionMass.push_back(2.3517000E-02); 
  elements.push_back("Ca");    fractionMass.push_back(8.2190500E-02); 
  elements.push_back("K");     fractionMass.push_back(3.3208000E-03);
  elements.push_back("P");     fractionMass.push_back(1.3092000E-03);
  elements.push_back("O");     fractionMass.push_back(4.4620000E-01);
  
  icelandicRF = nistMan->ConstructNewMaterial("IcelandicRF", elements, fractionMass, density);
  
  elements.clear();
  fractionMass.clear();

  //////////////////////////////////
  //  Gabbro -- WD composition  //
  //////////////////////////////////

  density = 3*g/cm3;

  elements.push_back("Si");    fractionMass.push_back(1.8696000E-01);    
  elements.push_back("Ti");    fractionMass.push_back(2.3380500E-02); 
  elements.push_back("Al");    fractionMass.push_back(4.6049100E-02);  
  elements.push_back("Fe");    fractionMass.push_back(1.2239500E-01);  
  elements.push_back("Mg");    fractionMass.push_back(8.3817000E-02); 
  elements.push_back("Ca");    fractionMass.push_back(1.0720500E-01); 
  elements.push_back("Na");    fractionMass.push_back(5.9352000E-03);  
  elements.push_back("K");     fractionMass.push_back(1.6604000E-03);
  elements.push_back("O");     fractionMass.push_back(4.2259780E-01);
  
  gabbroWD = nistMan->ConstructNewMaterial("GabbroWD", elements, fractionMass, density);
  
  elements.clear();
  fractionMass.clear();

  //////////////////////////////////
  //  Gabbro -- RF composition  //
  //////////////////////////////////

  density = 3*g/cm3;


  elements.push_back("Si");    fractionMass.push_back(1.6826400E-01);    
  elements.push_back("Ti");    fractionMass.push_back(2.2781000E-02); 
  elements.push_back("Al");    fractionMass.push_back(5.8223000E-02);  
  elements.push_back("Fe");    fractionMass.push_back(1.2729080E-01);  
  elements.push_back("Mn");    fractionMass.push_back(1.5490000E-03);
  elements.push_back("Mg");    fractionMass.push_back(8.3817000E-02); 
  elements.push_back("Ca");    fractionMass.push_back(1.1721080E-01); 
  elements.push_back("Na");    fractionMass.push_back(0.0000000E+00);  
  elements.push_back("K");     fractionMass.push_back(1.6604000E-03);
  elements.push_back("P");     fractionMass.push_back(1.7456000E-03);
  elements.push_back("O");     fractionMass.push_back(4.1845840E-01);
  
  gabbroRF = nistMan->ConstructNewMaterial("GabbroRF", elements, fractionMass, density);
  
  elements.clear();
  fractionMass.clear();


  ///////////////////////
  //     Anorthosite   //
  ///////////////////////


  density = 2.8*g/cm3;

  elements.push_back("Fe");    fractionMass.push_back(0.095283);    
  elements.push_back("Mn");    fractionMass.push_back(0.00137086);  
  elements.push_back("Ni");    fractionMass.push_back(5e-5);        
  elements.push_back("Cu");    fractionMass.push_back(5.2e-4);      
  elements.push_back("Na");    fractionMass.push_back(0.017635);    
  elements.push_back("Mg");    fractionMass.push_back(0.0245361);   
  elements.push_back("Al");    fractionMass.push_back(0.0800355);   
  elements.push_back("Si");    fractionMass.push_back(0.232204);    
  elements.push_back("Ca");    fractionMass.push_back(0.0635368);   
  elements.push_back("K");     fractionMass.push_back(0.00464912);  
  elements.push_back("C");     fractionMass.push_back(0.000837803); 
  elements.push_back("P");     fractionMass.push_back(0.00176742);  
  elements.push_back("Ti");    fractionMass.push_back(0.0240879);   
  elements.push_back("Cl");    fractionMass.push_back(0.00014);     
  elements.push_back("Pd");    fractionMass.push_back(0.00001);     
  elements.push_back("Cd");    fractionMass.push_back(0.00018);     
  elements.push_back("Ag");    fractionMass.push_back(0.00048);     
  elements.push_back("S");     fractionMass.push_back(0.00144);     
  elements.push_back("V");     fractionMass.push_back(0.00228);     
  elements.push_back("Ba");    fractionMass.push_back(0.00151);     
  elements.push_back("O");     fractionMass.push_back(0.447026);    

  anorthosite = nistMan->ConstructNewMaterial("Anorthosite", elements, fractionMass, density);

  elements.clear();
  fractionMass.clear();

  ////////////////////////////////////////
  //     Gabbro       0059.PP.0048      //
  ////////////////////////////////////////


  density = 3.0*g/cm3;

  elements.push_back("Si");    fractionMass.push_back(1.8284688E-01);  
  elements.push_back("Ti");    fractionMass.push_back(2.2601150E-02);  
  elements.push_back("Al");    fractionMass.push_back(4.4831710E-02);  
  elements.push_back("Fe");    fractionMass.push_back(1.2578402E-01);  
  elements.push_back("Mn");    fractionMass.push_back(1.3166500E-03);  
  elements.push_back("Mg");    fractionMass.push_back(8.1706500E-02);  
  elements.push_back("Ca");    fractionMass.push_back(1.0506090E-01);  
  elements.push_back("Na");    fractionMass.push_back(5.4900600E-03);  
  elements.push_back("K");     fractionMass.push_back(1.4943600E-03);  
  elements.push_back("P");     fractionMass.push_back(3.4912000E-04);  
  elements.push_back("O");     fractionMass.push_back(4.0651865E-01);

  gabbro = nistMan->ConstructNewMaterial("Gabbro", elements, fractionMass, density);

  elements.clear();
  fractionMass.clear();

  //define gallium arsenide

  elements.push_back("Ga");     natoms.push_back(1);  
  elements.push_back("As");     natoms.push_back(1);   

  density = 5.32 * g/cm3;
  GaAs = nistMan->ConstructNewMaterial("gallium arsenide", elements, natoms, density);

  elements.clear();
  natoms.clear();

  /*
  // define germanium
  
  density = 5.32 * g/cm3;
 
  elements.push_back("Ge");     natoms.push_back(1); 

  G4cout << elements[1] <<", "<<natoms[1] <<", " << elements.size() << ", " << natoms.size() << G4endl;
  

  HPGe = nistMan->ConstructNewMaterial("High Purity Germanium",elements, natoms, density);

  elements.clear();
  natoms.clear();
  */
  //define scintillator

  elements.push_back("C");     natoms.push_back(9);  
  elements.push_back("H");     natoms.push_back(10);   

  density = 1.032*g/cm3;
  Sci = nistMan->ConstructNewMaterial("Scintillator", elements, natoms, density);

  elements.clear();
  natoms.clear();
  
  //define vacuum
  
  density     = universe_mean_density;    //from PhysicalConstants.h
  G4double pressure    = 3.e-18*pascal;
  G4double temperature = 2.73*kelvin;
  Vacuum = new G4Material("Galactic", 1., 1.01*g/mole, density,
				       kStateGas,temperature,pressure);

  elements.clear();
  natoms.clear();

  //define basalt
  density = 3.*g/cm3;

  elements.push_back("Fe");     fractionMass.push_back(0.1200); 
  elements.push_back("Ti");     fractionMass.push_back(0.0160);   
  elements.push_back("Ca");     fractionMass.push_back(0.0750); 
  elements.push_back("Si");     fractionMass.push_back(0.2160);   
  elements.push_back("Al");     fractionMass.push_back(0.0710); 
  elements.push_back("Mg");     fractionMass.push_back(0.0590);   
  elements.push_back("O");      fractionMass.push_back(0.4430); 
 
  basalt = nistMan->ConstructNewMaterial("Basalt", elements, fractionMass, density);

  elements.clear();
  fractionMass.clear();

  
  // define silicon

  density = 2330*kg/m3;

  // workaround for a problem in nistMan: it doesn't like material with a single element.

  elements.push_back("Si");       natoms.push_back(1); 
  elements.push_back("Si");       natoms.push_back(1); 
  
  SiLi = nistMan->ConstructNewMaterial("SiLi",elements, natoms, density);

  elements.clear();
  natoms.clear();
 

  // define copper

  density = 8920*kg/m3;

  // workaround for a problem in nistMan: it doesn't like material with a single element.
  elements.push_back("Cu");       natoms.push_back(1); 
  elements.push_back("Cu");       natoms.push_back(1); 

  copper = nistMan->ConstructNewMaterial("Cu",elements, natoms, density);

  elements.clear();
  natoms.clear();
  /*
  // define nikel

  density = 8908*kg/m3;

  // workaround for a problem in nistMan: it doesn't like material with a single element.
  elements.push_back("Ni");       natoms.push_back(1); 
  elements.push_back("Ni");       natoms.push_back(1); 

  nickel = nistMan->ConstructNewMaterial("Nickel",elements, natoms, density);
  */
  elements.clear();
  natoms.clear();



  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

