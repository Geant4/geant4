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
//
//
// $Id: XrayFluoDetectorConstruction.hh
// GEANT4 tag $Name: xray_fluo-V03-02-00
//
// Author: Alfonso Mantero (Alfonso.Mantero@ge.infn.it)
//
// History:
// -----------
//  20 Aug 2001  Alfonso Mantero   Created
//
// -------------------------------------------------------------------

#include "XrayFluoMaterials.hh"

XrayFluoMaterials::XrayFluoMaterials()
{ CreateMaterials();}

XrayFluoMaterials* XrayFluoMaterials::instance = 0;

XrayFluoMaterials* XrayFluoMaterials::GetInstance()
{
  if (instance == 0)
    {
      instance = new XrayFluoMaterials;
     
    }
  return instance;
}

G4Material* XrayFluoMaterials::GetMaterial(G4String material)  
{
  G4Material* pttoMaterial = G4Material::GetMaterial(material);
  return pttoMaterial;
}


void XrayFluoMaterials::CreateMaterials()
{

  //define elements
  G4String name, symbol;             //a=mass of a mole;
  G4double a, z, density;            //z=mean number of protons;  
  G4int  natoms,ncomponents;
  G4double temperature, pressure;  
  G4double fractionmass;

  // Elements Definitions

  //define Niobium
  
  a = 92.906*g/mole;
  G4Element* Nb  = new G4Element(name="Niobium"  ,symbol="Nb" , z= 41., a);

  //define Zirconium
  
  a = 91.22*g/mole;
  G4Element* Zr  = new G4Element(name="Zirconium"  ,symbol="Zr" , z= 40., a);

  //define Yttrium
  
  a = 88.905*g/mole;
  G4Element* Y  = new G4Element(name="Yttrium"  ,symbol="Y" , z= 39., a);

  //define Stronzium
  
  a = 87.62*g/mole;
  G4Element* Sr  = new G4Element(name="Stronzium"  ,symbol="Sr" , z= 38., a);

  //define Rubidium
  
  a = 85.47*g/mole;
  G4Element* Rb  = new G4Element(name="Rubidium"  ,symbol="Rb" , z= 37., a);

  //define Zinc
  
  a = 65.37*g/mole;
  G4Element* Zn  = new G4Element(name="Zinc"  ,symbol="Zn" , z= 30., a);

  //define Nichel
  
  a = 58.71*g/mole;
  G4Element* Ni  = new G4Element(name="Nichel"  ,symbol="Ni" , z= 28., a);

  //define Scandio
  
  a = 44.956*g/mole;
  G4Element* Sc  = new G4Element(name="Scandium"  ,symbol="Sc" , z= 21., a);

  //define Vanadium
  
  a = 50.942*g/mole;
  G4Element* V  = new G4Element(name="Vanadium"  ,symbol="V" , z= 39., a);

  //define Vanadium
  
  a = 183.84*g/mole;
  G4Element* W  = new G4Element(name="Tungsten"  ,symbol="W" , z= 74., a);


  //define Cromium
 
  a = 51.996*g/mole;
  G4Element* Cr  = new G4Element(name="Cromium"  ,symbol="Cr" , z= 24., a);

  //define Cobalt
  
  a = 58.933*g/mole;
  G4Element* Co  = new G4Element(name="Cobalt"  ,symbol="Co" , z= 27., a);

  //define Copper
  
  a = 63.54*g/mole;
  G4Element* elCu  = new G4Element(name="Copper"  ,symbol="Cu" , z= 29., a);

  //define Barium
  
  a = 137.34*g/mole;
  G4Element* Ba  = new G4Element(name="Barium"  ,symbol="Ba" , z= 56., a);

  //define Cerium
  
  a = 140.12*g/mole;
  G4Element* Ce  = new G4Element(name="Cerium"  ,symbol="Ce" , z= 58., a);

  //define Neodimuim
  
  a = 144.24*g/mole;
  G4Element* Nd  = new G4Element(name="Neodimuim"  ,symbol="Nd" , z= 60., a);


  //Define Zolfo

  a = 32.064*g/mole;
  G4Element* elS  = new G4Element(name="Sulphur"  ,symbol="S" , z= 16., a);


  //define carbon
  
  a = 12.0107*g/mole;
  G4Element* C  = new G4Element(name="Carbon"  ,symbol="C" , z= 6., a);

  //define Nitrogen

  a = 14.01*g/mole;
  G4Element* N  = new G4Element(name="Nitrogen",symbol="N" , z= 7., a);

  // define Oxigen
  a = 15.9994*g/mole;
  G4Element* O  = new G4Element(name="Oxygen"  ,symbol="O" , z= 8., a);

  //define Arsenic
  a = 74.9216 * g/mole;
  G4Element * As = new G4Element( name="arsenic",symbol="As",z= 33.,a);

  //Define Gallium  
  a = 69.72* g/mole;
  G4Element * Ga = new G4Element(name="gallium",symbol="Ga",z= 31.,a);

  //define Iron  
  a = 55.847*g/mole;
  G4Element* Fe = new G4Element(name="Iron"  ,symbol="Fe", z=26., a);
  //define hydrogen
  
  a = 1.01*g/mole;
  G4Element* H  = new G4Element(name="Hydrogen",symbol="H" , z= 1., a);
  
  //define germanium
  a = 72.61*g/mole;
  G4Element* Ge = new G4Element(name="Germanium",symbol="Ge", z= 32.,a);

  //define phosporus
  a = 30.97*g/mole;
  G4Element* P =  new G4Element(name="Phosporus",symbol="P", z= 15., a);

  // define Titanium
  a = 47.88*g/mole;
  G4Element* elTi = new G4Element(name="Titanium",symbol="Ti" , z= 22., a);

  // define Calcium
  a = 40.078*g/mole;
  G4Element* Ca = new G4Element(name="Calcium",symbol="Ca" , z= 20., a);

  // define silicon
  a = 28.0855*g/mole;
  G4Element* elSi = new G4Element(name="Silicon",symbol="Si" , z= 14., a);

  // define Aluminium
  a = 26.98154*g/mole;
  G4Element* elAl = new G4Element(name="Aluminium",symbol="Al" , z= 13., a);

  // Define Magnesium
  a = 24.305*g/mole;
  G4Element* Mg = new G4Element(name="Magnesium",symbol="Mg" , z= 12., a);

  // Define Manganese
  a = 54.938*g/mole;
  G4Element* Mn = new G4Element(name="Manganese",symbol="Mn" , z= 25., a);

  // Define Sodium
  a = 22.989*g/mole;
  G4Element* Na = new G4Element(name="Sodium",symbol="Na" , z= 11., a);

  // Define Potassium
  a = 39.10*g/mole;
  G4Element* K = new G4Element(name="Potassium",symbol="K" , z= 19., a); 
  // Define lead

  a=207.19*g/mole;
  G4Element* elPb = new G4Element(name="Lead",symbol="pb", z=82.,a);

  //define Uranium 

  a =  238.02891*g/mole;
  G4Element* elU  = new G4Element(name="Uranium",symbol="U", z=92.,a);


  // define Palladium
  a= 106.4*g/mole;
  G4Element* Pd = new   G4Element(name="Palladium",symbol="Pd",z=46.,a);

  // define cadmium
  a = 112.4 *g/mole;
  G4Element* Cd = new   G4Element(name="Cadmium",symbol="Cd",z=48.,a);

  // define Silver
  a = 107.87 *g/mole;
  G4Element* Ag = new   G4Element(name="Silver",symbol="Ag",z=47.,a);

  // define Clorine
  a = 35.453 * g/mole;
  G4Element * Cl = new G4Element( name="Chlorine",symbol="Cl",z= 17.,a);

  // define Lantanium
  a = 138.91 * g/mole;
  G4Element * La = new G4Element( name="Lantanium",symbol="La",z= 57.,a);

  // define Molibdenum
  a = 95.94 * g/mole;
  G4Element * Mo = new G4Element( name="Molibdenum",symbol="Mo",z= 42.,a);

  // define Thorium
  a = 232.0381*g/mole;
  G4Element* Th = new G4Element(name="Thorium",symbol="Th" , z= 90., a);

  // define Samarium
  a = 150.36*g/mole;
  G4Element* Sm = new G4Element(name="Samarium",symbol="Sm" , z= 62., a);

  // define Europium
  a = 151.964*g/mole;
  G4Element* Eu = new G4Element(name="Europium",symbol="Eu" , z= 63., a);

  // define Gadolinium
  a = 157.25*g/mole;
  G4Element* Gd = new G4Element(name="Gadolinium",symbol="Gd" , z= 64., a);

  // define Terbium
  a = 158.92534*g/mole;
  G4Element* Tb = new G4Element(name="Terbium",symbol="Tb" , z= 65., a);

  // define Ytterbium
  a = 173.04*g/mole;
  G4Element* Yb = new G4Element(name="Ytterbium",symbol="Yb" , z= 70., a);

  // define Lutetium
  a = 174.967*g/mole;
  G4Element* Lu = new G4Element(name="Lutetium",symbol="Lu" , z= 71., a);

  // define Tantalum
  a = 180.9479*g/mole;
  G4Element* Ta = new G4Element(name="Tantalum",symbol="Ta" , z= 73., a);

  // define Hafnium
  a = 178.49*g/mole;
  G4Element* Hf = new G4Element(name="Hafnium",symbol="Hf" , z= 73., a);

  G4cout << "Elements created" << G4endl;


  // Materials Definitions

  // Define Madagascar Basalt main components  0054.PP.0044 sample
  density = 3*g/cm3;
  G4Material* madaBasaltMain= new G4Material(name="MadaBasaltMain", density, ncomponents=11);
  madaBasaltMain->AddElement(elSi,fractionmass=0.1992);  // 0.007093 mol/g(mat)
  madaBasaltMain->AddElement(elTi,fractionmass=0.02027); // 4.235e-4
  madaBasaltMain->AddElement(elAl,fractionmass=0.04758); // 0.001763
  madaBasaltMain->AddElement(Fe,  fractionmass=0.1303);  // 0.002333
  madaBasaltMain->AddElement(Mn,  fractionmass=0.001549);// 2.820e-5
  madaBasaltMain->AddElement(Mg,  fractionmass=0.08141); // 0.003350
  madaBasaltMain->AddElement(Ca,  fractionmass=0.06468); // 0.001614
  madaBasaltMain->AddElement(Na,  fractionmass=0.01692); // 7.360e-4
  madaBasaltMain->AddElement(K,   fractionmass=0.008576);// 2.193e-4
  madaBasaltMain->AddElement(P,   fractionmass=0.001977);// 6.383e-5
  madaBasaltMain->AddElement(O,   fractionmass=0.427538);// 0.02672
  // sum is 0.04434383 total number of moles of atoms in one gram of material
  // 2.248766e8 g per 10.000.000 moles.

// Define Madagascar Basalt traces components  0054.PP.0044 sample
  density = 3*g/cm3;
  G4Material* madaBasaltTraces= new G4Material(name="MadaBasaltTraces", density, ncomponents=24);
  
  madaBasaltTraces->AddElement(elTi,natoms=33); 
  madaBasaltTraces->AddElement(Ba  ,natoms=4131);
  madaBasaltTraces->AddElement(Ce  ,natoms=694);
  madaBasaltTraces->AddElement(Co  ,natoms=965);
  madaBasaltTraces->AddElement(Cr  ,natoms=5584);
  madaBasaltTraces->AddElement(La  ,natoms=269);
  madaBasaltTraces->AddElement(Nb  ,natoms=259);
  madaBasaltTraces->AddElement(Nd  ,natoms=410);
  madaBasaltTraces->AddElement(Ni  ,natoms=389);
  madaBasaltTraces->AddElement(Rb  ,natoms=227);
  madaBasaltTraces->AddElement(Sc  ,natoms=212);
  madaBasaltTraces->AddElement(Sr  ,natoms=8686);
  madaBasaltTraces->AddElement(V   ,natoms=4203);
  madaBasaltTraces->AddElement(Y   ,natoms=272);
  madaBasaltTraces->AddElement(Zn  ,natoms=1440);
  madaBasaltTraces->AddElement(Th  ,natoms=19);
  madaBasaltTraces->AddElement(Sm  ,natoms=93);
  madaBasaltTraces->AddElement(Eu  ,natoms=32);
  madaBasaltTraces->AddElement(Gd  ,natoms=89);
  madaBasaltTraces->AddElement(Tb  ,natoms=13);
  madaBasaltTraces->AddElement(Yb  ,natoms=15);
  madaBasaltTraces->AddElement(Lu  ,natoms=2);
  madaBasaltTraces->AddElement(Ta  ,natoms=15);
  madaBasaltTraces->AddElement(Hf  ,natoms=62); //tot 28114/10e7  weight: 2335253.28 g per 10e6 moles

  // Define Madacagascar Basalt complete material  0054.PP.0044 sample
  density = 3*g/cm3;
  G4Material* madaBasalt= new G4Material(name="MadaBasalt", density, ncomponents=2);
  madaBasalt->AddMaterial(madaBasaltMain,    fractionmass=0.9897);
  madaBasalt->AddMaterial(madaBasaltTraces,  fractionmass=0.0103);

  // Define Icelandic Basalt main components  0029.PP.0035 sample
  density = 3*g/cm3;
  G4Material* icelandicBasaltMain= new G4Material(name="IceBasaltMain", density, ncomponents=12);
  icelandicBasaltMain->AddElement(elSi,fractionmass=0.2313);
  icelandicBasaltMain->AddElement(elTi,fractionmass=0.0127);
  icelandicBasaltMain->AddElement(elAl,fractionmass=0.0702);
  icelandicBasaltMain->AddElement(Fe,  fractionmass=0.1134);
  icelandicBasaltMain->AddElement(Mn,  fractionmass=0.0019);
  icelandicBasaltMain->AddElement(Mg,  fractionmass=0.0349);
  icelandicBasaltMain->AddElement(Ca,  fractionmass=0.0756);
  icelandicBasaltMain->AddElement(Na,  fractionmass=0.0892);
  icelandicBasaltMain->AddElement(K,   fractionmass=0.0032);
  icelandicBasaltMain->AddElement(P,   fractionmass=0.00096);
  icelandicBasaltMain->AddElement(elS, fractionmass=0.0004);
  icelandicBasaltMain->AddElement(O,   fractionmass=0.36624);

  // Define Icelandic Basalt traces components  0029.PP.0035 sample
  density = 3*g/cm3;
  G4Material* icelandicBasaltTraces= new G4Material(name="IceBasaltTraces", density, ncomponents=20);
  icelandicBasaltTraces->AddElement(Ba,  natoms=756);
  icelandicBasaltTraces->AddElement(Ce  ,natoms=328);
  icelandicBasaltTraces->AddElement(Co  ,natoms=643);
  icelandicBasaltTraces->AddElement(Cr  ,natoms=1000);
  icelandicBasaltTraces->AddElement(elCu,natoms=1396);
  icelandicBasaltTraces->AddElement(Ga  ,natoms=190);
  icelandicBasaltTraces->AddElement(La  ,natoms=103);
  icelandicBasaltTraces->AddElement(Mo  ,natoms=9);
  icelandicBasaltTraces->AddElement(Nb  ,natoms=114);
  icelandicBasaltTraces->AddElement(Nd  ,natoms=104);
  icelandicBasaltTraces->AddElement(Ni  ,natoms=544);
  icelandicBasaltTraces->AddElement(Rb  ,natoms=78);
  icelandicBasaltTraces->AddElement(elS ,natoms=5550);
  icelandicBasaltTraces->AddElement(Sc  ,natoms=531);
  icelandicBasaltTraces->AddElement(Sr  ,natoms=1353);
  icelandicBasaltTraces->AddElement(elU ,natoms=22);
  icelandicBasaltTraces->AddElement(V   ,natoms=4533);
  icelandicBasaltTraces->AddElement(Y   ,natoms=408);
  icelandicBasaltTraces->AddElement(Zn  ,natoms=1259);
  icelandicBasaltTraces->AddElement(Zr  ,natoms=1274);


  // Define Icelandic Basalt complete material  0029.PP.0035 sample
  density = 3*g/cm3;
  G4Material* icelandicBasalt= new G4Material(name="IceBasalt", density, ncomponents=2);
  icelandicBasalt->AddMaterial(icelandicBasaltMain,    fractionmass=0.9978);
  icelandicBasalt->AddMaterial(icelandicBasaltTraces,  fractionmass=0.0022);


  // Define dolorite main components 0055.PP.0038 sample

  density = 3*g/cm3;
  G4Material* diorite = new G4Material(name="Diorite", density, ncomponents=11);
  diorite->AddElement(Fe,   fractionmass=0.1750);
  diorite->AddElement(elTi, fractionmass=0.0082);
  diorite->AddElement(Ca,   fractionmass=0.0753);
  diorite->AddElement(elSi, fractionmass=0.2188);
  diorite->AddElement(elAl, fractionmass=0.0676);
  diorite->AddElement(Mg,   fractionmass=0.0008);
  diorite->AddElement(O ,   fractionmass=0.4377);
  diorite->AddElement(Mn ,  fractionmass=0.0015);
  diorite->AddElement(Na ,  fractionmass=0.0134);
  diorite->AddElement(K ,   fractionmass=0.0011);
  diorite->AddElement(P ,   fractionmass=0.0006);

  // define traces in dolorite 0055.PP.0038 sample

  density = 3*g/cm3;
  G4Material* tracesOfDolorite = new G4Material(name="TracesOfDolorite", density, ncomponents=16);
  tracesOfDolorite->AddElement(Nb,   natoms=5);
  tracesOfDolorite->AddElement(Zr,   natoms=91);
  tracesOfDolorite->AddElement(Y,    natoms=29);
  tracesOfDolorite->AddElement(Sr,   natoms=140);
  tracesOfDolorite->AddElement(Rb,   natoms=3);
  tracesOfDolorite->AddElement(Ga,   natoms=20);
  tracesOfDolorite->AddElement(Zn,   natoms=99);
  tracesOfDolorite->AddElement(Ni,   natoms=77);
  tracesOfDolorite->AddElement(Sc,   natoms=32);
  tracesOfDolorite->AddElement(V,    natoms=314);
  tracesOfDolorite->AddElement(Cr,   natoms=130);
  tracesOfDolorite->AddElement(Co,   natoms=56);
  tracesOfDolorite->AddElement(elCu,   natoms=119);
  tracesOfDolorite->AddElement(Ba,   natoms=38);
  tracesOfDolorite->AddElement(Ce,   natoms=15);
  tracesOfDolorite->AddElement(Nd,   natoms=9);

  // define dolorite (full) --  0055.PP.0038 sample

  density = 3*g/cm3;
  dolorite = new G4Material(name="Dolorite", density, ncomponents=2);
  dolorite->AddMaterial(tracesOfDolorite, fractionmass=0.0027842352);
  dolorite->AddMaterial(diorite, fractionmass=0.9972157648);

  // define mars1 --  01.PP.0030 sample

  density = 3*g/cm3;
  G4Material* mars1Main = new G4Material(name="Mars1 Main components", density, ncomponents=11);
  mars1Main->AddElement(Fe,   fractionmass=0.100916);
  mars1Main->AddElement(elTi, fractionmass=0.0186804);
  mars1Main->AddElement(Ca,   fractionmass=0.0404091);
  mars1Main->AddElement(elSi, fractionmass=0.196378);
  mars1Main->AddElement(elAl, fractionmass=0.103282);
  mars1Main->AddElement(Mg,   fractionmass=0.0241622);
  mars1Main->AddElement(Mn ,  fractionmass=0.00184331);
  mars1Main->AddElement(Na ,  fractionmass=0.0177908);
  mars1Main->AddElement(K ,   fractionmass=0.00574498);
  mars1Main->AddElement(P ,   fractionmass=0.00280169);
  mars1Main->AddElement(O ,   fractionmass=0.48799152);

  density = 3*g/cm3;
  G4Material* tracesOfMars1 = new G4Material(name="TracesOfMars1", density, ncomponents=17);
  tracesOfMars1->AddElement(Nb,   natoms=55);
  tracesOfMars1->AddElement(Zr,   natoms=433);
  tracesOfMars1->AddElement(Y,    natoms=58);
  tracesOfMars1->AddElement(Sr,   natoms=968);
  tracesOfMars1->AddElement(Rb,   natoms=16);
  tracesOfMars1->AddElement(Ga,   natoms=24);
  tracesOfMars1->AddElement(Zn,   natoms=109);
  tracesOfMars1->AddElement(Ni,   natoms=70);
  tracesOfMars1->AddElement(Sc,   natoms=21);
  tracesOfMars1->AddElement(V,    natoms=134);
  tracesOfMars1->AddElement(Cr,   natoms=141);
  tracesOfMars1->AddElement(Co,   natoms=30);
  tracesOfMars1->AddElement(elCu, natoms=19);
  tracesOfMars1->AddElement(Ba,   natoms=580);
  tracesOfMars1->AddElement(elPb,   natoms=4);
  tracesOfMars1->AddElement(elS,  natoms=444);
  tracesOfMars1->AddElement(elU,    natoms=2);


  density = 3*g/cm3;
            mars1 = new G4Material(name="Mars1", density, ncomponents=2);
  mars1->AddMaterial(tracesOfMars1, fractionmass=0.0044963163);
  mars1->AddMaterial(mars1Main, fractionmass=0.9955036837);

  // define anorthosite

  density = 2.8*g/cm3;
  anorthosite = new G4Material(name="Anorthosite", density, ncomponents=21);
  anorthosite->AddElement(Fe,   fractionmass=0.095283);
  anorthosite->AddElement(Mn,   fractionmass=0.00137086);
  anorthosite->AddElement(Ni,   fractionmass=5e-5);
  anorthosite->AddElement(elCu, fractionmass=5.2e-4);
  anorthosite->AddElement(Na,   fractionmass=0.017635);
  anorthosite->AddElement(Mg,   fractionmass=0.0245361);
  anorthosite->AddElement(elAl, fractionmass=0.0800355);
  anorthosite->AddElement(elSi, fractionmass=0.232204);
  anorthosite->AddElement(Ca,   fractionmass=0.0635368);
  anorthosite->AddElement(K,    fractionmass=0.00464912);
  anorthosite->AddElement(C,    fractionmass=0.000837803);
  anorthosite->AddElement(P,    fractionmass=0.00176742);
  anorthosite->AddElement(elTi, fractionmass=0.0240879);
  anorthosite->AddElement(Cl,   fractionmass=0.00014);
  anorthosite->AddElement(Pd,   fractionmass=0.00001);
  anorthosite->AddElement(Cd,   fractionmass=0.00018);
  anorthosite->AddElement(Ag,   fractionmass=0.00048);
  anorthosite->AddElement(elS,  fractionmass=0.00144);
  anorthosite->AddElement(V,    fractionmass=0.00228);
  anorthosite->AddElement(Ba,   fractionmass=0.00151);
  anorthosite->AddElement(O,    fractionmass=0.447026);


  //define Neodimuim
  
  density =  6800*kg/m3;
              materialNd  = new G4Material(name="Neodimuim"  ,density , ncomponents=1);
  materialNd ->AddElement(Nd,natoms=1);

  // define Berillium
  density = 1848 * kg/m3;
  a = 9.012182 * g / mole;
              Be = new G4Material(name="Beryllium",z=4., a,density);


  // Define Magnesium
  density = 1738 * kg/m3;
              materialMg = new G4Material(name="Magnesium",density , ncomponents=1);
  materialMg->AddElement(Mg,natoms=1);


  //define Tungsten
  density = 19250 * kg/m3;
              materialW = new G4Material(name="Tungsten",density , ncomponents=1);
  materialW->AddElement(W,natoms=1);


  //define iron 

  density = 7.86 * g/cm3;
              FeMaterial = new G4Material(name="Iron",density,ncomponents=1);
  FeMaterial->AddElement(Fe,natoms=1);



  //define gallium arsenide
  
  density = 5.32 * g/cm3;
  G4Material * GaAs = new G4Material(name ="gallium arsenide",density,ncomponents=2);
  GaAs->AddElement(Ga,natoms=1);
  GaAs->AddElement(As,natoms=1);


  // define germanium
  
  density = 5.32 * g/cm3;
              HPGe = new G4Material(name="HPGe",density,ncomponents=1);
  HPGe ->AddElement(Ge,natoms=1);
  
  //define silicon
  
  density = 2.333*g/cm3;
  a = 28.0855*g/mole;
	      Si = new G4Material(name="Silicon",z=14., a,density);
  
  //define copper
  
  density = 8.960*g/cm3;
  a = 63.55*g/mole;
              Cu = new G4Material(name="Copper"   , z=29., a, density);


  
	      ////define Oxigen
	      //density = 1*g/cm3;
	      //a=16*g/mole;
	      //G4Material* matOx = new G4Material(name="Oxigen", z=8., a, density);
  
  //define aluminium
  
  density = 2.700*g/cm3;
  a = 26.98*g/mole;
	      Al = new G4Material(name="Aluminium", z=13., a, density);

  //define titanium 
  density = 4.54 *g/cm3;
  a = 47.867*g/mole;
              Ti  = new G4Material(name="Titanium",z=22.,a,density);


 //define Uranium 
 density = 19050*kg/m3;
 a =  238.02891*g/mole;
              U  = new G4Material(name="Uranium",z=92.,a,density);

  //define Tin
  density = 7310*kg/m3;
  a =  118.710*g/mole;
              Sn  = new G4Material(name="Tin",z=50.,a,density);

  //define lead
  
  density = 11.35*g/cm3;
  a=207.19*g/mole;
              Pb = new G4Material(name="Lead",z=82.,a,density);

  // define Silver material
  density =   10490*kg/m3;
              materialAg  = new G4Material(name="Silver"  ,density , ncomponents=1);
  materialAg ->AddElement(Ag,natoms=1);

  //define gold
  
  density = 19300*kg/m3;
  a= 196.96655*g/mole;
              Au = new G4Material(name="Gold",z=79.,a,density);

  //define caesium
  
  density = 1879*kg/m3;
  a= 132.90545*g/mole;
              Cs = new G4Material(name="Caesium",z=55.,a,density);

  // define Potassium material
  density =   1879*kg/m3;
              materialK  = new G4Material(name="Potassium"  ,density , ncomponents=1);
  materialK ->AddElement(K,natoms=1);

  // define Manganese material
  density =   7470*kg/m3;
              materialMn  = new G4Material(name="Manganese"  ,density , ncomponents=1);
  materialMn ->AddElement(Mn,natoms=1);

  // define Phosphorus material
  density =   1823*kg/m3;
              materialP  = new G4Material(name="Phosphorus"  ,density , ncomponents=1);
  materialP ->AddElement(P,natoms=1);

  // define Sufur material
  density =   1960*kg/m3;
              materialS  = new G4Material(name="Sulphur"  ,density , ncomponents=1);
  materialS ->AddElement(elS,natoms=1);

  // define Calcium material
  density =   1550*kg/m3;
              materialCa  = new G4Material(name="Calcium"  ,density , ncomponents=1);
  materialCa ->AddElement(Ca,natoms=1);

  // define Sodium material
  density =   968*kg/m3;
              materialNa  = new G4Material(name="Sodium"  ,density , ncomponents=1);
  materialNa ->AddElement(Na,natoms=1);







  //define scintillator
  
  density = 1.032*g/cm3;
  Sci = new G4Material(name="Scintillator", density, ncomponents=2);
  Sci->AddElement(C, natoms=9);
  Sci->AddElement(H, natoms=10);  

  
  //define air

  density = 1.290*mg/cm3;
  Air = new G4Material(name="Air"  , density, ncomponents=2);
  Air->AddElement(N, fractionmass=0.7);
  Air->AddElement(O, fractionmass=0.3);
  
  //define vacuum
  
  density     = universe_mean_density;    //from PhysicalConstants.h
  pressure    = 3.e-18*pascal;
  temperature = 2.73*kelvin;
  Vacuum = new G4Material(name="Galactic", z=1., a=1.01*g/mole, density,
				       kStateGas,temperature,pressure);

  //define basalt
  density = 3.*g/cm3; 
  basalt = new G4Material(name="Basalt", density, ncomponents=7);
  basalt->AddElement(Fe, fractionmass=0.1200);
  basalt->AddElement(elTi, fractionmass=0.0160);
  basalt->AddElement(Ca, fractionmass=0.0750);
  basalt->AddElement(elSi, fractionmass=0.2160);
  basalt->AddElement(elAl, fractionmass=0.0710);
  basalt->AddElement(Mg, fractionmass=0.0590);
  basalt->AddElement(O , fractionmass=0.4430);


  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

