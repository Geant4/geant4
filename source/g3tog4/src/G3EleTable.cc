// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3EleTable.cc,v 1.1 1999-05-06 04:22:22 lockman Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4strstreambuf.hh"
#include "G4ios.hh"
#include "G3EleTable.hh"

G3EleTable::G3EleTable() :_MaxEle(109){
  _EleNames = new char*[_MaxEle];
  _Ele = new G4Element*[_MaxEle];
  LoadUp();
}

G3EleTable::~G3EleTable(){
  G4cout << "Destructing G3EleTable." << endl;
  delete [] _EleNames;
  delete [] _Ele;
};

G4Element* 
G3EleTable::GetEle(G4double Z){
  G4double A;
  G4Element* E=0;
  char name[20], sym[3];
  if (!parse(Z, name, sym, A)) {
    G4String nm(name);
    G4String sy(sym);
    G4int z = Z-1;
    if (_Ele[z] == 0) {
      _Ele[z] = new G4Element(nm, sy, Z, A*g/mole);
    }
    E = _Ele[z];
  }
  return E;
}

int 
G3EleTable::parse(G4double& Z, char* name, char* sym, G4double& A){ 
 int rc = 0;
  if (Z>0 && Z <=_MaxEle){
    G4int z = Z-1;
    istrstream in(_EleNames[z]);
    in >> name >> sym >> A;
  } else {
    rc = -1;
  }
  return rc;
};

void
G3EleTable::LoadUp(){
  int i=0;
  _EleNames[i]="Hydrogen H 1.00794"; i++;
  _EleNames[i]="Helium He 4.0026"; i++;
  _EleNames[i]="Lithium Li 6.941"; i++;
  _EleNames[i]="Beryllium Be 9.012182"; i++;
  _EleNames[i]="Boron B 10.811"; i++;
  _EleNames[i]="Carbon C 12.011"; i++;
  _EleNames[i]="Nitrogen N 14.00674"; i++;
  _EleNames[i]="Oxygen O 15.9994"; i++;
  _EleNames[i]="Fluorine F 18.9984032"; i++;
  _EleNames[i]="Neon Ne 20.1797"; i++;

  _EleNames[i]="Sodium Na 22.989768"; i++;
  _EleNames[i]="Magnesium Mg 24.3050"; i++;
  _EleNames[i]="Aluminum Al 26.981539"; i++;
  _EleNames[i]="Silicon Si 28.0855"; i++;
  _EleNames[i]="Phosphorus P 30.973762"; i++;
  _EleNames[i]="Sulfur S 32.066"; i++;
  _EleNames[i]="Chlorine Cl 35.4527"; i++;
  _EleNames[i]="Argon Ar 39.948"; i++;
  _EleNames[i]="Potassium K 39.0983"; i++;
  _EleNames[i]="Calcium Ca 40.078"; i++;

  _EleNames[i]="Scandium Sc 44.955910"; i++;
  _EleNames[i]="Titanium Ti 47.867"; i++;
  _EleNames[i]="Vanadium V 50.9415"; i++;
  _EleNames[i]="Chromium Cr 51.9961"; i++;
  _EleNames[i]="Manganese Mn 54.93805"; i++;
  _EleNames[i]="Iron Fe 55.845"; i++;
  _EleNames[i]="Cobalt Co 58.93320"; i++;
  _EleNames[i]="Nickel Ni 58.6934"; i++;
  _EleNames[i]="Copper Cu 63.546"; i++;
  _EleNames[i]="Zinc Zn 65.39"; i++;

  _EleNames[i]="Gallium Ga 69.723"; i++;
  _EleNames[i]="Germanium Ge 72.61"; i++;
  _EleNames[i]="Arsenic As 74.92159"; i++;
  _EleNames[i]="Selenium Se 78.96"; i++;
  _EleNames[i]="Bromine Br 79.904"; i++;
  _EleNames[i]="Krypton Kr 83.80"; i++;
  _EleNames[i]="Rubidium Rb 85.4678"; i++;
  _EleNames[i]="Strontium Sr 87.62"; i++;
  _EleNames[i]="Yttrium Y 88.90585"; i++;
  _EleNames[i]="Zirconium Zr 91.224"; i++;

  _EleNames[i]="Niobium Nb 92.90638"; i++;
  _EleNames[i]="Molybdenum Mo 95.94"; i++;
  _EleNames[i]="Technetium Tc 97.907215"; i++;
  _EleNames[i]="Ruthenium Ru 101.07"; i++;
  _EleNames[i]="Rhodium Rh 102.90550"; i++;
  _EleNames[i]="Palladium Pd 106.42"; i++;
  _EleNames[i]="Silver Ag 107.8682"; i++;
  _EleNames[i]="Cadmium Cd 112.41"; i++;
  _EleNames[i]="Indium In 114.818"; i++;
  _EleNames[i]="Tin Sn 118.710"; i++;

  _EleNames[i]="Antimony Sb 121.760"; i++;
  _EleNames[i]="Tellurium Te 127.60"; i++;
  _EleNames[i]="Iodine I 126.90447"; i++;
  _EleNames[i]="Xenon Xe 131.29"; i++;
  _EleNames[i]="Cesium Cs 132.90543"; i++;
  _EleNames[i]="Barium Ba 137.27"; i++;
  _EleNames[i]="Lanthanum La 138.9055"; i++;
  _EleNames[i]="Cerium Ce 140.115"; i++;
  _EleNames[i]="Praeseodymium Pr 140.90765"; i++;
  _EleNames[i]="NeoDymium Nd 144.24"; i++;
  
  _EleNames[i]="Promethium Pm 144.912745"; i++;
  _EleNames[i]="Samarium Sm 150.36"; i++;
  _EleNames[i]="Europium Eu 151.965"; i++;
  _EleNames[i]="Gadolinium Gd 157.25"; i++;
  _EleNames[i]="Terbium Tb 158.92534"; i++;
  _EleNames[i]="Dysprosium Dy 162.50"; i++;
  _EleNames[i]="Holmium Ho 164.93032"; i++;
  _EleNames[i]="Erbium Er 167.26"; i++;
  _EleNames[i]="Thulium Tm 168.93421"; i++;
  _EleNames[i]="Ytterbium Yb 173.04"; i++;

  _EleNames[i]="Lutetium Lu 174.967"; i++;
  _EleNames[i]="Hafnium Hf 178.49"; i++;
  _EleNames[i]="Tantalum Ta 180.9479"; i++;
  _EleNames[i]="Tungsten W 183.84"; i++;
  _EleNames[i]="Rhenium Re 186.207"; i++;
  _EleNames[i]="Osmium Os 190.23"; i++;
  _EleNames[i]="Iridium Ir 192.217"; i++;
  _EleNames[i]="Platinum Pt 195.08"; i++;
  _EleNames[i]="Gold Au 196.96654"; i++;
  _EleNames[i]="Mercury Hg 200.59"; i++;

  _EleNames[i]="Thallium Tl 204.3833"; i++;
  _EleNames[i]="Lead Pb 207.2"; i++;
  _EleNames[i]="Bismuth Bi 208.98037"; i++;
  _EleNames[i]="Polonium Po 208.982415"; i++;
  _EleNames[i]="Astatine At 209.987131"; i++;
  _EleNames[i]="Radon Rn 222.017570"; i++;
  _EleNames[i]="Francium Fr 223.019731"; i++;
  _EleNames[i]="Radium Ra 226.025402"; i++;
  _EleNames[i]="Actinium Ac 227.027747"; i++;
  _EleNames[i]="Thorium Th 232.0381"; i++;

  _EleNames[i]="Protactinium Pa 231.03588"; i++;
  _EleNames[i]="Uranium U 238.0289"; i++;
  _EleNames[i]="Neptunium Np 237.048166"; i++;
  _EleNames[i]="Plutonium Pu 244.064197"; i++;
  _EleNames[i]="Americium Am 243.061372"; i++;
  _EleNames[i]="Curium Cm 247.070346"; i++;
  _EleNames[i]="Berkelium Bk 247.070298"; i++;
  _EleNames[i]="Californium Cf 251.079579"; i++;
  _EleNames[i]="Einsteinium Es 252.08297"; i++;
  _EleNames[i]="Fermium Fm 257.095096"; i++;

  _EleNames[i]="Mendelevium Md 258.098427"; i++;
  _EleNames[i]="Nobelium No 259.1011"; i++;
  _EleNames[i]="Lawrencium Lr 262.1098"; i++;
  _EleNames[i]="Rutherfordium Rf 261.1089";  i++;
  _EleNames[i]="Hahnium Ha 262.1144";  i++;
  _EleNames[i]="Seaborgium Sg 263.1186";  i++;
  _EleNames[i]="Nielsborium Ns 262.1231";  i++;
  _EleNames[i]="Hassium Hs 265.1306";  i++;
  _EleNames[i]="Meitnerium Mt 266.1378";  i++;

  // initialize element pointers to 0
  for (int j=0; j<i; j++) {
    _Ele[j]=0;
  }
}


