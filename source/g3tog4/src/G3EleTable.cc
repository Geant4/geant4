// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3EleTable.cc,v 1.8 1999-12-15 14:49:43 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4strstreambuf.hh"
#include "G4ios.hh"
#include "G3EleTable.hh"

G3EleTable::G3EleTable() :_MaxEle(109){
  _EleNames = new char*[_MaxEle];
  // create an array of pointers to elements
  _Ele = new G4Element*[_MaxEle];
  LoadUp();
}

G3EleTable::~G3EleTable(){
  delete [] _EleNames;
  delete [] _Ele;
};

G4Element* 
G3EleTable::GetEle(G4double Z){
  G4double A;
  char name[20], sym[3];
  G4int index = (G4int) Z-1;
  if (!parse(Z, name, sym, A)) {
    G4String nm(name);
    G4String sy(sym);
    if (_Ele[index] == 0) {
      // add an element to the element table here
      _Ele[index] = new G4Element(nm, sy, Z, A*g/mole);
    }
  }
  return _Ele[index];
}

int 
G3EleTable::parse(G4double& Z, char* name, char* sym, G4double& A){ 
 int rc = 0;
  if (Z>0 && Z <=_MaxEle){
    G4int z = (G4int) Z-1;
    G4std::istrstream in(_EleNames[z]);
    in >> name >> sym >> A;
  } else {
    rc = -1;
  }
  return rc;
};

void
G3EleTable::LoadUp(){
  int i=0;
  _EleNames[i]=(char *)"Hydrogen H 1.00794"; i++;
  _EleNames[i]=(char *)"Helium He 4.0026"; i++;
  _EleNames[i]=(char *)"Lithium Li 6.941"; i++;
  _EleNames[i]=(char *)"Beryllium Be 9.012182"; i++;
  _EleNames[i]=(char *)"Boron B 10.811"; i++;
  _EleNames[i]=(char *)"Carbon C 12.011"; i++;
  _EleNames[i]=(char *)"Nitrogen N 14.00674"; i++;
  _EleNames[i]=(char *)"Oxygen O 15.9994"; i++;
  _EleNames[i]=(char *)"Fluorine F 18.9984032"; i++;
  _EleNames[i]=(char *)"Neon Ne 20.1797"; i++;

  _EleNames[i]=(char *)"Sodium Na 22.989768"; i++;
  _EleNames[i]=(char *)"Magnesium Mg 24.3050"; i++;
  _EleNames[i]=(char *)"Aluminum Al 26.981539"; i++;
  _EleNames[i]=(char *)"Silicon Si 28.0855"; i++;
  _EleNames[i]=(char *)"Phosphorus P 30.973762"; i++;
  _EleNames[i]=(char *)"Sulfur S 32.066"; i++;
  _EleNames[i]=(char *)"Chlorine Cl 35.4527"; i++;
  _EleNames[i]=(char *)"Argon Ar 39.948"; i++;
  _EleNames[i]=(char *)"Potassium K 39.0983"; i++;
  _EleNames[i]=(char *)"Calcium Ca 40.078"; i++;

  _EleNames[i]=(char *)"Scandium Sc 44.955910"; i++;
  _EleNames[i]=(char *)"Titanium Ti 47.867"; i++;
  _EleNames[i]=(char *)"Vanadium V 50.9415"; i++;
  _EleNames[i]=(char *)"Chromium Cr 51.9961"; i++;
  _EleNames[i]=(char *)"Manganese Mn 54.93805"; i++;
  _EleNames[i]=(char *)"Iron Fe 55.845"; i++;
  _EleNames[i]=(char *)"Cobalt Co 58.93320"; i++;
  _EleNames[i]=(char *)"Nickel Ni 58.6934"; i++;
  _EleNames[i]=(char *)"Copper Cu 63.546"; i++;
  _EleNames[i]=(char *)"Zinc Zn 65.39"; i++;

  _EleNames[i]=(char *)"Gallium Ga 69.723"; i++;
  _EleNames[i]=(char *)"Germanium Ge 72.61"; i++;
  _EleNames[i]=(char *)"Arsenic As 74.92159"; i++;
  _EleNames[i]=(char *)"Selenium Se 78.96"; i++;
  _EleNames[i]=(char *)"Bromine Br 79.904"; i++;
  _EleNames[i]=(char *)"Krypton Kr 83.80"; i++;
  _EleNames[i]=(char *)"Rubidium Rb 85.4678"; i++;
  _EleNames[i]=(char *)"Strontium Sr 87.62"; i++;
  _EleNames[i]=(char *)"Yttrium Y 88.90585"; i++;
  _EleNames[i]=(char *)"Zirconium Zr 91.224"; i++;

  _EleNames[i]=(char *)"Niobium Nb 92.90638"; i++;
  _EleNames[i]=(char *)"Molybdenum Mo 95.94"; i++;
  _EleNames[i]=(char *)"Technetium Tc 97.907215"; i++;
  _EleNames[i]=(char *)"Ruthenium Ru 101.07"; i++;
  _EleNames[i]=(char *)"Rhodium Rh 102.90550"; i++;
  _EleNames[i]=(char *)"Palladium Pd 106.42"; i++;
  _EleNames[i]=(char *)"Silver Ag 107.8682"; i++;
  _EleNames[i]=(char *)"Cadmium Cd 112.41"; i++;
  _EleNames[i]=(char *)"Indium In 114.818"; i++;
  _EleNames[i]=(char *)"Tin Sn 118.710"; i++;

  _EleNames[i]=(char *)"Antimony Sb 121.760"; i++;
  _EleNames[i]=(char *)"Tellurium Te 127.60"; i++;
  _EleNames[i]=(char *)"Iodine I 126.90447"; i++;
  _EleNames[i]=(char *)"Xenon Xe 131.29"; i++;
  _EleNames[i]=(char *)"Cesium Cs 132.90543"; i++;
  _EleNames[i]=(char *)"Barium Ba 137.27"; i++;
  _EleNames[i]=(char *)"Lanthanum La 138.9055"; i++;
  _EleNames[i]=(char *)"Cerium Ce 140.115"; i++;
  _EleNames[i]=(char *)"Praeseodymium Pr 140.90765"; i++;
  _EleNames[i]=(char *)"NeoDymium Nd 144.24"; i++;
  
  _EleNames[i]=(char *)"Promethium Pm 144.912745"; i++;
  _EleNames[i]=(char *)"Samarium Sm 150.36"; i++;
  _EleNames[i]=(char *)"Europium Eu 151.965"; i++;
  _EleNames[i]=(char *)"Gadolinium Gd 157.25"; i++;
  _EleNames[i]=(char *)"Terbium Tb 158.92534"; i++;
  _EleNames[i]=(char *)"Dysprosium Dy 162.50"; i++;
  _EleNames[i]=(char *)"Holmium Ho 164.93032"; i++;
  _EleNames[i]=(char *)"Erbium Er 167.26"; i++;
  _EleNames[i]=(char *)"Thulium Tm 168.93421"; i++;
  _EleNames[i]=(char *)"Ytterbium Yb 173.04"; i++;

  _EleNames[i]=(char *)"Lutetium Lu 174.967"; i++;
  _EleNames[i]=(char *)"Hafnium Hf 178.49"; i++;
  _EleNames[i]=(char *)"Tantalum Ta 180.9479"; i++;
  _EleNames[i]=(char *)"Tungsten W 183.84"; i++;
  _EleNames[i]=(char *)"Rhenium Re 186.207"; i++;
  _EleNames[i]=(char *)"Osmium Os 190.23"; i++;
  _EleNames[i]=(char *)"Iridium Ir 192.217"; i++;
  _EleNames[i]=(char *)"Platinum Pt 195.08"; i++;
  _EleNames[i]=(char *)"Gold Au 196.96654"; i++;
  _EleNames[i]=(char *)"Mercury Hg 200.59"; i++;

  _EleNames[i]=(char *)"Thallium Tl 204.3833"; i++;
  _EleNames[i]=(char *)"Lead Pb 207.2"; i++;
  _EleNames[i]=(char *)"Bismuth Bi 208.98037"; i++;
  _EleNames[i]=(char *)"Polonium Po 208.982415"; i++;
  _EleNames[i]=(char *)"Astatine At 209.987131"; i++;
  _EleNames[i]=(char *)"Radon Rn 222.017570"; i++;
  _EleNames[i]=(char *)"Francium Fr 223.019731"; i++;
  _EleNames[i]=(char *)"Radium Ra 226.025402"; i++;
  _EleNames[i]=(char *)"Actinium Ac 227.027747"; i++;
  _EleNames[i]=(char *)"Thorium Th 232.0381"; i++;

  _EleNames[i]=(char *)"Protactinium Pa 231.03588"; i++;
  _EleNames[i]=(char *)"Uranium U 238.0289"; i++;
  _EleNames[i]=(char *)"Neptunium Np 237.048166"; i++;
  _EleNames[i]=(char *)"Plutonium Pu 244.064197"; i++;
  _EleNames[i]=(char *)"Americium Am 243.061372"; i++;
  _EleNames[i]=(char *)"Curium Cm 247.070346"; i++;
  _EleNames[i]=(char *)"Berkelium Bk 247.070298"; i++;
  _EleNames[i]=(char *)"Californium Cf 251.079579"; i++;
  _EleNames[i]=(char *)"Einsteinium Es 252.08297"; i++;
  _EleNames[i]=(char *)"Fermium Fm 257.095096"; i++;

  _EleNames[i]=(char *)"Mendelevium Md 258.098427"; i++;
  _EleNames[i]=(char *)"Nobelium No 259.1011"; i++;
  _EleNames[i]=(char *)"Lawrencium Lr 262.1098"; i++;
  _EleNames[i]=(char *)"Rutherfordium Rf 261.1089";  i++;
  _EleNames[i]=(char *)"Hahnium Ha 262.1144";  i++;
  _EleNames[i]=(char *)"Seaborgium Sg 263.1186";  i++;
  _EleNames[i]=(char *)"Nielsborium Ns 262.1231";  i++;
  _EleNames[i]=(char *)"Hassium Hs 265.1306";  i++;
  _EleNames[i]=(char *)"Meitnerium Mt 266.1378";  i++;

  // initialize element pointers to 0
  for (int j=0; j<i; j++) {
    _Ele[j]=0;
  }
}


