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
// $Id: G4UnitsTable.cc 99523 2016-09-26 10:25:48Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
// 17-05-98: first version, M.Maire
// 05-08-98: angstrom,microbarn,picobarn,petaelectronvolt, M.Maire
// 13-10-98: units and symbols printed in fixed length, M.Maire
// 01-03-01: parsec, M.Maire
// 06-03-01: migration to STL vectors, G.Cosmo
// 06-05-02: BestUnit operator<<  flux instead of G4cout (mma)
// 12-08-05: cm2/g ("Surface/Mass")  (mma)
// 30-06-05: um for micrometer (mma)
// 07-02-06: GeV/cm MeV/cm keV/cm eV/cm ("Energy/Length")  (mma)
// 15-02-06: g/cm2 ("Mass/Surface")
//           MeV*cm2/g ..etc.. ("Energy*Surface/Mass")
// 18-08-06: remove symbol mum (mma)
// 06-05-08: V/m ("Electric field")  (mma)
// 09-08-10: new category "Solid angle"  (mma)
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include <iomanip>
#include <sstream>

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4Threading.hh"

G4ThreadLocal G4UnitsTable* G4UnitDefinition::pUnitsTable = 0;
G4ThreadLocal G4bool G4UnitDefinition::unitsTableDestroyed = false;

#ifdef G4MULTITHREADED
G4UnitsTable* G4UnitDefinition::pUnitsTableShadow = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4UnitsTable::G4UnitsTable() {;}
G4UnitsTable::~G4UnitsTable() 
{
  G4UnitsTable::iterator itr = begin();
  for(;itr!=end();itr++)
  { delete *itr; }
  clear();
}

#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4UnitDefinition::G4UnitDefinition(const G4String& name,
                                   const G4String& symbol,
                                   const G4String& category, G4double value)
  : Name(name),SymbolName(symbol),Value(value)
{
    if (!pUnitsTable)
    {
      if(unitsTableDestroyed)
      {
        G4Exception("G4UnitDefinition::G4UnitDefinition","UnitsTable0000",
                    FatalException,"G4UnitsTable had already deleted.");
      }
      pUnitsTable = new G4UnitsTable;
#ifdef G4MULTITHREADED
      if(G4Threading::IsMasterThread())
      { pUnitsTableShadow = pUnitsTable; }
#endif
    }

    // Does the Category objet already exist ?
    //
    size_t nbCat = pUnitsTable->size();
    size_t i = 0;
    while ((i<nbCat)&&((*pUnitsTable)[i]->GetName()!=category))  { i++; }
    if (i == nbCat)
      { pUnitsTable->push_back( new G4UnitsCategory(category)); }
    CategoryIndex = i;

    // Insert this Unit in the Units table
    //
    ((*pUnitsTable)[CategoryIndex]->GetUnitsList()).push_back(this);
    
    // Update string max length for name and symbol
    //
    (*pUnitsTable)[i]->UpdateNameMxLen((G4int)name.length());
    (*pUnitsTable)[i]->UpdateSymbMxLen((G4int)symbol.length());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4UnitDefinition::~G4UnitDefinition()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4UnitDefinition::G4UnitDefinition(const G4UnitDefinition& right)
{
    *this = right;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4UnitDefinition& G4UnitDefinition::operator=(const G4UnitDefinition& right)
{
  if (this != &right)
    {
      Name          = right.Name;
      SymbolName    = right.SymbolName;
      Value         = right.Value;
      CategoryIndex = right.CategoryIndex;
    }
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4int G4UnitDefinition::operator==(const G4UnitDefinition& right) const
{
  return (this == (G4UnitDefinition *) &right);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4int G4UnitDefinition::operator!=(const G4UnitDefinition &right) const
{
  return (this != (G4UnitDefinition *) &right);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4UnitsTable& G4UnitDefinition::GetUnitsTable()
{
  if (!pUnitsTable)  { pUnitsTable = new G4UnitsTable; }
  if(pUnitsTable->size()==0)  { BuildUnitsTable(); }
#ifdef G4MULTITHREADED
  if(G4Threading::IsMasterThread() && !pUnitsTableShadow)
  { pUnitsTableShadow = pUnitsTable; }
#endif
  return *pUnitsTable;
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4bool G4UnitDefinition::IsUnitDefined(const G4String& str)
{
  G4String name, symbol;
  for (size_t i=0;i<(GetUnitsTable()).size();++i)
  {
    G4UnitsContainer& units = (*pUnitsTable)[i]->GetUnitsList();
    for (size_t j=0;j<units.size();++j)
    {
      name=units[j]->GetName(); symbol=units[j]->GetSymbol();
      if(str==name||str==symbol) { return true; }
    }
  }
  return false;             
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4double G4UnitDefinition::GetValueOf(const G4String& str)
{
  G4String name, symbol;
  for (size_t i=0;i<(GetUnitsTable()).size();++i)
     {
       G4UnitsContainer& units = (*pUnitsTable)[i]->GetUnitsList();
       for (size_t j=0;j<units.size();++j)
          {
            name=units[j]->GetName(); symbol=units[j]->GetSymbol();
            if(str==name||str==symbol) { return units[j]->GetValue(); }
          }
     }
  std::ostringstream message;
  message << "The unit '" << str << "' does not exist in the Units Table!";
  G4Exception("G4UnitDefinition::GetValueOf()", "InvalidUnit",
              FatalException, message);
  return 0.;             
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
G4String G4UnitDefinition::GetCategory(const G4String& str)
{
  G4String name, symbol;
  for (size_t i=0;i<(GetUnitsTable()).size();++i)
     {
       G4UnitsContainer& units = (*pUnitsTable)[i]->GetUnitsList();
       for (size_t j=0;j<units.size();++j)
          {
            name=units[j]->GetName(); symbol=units[j]->GetSymbol();
            if(str==name||str==symbol) { return (*pUnitsTable)[i]->GetName(); }
          }
     }
  std::ostringstream message;
  message << "The unit '" << str << "' does not exist in the Units Table!";
  G4Exception("G4UnitDefinition::GetCategory()", "InvalidUnit",
              FatalException, message);
  name = "None";     
  return name;             
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void G4UnitDefinition::PrintDefinition()
{
  G4int nameL = (*pUnitsTable)[CategoryIndex]->GetNameMxLen();
  G4int symbL = (*pUnitsTable)[CategoryIndex]->GetSymbMxLen();
  G4cout << std::setw(nameL) << Name << " (" 
         << std::setw(symbL) << SymbolName << ") = " << Value << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void G4UnitDefinition::BuildUnitsTable()
{
 //Length
 new G4UnitDefinition(    "parsec","pc"      ,"Length",parsec); 
 new G4UnitDefinition( "kilometer","km"      ,"Length",kilometer);
 new G4UnitDefinition(     "meter","m"       ,"Length",meter);
 new G4UnitDefinition("centimeter","cm"      ,"Length",centimeter); 
 new G4UnitDefinition("millimeter","mm"      ,"Length",millimeter);
 new G4UnitDefinition("micrometer","um"      ,"Length",micrometer);
 new G4UnitDefinition( "nanometer","nm"      ,"Length",nanometer);
 new G4UnitDefinition(  "angstrom","Ang"     ,"Length",angstrom);    
 new G4UnitDefinition(     "fermi","fm"      ,"Length",fermi);
 
 //Surface
 new G4UnitDefinition( "kilometer2","km2"    ,"Surface",kilometer2);
 new G4UnitDefinition(     "meter2","m2"     ,"Surface",meter2);
 new G4UnitDefinition("centimeter2","cm2"    ,"Surface",centimeter2); 
 new G4UnitDefinition("millimeter2","mm2"    ,"Surface",millimeter2);
 new G4UnitDefinition(       "barn","barn"   ,"Surface",barn);
 new G4UnitDefinition(  "millibarn","mbarn"  ,"Surface",millibarn);   
 new G4UnitDefinition(  "microbarn","mubarn" ,"Surface",microbarn);
 new G4UnitDefinition(   "nanobarn","nbarn"  ,"Surface",nanobarn);
 new G4UnitDefinition(   "picobarn","pbarn"  ,"Surface",picobarn);
 
 //Volume
 new G4UnitDefinition( "kilometer3","km3"    ,"Volume",kilometer3);
 new G4UnitDefinition(     "meter3","m3"     ,"Volume",meter3);
 new G4UnitDefinition("centimeter3","cm3"    ,"Volume",centimeter3); 
 new G4UnitDefinition("millimeter3","mm3"    ,"Volume",millimeter3);
 
 new G4UnitDefinition(      "liter","L"      ,"Volume",liter);
 new G4UnitDefinition(        "dL","dL"      ,"Volume",dL);
 new G4UnitDefinition(        "cL","cL"      ,"Volume",cL);
 new G4UnitDefinition(        "mL","mL"      ,"Volume",mL);    

 //Angle
 new G4UnitDefinition(     "radian","rad"    ,"Angle",radian);
 new G4UnitDefinition("milliradian","mrad"   ,"Angle",milliradian); 
 new G4UnitDefinition(     "degree","deg"    ,"Angle",degree);
 
 //Solid angle
 new G4UnitDefinition(  "steradian","sr"     ,"Solid angle",steradian);
 new G4UnitDefinition("millisteradian","msr" ,"Solid angle",steradian*0.001);
   
 //Time
 new G4UnitDefinition(     "second","s"      ,"Time",second);
 new G4UnitDefinition("millisecond","ms"     ,"Time",millisecond);
 new G4UnitDefinition("microsecond","us"     ,"Time",microsecond);
 new G4UnitDefinition( "nanosecond","ns"     ,"Time",nanosecond);
 new G4UnitDefinition( "picosecond","ps"     ,"Time",picosecond);
 
 //Frequency
 new G4UnitDefinition(    "hertz","Hz"       ,"Frequency",hertz);
 new G4UnitDefinition("kilohertz","kHz"      ,"Frequency",kilohertz);
 new G4UnitDefinition("megahertz","MHz"      ,"Frequency",megahertz);
 
 //Electric charge
 new G4UnitDefinition(  "eplus","e+"         ,"Electric charge",eplus);
 new G4UnitDefinition("coulomb","C"          ,"Electric charge",coulomb); 
 
 //Energy
 new G4UnitDefinition(    "electronvolt","eV" ,"Energy",electronvolt);
 new G4UnitDefinition("kiloelectronvolt","keV","Energy",kiloelectronvolt);
 new G4UnitDefinition("megaelectronvolt","MeV","Energy",megaelectronvolt);
 new G4UnitDefinition("gigaelectronvolt","GeV","Energy",gigaelectronvolt);
 new G4UnitDefinition("teraelectronvolt","TeV","Energy",teraelectronvolt);
 new G4UnitDefinition("petaelectronvolt","PeV","Energy",petaelectronvolt);
 new G4UnitDefinition(           "joule","J"  ,"Energy",joule);
 
 // Energy/Length
 new G4UnitDefinition( "GeV/cm", "GeV/cm","Energy/Length", GeV/cm);
 new G4UnitDefinition( "MeV/cm", "MeV/cm","Energy/Length", MeV/cm);
 new G4UnitDefinition( "keV/cm", "keV/cm","Energy/Length", keV/cm);
 new G4UnitDefinition(  "eV/cm",  "eV/cm","Energy/Length",  eV/cm); 
  
 //Mass
 new G4UnitDefinition("milligram","mg","Mass",milligram);
 new G4UnitDefinition(     "gram","g" ,"Mass",gram);
 new G4UnitDefinition( "kilogram","kg","Mass",kilogram);
 
 //Volumic Mass
 new G4UnitDefinition( "g/cm3", "g/cm3","Volumic Mass", g/cm3);
 new G4UnitDefinition("mg/cm3","mg/cm3","Volumic Mass",mg/cm3);
 new G4UnitDefinition("kg/m3", "kg/m3", "Volumic Mass",kg/m3);
 
 // Mass/Surface
 new G4UnitDefinition(  "g/cm2",  "g/cm2","Mass/Surface",  g/cm2);
 new G4UnitDefinition( "mg/cm2", "mg/cm2","Mass/Surface", mg/cm2);
 new G4UnitDefinition( "kg/cm2", "kg/cm2","Mass/Surface", kg/cm2);
   
 // Surface/Mass
 new G4UnitDefinition( "cm2/g", "cm2/g","Surface/Mass", cm2/g);
 
 // Energy.Surface/Mass
 new G4UnitDefinition( "eV*cm2/g", " eV*cm2/g","Energy*Surface/Mass", eV*cm2/g);
 new G4UnitDefinition("keV*cm2/g", "keV*cm2/g","Energy*Surface/Mass",keV*cm2/g);
 new G4UnitDefinition("MeV*cm2/g", "MeV*cm2/g","Energy*Surface/Mass",MeV*cm2/g);
 new G4UnitDefinition("GeV*cm2/g", "GeV*cm2/g","Energy*Surface/Mass",GeV*cm2/g);
     
 //Power
 new G4UnitDefinition("watt","W","Power",watt);
 
 //Force
 new G4UnitDefinition("newton","N","Force",newton);
 
 //Pressure
 new G4UnitDefinition(    "pascal","Pa" ,"Pressure",hep_pascal);
 new G4UnitDefinition(       "bar","bar","Pressure",bar); 
 new G4UnitDefinition("atmosphere","atm","Pressure",atmosphere);
 
 //Electric current
 new G4UnitDefinition(     "ampere","A"  ,"Electric current",ampere);
 new G4UnitDefinition("milliampere","mA" ,"Electric current",milliampere);
 new G4UnitDefinition("microampere","muA","Electric current",microampere);
 new G4UnitDefinition( "nanoampere","nA" ,"Electric current",nanoampere);   
 
 //Electric potential
 new G4UnitDefinition(    "volt","V" ,"Electric potential",volt); 
 new G4UnitDefinition("kilovolt","kV","Electric potential",kilovolt);
 new G4UnitDefinition("megavolt","MV","Electric potential",megavolt);
 
 //Electric field
 new G4UnitDefinition(    "volt/m","V/m" ,"Electric field",volt/m);
 new G4UnitDefinition("kilovolt/m","kV/m","Electric field",kilovolt/m);
 new G4UnitDefinition("megavolt/m","MV/m","Electric field",megavolt/m);

 //Magnetic flux
 new G4UnitDefinition("weber","Wb","Magnetic flux",weber);
 
 //Magnetic flux density
 new G4UnitDefinition(    "tesla","T" ,"Magnetic flux density",tesla);
 new G4UnitDefinition("kilogauss","kG","Magnetic flux density",kilogauss);
 new G4UnitDefinition(    "gauss","G" ,"Magnetic flux density",gauss);
 
 //Temperature
 new G4UnitDefinition("kelvin","K","Temperature",kelvin);
 
 //Amount of substance
 new G4UnitDefinition("mole","mol","Amount of substance",mole);
 new G4UnitDefinition("g/mole","g/mol","Molar mass",g/mole);
 
 //Activity
 new G4UnitDefinition("becquerel","Bq","Activity",becquerel);
 new G4UnitDefinition(    "curie","Ci","Activity",curie);
 
 //Dose
 new G4UnitDefinition("gray","Gy","Dose",gray);                          
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void G4UnitDefinition::PrintUnitsTable()
{
  G4cout << "\n          ----- The Table of Units ----- \n";
  if (!pUnitsTable)  { pUnitsTable = new G4UnitsTable; }
  for(size_t i=0;i<pUnitsTable->size();i++)
  {
    (*pUnitsTable)[i]->PrintCategory();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4UnitDefinition::ClearUnitsTable()
{
#ifdef G4MULTITHREADED
  delete pUnitsTable;
  pUnitsTable = nullptr;
  if(G4Threading::IsMasterThread())
  { pUnitsTableShadow = nullptr; }
#else
  for (size_t i=0;i<pUnitsTable->size();i++)
  {
    delete (*pUnitsTable)[i];
  }
  pUnitsTable->clear();
#endif
  unitsTableDestroyed = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
   
G4UnitsCategory::G4UnitsCategory(const G4String& name)
  : Name(name),UnitsList(),NameMxLen(0),SymbMxLen(0)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4UnitsCategory::~G4UnitsCategory()
{
  for(size_t i=0;i<UnitsList.size();i++)
  {
    delete UnitsList[i];
  }
  UnitsList.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4UnitsCategory::G4UnitsCategory(const G4UnitsCategory& right)
{
  *this = right;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4UnitsCategory& G4UnitsCategory::operator=(const G4UnitsCategory& right)
{
  if (this != &right)
    {
      Name      = right.Name;
      UnitsList = right.UnitsList;
      NameMxLen = right.NameMxLen;
      SymbMxLen = right.SymbMxLen;
    }
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4int G4UnitsCategory::operator==(const G4UnitsCategory& right) const
{
  return (this == (G4UnitsCategory *) &right);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4int G4UnitsCategory::operator!=(const G4UnitsCategory &right) const
{
  return (this != (G4UnitsCategory *) &right);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void G4UnitsCategory::PrintCategory()
{
  G4cout << "\n  category: " << Name << G4endl;
  for(size_t i=0;i<UnitsList.size();i++)
    { UnitsList[i]->PrintDefinition(); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
       
G4BestUnit::G4BestUnit(G4double value, const G4String& category)
  : nbOfVals(1)
{
 // find the category
    G4UnitsTable& theUnitsTable = G4UnitDefinition::GetUnitsTable();
    size_t nbCat = theUnitsTable.size();
    size_t i = 0;
    while ((i<nbCat)&&(theUnitsTable[i]->GetName()!=category)) { i++; }
    if (i == nbCat) 
       {
         G4cout << " G4BestUnit: the category " << category 
                << " does not exist !!" << G4endl;
         G4Exception("G4BestUnit::G4BestUnit()", "InvalidCall",
                     FatalException, "Missing unit category !") ;
       }  
  //
    Value[0] = value;
    Value[1] = 0.;
    Value[2] = 0.;
    IndexOfCategory = i;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
       
G4BestUnit::G4BestUnit(const G4ThreeVector& value, const G4String& category)
  : nbOfVals(3)
{
 // find the category
    G4UnitsTable& theUnitsTable = G4UnitDefinition::GetUnitsTable();
    size_t nbCat = theUnitsTable.size();
    size_t i = 0;
    while ((i<nbCat)&&(theUnitsTable[i]->GetName()!=category)) { i++; }
    if (i == nbCat) 
       {
         G4cerr << " G4BestUnit: the category " << category 
                << " does not exist." << G4endl;
         G4Exception("G4BestUnit::G4BestUnit()", "InvalidCall",
                     FatalException, "Missing unit category !") ;
       }  
  //
    Value[0] = value.x();
    Value[1] = value.y();
    Value[2] = value.z();
    IndexOfCategory = i;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4BestUnit::~G4BestUnit()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4BestUnit::operator G4String () const
{
  std::ostringstream oss;
  oss << *this;
  return oss.str();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
std::ostream& operator<<(std::ostream& flux, G4BestUnit a)
{
  G4UnitsTable& theUnitsTable = G4UnitDefinition::GetUnitsTable();
  G4UnitsContainer& List = theUnitsTable[a.IndexOfCategory]
                           ->GetUnitsList();
  G4int len = theUnitsTable[a.IndexOfCategory]->GetSymbMxLen();
                           
  G4int    ksup(-1), kinf(-1);
  G4double umax(0.), umin(DBL_MAX);
  G4double rsup(DBL_MAX), rinf(0.);

  //for a ThreeVector, choose the best unit for the biggest value 
  G4double value = std::max(std::max(std::fabs(a.Value[0]),
                                     std::fabs(a.Value[1])),
                            std::fabs(a.Value[2]));

  for (size_t k=0; k<List.size(); k++)
     {
       G4double unit = List[k]->GetValue();
       if (!(value!=DBL_MAX))
         {if(unit>umax) {umax=unit; ksup=k;}}
       else if (value<=DBL_MIN)
         {if(unit<umin) {umin=unit; kinf=k;}}
       else
       {
         G4double ratio = value/unit;
         if ((ratio>=1.)&&(ratio<rsup)) {rsup=ratio; ksup=k;}
         if ((ratio< 1.)&&(ratio>rinf)) {rinf=ratio; kinf=k;}
       } 
     }
 
  G4int index=ksup;
  if(index==-1) { index=kinf; }
  if(index==-1) { index=0; }
  
  for (G4int j=0; j<a.nbOfVals; j++) 
     { flux << a.Value[j]/(List[index]->GetValue()) << " "; }

  std::ios::fmtflags oldform = flux.flags();

  flux.setf(std::ios::left,std::ios::adjustfield);
  flux << std::setw(len) << List[index]->GetSymbol();       
  flux.flags(oldform);

  return flux;
}       

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#ifdef G4MULTITHREADED

void G4UnitsTable::Synchronize()
{
  G4UnitsTable* orig = &(G4UnitDefinition::GetUnitsTableShadow());
  if(this==orig) return;

  G4UnitsTable::iterator utItr = orig->begin();
  for(;utItr!=orig->end();utItr++)
  {
    G4UnitsCategory* category = *utItr;
    G4String catName = category->GetName();
    G4UnitsContainer* units = &(category->GetUnitsList());
    G4UnitsContainer::iterator ucItr = units->begin();
    for(;ucItr!=units->end();ucItr++)
    {
      G4UnitDefinition* unit = *ucItr;
      if(!Contains(unit,catName))
      {
        new G4UnitDefinition(unit->GetName(),unit->GetSymbol(),
                               catName,unit->GetValue());
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4UnitsTable::Contains(const G4UnitDefinition* unit,
                              const G4String& categoryName)
{
  G4UnitsTable::iterator utItr = begin();
  for(;utItr!=end();utItr++)
  {
    G4UnitsCategory* category = *utItr;
    G4String catName = category->GetName();
    if(catName!=categoryName) continue;
    G4UnitsContainer* units = &(category->GetUnitsList());
    G4UnitsContainer::iterator ucItr = units->begin();
    for(;ucItr!=units->end();ucItr++)
    {
      if((*ucItr)->GetName()==unit->GetName() &&
         (*ucItr)->GetSymbol()==unit->GetSymbol())
      { return true; }
    }
  }
  return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
 
