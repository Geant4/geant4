
#include "UIcmdWithNucleusAndUnit.hh"
#ifdef WIN32
#  include "g4std/Strstream"
#else
#  include "g4std/strstream"
#endif
////////////////////////////////////////////////////////////////////////////////
//
UIcmdWithNucleusAndUnit::UIcmdWithNucleusAndUnit
(const char * theCommandPath,G4UImessenger * theMessenger)
:G4UIcommand(theCommandPath,theMessenger)
{
  G4UIparameter * intParamA = new G4UIparameter('i');
  SetParameter(intParamA);
  G4UIparameter * intParamZ = new G4UIparameter('i');
  SetParameter(intParamZ);
  G4UIparameter * dblParamE = new G4UIparameter('d');
  SetParameter(dblParamE);
  G4UIparameter * untParam = new G4UIparameter('s');
  SetParameter(untParam);
  untParam->SetParameterName("Unit");
}
////////////////////////////////////////////////////////////////////////////////
//
UIcmdWithNucleusAndUnit::~UIcmdWithNucleusAndUnit()
{
  ;
}
////////////////////////////////////////////////////////////////////////////////
//
Nucleus UIcmdWithNucleusAndUnit::GetNewNucleusValue(G4String paramString)
{
  G4int a;
  G4int z;
  G4double e;
  char unts[30];

  const char* t = paramString;
  G4std::istrstream is((char*)t);
  is >> a >> z >> e >>unts;
  G4String unt = unts;

  return Nucleus(a,z,e*ValueOf(unt));
}

G4double UIcmdWithNucleusAndUnit::GetNewUnitValue(G4String paramString)
{
  G4int a;
  G4int z;
  G4double e;  

  char unts[30];
  
  const char* t = paramString;
  G4std::istrstream is((char*)t);
  is >> a >> z >> e  >> unts;

  G4String unt = unts;
  
  return ValueOf(unt);
}

////////////////////////////////////////////////////////////////////////////////
//
G4String UIcmdWithNucleusAndUnit::ConvertToString(Nucleus def, 
						    const char *unitName)
{
  G4String unt = unitName;
  G4double uv = ValueOf(unitName);

  char st[100];
  G4std::ostrstream os(st,100);
  os << def.GetA() << " " << def.GetZ()
     << " "<< def.GetE()/uv<<" "<< unitName <<  '\0';
  G4String vl = st;
  return vl;
}                         
////////////////////////////////////////////////////////////////////////////////
//
void UIcmdWithNucleusAndUnit::SetParameterName
(const char * theNameA, const char * theNameZ,
const char * theNameE,G4bool omittable,G4bool currentAsDefault)
{
  G4UIparameter * theParamA = GetParameter(0);
  theParamA->SetParameterName(theNameA);
  theParamA->SetOmittable(omittable);
  theParamA->SetCurrentAsDefault(currentAsDefault);
  G4UIparameter * theParamZ = GetParameter(1);
  theParamZ->SetParameterName(theNameZ);
  theParamZ->SetOmittable(omittable);
  theParamZ->SetCurrentAsDefault(currentAsDefault);
  G4UIparameter * theParamE = GetParameter(2);
  theParamE->SetParameterName(theNameE);
  theParamE->SetOmittable(omittable);
  theParamE->SetCurrentAsDefault(currentAsDefault);
}
////////////////////////////////////////////////////////////////////////////////
//
void UIcmdWithNucleusAndUnit::SetDefaultValue(Nucleus def)
{
  G4UIparameter * theParamA = GetParameter(0);
  theParamA->SetDefaultValue(def.GetA());
  G4UIparameter * theParamZ = GetParameter(1);
  theParamZ->SetDefaultValue(def.GetZ());
  G4UIparameter * theParamE = GetParameter(2);
  theParamE->SetDefaultValue(def.GetE());
}


void UIcmdWithNucleusAndUnit::SetUnitCandidates(const char * candidateList)
{
  G4UIparameter * untParam = GetParameter(3);
  G4String canList = candidateList;
  untParam->SetParameterCandidates(canList);
}

void UIcmdWithNucleusAndUnit::SetDefaultUnit(const char * defUnit)
{
  G4UIparameter * untParam = GetParameter(3);
  untParam->SetOmittable(true);
  untParam->SetDefaultValue(defUnit);
}












