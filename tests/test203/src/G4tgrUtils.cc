#include "G4tgrUtils.hh"
#include "G4tgrParameterMgr.hh"
#include "G4tgrMessenger.hh"
#include "G4UnitsTable.hh"
#include <iomanip>

using namespace CLHEP;

HepTool::Evaluator* G4tgrUtils::theCLHEPevaluator = new HepTool::Evaluator;

//-------------------------------------------------------------
bool G4tgrUtils::IsSeparator( const char ch)
{
  char nonCharacters[8] = {".()+-*/"};
  for( uint ii = 0; ii < 7; ii++ ){
    //-    G4cout << ii << " IsSeparator " << ch << " =?= " << nonCharacters[ii] << G4endl;
    if( ch == nonCharacters[ii] ){
      return TRUE;
    }
  }

  return FALSE;

}

//-------------------------------------------------------------
bool G4tgrUtils::IsNumber( const G4String& str)
{
  int isnum = 1;
  int numE = 0;
  for(uint ii=0; ii<str.length(); ii++){
    if(!isdigit(str[ii]) && str[ii]!='.' && str[ii]!='-' && str[ii]!='+') {
      //--- check for E(xponential)
      if(str[ii] == 'E' || str[ii] == 'e' ) {
        if(numE != 0 || ii == str.length()-1)  {
	  isnum = 0;
	  break;
	}
	numE++;
      } else {
	isnum = 0; 
	break;
      }
    }
  }
 
  return isnum;
}


//-------------------------------------------------------------
void G4tgrUtils::Dump3v( const Hep3Vector& vec, const char* msg) 
{
  G4cout << msg << setprecision(8) << vec << G4endl;
}


//-------------------------------------------------------------
void G4tgrUtils::Dumprm( const HepRotation& rm, const char* msg) {

  G4cout << msg << G4endl
       << " xx=" << rm.xx() << " yx=" << rm.yx() << " zx=" << rm.zx() << G4endl
       << " xy=" << rm.xy() << " yy=" << rm.yy() << " zy=" << rm.zy() << G4endl
       << " xz=" << rm.xz() << " yz=" << rm.yz() << " zz=" << rm.zz() << G4endl;

}


//-------------------------------------------------------------
void G4tgrUtils::DumpVS( const vector<G4String>& wl , const char* msg, ostream& outs) 
{
  outs << msg << G4endl;
  /*  ostream_iterator<G4String> os(outs," ");
  copy(wl.begin(), wl.end(), os);
  outs << G4endl;*/
  vector<G4String>::const_iterator ite;
  for( ite = wl.begin(); ite != wl.end(); ite++ ){
    outs << *ite << " ";
  }
  outs << G4endl;
}


//-------------------------------------------------------------
void G4tgrUtils::DumpVS( const vector<G4String>& wl , const char* msg) 
{
  DumpVS( wl, msg, G4cout);
}


//-------------------------------------------------------------
G4String G4tgrUtils::SubColon( const G4String& str ) 
{
  
  if( str.find(':') != 0 ) { 
    cerr << "!!!EXITING trying to substract leading colon from a word that has no leading colon: " << str << G4endl;
    exit(1);
  }

  //  str = str.strip(G4String::leading, '\:');
  G4String strt = str.substr(1,str.size()-1);
  return strt;

}

//-------------------------------------------------------------
G4String G4tgrUtils::SubQuotes( const G4String& str ) 
{
  return str;

  //---------- Take out leading and trailing '"'
  if( str.find('"') != 0 || str.rfind('"') != str.length()-1 ) { 
    cerr << "!!!EXITING tryin g to substract quotes from a word that has no quotes: " << str << G4endl;
    exit(1);
  }

  //  str = str.strip(G4String::both, '\"');
  //---------- Take out leading and trallling '"'
  G4String strt = str.substr(1,str.size()-2);

#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
     G4cout << " G4tgrUtils subquotes " << str << G4endl;
#endif
  //---------- Look for leading spaces
  while( strt.c_str()[0] == ' ' ) {
   strt = strt.substr(1,strt.size()-1);
  }

  //---------- Look for trailing spaces
  while( strt[strt.size()-1] == ' ' ) {
   strt = strt.substr(0,strt.size()-1);
  }

  return strt;

}

//-------------------------------------------------------------
/*G4String G4tgrUtils::GetFloatString( const G4String& str ) 
{
  fl = atof( str.c_str() );
  return ftoa( fl );
}*/

//-------------------------------------------------------------
double G4tgrUtils::GetFloat( const G4String& str, G4double unitval ) 
{
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 3 ) 
    G4cout << "GetFloat processing " << str << " unitval " << unitval << G4endl;
#endif

  G4String strnew = str;
 
  G4int strlen = strnew.length();
 
  //----- Check if there is a unit defined, else take unitval
  G4int iast = strnew.rfind('*');
  if( iast != -1 ) {
    G4String unitstr = strnew.substr( iast+1, strlen );
    if(G4UnitDefinition::GetCategory(unitstr) != "None" ){
      unitstr = ftoa( G4UnitDefinition::GetValueOf(unitstr) );
      strnew = strnew.substr( 0, iast+1 ) + unitstr;
      strlen = strnew.length();
#ifdef G4VERBOSE
      if( G4tgrMessenger::GetVerboseLevel() >= 3 ) 
	G4cout << " unit string found " << unitstr << " in string " << strnew << G4endl;
#endif

      unitval = 1.; // do not use unitval value, dimension is written in string
    }
  }

  //----- Replace all parameters 
  G4int paramStart = -1;
  G4String strnew2;
  G4tgrParameterMgr* parmgr = G4tgrParameterMgr::GetInstance();
  for(uint ii=0; ii<strlen; ii++) {
    if( strnew.c_str()[ii] == '$' ) {
      paramStart = ii;
#ifdef G4VERBOSE
      if( G4tgrMessenger::GetVerboseLevel() >= 3 ) 
	G4cout << "GetFloat paramStart " << paramStart << G4endl;
#endif
    } else {
      if( paramStart != -1 ) {
	if( IsSeparator( strnew.c_str()[ii] ) | ii == strlen-1 ) { 
	  G4int paramLength = ii-paramStart-1; 
	  if( ii == strlen-1 ) { 
	    paramLength++;
	  }
	  
	  G4String parstr = parmgr->FindParameter( strnew.substr(paramStart+1, paramLength ));
	  if( parstr.substr(0,1) == "-" ){
	    strnew2 += "(";
	  }
	  strnew2 += parstr;
	  if( parstr.substr(0,1) == "-" ){
	    strnew2 += ")";
	  }
#ifdef G4VERBOSE
	  if( G4tgrMessenger::GetVerboseLevel() >= 3 ) 
	    G4cout << " param found " <<strnew.substr(paramStart+1, paramLength )<< " in string " << strnew << G4endl;
#endif	  
	  if( ii != strlen-1 ) { 
	    strnew2 += strnew.substr(ii,1);
	  }
	  paramStart = -1;
	}
      } else {
	strnew2 += strnew.substr(ii,1);
      }

    }
  }

  G4double val = theCLHEPevaluator->evaluate( strnew2 );
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 3 ) 
    G4cout << " GetFloat:: val " << val << " * " << unitval << " from string " << str << " converted to " << strnew2 << G4endl;
#endif
  return val * unitval;

}


//-------------------------------------------------------------
int G4tgrUtils::GetInt( const G4String& str ) 
{
  double val;
  //----------- first check if it is parameter
  if( str.c_str()[0] == '$' ) {
    val = GetFloat(G4tgrParameterMgr::GetInstance()->FindParameter( str.substr(1,str.size()), 1) );
  } else {
    //----------- first check that it is an integer
    if(!IsNumber(str) ) {
      //----- Check that it is a number 
      cerr << "!!!! EXITING: trying to get the integer from a G4String that is not a number " << str << G4endl;
      exit(1);
    } else {
      //----- Check that it is not a float, no decimal or E-n
      uint ii;
      
      bool isFloat = 0;
      int ch = str.find('.');
      if(ch != -1 ) {
	for(ii = ch+1; ii < str.size(); ii++) {
	if( str[ii] != '0' ) isFloat = 1;
	}
      }
      
      ch = str.find('E');
      if(ch != -1 ) ch = str.find('e');
      if(ch != -1 ) {
	if(str.c_str()[ch+1] == '-') isFloat = 1;
      }
      
      if(isFloat) {
	cerr << "!!!! EXITING: trying to get the integer from a G4String that is a float: " << str << G4endl;
	cerr << ii << " ii "  << ch <<G4endl;
	exit(1);
      }
    }
    val = atof( str.c_str() );
  }
  return int( val );
}



//-------------------------------------------------------------
bool G4tgrUtils::GetBool( const G4String& str ) 
{
   bool val;
  
 //t str = upper( str );
  //----------- first check that it is a not number
  if( str == "ON" || str == "TRUE"  ) {
    val = true;
  } else if( str == "OFF" || str == "FALSE" ) {
    val = false;
  } else {
    cerr << "!!!! EXITING: trying to get the float from a G4String that is not 'ON'/'OFF'/'TRUE'/'FALSE' " << str << G4endl;
    exit(1);
  }

  return val;
}

//-------------------------------------------------------------
G4String G4tgrUtils::ftoa( const G4double dou ) 
{
  char chartmp[20];
  gcvt( dou, 10, chartmp );

  G4String str = G4String(chartmp);
  return str;

}


//-------------------------------------------------------------
void G4tgrUtils::CheckWLsize( const vector<G4String>& wl, uint noWords, WLSIZEtype st, const G4String& methodName )
{
  G4String outstr = G4String("!!!! EXITING: ") + methodName + G4String(".  Line read with number of words ");
  bool isOK = 1;

  uint wlsize = wl.size();
  switch (st) {
  case WLSIZE_EQ:
    if( wlsize != noWords ) {
      isOK = 0;
      outstr += G4String("not equal than ");
    }
    break;
  case WLSIZE_NE:
    if( wlsize = noWords ) {
      isOK = 0;
      outstr += G4String("equal than ");
    }
    break;
  case WLSIZE_LE:
    if( wlsize > noWords ) {
      isOK = 0;
      outstr += G4String("greater than ");
    }
    break;
  case WLSIZE_LT:
    if( wlsize >= noWords ) {
      isOK = 0;
      outstr += G4String("greater or equal than ");
    }
    break;
  case WLSIZE_GE:
    if( wlsize < noWords ) {
      isOK = 0;
      outstr += G4String("less than ");
    }
    break;
  case WLSIZE_GT:
    if( wlsize <= noWords ) {
      isOK = 0;
      outstr += G4String("less or equal than ");
    }
    break;
  default:
    cerr << " !!!! EXITING: G4tgrUtils::checkWLsize. type of WLSIZE type not found " << st << G4endl;
    break;
  }

  if( !isOK ){ 
    char chartmp[20];
    gcvt( noWords, 10, chartmp );
    outstr += (G4String(chartmp) + G4String(" words") );
    DumpVS( wl, outstr.c_str() );
    //    cerr << noWords << " words " << G4endl;
    G4Exception("");
  }

}
