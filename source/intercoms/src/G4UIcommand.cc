// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIcommand.cc,v 1.3 1999-11-11 15:36:06 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "G4UIcommand.hh"
#include "G4UImessenger.hh"
#include "G4UImanager.hh"
#include "G4UIcommandStatus.hh"
#include "G4StateManager.hh"
#include "G4UnitsTable.hh"
#include "g4rw/ctoken.h"
#include "G4ios.hh"

G4UIcommand::G4UIcommand():paramERR(0) { }

G4UIcommand::G4UIcommand(const char * theCommandPath,
			 G4UImessenger * theMessenger)
:messenger(theMessenger), paramERR(0)
{
  G4String comStr = theCommandPath;
  if(!theMessenger)
  { // this must be a directory
    if(comStr(comStr.length()-1)!='/')
    {
      G4cerr << "G4UIcommand Warning : " << endl;
      G4cerr << "  <" << theCommandPath << "> must be a directory." << endl;
      G4cerr << "  '/' is appended." << endl;
      comStr += "/";
    }
  }
  G4UIcommandCommonConstructorCode (comStr);
  G4String nullString;
  availabelStateList.clear();
  availabelStateList.insert(PreInit);
  availabelStateList.insert(Init);
  availabelStateList.insert(Idle);
  availabelStateList.insert(GeomClosed);
  availabelStateList.insert(EventProc);
}

void G4UIcommand::G4UIcommandCommonConstructorCode
(const char * theCommandPath)
{ 
  commandPath = theCommandPath;
  commandName = theCommandPath;
  int commandNameIndex = commandName.last('/');
  commandName.remove(0,commandNameIndex+1);

  G4UImanager::GetUIpointer()->AddNewCommand(this);
}

G4UIcommand::~G4UIcommand()
{
  G4UImanager* fUImanager = G4UImanager::GetUIpointer();
  if(fUImanager) fUImanager->RemoveCommand(this);
  
  int n_parameterEntry = parameter.entries();
  for( int i_thParameter=0; i_thParameter < n_parameterEntry; i_thParameter++ )
  { delete parameter[i_thParameter]; }
}

int G4UIcommand::operator==(const G4UIcommand &right) const
{
  return ( commandPath == right.GetCommandPath() );
}

int G4UIcommand::operator!=(const G4UIcommand &right) const
{
  return ( commandPath != right.GetCommandPath() );
}

G4int G4UIcommand::DoIt(G4String parameterList)
{
  G4String correctParameters;
  int n_parameterEntry = parameter.entries();
  if( n_parameterEntry != 0 )
  {
    G4String aToken;
    G4String correctToken;
    RWCTokenizer parameterToken( parameterList );
    for( int i_thParameter=0; i_thParameter<n_parameterEntry; i_thParameter++ )
    {
      if(i_thParameter > 0)
      {
        correctParameters.append(" ");
      }
      aToken = parameterToken();
      if( aToken.length()>0 && aToken(0)=='"' )
      {
        while( aToken(aToken.length()-1) != '"' )
        {
          G4String additionalToken = parameterToken();
          if( additionalToken.isNull() )
          { return fParameterUnreadable+i_thParameter; }
          aToken += " ";
          aToken += additionalToken;
        }
        // aToken.strip(G4String::both,'"');
      }
      if( aToken.isNull() || aToken == "!" )
      {
        if(parameter[i_thParameter]->IsOmittable())
        { 
          if(parameter[i_thParameter]->GetCurrentAsDefault())
          {
            RWCTokenizer cvt(messenger->GetCurrentValue(this));
            G4String parVal;
            for(int ii=0;ii<i_thParameter;ii++)
            { parVal = cvt(); }
	    G4String aCVToken = cvt();
	    if (aCVToken(0)=='"')
	    {
	      while( aCVToken(aCVToken.length()-1) != '"' )
	      {
		G4String additionalToken = cvt();
		if( additionalToken.isNull() )
		{ return fParameterUnreadable+i_thParameter; }
		aCVToken += " ";
		aCVToken += additionalToken;
	      }
	      // aCVToken.strip(G4String::both,'"');
	    }
            correctParameters.append(aCVToken);
          }
          else
          { correctParameters.append(parameter[i_thParameter]->GetDefaultValue()); }
        }
        else
        { return fParameterUnreadable+i_thParameter; }
      }
      else
      {
        int stat = parameter[i_thParameter]->CheckNewValue( aToken );
        if(stat) return stat+i_thParameter;
        correctParameters.append(aToken);
      }
    }
  }

  if(CheckNewValue( correctParameters ))
  { return fParameterOutOfRange+99; }

  messenger->SetNewValue( this, correctParameters );
  return 0;
}

G4String G4UIcommand::GetCurrentValue()
{
  return messenger->GetCurrentValue(this);
}

void G4UIcommand::AvailableForStates(G4ApplicationState s1)
{
  availabelStateList.clear();
  availabelStateList.insert(s1);
}

void G4UIcommand::AvailableForStates(G4ApplicationState s1,
                                     G4ApplicationState s2)
{
  availabelStateList.clear();
  availabelStateList.insert(s1);
  availabelStateList.insert(s2);
}

void G4UIcommand::AvailableForStates(G4ApplicationState s1,
                                     G4ApplicationState s2,
                                     G4ApplicationState s3)
{
  availabelStateList.clear();
  availabelStateList.insert(s1);
  availabelStateList.insert(s2);
  availabelStateList.insert(s3);
}

void G4UIcommand::AvailableForStates(G4ApplicationState s1,
                                     G4ApplicationState s2,
                                     G4ApplicationState s3,
                                     G4ApplicationState s4)
{
  availabelStateList.clear();
  availabelStateList.insert(s1);
  availabelStateList.insert(s2);
  availabelStateList.insert(s3);
  availabelStateList.insert(s4);
}

void G4UIcommand::AvailableForStates(G4ApplicationState s1,
                                     G4ApplicationState s2,
                                     G4ApplicationState s3,
                                     G4ApplicationState s4,
                                     G4ApplicationState s5)
{
  availabelStateList.clear();
  availabelStateList.insert(s1);
  availabelStateList.insert(s2);
  availabelStateList.insert(s3);
  availabelStateList.insert(s4);
  availabelStateList.insert(s5);
}

G4bool G4UIcommand::IsAvailable()
{
  G4bool av = false;
  G4ApplicationState currentState 
   = G4StateManager::GetStateManager()->GetCurrentState();
   
  int nState = availabelStateList.entries();
  for(int i=0;i<nState;i++)
  {
  	if(availabelStateList[i]==currentState)
  	{
  	  av = true;
  	  break;
  	}
  }

  return av;
}

G4double G4UIcommand::ValueOf(G4String unitName)
{
   G4double value = 0.;
   value = G4UnitDefinition::GetValueOf(unitName);
   return value;              
}

G4String G4UIcommand::CategoryOf(G4String unitName)
{
   return G4UnitDefinition::GetCategory(unitName);
}

G4String G4UIcommand::UnitsList(G4String unitCategory)
{
  G4String retStr;
  G4UnitsTable& UTbl = G4UnitDefinition::GetUnitsTable();
  G4int i;
  for(i=0;i<UTbl.entries();i++)
  { if(UTbl[i]->GetName()==unitCategory) break; }
  if(i==UTbl.entries())
  { 
    G4cerr << "Unit category <" << unitCategory << "> is not defined." << endl;
    return retStr;
  }
  G4UnitsContainer& UCnt = UTbl[i]->GetUnitsList();
  retStr = UCnt[0]->GetSymbol();
  G4int je = UCnt.entries();
  for(int j=1;j<je;j++)
  {
    retStr += " ";
    retStr += UCnt[j]->GetSymbol();
  }
  for(int k=0;k<je;k++)
  {
    retStr += " ";
    retStr += UCnt[k]->GetName();
  }
  return retStr;
}
  
void G4UIcommand::List()
{
  G4cout << endl;
  G4cout << endl;
  if(commandPath(commandPath.length()-1)!='/')
  { G4cout << "Command " << commandPath << endl; }
  G4cout << "Guidance :" << endl;
  int n_guidanceEntry = commandGuidance.entries();
  for( int i_thGuidance=0; i_thGuidance < n_guidanceEntry; i_thGuidance++ )
  { G4cout << commandGuidance[i_thGuidance] << endl; }
  if( ! rangeString.isNull() )
  { G4cout << " Range of parameters : " << rangeString << endl; }
  int n_parameterEntry = parameter.entries();
  if( n_parameterEntry > 0 )
  {
    for( int i_thParameter=0; i_thParameter<n_parameterEntry; i_thParameter++ )
    { parameter[i_thParameter]->List(); }
  }
  G4cout << endl;
}

// ----- the following is used by CheckNewValue()  ------------


#include <ctype.h>          // isalpha(), toupper()
#ifdef WIN32
#include <strstrea.h>
#else
#include <strstream.h>
#endif

//#include "checkNewValue_debug.icc"
//#define DEBUG 1

int G4UIcommand::
CheckNewValue(G4String newValue)
{
   yystype result;
   // if( TypeCheck(newValue) == 0 ) return 1;
   if( ! rangeString.isNull() )
   { if( RangeCheck(newValue) == 0 ) return fParameterOutOfRange; }
   return 0;  // succeeded
}

// ------------------ type check routines -------------------

int G4UIcommand::
TypeCheck(G4String newValues)
{
    G4String aNewValue;
    char type;
    const char* t = newValues;
    istrstream is((char*)t);
    for (unsigned i=0; i< parameter.entries(); i++) {
        is >> aNewValue;
        type = toupper(parameter(i)->GetParameterType());
        switch ( type ) {
            case 'D':
                if( IsDouble(aNewValue)==0 ){
                    G4cerr << aNewValue << ": double value expected."
                         << endl;
                    return 0;
                } break;
            case 'I':
                if( IsInt(aNewValue,20)==0 ){
                    G4cerr <<aNewValue<<": integer expected."
                         <<endl;
                     return 0;
                } break;
            case 'S':
                break;
            case 'B':
                aNewValue.toUpper(); 
                if (aNewValue == "Y"   || aNewValue == "N"
                  ||aNewValue == "YES" || aNewValue == "NO"
                  ||aNewValue == "1"   || aNewValue == "0"
                  ||aNewValue == "T"   || aNewValue == "F"
                  ||aNewValue == "TRUE" || aNewValue == "FALSE")
                return 1;
                else return 0; 
                break;
            default:  ;
        }
    }
    return 1;
}


int G4UIcommand::
IsInt(const char* buf, short maxDigits)
{
    const char* p= buf;
    int length=0;
    if( *p == '+' || *p == '-') { ++p; }
    if( isdigit( (int)(*p) )) {
        while( isdigit( (int)(*p) )) { ++p;  ++length; }
        if( *p == '\0' ) {
            if( length > maxDigits) {
                G4cerr <<"digit length exceeds"<<endl;
                return 0;
            }
            return 1;
        } else {
            // G4cerr <<"illegal character after int:"<<buf<<endl;
        }
    } else {
        // G4cerr <<"illegal int:"<<buf<<endl;
    }
    return 0;
}


int G4UIcommand::
ExpectExponent(const char* str)   // used only by IsDouble()
{
    int maxExplength;
    if( IsInt( str, maxExplength=7 )) return 1;
    else return 0;
}


int G4UIcommand::
IsDouble(const char* buf)
{
    const char* p= buf;
    switch( *p) {
        case '+':  case '-': ++p;
            if( isdigit(*p) ) {
                 while( isdigit( (int)(*p) )) { ++p; }
                 switch ( *p ) {
                     case '\0':    return 1;
		        // break;
                     case 'E':  case 'e':
                         return ExpectExponent(++p );
			// break;
                     case '.':  ++p;
                         if( *p == '\0' )  return 1;
                         if( *p == 'e' || *p =='E' ) return ExpectExponent(++p );
                         if( isdigit(*p) ) {
                             while( isdigit( (int)(*p) )) { ++p; }
                             if( *p == '\0' )  return 1;
                             if( *p == 'e' || *p =='E') return ExpectExponent(++p);
                         } else return 0;   break;
                     default: return 0;
                 }
            }
            if( *p == '.' ) { ++p;
                 if( isdigit(*p) ) {
                     while( isdigit( (int)(*p) )) { ++p; }
                     if( *p == '\0' )  return 1;
                     if( *p == 'e' || *p =='E')  return ExpectExponent(++p);
                 }
            }
            break;
        case '.':  ++p;
            if( isdigit(*p) ) {
                 while( isdigit( (int)(*p) )) { ++p; }
                 if( *p == '\0' )  return 1;
                 if( *p == 'e' || *p =='E' )  return ExpectExponent(++p);
            }    break;
        default: // digit is expected
            if( isdigit(*p) ) {
                 while( isdigit( (int)(*p) )) { ++p; }
                 if( *p == '\0' )  return 1;
                 if( *p == 'e' || *p =='E')  return ExpectExponent(++p);
                 if( *p == '.' ) { ++p;
                      if( *p == '\0' )  return 1;
                      if( *p == 'e' || *p =='E')  return ExpectExponent(++p);
                      if( isdigit(*p) ) {
                          while( isdigit( (int)(*p) )) { ++p; }
                          if( *p == '\0' )  return 1;
                          if( *p == 'e' || *p =='E') return ExpectExponent(++p);
                      }
                 }
            }
     }
     return 0;
}


// ------------------ range Check routines -------------------
int G4UIcommand::
RangeCheck(G4String newValue) {
    yystype result;
    char type;
    bp = 0;                 // reset buffer pointer for G4UIpGetc()
    const char* t = newValue;
    istrstream is((char*)t);
    for (unsigned i=0; i< parameter.entries(); i++) {
        type= toupper(parameter(i)->GetParameterType());
        switch ( type ) {
            case 'D':  is >> newVal(i).D;  break;
            case 'I':  is >> newVal(i).I;  break;
            case 'S':
            case 'B':
            default:  ;
        }
   }
   // PrintToken();          // Print tokens (consumes all tokens)
   token= Yylex();
   result = Expression();

   if( paramERR == 1 ) return 0;
   if( result.type != CONSTINT) {
      G4cerr << "Illegal Expression in parameter range." << endl;
      return 0;
   }
   if ( result.I ) return 1;
   G4cerr << "parameter out of range: "<< rangeString << endl;
   return 0;
}

// ------------------ syntax node functions  ------------------
yystype G4UIcommand:: 
Expression(void)
{
    yystype result;
    #ifdef DEBUG
        G4cerr << " Expression()" << endl;
    #endif
    result = LogicalORExpression();
    return result;
}

yystype G4UIcommand:: 
LogicalORExpression(void)
{
    yystype result;
    yystype p;
    p = LogicalANDExpression();
    if( token != LOGICALOR)  return p;
    if( p.type == CONSTSTRING || p.type == IDENTIFIER ) {
        G4cerr << "Parameter range: illegal type at '||'" << endl;
        paramERR = 1;
    }
    result.I = p.I;
    while (token == LOGICALOR) 
    {  
        token = Yylex();
        p = LogicalANDExpression();
        if( p.type == CONSTSTRING || p.type == IDENTIFIER ) {
            G4cerr << "Parameter range: illegal type at '||'" <<endl;
            paramERR = 1;
        }
        switch (p.type) {
            case CONSTINT: 
                result.I  += p.I; 
                result.type = CONSTINT;      break;
            case CONSTDOUBLE:
                result.I += (p.D != 0.0); 
                result.type = CONSTINT;      break;
            default: 
                G4cerr << "Parameter range: unknown type"<<endl; 
                paramERR = 1;
        } 
    }
    return result;
}

yystype G4UIcommand:: 
LogicalANDExpression(void)
{
    yystype result;
    yystype p;
    p = EqualityExpression();
    if( token != LOGICALAND)  return p;
    if( p.type == CONSTSTRING || p.type == IDENTIFIER ) {
        G4cerr << "Parameter range: illegal type at '&&'" << endl;
        paramERR = 1;
    }
    result.I = p.I;
    while (token == LOGICALAND)
    {
        token = Yylex();
        p = EqualityExpression();
        if( p.type == CONSTSTRING || p.type == IDENTIFIER ) {
            G4cerr << "Parameter range: illegal type at '&&'" << endl;
            paramERR = 1;
        }
        switch (p.type) {
            case CONSTINT:
                result.I  *= p.I;
                result.type = CONSTINT;      break;
            case CONSTDOUBLE:
                result.I *= (p.D != 0.0);
                result.type = CONSTINT;      break;
            default:
                G4cerr << "Parameter range: unknown type."<< endl;
                paramERR = 1;
        } 
    }
    return result;
}


yystype G4UIcommand:: 
EqualityExpression(void)
{ 
    yystype  arg1, arg2;
    int operat;
    yystype result;
    #ifdef DEBUG
        G4cerr << " EqualityExpression()" <<endl;
    #endif
    result = RelationalExpression();
    if( token==EQ || token==NE ) {
        operat = token;
        token =  Yylex();
        arg1 = result;
        arg2 = RelationalExpression();
        result.I = Eval2( arg1, operat, arg2 );   // semantic action
        result.type = CONSTINT;
        #ifdef DEBUG
            G4cerr << " return code of Eval2(): " << result.I <<endl;
        #endif
    } else {
        if (result.type != CONSTINT && result.type != CONSTDOUBLE) {  
            G4cerr << "Parameter range: error at EqualityExpression"
                 << endl;
            paramERR = 1;
        }
    }
    return  result;
}


yystype G4UIcommand:: 
RelationalExpression(void)
{ 
    yystype  arg1, arg2;
    int operat;
    yystype result;
    #ifdef DEBUG
        G4cerr << " RelationalExpression()" <<endl;
    #endif

    arg1 = AdditiveExpression();
    if( token==GT || token==GE || token==LT || token==LE  ) {
        operat = token;
        token =  Yylex();
        arg2 = AdditiveExpression();
        result.I = Eval2( arg1, operat, arg2 );    // semantic action
        result.type = CONSTINT;
        #ifdef DEBUG
            G4cerr << " return code of Eval2(): " << result.I << endl;
        #endif
    } else {
              result = arg1;
    }
    #ifdef DEBUG
       G4cerr <<" return RelationalExpression()"<< endl;
    #endif
    return  result;
}

yystype G4UIcommand::
AdditiveExpression(void)
{   yystype result;
    result = MultiplicativeExpression();
    if( token != '+' && token != '-' )  return result;
    G4cerr << "Parameter range: operator " 
         << (char)token 
         << " is not supported." << endl;
    paramERR = 1;
    return  result;
}

yystype G4UIcommand::
MultiplicativeExpression(void)
{   yystype result;
    result = UnaryExpression();
    if( token != '*' && token != '/' && token != '%' ) return result;
    G4cerr << "Parameter range: operator "
         << (char)token
         << " is not supported." << endl;
    paramERR = 1;
    return  result;
}

yystype G4UIcommand::
UnaryExpression(void)
{
    yystype result;
    yystype p;
    #ifdef DEBUG
        G4cerr <<" UnaryExpression"<< endl;
    #endif
    switch(token) {
        case '-':
            token = Yylex();
            p = UnaryExpression();
            if (p.type == CONSTINT) {
                result.I = - p.I;
                result.type = CONSTINT;
            }
            if (p.type == CONSTDOUBLE) {
                result.D = - p.D;
                result.type = CONSTDOUBLE;
            }                              break;
        case '+':
            token = Yylex();
            result = UnaryExpression();   break;
        case '!':
            token = Yylex();
            G4cerr << "Parameter range error: "
                 << "operator '!' is not supported (sorry)."
                 << endl;
            paramERR = 1;
            result = UnaryExpression();   break;
        default:
            result = PrimaryExpression();
    }
    return result;
}


yystype G4UIcommand:: 
PrimaryExpression(void)
{
     yystype result;
     #ifdef DEBUG
         G4cerr <<" primary_exp"<<endl;
     #endif
     switch (token) {
         case IDENTIFIER:
              result.S = yylval.S;
              result.type =  token;
              token = Yylex();           break;
         case CONSTINT:
              result.I = yylval.I;
              result.type =  token;
              token= Yylex();            break;
         case CONSTDOUBLE:
              result.D = yylval.D;
              result.type =  token;
              token = Yylex();           break;
         case '(' :
              token= Yylex();
              result = Expression();
              if( token !=  ')'  ) {
                  G4cerr << " ')' expected" << endl;
                  paramERR = 1;
              }
              token = Yylex();
                                         break;
         default:
         return result;
    }
    return result; // never executed
}

//---------------- semantic routines ---------------------------------

int G4UIcommand::
Eval2(yystype arg1, int op, yystype arg2)
{
    char newValtype;
    if( (arg1.type != IDENTIFIER) && (arg2.type != IDENTIFIER)) {
        G4cerr << commandName
             << ": meaningless comparison"
             << endl;
        paramERR = 1;
    }

    if( arg1.type == IDENTIFIER) {
        unsigned i = IndexOf( arg1.S );
        newValtype = toupper(parameter(i)->GetParameterType());
        switch ( newValtype ) {
            case 'I': 
                if( arg2.type == CONSTINT ) {
                    return CompareInt( newVal(i).I, op, arg2.I );
                } else {
                    G4cerr << "integer operand expected for "
                         <<  rangeString 
                         << '.' << endl;
                } break;
            case 'D':
                if( arg2.type == CONSTDOUBLE ) {
                    return CompareDouble( newVal(i).D, op, arg2.D );
                } else
                if ( arg2.type == CONSTINT ) {  // integral promotion
                    return CompareDouble( newVal(i).D, op, arg2.I );
                } break;
            default: ;
        }
    }
    if( arg2.type == IDENTIFIER) {
        unsigned i = IndexOf( arg2.S );
        newValtype = toupper(parameter(i)->GetParameterType());
        switch ( newValtype ) {
            case 'I': 
                if( arg1.type == CONSTINT ) {
                    return CompareInt( arg1.I, op, newVal(i).I );
                } else {
                    G4cerr << "integer operand expected for "
                         <<  rangeString
                         << '.' << endl;
                } break;
            case 'D':
                if( arg1.type == CONSTDOUBLE ) {
                    return CompareDouble( arg1.D, op, newVal(i).D );
                } else
                if ( arg1.type == CONSTINT ) {  // integral promotion
                    return CompareDouble( arg1.I, op, newVal(i).D );
                } break;
            default: ;
        }
    }
    return 0;
}

int G4UIcommand::
CompareInt(int arg1, int op, int arg2)
{   
    int result;
    G4String opr;
    switch (op) {
       case GT:  result = ( arg1 >  arg2); opr= ">" ;  break;
       case GE:  result = ( arg1 >= arg2); opr= ">=";  break;
       case LT:  result = ( arg1 <  arg2); opr= "<" ;  break;
       case LE:  result = ( arg1 <= arg2); opr= "<=";  break;
       case EQ:  result = ( arg1 == arg2); opr= "==";  break;
       case NE:  result = ( arg1 != arg2); opr= "!=";  break;
       default: 
           G4cerr << "Parameter range: error at CompareInt" << endl;
           paramERR = 1;
    }
    #ifdef DEBUG
        G4cerr << "CompareInt "
             << arg1 << " " << opr << arg2 
             << " result: " << result
             << endl;
    #endif
    return result;
}

int G4UIcommand::
CompareDouble(double arg1, int op, double arg2)
{   
    int result;
    G4String opr;
    switch (op) {
        case GT:  result = ( arg1 >  arg2); opr= ">";   break;
        case GE:  result = ( arg1 >= arg2); opr= ">=";  break;
        case LT:  result = ( arg1 <  arg2); opr= "<";   break;
        case LE:  result = ( arg1 <= arg2); opr= "<=";  break;
        case EQ:  result = ( arg1 == arg2); opr= "==";  break;
        case NE:  result = ( arg1 != arg2); opr= "!=";  break;
        default:
           G4cerr << "Parameter range: error at CompareDouble"
                << endl;
           paramERR = 1;
    }
    #ifdef DEBUG
        G4cerr << "CompareDouble " 
             << arg1 <<" " << opr << " "<< arg2
             << " result: " << result
             << endl;
    #endif
    return result;
}

unsigned G4UIcommand::
IndexOf( G4String nam)
{
    unsigned i;
    G4String pname;
    for( i=0;  i<parameter.entries(); i++)
    {
        pname = parameter(i)-> GetParameterName();
        if( pname == nam ) {
            return i;
        }
    }
    paramERR = 1;
    G4cerr << "parameter name:"<<nam<<" not found."<< endl;
    return 0;
}


unsigned G4UIcommand::
IsParameter(G4String nam)
{
    G4String pname;
    for(unsigned i=0;  i<parameter.entries(); i++)  
    {
        pname = parameter(i)-> GetParameterName();
        if( pname == nam ) return 1;
    }
    return 0;
}


// --------------------- utility functions --------------------------

tokenNum G4UIcommand::
Yylex()         // reads input and returns token number, KR486
{               // (returns EOF)
    int c;             
    G4String buf;

    while(( c= G4UIpGetc())==' '|| c=='\t' || c== '\n' )
        ;
    if (c== EOF)
        return (tokenNum)EOF;            // KR488 
    buf= "";
    if (isdigit(c) || c== '.') {         // I or D
        do {
             buf += G4String((unsigned char)c);
             c=G4UIpGetc();
         }  while (c=='.' || isdigit(c) || 
                   c=='e' || c=='E' || c=='+' || c=='-');
         G4UIpUngetc(c);
         const char* t = buf;
	 istrstream is((char*)t);
         if ( IsInt(buf.data(),20) ) {
             is >> yylval.I;
             return  CONSTINT;
         } else 
         if ( IsDouble(buf.data()) ) {
             is >> yylval.D;
             return  CONSTDOUBLE;
         } else {
             G4cerr << buf<<": numeric format error."<<endl;
         }
    }
    buf="";
    if (isalpha(c)|| c=='_') {           // IDENTIFIER
        do {
            buf += G4String((unsigned char)c); 
        } while ((c=G4UIpGetc()) != EOF && (isalnum(c) || c=='_'));
        G4UIpUngetc(c);
        if( IsParameter(buf) ) {
            yylval.S =buf;
            return IDENTIFIER;
        } else {
            G4cerr << buf << " is not a parameter name."<< endl;
            paramERR = 1;
       }
    }
    switch (c) {
      case '>':   return (tokenNum) Follow('=', GE,        GT);
      case '<':   return (tokenNum) Follow('=', LE,        LT);
      case '=':   return (tokenNum) Follow('=', EQ,        '=');
      case '!':   return (tokenNum) Follow('=', NE,        '!');
      case '|':   return (tokenNum) Follow('|', LOGICALOR, '|');
      case '&':   return (tokenNum) Follow('&', LOGICALAND, '&');
      default:
          return (tokenNum) c;
    }
}


int G4UIcommand::
Follow(int expect, int ifyes, int ifno)
{
    int c = G4UIpGetc();
    if ( c== expect)
          return ifyes;
    G4UIpUngetc(c);
    return ifno;
}

//------------------ low level routines -----------------------------
int G4UIcommand::
G4UIpGetc() {                        // emulation of getc() 
    int length = rangeString.length();
    if( bp < length)
        return  rangeString(bp++);
    else 
        return EOF;
}
int G4UIcommand::
G4UIpUngetc(int c) {                 // emulation of ungetc() 
    if (c<0) return -1;
    if (bp >0 && c == rangeString(bp-1)) {
         --bp;
    } else {
         G4cerr << "G4UIpUngetc() failed." << endl;
         G4cerr << "bp="<<bp <<" c="<<c
              << " pR(bp-1)=" << rangeString(bp-1)
              << endl;
         paramERR = 1;
         return -1;
    }
    return 0;
}
