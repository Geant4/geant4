// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIcommand.hh,v 1.2 1999-05-07 10:50:42 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#ifndef G4UIcommand_h
#define G4UIcommand_h 1

#include "G4UIparameter.hh"
class G4UImessenger;
#include "globals.hh"
#include "G4ApplicationState.hh"
#include <rw/tpordvec.h>
#include <rw/tvordvec.h>
#include "G4UItokenNum.hh"

class G4UIcommand 
{
  public:
      G4UIcommand();
      G4UIcommand(const char * theCommandPath, G4UImessenger * theMessenger);
      virtual ~G4UIcommand();

      int operator==(const G4UIcommand &right) const;
      int operator!=(const G4UIcommand &right) const;

      G4int DoIt(G4String parameterList);
      G4String GetCurrentValue();
      void AvailableForStates(G4ApplicationState s1);
      void AvailableForStates(G4ApplicationState s1,G4ApplicationState s2);
      void AvailableForStates(G4ApplicationState s1,G4ApplicationState s2,
                              G4ApplicationState s3);
      void AvailableForStates(G4ApplicationState s1,G4ApplicationState s2,
                              G4ApplicationState s3,G4ApplicationState s4);
      void AvailableForStates(G4ApplicationState s1,G4ApplicationState s2,
                              G4ApplicationState s3,G4ApplicationState s4,
                              G4ApplicationState s5);
      G4bool IsAvailable();
      virtual void List();

  public:
      static G4double ValueOf(G4String unitName);
      static G4String CategoryOf(G4String unitName);
      static G4String UnitsList(G4String unitCategory);

  private:  
      void G4UIcommandCommonConstructorCode (const char * theCommandPath);
      G4UImessenger *messenger;
      G4String commandPath;
      G4String commandName;
      G4String rangeString;
      RWTPtrOrderedVector<G4UIparameter> parameter;
      RWTValOrderedVector<G4String> commandGuidance;
      RWTValOrderedVector<G4ApplicationState> availabelStateList;

  public:
      inline void SetRange(const char* rs)
      { rangeString = rs; }
      inline const G4String GetRange() const
      { return rangeString; };
      inline G4int GetGuidanceEntries() const
      { return commandGuidance.entries(); }
      inline const G4String GetGuidanceLine(int i) const
      { return commandGuidance[i]; }
      inline const G4String GetCommandPath() const
      { return commandPath; }
      inline const G4String GetCommandName() const
      { return commandName; }
      inline G4int GetParameterEntries() const
      { return parameter.entries(); }
      inline G4UIparameter * GetParameter(int i) const
      { return parameter[i]; }
      inline void SetParameter(G4UIparameter *const newParameter)
      {
	parameter.insert( newParameter );
	newVal.resize( parameter.entries() );
      }
      inline void SetGuidance(const char * aGuidance)
      { 
        G4String * theGuidance = new G4String( aGuidance );
        commandGuidance.insert( *theGuidance ); 
      }
      inline const G4String GetTitle() const
      {
	    if(commandGuidance.entries() == 0)
	    { return G4String("...Title not available..."); }
 	    else
	    { return commandGuidance(0); }
      }

  protected:
    int CheckNewValue(G4String newValue);

    // --- the following is used by CheckNewValue() --------
  private:
    int TypeCheck(G4String newValue);
    int RangeCheck(G4String newValue);
    int IsInt(const char* str, short maxLength);
    int IsDouble(const char* str);
    int ExpectExponent(const char* str);
    //  syntax nodes
    yystype Expression( void );
    yystype LogicalORExpression( void );
    yystype LogicalANDExpression( void );
    yystype EqualityExpression ( void );
    yystype RelationalExpression( void );
    yystype AdditiveExpression( void );
    yystype MultiplicativeExpression( void );
    yystype UnaryExpression( void );
    yystype PrimaryExpression( void );
    //  semantics routines
    int Eval2( yystype arg1, int op, yystype arg2 );
    int CompareInt( int arg1, int op, int arg2);
    int CompareDouble( double arg1, int op, double arg2);
    //  utility 
    tokenNum Yylex( void );      // returns next token
    unsigned IndexOf( G4String ); // returns the index of the var name
    unsigned IsParameter( G4String ); // returns 1 or 0
    int G4UIpGetc( void );      // read one char from rangeBuf
    int G4UIpUngetc( int c );   // put back  
    int Backslash( int c );
    int Follow( int expect, int ifyes, int ifno );
    G4String TokenToStr(int token);
    void PrintToken(void);      // for debug
    //  data
    G4String rangeBuf;
    int bp;                      // buffer pointer for rangeBuf
    tokenNum token;
    yystype yylval;
    RWTValOrderedVector<yystype>  newVal;
    int paramERR;
};

#endif

