// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIcommand.hh,v 1.6 2001-02-08 06:07:18 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#ifndef G4UIcommand_h
#define G4UIcommand_h 1

#include "G4UIparameter.hh"
class G4UImessenger;
#include "globals.hh"
#include "G4ApplicationState.hh"
//#include "g4rw/tpordvec.h"
//#include "g4rw/tvordvec.h"
#include "g4std/vector"
#include "G4UItokenNum.hh"

// class description:
//  
//  This G4UIcommand is the "concrete" base class which represents a command
// used by Geant4 (G)UI. The user can use this class in case the parameter
// arguments of a command are not suitable with respect to the derived command
// classes.
//  Some methods defined in this base class are used by the derived classes.
//

class G4UIcommand 
{
  public: 
      G4UIcommand();
  public: // with description
      G4UIcommand(const char * theCommandPath, G4UImessenger * theMessenger);
      //  Constructor. The command string with full path directory
      // and the pointer to the messenger must be given.
  public: 
      virtual ~G4UIcommand();

      int operator==(const G4UIcommand &right) const;
      int operator!=(const G4UIcommand &right) const;

      G4int DoIt(G4String parameterList);
      G4String GetCurrentValue();
  public: // with description
      void AvailableForStates(G4ApplicationState s1);
      void AvailableForStates(G4ApplicationState s1,G4ApplicationState s2);
      void AvailableForStates(G4ApplicationState s1,G4ApplicationState s2,
                              G4ApplicationState s3);
      void AvailableForStates(G4ApplicationState s1,G4ApplicationState s2,
                              G4ApplicationState s3,G4ApplicationState s4);
      void AvailableForStates(G4ApplicationState s1,G4ApplicationState s2,
                              G4ApplicationState s3,G4ApplicationState s4,
                              G4ApplicationState s5);
      //  These methods define the states where the command is available.
      // Once one of these commands is invoked, the command application will
      // be denied when Geant4 is NOT in the assigned states.
  public:
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
      G4std::vector<G4UIparameter*> parameter;
      G4std::vector<G4String> commandGuidance;
      G4std::vector<G4ApplicationState> availabelStateList;

  public: // with description
      inline void SetRange(const char* rs)
      { rangeString = rs; }
      //  Defines the range the command parameter(s) can take.
      //  The variable name(s) appear in the range expression must be same
      // as the name(s) of the parameter(s).
      //  All the C++ syntax of relational operators are allowed for the
      // range expression.
  public:
      inline const G4String GetRange() const
      { return rangeString; };
      inline G4int GetGuidanceEntries() const
      { return commandGuidance.size(); }
      inline const G4String GetGuidanceLine(int i) const
      { return commandGuidance[i]; }
      inline const G4String GetCommandPath() const
      { return commandPath; }
      inline const G4String GetCommandName() const
      { return commandName; }
      inline G4int GetParameterEntries() const
      { return parameter.size(); }
      inline G4UIparameter * GetParameter(int i) const
      { return parameter[i]; }
  public: // with description
      inline void SetParameter(G4UIparameter *const newParameter)
      {
	parameter.push_back( newParameter );
	newVal.resize( parameter.size() );
      }
      //  Defines a parameter. This method is used by the derived command classes
      // but the user can directly use this command when he/she defines a command
      // by hem(her)self without using the derived class. For this case, the order
      // of the parameters is the order of invoking this method.
      inline void SetGuidance(const char * aGuidance)
      { 
        G4String * theGuidance = new G4String( aGuidance );
        commandGuidance.push_back( *theGuidance ); 
      }
      //  Adds a guidance line. Unlimitted number of invokation of this method is
      // allowed. The given lines of guidance will appear for the help. The first
      // line of the guidance will be used as the title of the command, i.e. one
      // line list of the commands.
  public:
      inline const G4String GetTitle() const
      {
	    if(commandGuidance.size() == 0)
	    { return G4String("...Title not available..."); }
 	    else
	    { return commandGuidance[0]; }
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
    G4std::vector<yystype>  newVal;
    int paramERR;
};

#endif

