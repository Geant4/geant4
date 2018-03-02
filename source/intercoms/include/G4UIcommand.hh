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
// $Id: G4UIcommand.hh 108488 2018-02-15 14:52:00Z gcosmo $
//
//

#ifndef G4UIcommand_h
#define G4UIcommand_h 1

#include "G4UIparameter.hh"
class G4UImessenger;
#include "globals.hh"
#include "G4ApplicationState.hh"
#include <vector>
#include "G4UItokenNum.hh"
#include "G4ThreeVector.hh"

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
      G4UIcommand(const char * theCommandPath, G4UImessenger * theMessenger,
                  G4bool tBB = true);
      //  Constructor. The command string with full path directory
      // and the pointer to the messenger must be given.
      // If tBB is set to false, this command won't be sent to worker threads.
      // This tBB parameter could be changed with SetToBeBroadcasted() method
      // except for G4UIdirectory.

  public:
      virtual ~G4UIcommand();

      G4int operator==(const G4UIcommand &right) const;
      G4int operator!=(const G4UIcommand &right) const;

      virtual G4int DoIt(G4String parameterList);
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

  public: // with description
      static G4String ConvertToString(G4bool boolVal);
      static G4String ConvertToString(G4int intValue);
      static G4String ConvertToString(G4double doubleValue);
      static G4String ConvertToString(G4double doubleValue,const char* unitName);
      static G4String ConvertToString(G4ThreeVector vec);
      static G4String ConvertToString(G4ThreeVector vec,const char* unitName);
      // Static methods for conversion from value(s) to a string. These methods are to be
      // used by GetCurrentValues() methods of concrete messengers.

      static G4bool ConvertToBool(const char* st);
      static G4int ConvertToInt(const char* st);
      static G4double ConvertToDouble(const char* st);
      static G4double ConvertToDimensionedDouble(const char* st);
      static G4ThreeVector ConvertTo3Vector(const char* st);
      static G4ThreeVector ConvertToDimensioned3Vector(const char* st);
      // Static methods for conversion from a string to a value of the returning type.
      // These methods are to be used directly by SetNewValues() methods of concrete
      // messengers, or GetNewXXXValue() of classes derived from this G4UIcommand class.

      static G4double ValueOf(const char* unitName);
      static G4String CategoryOf(const char* unitName);
      static G4String UnitsList(const char* unitCategory);
      // Static methods for unit and its category.

  private:
      void G4UIcommandCommonConstructorCode (const char * theCommandPath);
      G4UImessenger *messenger;
      G4String commandPath;
      G4String commandName;
      G4String rangeString;
      std::vector<G4UIparameter*> parameter;
      std::vector<G4String> commandGuidance;
      std::vector<G4ApplicationState> availabelStateList;

  public: // with description
      inline void SetRange(const char* rs)
      { rangeString = rs; }
      //  Defines the range the command parameter(s) can take.
      //  The variable name(s) appear in the range expression must be same
      // as the name(s) of the parameter(s).
      //  All the C++ syntax of relational operators are allowed for the
      // range expression.
  public:
      inline const G4String & GetRange() const
      { return rangeString; };
      inline G4int GetGuidanceEntries() const
      { return commandGuidance.size(); }
      inline const G4String & GetGuidanceLine(G4int i) const
      { return commandGuidance[i]; }
      inline const G4String & GetCommandPath() const
      { return commandPath; }
      inline const G4String & GetCommandName() const
      { return commandName; }
      inline G4int GetParameterEntries() const
      { return parameter.size(); }
      inline G4UIparameter * GetParameter(G4int i) const
      { return parameter[i]; }
      inline std::vector<G4ApplicationState>* GetStateList()
      { return &availabelStateList; }
      inline G4UImessenger * GetMessenger() const
      { return messenger; }
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
        commandGuidance.push_back( G4String( aGuidance ) );
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
    G4bool toBeBroadcasted;
    G4bool toBeFlushed;
    G4bool workerThreadOnly;

  public:
    inline void SetToBeBroadcasted(G4bool val)
    { toBeBroadcasted = val; }
    inline G4bool ToBeBroadcasted() const
    { return toBeBroadcasted; }
    inline void SetToBeFlushed(G4bool val)
    { toBeFlushed = val; }
    inline G4bool ToBeFlushed() const
    { return toBeFlushed; }
    inline void SetWorkerThreadOnly(G4bool val=true)
    { workerThreadOnly = val; }
    inline G4bool IsWorkerThreadOnly() const
    { return workerThreadOnly; }

  protected:
    G4int commandFailureCode;
    G4String failureDescription;

  public:
    inline void CommandFailed( G4int errCode, G4ExceptionDescription& ed )
    { commandFailureCode = errCode; failureDescription = ed.str(); }
    inline void CommandFailed( G4ExceptionDescription& ed )
    { commandFailureCode = 1; failureDescription = ed.str(); }
    inline G4int IfCommandFailed() 
    { return commandFailureCode; }
    inline G4String& GetFailureDescription()
    { return failureDescription; }
    inline void ResetFailure()
    { commandFailureCode = 0; failureDescription = ""; }
    
  protected:
    G4int CheckNewValue(const char* newValue);

    // --- the following is used by CheckNewValue() --------
  private:
    G4int TypeCheck(const char* t);
    G4int RangeCheck(const char* t);
    G4int IsInt(const char* str, short maxLength);
    G4int IsDouble(const char* str);
    G4int ExpectExponent(const char* str);
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
    G4int Eval2( yystype arg1, G4int op, yystype arg2 );
    G4int CompareInt( G4int arg1, G4int op, G4int arg2);
    G4int CompareDouble( G4double arg1, G4int op, G4double arg2);
    //  utility
    tokenNum Yylex( void );      // returns next token
    unsigned IndexOf( const char* ); // returns the index of the var name
    unsigned IsParameter( const char* ); // returns 1 or 0
    G4int G4UIpGetc( void );      // read one char from rangeBuf
    G4int G4UIpUngetc( G4int c );   // put back
    G4int Backslash( G4int c );
    G4int Follow( G4int expect, G4int ifyes, G4int ifno );
    G4String TokenToStr(G4int token);
    void PrintToken(void);      // for debug
    //  data
    G4String rangeBuf;
    G4int bp;                      // buffer pointer for rangeBuf
    tokenNum token;
    yystype yylval;
    std::vector<yystype>  newVal;
    G4int paramERR;
};

#endif

