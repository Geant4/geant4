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
// G4UIcommand
//
// Class description:
//
// This G4UIcommand is the "concrete" base class which represents a command
// used by Geant4 (G)UI. The user can use this class in case the parameter
// arguments of a command are not suitable with respect to the derived command
// classes.
// Some methods defined in this base class are used by the derived classes

// Author: Makoto Asai (SLAC), 1998
// --------------------------------------------------------------------
#ifndef G4UIcommand_hh
#define G4UIcommand_hh 1

#include <vector>

#include "G4UIparameter.hh"
#include "globals.hh"
#include "G4ApplicationState.hh"
#include "G4UItokenNum.hh"
#include "G4ThreeVector.hh"

class G4UImessenger;

class G4UIcommand
{
  public:

    G4UIcommand() = default;
      // Dummy default constructor

    G4UIcommand(const char* theCommandPath, G4UImessenger* theMessenger,
                G4bool tBB = true);
      // Constructor. The command string with full path directory
      // and the pointer to the messenger must be given.
      // If tBB is set to false, this command won't be sent to worker threads.
      // This tBB parameter could be changed with SetToBeBroadcasted() method
      // except for G4UIdirectory

    virtual ~G4UIcommand();

    G4bool operator==(const G4UIcommand& right) const;
    G4bool operator!=(const G4UIcommand& right) const;

    virtual G4int DoIt(G4String parameterList);

    G4String GetCurrentValue();

    void AvailableForStates(G4ApplicationState s1);
    void AvailableForStates(G4ApplicationState s1, G4ApplicationState s2);
    void AvailableForStates(G4ApplicationState s1, G4ApplicationState s2,
                            G4ApplicationState s3);
    void AvailableForStates(G4ApplicationState s1, G4ApplicationState s2,
                            G4ApplicationState s3, G4ApplicationState s4);
    void AvailableForStates(G4ApplicationState s1, G4ApplicationState s2,
                            G4ApplicationState s3, G4ApplicationState s4,
                            G4ApplicationState s5);
      // These methods define the states where the command is available.
      // Once one of these commands is invoked, the command application will
      // be denied when Geant4 is NOT in the assigned states

    G4bool IsAvailable();

    virtual void List();

    static G4String ConvertToString(G4bool boolVal);
    static G4String ConvertToString(G4int intValue);
    static G4String ConvertToString(G4long longValue);
    static G4String ConvertToString(G4double doubleValue);
    static G4String ConvertToString(G4double doubleValue, const char* unitName);
    static G4String ConvertToString(const G4ThreeVector& vec);
    static G4String ConvertToString(const G4ThreeVector& vec,
                                    const char* unitName);
    // Static methods for conversion from value(s) to a string.
    // These methods are to be used by GetCurrentValues() methods
    // of concrete messengers

    static G4bool ConvertToBool(const char* st);
    static G4int ConvertToInt(const char* st);
    static G4long ConvertToLongInt(const char* st);
    static G4double ConvertToDouble(const char* st);
    static G4double ConvertToDimensionedDouble(const char* st);
    static G4ThreeVector ConvertTo3Vector(const char* st);
    static G4ThreeVector ConvertToDimensioned3Vector(const char* st);
      // Static methods for conversion from a string to a value of the returning
      // type. These methods are to be used directly by SetNewValues() methods
      // of concrete messengers, or GetNewXXXValue() of classes derived from
      // this G4UIcommand class

    static G4double ValueOf(const char* unitName);
    static G4String CategoryOf(const char* unitName);
    static G4String UnitsList(const char* unitCategory);
      // Static methods for unit and its category

    inline void SetRange(const char* rs) { rangeString = rs; }
      // Defines the range the command parameter(s) can take.
      // The variable name(s) appear in the range expression must be the same
      // as the name(s) of the parameter(s).
      // All the C++ syntax of relational operators are allowed for the
      // range expression

    inline const G4String& GetRange() const { return rangeString; }
    inline std::size_t GetGuidanceEntries() const
    {
      return commandGuidance.size();
    }
    inline const G4String& GetGuidanceLine(G4int i) const
    {
      return commandGuidance[i];
    }
    inline const G4String& GetCommandPath() const { return commandPath; }
    inline const G4String& GetCommandName() const { return commandName; }
    inline std::size_t GetParameterEntries() const { return parameter.size(); }
    inline G4UIparameter* GetParameter(G4int i) const { return parameter[i]; }
    inline std::vector<G4ApplicationState>* GetStateList()
    {
      return &availabelStateList;
    }
    inline G4UImessenger* GetMessenger() const { return messenger; }

    inline void SetParameter(G4UIparameter* const newParameter)
      // Defines a parameter. This method is used by the derived command
      // classes but the user can directly use this command when defining
      // a command, without using the derived class. For this case, the order
      // of the parameters is the order of invoking this method
    {
      parameter.push_back(newParameter);
      newVal.resize(parameter.size());
    }

    inline void SetGuidance(const char* aGuidance)
      // Adds a guidance line. Unlimited times of invokation of this method is
      // allowed. The given lines of guidance will appear for the help.
      // The first line of the guidance will be used as the title of the
      // command, i.e. one line list of the commands
    {
      commandGuidance.emplace_back(aGuidance);
    }

    inline const G4String GetTitle() const
    {
      return (commandGuidance.empty()) ? G4String("...Title not available...")
                                       : commandGuidance[0];
    }

    inline void SetToBeBroadcasted(G4bool val) { toBeBroadcasted = val; }
    inline G4bool ToBeBroadcasted() const { return toBeBroadcasted; }
    inline void SetToBeFlushed(G4bool val) { toBeFlushed = val; }
    inline G4bool ToBeFlushed() const { return toBeFlushed; }
    inline void SetWorkerThreadOnly(G4bool val = true) { workerThreadOnly=val; }
    inline G4bool IsWorkerThreadOnly() const { return workerThreadOnly; }

    inline void CommandFailed(G4int errCode, G4ExceptionDescription& ed)
    {
      commandFailureCode = errCode;
      failureDescription = ed.str();
    }
    inline void CommandFailed(G4ExceptionDescription& ed)
    {
      commandFailureCode = 1;
      failureDescription = ed.str();
    }
    inline G4int IfCommandFailed() { return commandFailureCode; }
    inline const G4String& GetFailureDescription() {return failureDescription;}
    inline void ResetFailure()
    {
      commandFailureCode = 0;
      failureDescription = "";
    }

  public:
    enum CommandType 
    { BaseClassCmd, WithoutParameterCmd,
      WithABoolCmd, WithAnIntegerCmd, WithALongIntCmd, 
      WithADoubleCmd, WithADoubleAndUnitCmd, With3VectorCmd, With3VectorAndUnitCmd,
      WithAStringCmd, CmdDirectory = -1 };
      
    inline CommandType GetCommandType() const
    { return commandType; }
    void SetCommandType(CommandType);

    inline void SetDefaultSortFlag(G4bool val)
    { ifSort = val; }

  protected:

    // --- the following is used by CheckNewValue() --------
    using yystype  = G4UItokenNum::yystype;
    using tokenNum = G4UItokenNum::tokenNum;

    G4int CheckNewValue(const char* newValue);

    G4bool toBeBroadcasted = false;
    G4bool toBeFlushed = false;
    G4bool workerThreadOnly = false;

    G4int commandFailureCode = 0;
    G4String failureDescription = "";

    G4bool ifSort = false;

  private:

    void G4UIcommandCommonConstructorCode(const char* theCommandPath);

    G4int TypeCheck(const char* t);
    G4int RangeCheck(const char* t);
    G4int IsInt(const char* str, short maxLength); // used for both int and long int
    G4int IsDouble(const char* str);
    G4int ExpectExponent(const char* str);
    //  syntax nodes
    yystype Expression();
    yystype LogicalORExpression();
    yystype LogicalANDExpression();
    yystype EqualityExpression();
    yystype RelationalExpression();
    yystype AdditiveExpression();
    yystype MultiplicativeExpression();
    yystype UnaryExpression();
    yystype PrimaryExpression();
    //  semantics routines
    G4int Eval2(const yystype& arg1, G4int op, const yystype& arg2);
    G4int CompareInt(G4int arg1, G4int op, G4int arg2);
    G4int CompareLong(G4long arg1, G4int op, G4long arg2);
    G4int CompareDouble(G4double arg1, G4int op, G4double arg2);
    //  utility
    tokenNum Yylex();                   // returns next token
    unsigned IndexOf(const char*);      // returns the index of the var name
    unsigned IsParameter(const char*);  // returns 1 or 0
    G4int G4UIpGetc();                  // read one char from rangeBuf
    G4int G4UIpUngetc(G4int c);         // put back
    G4int Backslash(G4int c);
    G4int Follow(G4int expect, G4int ifyes, G4int ifno);
    //G4String TokenToStr(G4int token);
    //void PrintToken(void);  // for debug

    // Data -----------------------------------------------------------

  private:
    CommandType commandType = BaseClassCmd;
    G4UImessenger* messenger = nullptr;
    G4String commandPath;
    G4String commandName;
    G4String rangeString;
    std::vector<G4UIparameter*> parameter;
    std::vector<G4String> commandGuidance;
    std::vector<G4ApplicationState> availabelStateList;

    G4String rangeBuf;
    G4int bp = 0;  // buffer pointer for rangeBuf
    tokenNum token = G4UItokenNum::IDENTIFIER;
    yystype yylval;
    std::vector<yystype> newVal;
    G4int paramERR = 0;
};

#endif
