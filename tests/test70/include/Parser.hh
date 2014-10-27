#ifndef PARSER_HH
#define PARSER_HH

#include "globals.hh"
#include <map>

enum CommandType
{
    WithOption,
    WithoutOption,
    OptionNotCompulsory
};

struct Command
{
    friend class CommandLineParser;
    CommandType fType;
    G4String fOption;
    G4bool fActive;
    G4String fDescription;
    G4String fOptionName;

    Command(CommandType, const G4String &description = "", const G4String &optionName = "optionName");
    ~Command(){;}

public :
    const G4String& GetOption()     {return fOption;}
    CommandType GetType()           {return fType;}
    G4bool IsActive()               {return fActive;}
    const G4String& GetDescription()               {return fDescription;}
    const G4String& GetOptionName()               {return fOptionName;}
};

class CommandLineParser
{
    static CommandLineParser* fpInstance;
    std::map<G4String, Command*> fCommandMap;
    G4bool fOptionsWereSetup;
    G4int fMaxMarkerLength;
    G4int fMaxOptionNameLength;
    G4int fVerbose;

public:
    static CommandLineParser* GetParser();
    CommandLineParser();
    ~CommandLineParser();
    static void DeleteInstance();
    int Parse(int& argc, char **argv);
    void PrintHelp();
    void CorrectRemainingOptions(int& argc, char **argv);
    void AddCommand(const G4String & marker,CommandType, const G4String& description = "", const G4String& optionName = "");
    Command* FindCommand(const G4String &marker);
    Command* GetCommandIfActive(const G4String &marker);
    G4bool WereOptionsSetup(){return fOptionsWereSetup;}
};

#endif // PARSER_HH
