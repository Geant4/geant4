// $Id: MyApplication.hh,v 1.1 2006-02-27 10:05:24 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   MyApplication.hh
//
//                                         2005 Q
// ====================================================================
#ifndef MY_APPLICATION_H
#define MY_APPLICATION_H

#include "globals.hh"

// ====================================================================
//
// class definition
//
// ====================================================================
class MyApplication {
protected:
  G4String version;
  G4String batchname;

public:
  MyApplication();
  MyApplication(int argc, char** argv, const G4String& ver="");
  ~MyApplication();

  void SetVersion(const G4String& aname);
  G4String GetVersion() const;

  void ParseArguments(int argc, char** argv);
  void ShowHelp() const;

  void StartSession();
  void ExecuteBatch();

  void Configure();
};

// ====================================================================
//   inline functions
// ====================================================================
inline void MyApplication::SetVersion(const G4String& aname)
{  version= aname; }

inline G4String MyApplication::GetVersion() const
{ return version; }

#endif
