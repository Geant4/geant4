#ifndef __ARGUMENTS__
#define __ARGUMENTS__

#ifdef IS_GCC
#pragma interface
#endif

#include "DList.hh"
#include "Boolean.h"
#include "Fallible.hh"

class Arguments;
class ostream;

class ArgumentEntryBase
{
public:
   virtual ~ArgumentEntryBase();
   virtual int getNumber() { return N; }
   virtual char* getKey() { return name ? name : ""; }
   virtual Boolean hasKey() { return (name != 0); }
   virtual void convertChar(char*,int = 0) = 0;
protected:
   ArgumentEntryBase(char*,int);
private:
   int N;
   char* name;
};

template<class t>
class ArgumentEntry : public ArgumentEntryBase,public Fallible<t>
{
public:
  ArgumentEntry(char*,const t&,int = 1);
  ArgumentEntry(char*,const t*,int);
  ~ArgumentEntry() {}
  inline void convertChar(char*,int = 0);
  operator t() { return getValue(); }
};

class SwitchArgument : public ArgumentEntry<Boolean>
{
public:
  SwitchArgument(char* key,Boolean flag=Boolean::False)
    : ArgumentEntry<Boolean>(key,flag,0) {}
};

template<class t>
class NoKeyArgument : public ArgumentEntry<t>
{
public:
   NoKeyArgument(const t& def) : ArgumentEntry<t>(0,def) {}
   NoKeyArgument(const t* def,int n = 1) 
     : ArgumentEntry<t>(0,def,n) {}
   Boolean hasKey() { return Boolean::False; }
};

// -------------------------------------------------------------------------

class Arguments
{
friend Arguments& operator>>(Arguments&,ArgumentEntryBase&);
public:
  Arguments(int argc,char* argv[]);
  ~Arguments();
  void Apply();

  class WrongOptionDeclarator : public Error
  {
    char ch;
    public:
      WrongOptionDeclarator(char c) : ch(c) {}
      void writeMessage(ostream&) const;
  };

  class IllegalOption : public Error
  {
    char* command,* option;
  public:
    IllegalOption(char* c,char* o);
    void writeMessage(ostream&) const;
  };

private:
  int argc;
  char** argv;
  DList< ArgumentEntryBase > args;
};

#ifndef IS_GCC
#include "Arguments.tcc"
#endif

#endif
