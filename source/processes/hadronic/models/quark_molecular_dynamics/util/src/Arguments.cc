#include <string.h>
#include "Error.hh"
#include "Conversions.hh"
#include "Arguments.hh"

Arguments& operator>>(Arguments& list,ArgumentEntryBase& x)
{
  if ( x.hasKey() && !list.args.isEmpty() ) {
    unsigned int len = strlen( x.getKey() );
    DList<ArgumentEntryBase>::Iterator k(list.args);
    k.rewind();
    ArgumentEntryBase* s = k.current();
    for (k.rewind(); k.more() && strlen( k.current()->getKey() ) > len; k.advance() );
    if (k.more() )
      list.args.insertBefore(k.current(),&x);
    else
      list.args.appendLast(&x);
  }
  else
    list.args.appendLast(&x);
  return list;
}

ArgumentEntryBase::ArgumentEntryBase(char* s,int n) 
  : N(n)
{
    name = s ? strdup(s) : 0;
}

ArgumentEntryBase::~ArgumentEntryBase()
{
  if (name)
    delete name;
}

Arguments::Arguments(int n,char** v) : argc(n),argv(new char*[n])
{
  for (int i=0; i<n; i++) {
    argv[i] = strdup(v[i]);
  }
}

Arguments::~Arguments()
{
  args.removeAll();
  for (int i=0; i<argc; i++)
    if ( argv[i] )
      delete argv[i];
  delete [] argv;
}

void Arguments::WrongOptionDeclarator::writeMessage(ostream& o) const
{
  o << "Illegal Option declarator '" << ch << "'.\n"; 
}

Arguments::IllegalOption::IllegalOption(char* c,char* o)
{
  command = strdup(c);
  option = strdup(o);
}

void Arguments::IllegalOption::writeMessage(ostream& o) const
{
  o << command << ": illegal option -- " << option << endl;
}

void Arguments::Apply()
{
    DList< ArgumentEntryBase >::Iterator k(args);
    ArgumentEntryBase* option = 0;
    char* arg = 0;
    int ArgumentsDetected = 1;  // Remember: program name counts as one argum.
    for (k.rewind(); k.more() && k.current()->hasKey(); k.advance() ) {
       int found = 0;
       option = k.current();
       for (int i=1; i<argc; i++) 
         if (argv[i]) {
	   arg = argv[i];
	   int len = strlen(option->getKey());
	   if ( !strncmp(option->getKey(),arg,len) ) {
	     arg += len;
	     found = i;
	     break;
	   }
	 }
       if (found) {
	 if ( option->getNumber() ) {
	   int n = 0;
	   if ( strlen(arg) > 0 ) {
	     option->convertChar(arg);
	     ++n;
	   }
	   for (int i=n; i<option->getNumber(); i++)
	     {
	       option->convertChar(argv[found+i+1-n],i);
	       delete argv[found+i+1-n];
	       argv[found+i+1-n] = 0;
	       ++ArgumentsDetected;
	     }
	   delete argv[found];
	   argv[found] = 0;
	   ++ArgumentsDetected;
	 }
	 else 
	   if ( strlen(arg) == 0 ) {
	     option->convertChar(0);
	     delete [] argv[found];
	     argv[found] = 0;
	     ++ArgumentsDetected;
	   }
	   else 
	     throw IllegalOption(argv[0],argv[found]);
       }
    }
    for (;k.more();k.advance() ) {
      option = k.current();
      for (int i=1; i<argc; i++) 
	if (argv[i]) {
	  option->convertChar(argv[i]);
	  delete [] argv[i];
	  argv[i] = 0;
	  ++ArgumentsDetected;
	  break;
	}
    }
    if (ArgumentsDetected != argc) { // illegal option
      int i;
      for (i=1; i<argc && !argv[i]; i++);
      throw IllegalOption(argv[0],argv[i]);
    }
}


