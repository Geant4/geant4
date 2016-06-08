#ifndef READ_FUNC_H
#define	READ_FUNC_H

#include <sdai.h>

#define MAX_COMMENT_LENGTH 512

// print Error information for debugging purposes
extern void PrintErrorState(ErrorDescriptor &err);

// print G4std::istream error information for debugging purposes
extern void IStreamState(G4std::istream &in);

extern int 
ReadInteger(SCLP23(Integer) &val, G4std::istream &in, ErrorDescriptor *err, 
	    char *tokenList);

extern int 
ReadInteger(SCLP23(Integer) &val, const char *s, ErrorDescriptor *err, 
	    char *tokenList);

extern Severity 
IntValidLevel (const char *attrValue, ErrorDescriptor *err,
	       int clearError, int optional, char *tokenList);

extern char * 
WriteReal(SCLP23(Real) val, SCLstring &s);

extern void
WriteReal(SCLP23(Real) val, G4std::ostream &out);

extern int
ReadReal(SCLP23(Real) &val, G4std::istream &in, ErrorDescriptor *err, 
	 char *tokenList);

extern int
ReadReal(SCLP23(Real) &val, const char *s, ErrorDescriptor *err, 
	 char *tokenList);

extern Severity 
RealValidLevel (const char *attrValue, ErrorDescriptor *err,
		int clearError, int optional, char *tokenList);

extern int
ReadNumber(SCLP23(Real) &val, G4std::istream &in, ErrorDescriptor *err, 
	   char *tokenList);

extern int
ReadNumber(SCLP23(Real) &val, const char *s, ErrorDescriptor *err, 
	   char *tokenList);

extern Severity 
NumberValidLevel (const char *attrValue, ErrorDescriptor *err,
		  int clearError, int optional, char *tokenList);


////////////////////

extern int   QuoteInString(G4std::istream& in);

extern void 
PushPastString (G4std::istream& in, SCLstring &s, ErrorDescriptor *err);

extern void 
PushPastImbedAggr (G4std::istream& in, SCLstring &s, ErrorDescriptor *err);

extern void 
PushPastAggr1Dim(G4std::istream& in, SCLstring &s, ErrorDescriptor *err);

////////////////////

extern Severity 
FindStartOfInstance(G4std::istream& in, SCLstring&  inst);

	//  used for instances that aren\'t valid - reads to next \';\'
extern Severity 
SkipInstance (G4std::istream& in, SCLstring & inst);

extern const char *
SkipSimpleRecord(G4std::istream &in, SCLstring &buf, ErrorDescriptor *err);

 // this includes entity names
extern const char *
ReadStdKeyword(G4std::istream& in, SCLstring &buf, int skipInitWS = 1);

extern const char* 
GetKeyword(G4std::istream& in, const char* delims, ErrorDescriptor &err);

extern int 
FoundEndSecKywd(G4std::istream& in, ErrorDescriptor &err);

extern const char *ReadComment(SCLstring &ss, const char *s);

extern const char *ReadComment(G4std::istream& in, SCLstring &s);

extern Severity    ReadPcd(G4std::istream& in);   //print control directive

extern void        ReadTokenSeparator(G4std::istream& in, SCLstring *comments = 0);

#endif
