// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: read_func.h,v 1.1 1999-01-07 16:08:09 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef READ_FUNC_H
#define	READ_FUNC_H

#include <sdai.h>

#define MAX_COMMENT_LENGTH 512

// Print Error information for debugging purposes
extern void PrintErrorState(ErrorDescriptor &err);

// Print istream error information for debugging purposes
extern void IStreamState(istream &in);

extern int 
ReadInteger(SdaiInteger &val, istream &in, ErrorDescriptor *err, 
	    char *tokenList);

extern int 
ReadInteger(SdaiInteger &val, const char *s, ErrorDescriptor *err, 
	    char *tokenList);

extern Severity 
IntValidLevel (const char *attrValue, ErrorDescriptor *err,
	       int clearError, int optional, char *tokenList);

extern int
ReadReal(SdaiReal &val, istream &in, ErrorDescriptor *err, 
	 char *tokenList);

extern int
ReadReal(SdaiReal &val, const char *s, ErrorDescriptor *err, 
	 char *tokenList);

extern Severity 
RealValidLevel (const char *attrValue, ErrorDescriptor *err,
		int clearError, int optional, char *tokenList);

extern int
ReadNumber(SdaiReal &val, istream &in, ErrorDescriptor *err, 
	   char *tokenList);

extern int
ReadNumber(SdaiReal &val, const char *s, ErrorDescriptor *err, 
	   char *tokenList);

extern Severity 
NumberValidLevel (const char *attrValue, ErrorDescriptor *err,
		  int clearError, int optional, char *tokenList);


////////////////////

extern int   QuoteInString(istream& in);

extern void 
PushPastString (istream& in, SCLstring &s, ErrorDescriptor *err);

extern void 
PushPastImbedAggr (istream& in, SCLstring &s, ErrorDescriptor *err);

extern void 
PushPastAggr1Dim(istream& in, SCLstring &s, ErrorDescriptor *err);

////////////////////

extern Severity 
FindStartOfInstance(istream& in, SCLstring&  inst);

	//  used for instances that aren\'t valid - reads to next \';\'
extern Severity 
SkipInstance (istream& in, SCLstring & inst);

extern const char *
SkipSimpleRecord(istream &in, SCLstring &buf, ErrorDescriptor *err);

 // this includes entity names
extern const char *
ReadStdKeyword(istream& in, SCLstring &buf, int skipInitWS = 1);

extern const char* 
GetKeyword(istream& in, const char* delims, ErrorDescriptor &err);

extern int 
FoundEndSecKywd(istream& in, ErrorDescriptor &err);

extern const char *ReadComment(SCLstring &ss, const char *s);

extern const char *ReadComment(istream& in, SCLstring &s);

extern Severity    ReadPcd(istream& in);   //Print control directive

extern void        ReadTokenSeparator(istream& in, SCLstring *comments = 0);

#endif
