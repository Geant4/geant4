

//



//
// $Id: dirobj.h,v 1.2 1999-05-21 20:20:38 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef dirobj_h
#define dirobj_h

/*
* NIST Utils Class Library
* clutils/dirobj.h
* May 1995
* David Sauder
* K. C. Morris

* Development of this software was funded by the United States Government,
* and is not subject to copyright.
*/

/*   */ 

#ifdef __O3DB__
#include <OpenOODB.h>
#endif


/*
 * DirObj - This object contains a List of files in a directory.
 *		Can be used for interpreting a tilde, checking paths etc.
 *	based on InterViews FBDirectory 
 *		Copyright (c) 1987, 1988, 1989 Stanford University from 
 *		/depot/iv/src/libInterViews/filebrowser.c <= notice
 *	I commented it, changed variable and function names, added a
 *	few things, and made it so we can actually use it.
 *		David Sauder
 */

//#include <InterViews/filebrowser.h>

///////////////////////////////////////////////////////////////////////////////
//
// These header files are a pain between compilers etc ... and 
// I agree they are in a mess.  You will have to CHECK them yourself
// depending on what you are using -- it works the way it is with g++
// but you need -I/depot/gnu/lib/g++-include
//
///////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#ifdef WIN32
#  include <sys/Stat.h>
#  ifndef MAXNAMLEN
#    define MAXNAMLEN 255
#  endif
#else
#  include <stat.h>
#  include <unistd.h>
#endif

#ifdef __OBJECTCENTER__
#include <sysent.h> // needed to declare getuid()
#include <sys/stat.h>
// extern        int getuid();
#endif
#ifdef __GNUC__
#include <unistd.h>  // needed to declare getuid() with new GNU 2.3.3
//extern _G_uid_t getuid _G_ARGS((void));
#include <stat.h> 
#endif
#ifdef __SUNCPLUSPLUS__
#include <sysent.h>
#include <sys/stat.h>
#endif

#ifndef WIN32
#  include <pwd.h>
#endif
//#include <InterViews/Std/os/auth.h>
//#include <InterViews/Std/os/fs.h>
	// the below is good for g++
	// if getting it from /usr/include change to __string_h
#ifndef _std_h
#include <string.h>
#endif
#ifdef WIN32
#  ifndef MAXPATHLEN
#    define MAXPATHLEN 1024
#  endif
   typedef struct __dirstream DIR;
#else
#  include <sys/param.h>
#endif
//#include <InterViews/Std/sys/dir.h>
//#include <InterViews/Std/sys/stat.h>
#include "dir.h"
#include "globals.hh"
//typedef unsigned boolean;

#ifndef nil
#define nil 0
#endif

/*****************************************************************************/

class DirObj {       
public:
    DirObj(const char* dirName);
    virtual ~DirObj();

    HepBoolean LoadDirectory(const char*);
    const char* Normalize(const char*);
    const char* ValidDirectories(const char*);

    int Index(const char*);
    const char* File(int index);
	// check for file in the currently loaded directory
    HepBoolean FileExists(const char* file) { return Index(file) ? true : false;};
    int Count();

    HepBoolean IsADirectory(const char*);
private:
    const char* Home(const char* x = nil);
    const char* ElimDot(const char*);
    const char* ElimDotDot(const char*);
    const char* InterpSlashSlash(const char*);
    const char* InterpTilde(const char*);
    const char* ExpandTilde(const char*, int);
    const char* RealPath(const char*);

    HepBoolean Reset(const char*);
    void ClearFileList();
    void CheckIndex(int index);
    void InsertFile(const char*, int index);
    void AppendFile(const char*);
    void RemoveFile(int index);
    virtual int Position(const char*);
private:
    char** fileList;
    int fileCount;
    int fileListSize;
};

//////////////////////////////// Count() //////////////////////////////////////
//
// Return the number of files in the loaded directory.
//
///////////////////////////////////////////////////////////////////////////////

inline int DirObj::Count () { return fileCount; }

//////////////////////////////// AppendFile() /////////////////////////////////
//
// Insert a new file into the fileList.
//
///////////////////////////////////////////////////////////////////////////////

inline void DirObj::AppendFile (const char* s) { InsertFile(s, fileCount); }

//////////////////////////////// File() ///////////////////////////////////////
//
// Return the file at the given index (starting at 0) in the fileList
//
///////////////////////////////////////////////////////////////////////////////

inline const char* DirObj::File (int index) { 
    return (0 <= index && index < fileCount) ? fileList[index] : nil;
}

//////////////////////////////// strdup() /////////////////////////////////////
//
// Duplicate a string.
//
///////////////////////////////////////////////////////////////////////////////

//inline char* strdup (const char* s) {
//    char* dup = new char[strlen(s) + 1];
//    strcpy(dup, s);
//    return dup;
//}

//////////////////////////////// DotSlash() ///////////////////////////////////
//
// Return 1 if the first char in 'path' is '.' and second char is
// '/' or '\0'; otherwise return 0
//
///////////////////////////////////////////////////////////////////////////////

inline HepBoolean DotSlash (const char* path) {
    return 
        path[0] != '\0' && path[0] == '.' &&
        (path[1] == '/' || path[1] == '\0')? true : false;
}

//////////////////////////////// DotDotSlash() ////////////////////////////////
//
// Return 1 if the first char in 'path' is '.', the second char is '.',
// and the third char is '/' or '\0'; otherwise return 0
//
///////////////////////////////////////////////////////////////////////////////

inline HepBoolean DotDotSlash (const char* path) {
    return 
        path[0] != '\0' && path[1] != '\0' &&
        path[0] == '.' && path[1] == '.' &&
        (path[2] == '/' || path[2] == '\0')? true : false;
}

#endif
