//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
#include "FreeHepTypes.h"

#include <iostream>
#include <string>

/**
 * A PrintWriter that keeps track of an indentation level
 * and indents the output appropriately.
 *
 * <b>Warning:</b> Only print and println methods taking strings have been overriden,
 * print, println methods taking other arguments may not be indented properly.
 *
 * @author Mark Donszelmann
 */
class IndentPrintWriter {

	public:
	    IndentPrintWriter(std::ostream* out, int level = 0);
        virtual ~IndentPrintWriter();

        void close();
        IndentPrintWriter& operator<< (const char *s);
        IndentPrintWriter& operator<< (std::ostream& (*pf)(std::ostream&));
	    void println(std::string s);
        void print(std::string s);
	    void println();
	    void indent();
	    void outdent();
	    int getIndent();
        void setIndent(int level);
        std::string getIndentString();
        void setIndentString(std::string indentString);

    private:
        void doIndent();

	std::ostream* out;
        int indentLevel;
	bool indented;
	std::string indentString;
};

#endif

