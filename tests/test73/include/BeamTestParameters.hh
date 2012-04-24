// $Id: Parameters.h,v 1.14 2010/09/29 07:59:26 paula Exp $
#ifndef BEAMTESTPARAMETERS_H 
#define BEAMTESTPARAMETERS_H 1

// Include files
#include<iostream>
#include<ctype.h>
#include<string>
#include<list>
#include<map>

/** @class Parameters Parameters.h
 *
 */

class Parameters {
	public: 
		/// Standard constructor
		Parameters( ); 
		int readCommandLine(int,char**);
		//int readConditions();
		bool exists(std::string);
		void help();

		virtual ~Parameters( ); ///< Destructor
		
		double pTransverse;
		double momentum;
		double zThickness;
		std::string filename;
        int numberOfChambers;
        double chamberSpacing;
};
#endif // BEAMTESTPARAMETERS_H
