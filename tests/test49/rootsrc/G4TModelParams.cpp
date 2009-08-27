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
//
// --------------------------------------------------------------------
// Class implementation
// --------------------------------------------------------------------

#include "G4TModelParams.h"

ClassImp(G4TModelParams)

using namespace std;
using namespace ROOT;
using namespace TMath;



//______________________________________________________________________________
void G4TModelParams::Load(const TString& pf)
{
	ifstream 		asciiFile;
	string   		asciiLine;
	ParamsData_t	result;

	// open the file
	asciiFile.open(pf.Data());
	if (!asciiFile) {
		cout << "Error opening file "<< pf.Data() << endl;
		return;
	}

	// process
	getline(asciiFile, asciiLine); //STD Call

	string buf; // Have a buffer string
	stringstream ss(asciiLine); // Insert the string into a stream

	vector<string> pFields;
	Tokenize(asciiLine, pFields, " ");

	result.temp 	= atof(pFields[0].data());
	result.ssse 	= atof(pFields[1].data());
	result.eepr 	= atof(pFields[2].data());
	result.nop 		= atof(pFields[3].data());
	result.momb 	= atof(pFields[4].data());
	result.enb 		= atof(pFields[5].data());
	result.pdgpr 	= atoi(pFields[6].data());
	result.pdgtg 	= atoi(pFields[7].data());
	result.nevnt 	= atoi(pFields[8].data());
	result.freeN 	= atof(pFields[9].data());
	result.freeD 	= atof(pFields[10].data());
	result.clustP 	= atof(pFields[11].data());
	result.rMed 	= atof(pFields[12].data());
	result.solA 	= atof(pFields[13].data());


	asciiFile.close();
	SetData(result);
}

//______________________________________________________________________________
void G4TModelParams::Save(const TString& pf )
{
	ofstream out(pf.Data(), ios::out);

	if(!out) {
		cout << "Cannot open output chips file.\n";
		return;
	}

	stringstream stream;
	stream << fData.temp 	<< ' ';
	stream << fData.ssse 	<< ' ';
	stream << fData.eepr 	<< ' ';
	stream << fData.nop 	<< ' ';
	stream << fData.momb 	<< ' ';
	stream << fData.enb 	<< ' ';
	stream << fData.pdgpr 	<< ' ';
	stream << fData.pdgtg 	<< ' ';
	stream << fData.nevnt 	<< ' ';
	stream << fData.freeN 	<< ' ';
	stream << fData.freeD 	<< ' ';
	stream << fData.clustP  << ' ';
	stream << fData.rMed 	<< ' ';
	stream << fData.solA 	<< ' ';
	out.write(stream.str().data(), stream.str().size());
	out.close();
}

//______________________________________________________________________________
inline void G4TModelParams::Tokenize(const string& str, vector<string>& tokens, const string& delimiters)
{
	  string::size_type lastPos = str.find_first_not_of(delimiters, 0);
	  string::size_type pos     = str.find_first_of(delimiters, lastPos);

	  while (string::npos != pos || string::npos != lastPos)
	  {
		  tokens.push_back(str.substr(lastPos, pos - lastPos));
		  lastPos = str.find_first_not_of(delimiters, pos);
		  pos = str.find_first_of(delimiters, lastPos);
	  }
}




