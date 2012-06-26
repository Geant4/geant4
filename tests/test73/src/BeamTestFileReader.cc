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
#include <sstream>
#include "BeamTestFileReader.hh"

FileReader::FileReader(const char* filename)
	: file(filename), line(), failed(false)
{}

FileReader::FileReader(const std::string& filename)
	: file(filename.c_str()), line(), failed(false)
{}	

bool FileReader::isValid() const
{
	if ( ! file )
		return false;
	else
		return true;
}

bool FileReader::nextLine()
{
	getline(file, line);
	if ( file.eof() )
		return false;
	else return true;
}

void FileReader::skip_fields(std::istringstream& ist, const int n)
{
	if ( n < 1 )
		return;
	std::string tmp;
	for(int i = 1; i <= n; ++i)
	{
		ist >> tmp;
	}
}

int FileReader::getFieldAsInt(const int n) 
{	
	failed = false;
	std::istringstream ist(line);
	this->skip_fields(ist, n-1);
	int rval;
	ist >> rval;
	if ( ist.fail() )
	{
		failed = true;
		return 0;
	}
	else
		return rval;
}

float FileReader::getFieldAsFloat(const int n)
{
	failed = false;
	std::istringstream ist(line);
	this->skip_fields(ist, n-1);
	float rval;
	ist >> rval;
	if ( ist.fail() )
	{
		failed = true;
		return 0.0;
	}
	else
		return rval;
}

double FileReader::getFieldAsDouble(const int n)
{
	failed = false;
	std::istringstream ist(line);
	this->skip_fields(ist, n-1);
	double rval;
	ist >> rval;
	if ( ist.fail() )
	{
		failed = true;
		return 0.0;
	}
	else
		return rval;
}

std::string FileReader::getFieldAsString(const int n)
{
	failed = false;
	std::istringstream ist(line);
	this->skip_fields(ist, n-1);
	std::string rval;
	ist >> rval;
	if ( ist.fail() )
	{
		failed = true;
		return std::string();
	}
	else
		return rval;
}

bool FileReader::inputFailed() const
{
	return failed;
}
