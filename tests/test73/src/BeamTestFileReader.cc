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
