#ifndef BEAMTESTFILEREADER_HH
#define BEAMTESTFILEREADER_HH

#include <fstream>
#include <string>

class FileReader {
	public:
		FileReader(const char* filename);
		FileReader(const std::string& filename);
		
		bool nextLine();
	
		int getFieldAsInt(const int n);
		float getFieldAsFloat(const int n);
		double getFieldAsDouble(const int n);
		std::string getFieldAsString(const int n);

		bool inputFailed() const;
		bool isValid() const;
	private:
		void skip_fields(std::istringstream& ist, const int n);
		std::ifstream file;
		std::string line;
		bool failed;
};

#endif
