#include "XrayFluoDataSet.hh"
#include "g4std/fstream"
#include "g4std/strstream"
#include "G4VDataSetAlgorithm.hh"

XrayFluoDataSet::XrayFluoDataSet(G4int Z,
			 G4DataVector* points, 
			 G4DataVector* values,
			 const G4VDataSetAlgorithm* interpolation,
			 G4double unitE, G4double unitData)
  :z(Z), energies(points), data(values), algorithm(interpolation)
{
  numberOfBins = energies->size();
  unit1 = unitE;
  unit2 = unitData;
}

XrayFluoDataSet:: XrayFluoDataSet(G4int Z, 
			  const G4String& dataFile,
			  const G4VDataSetAlgorithm* interpolation,
			  G4double unitE, G4double unitData)
  :z(Z), algorithm(interpolation)
{
  energies = new G4DataVector;
  data = new G4DataVector;
  unit1 = unitE;
  unit2 = unitData;  
  LoadData(dataFile);
  numberOfBins = energies->size();
}


// Destructor

XrayFluoDataSet::~XrayFluoDataSet()
{ 
  delete energies;
  delete data;
}


G4double XrayFluoDataSet::FindValue(G4double e, G4int id) const
{
  G4double value;
  G4double e0 = (*energies)[0];
  // Protections
  size_t bin = FindBinLocation(e);
  if (bin == numberOfBins)
    {
     
      value = (*data)[bin];
    }
  else if (e <= e0)
    {
    
      value = (*data)[0];
    }
  else
    {
      value = algorithm->Calculate(e,bin,*energies,*data);
    }
  
  return value;
}

G4int XrayFluoDataSet::FindBinLocation(G4double energy) const
{
  // Protection against call outside allowed range
  G4double e0 = (*energies)[0];
  if (energy < e0)
    {
    
      energy = e0;
    }

  size_t lowerBound = 0;
  size_t upperBound = numberOfBins - 1;
  
  // Binary search
  while (lowerBound <= upperBound) 
    {
      size_t midBin = (lowerBound + upperBound)/2;
      if ( energy < (*energies)[midBin] ) upperBound = midBin-1;
      else lowerBound = midBin+1;
  }
 
  return upperBound;
}


void XrayFluoDataSet::LoadData(const G4String& fileName)
{
  // Build the complete string identifying the file with the data set
  
  char nameChar[100] = {""};
  G4std::ostrstream ost(nameChar, 100, G4std::ios::out);
  
  ost << fileName <<".dat";
  
  G4String name(nameChar);
  
  char* path = "/mnt/home/guardi/workdir/diffuseFluo/XrayFluo";
 
  G4String pathString(path);
  G4String dirFile = pathString + "/" + name;
  G4std::ifstream file(dirFile);
  G4std::filebuf* lsdp = file.rdbuf();
  
  if (! (lsdp->is_open()) )
	{
	  G4String excep = "XrayFluoDataSet - data file: " + dirFile + " not found";
	  G4Exception(excep);
	}
  G4double a = 0;
  G4int k = 1;

  do
    {
      file >> a;
      G4int nColumns = 2;
      // The file is organized into two columns:
      // 1st column is the energy
      // 2nd column is the corresponding value
      // The file terminates with the pattern: -1   -1
      //                                       -2   -2
      if (a == -1 || a == -2)
	{
	 
	}
      else
	{
	  if (k%nColumns != 0)
	    {	
	      G4double e = a * unit1;
	      energies->push_back(e);
	    
	      k++;

	    }
	  else if (k%nColumns == 0)
	    {
	      G4double value = a * unit2;
	      data->push_back(value);
	     
	      k = 1;
	    }
	}
      
    } while (a != -2); // end of file
  
  file.close();
}
void XrayFluoDataSet::PrintData() const
{
  size_t size = numberOfBins;
  for (size_t i=0; i<size; i++)
    {
      G4double e = (*energies)[i]  / unit1;
      G4double sigma = (*data)[i] / unit2 ;
      G4cout << "Point: "
	     << e
	     << " - Data value : "
	     << sigma 
	     << G4endl; 
    }
}





