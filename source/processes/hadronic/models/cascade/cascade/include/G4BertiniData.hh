//This class reads HETC cross section data from disk and splits it up to 
//corresponding cross section vectors. It can also interpolate cross section for a given energy
//and interaction type. This class is implemented as a SINGLETON class.

#ifndef G4BERTINIDATA
#define G4BERTINIDATA
#include <fstream.h>
#include <math.h>
#include <vector>
#include "G4Types.hh"

//typedef int G4int;
//typedef float G4double;
typedef std::vector<G4double>::const_iterator iterator; //iterator for traversing vectors

//possible interaction types
enum G4Interaction {G4ProtonProtonSingleProd=3395, G4NeutronProtonSingleProd=3683,
		    G4PionPlusProtonSingleProd=4409, G4PionMinusProtonSingleProd=4643,
		    G4PionZeroProtonSingleProd=4526, G4PionMinusNeutronSingleProd=4760,
		    G4ProtonProtonDoubleProd=3553, G4NeutronProtonDoubleProd=3841,
		    G4ProtonProtonElastic=6193, G4NeutronProtonElastic=6369,
		    G4PionPlusProtonElastic=6067, G4PionMinusProtonElastic=5941,
		    G4PionZeroProtonElastic=5563, G4PionMinusNeutronElastic=5689,
		    G4PionMinusProtonExchange=5815,
                    G4ProtonProtonTotal, G4NeutronProtonTotal, G4PionPlusProtonTotal,
                    G4PionMinusProtonTotal, G4PionZeroProtonTotal, G4PionMinusNeutronTotal};

//sizes of data vectors 
enum { DATASIZE = 28950, NUCLEON_TOTAL_SIZE = 176, PION_TOTAL_SIZE = 126,
       NUCLEON_SINGLE_SIZE = 158, PION_SINGLE_SIZE = 117, DOUBLE_SIZE = 130,
       LOW_ENERGY_SIZE = 9}; 
 
class G4BertiniData
{
protected:
    
  G4BertiniData(); 

public:
  ~G4BertiniData();

  static G4BertiniData* Instance(){
    if(verboseLevel > 0){
      if(theInstance!=0) cout << "Warning, using existing instance!" << endl;
    }
    if (!theInstance) theInstance = new G4BertiniData();

    return theInstance;
  }
  
  //interpolates the cross section for a given energy and interaction
  G4double GetCrossSection(G4Interaction interaction,
			   G4double particleEnergy);

  void SetVerboseLevel(G4int level);

  //variables:
  vector<G4double> ProtonProtonSingleProdXSec;
  vector<G4double> NeutronProtonSingleProdXSec;
  vector<G4double> PionPlusProtonSingleProdXSec;
  vector<G4double> PionMinusProtonSingleProdXSec;
  vector<G4double> PionZeroProtonSingleProdXSec;
  vector<G4double> PionMinusNeutronSingleProdXSec;

  vector<G4double> ProtonProtonDoubleProdXSec;
  vector<G4double> NeutronProtonDoubleProdXSec;

  vector<G4double> ProtonProtonElasticXSec;
  vector<G4double> NeutronProtonElasticXSec;
  vector<G4double> PionPlusProtonElasticXSec;
  vector<G4double> PionMinusProtonElasticXSec;
  vector<G4double> PionZeroProtonElasticXSec;
  vector<G4double> PionMinusNeutronElasticXSec;

  vector<G4double> PionMinusProtonExchangeXSec;

  vector<G4double> ProtonProtonTotalXSec;
  vector<G4double> NeutronProtonTotalXSec;
  vector<G4double> PionPlusProtonTotalXSec;
  vector<G4double> PionMinusProtonTotalXSec;
  vector<G4double> PionZeroProtonTotalXSec;
  vector<G4double> PionMinusNeutronTotalXSec;

  vector<G4double> ProtonProtonLowEnergyElasticXSec;
  vector<G4double> NeutronProtonLowEnergyElasticXSec;
  vector<G4double> LowEnergy;
 
private:
  
//constants:
  
//starting energies for single and double production cross sections
  static const G4double NUCLEON_SINGLE_START_ENERGY = 360.0; 
  static const G4double PION_SINGLE_START_ENERGY = 180.0; 
  static const G4double NUCLEON_DOUBLE_START_ENERGY = 920.0; 
  static const G4double ENERGY_SPACING = 20.0; //difference of discrete energy points used in cross section calculation
 
//variables:
  vector<G4double> XSec;  //whole cross section data
  static G4BertiniData* theInstance;
  static G4int verboseLevel;
  G4double LowEnergyTable[LOW_ENERGY_SIZE];
  G4double PPLowEnergyElasticXSec[LOW_ENERGY_SIZE]; //energy < 20 MeV
  G4double NPLowEnergyElasticXSec[LOW_ENERGY_SIZE]; //energy < 20 MeV

//methods:

//returns the index which is used in cross section tabulation
//(index to upper tabulated value for 'particleEnergy')
G4int VectorIndex(G4Interaction interaction,
		  G4double particleEnergy);
 
//returns the tabeled value of energy corresponding to index 'index'
// (lower tabulated value for 'particleEnergy')  
G4double Energy(G4Interaction interaction,
		G4int index);

};
#endif






