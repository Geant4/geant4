#include "../include/G4BertiniData.hh"

// Initialize static pointer for singleton instance
G4BertiniData* G4BertiniData::theInstance = 0;

//Initialize verboseLevel
G4int G4BertiniData::verboseLevel = 0;

G4BertiniData::G4BertiniData(){

  ifstream dataFile;
  dataFile.open("chetc.dat");

  G4double temp;
  G4int i;

  //let's read whole data from file
  for (i = 0; i < DATASIZE; i++){
    dataFile >> temp;
    XSec.push_back(temp);
  }
 
  //let's split up the data to corresponding vectors 
  for(i=0; i < DOUBLE_SIZE; i++)
    ProtonProtonDoubleProdXSec.push_back(XSec[G4ProtonProtonDoubleProd + i]);

  for(i=0; i < DOUBLE_SIZE; i++)
    NeutronProtonDoubleProdXSec.push_back(XSec[G4NeutronProtonDoubleProd + i]);

   for(i=0; i < NUCLEON_SINGLE_SIZE; i++)
    ProtonProtonSingleProdXSec.push_back(XSec[G4ProtonProtonSingleProd + i]);

  for(i=0; i < NUCLEON_SINGLE_SIZE; i++)
    NeutronProtonSingleProdXSec.push_back(XSec[G4NeutronProtonSingleProd + i]);

  for(i=0; i < PION_SINGLE_SIZE; i++)
    PionPlusProtonSingleProdXSec.push_back(XSec[G4PionPlusProtonSingleProd + i]); 
 
  for(i=0; i < PION_SINGLE_SIZE; i++)
    PionMinusProtonSingleProdXSec.push_back(XSec[G4PionMinusProtonSingleProd + i]); 
  
  for(i=0; i < PION_SINGLE_SIZE ; i++)
    PionZeroProtonSingleProdXSec.push_back(XSec[G4PionZeroProtonSingleProd + i]); 
 
  for(i=0; i < PION_SINGLE_SIZE; i++)
    PionMinusNeutronSingleProdXSec.push_back(XSec[G4PionMinusNeutronSingleProd + i]);

  for(i=0; i < NUCLEON_TOTAL_SIZE; i++)
    ProtonProtonElasticXSec.push_back(XSec[G4ProtonProtonElastic + i]);

  for(i=0; i < NUCLEON_TOTAL_SIZE; i++)
    NeutronProtonElasticXSec.push_back(XSec[G4NeutronProtonElastic + i]); 

  for(i=0; i < PION_TOTAL_SIZE; i++)
    PionPlusProtonElasticXSec.push_back(XSec[G4PionPlusProtonElastic + i]); 
 
  for(i=0; i < PION_TOTAL_SIZE; i++)
    PionMinusProtonElasticXSec.push_back(XSec[G4PionMinusProtonElastic + i]); 
  
  for(i=0; i < PION_TOTAL_SIZE; i++)
    PionZeroProtonElasticXSec.push_back(XSec[G4PionZeroProtonElastic + i]); 
 
  for(i=0; i < PION_TOTAL_SIZE; i++)
    PionMinusNeutronElasticXSec.push_back(XSec[G4PionMinusNeutronElastic + i]); 

  for(i=0; i < PION_TOTAL_SIZE; i++)
    PionMinusProtonExchangeXSec.push_back(XSec[G4PionMinusProtonExchange + i]); 
  
  //let's initialize private data tables
  G4double LowEnergyTable[LOW_ENERGY_SIZE] = 
  { 
    1.0,  	3.0,  	5.0,  
    8.0,  	10.0, 	12.5,
    15.0, 	17.5,	20.0 
  };
  
  G4double PPLowEnergyElasticXSec[LOW_ENERGY_SIZE] =
  {
  1.20,	0.620, 0.310,	       	
  0.208,0.155, 2.28,	    
  1.14, 0.765, 0.555	     
  }; 

  G4double NPLowEnergyElasticXSec[LOW_ENERGY_SIZE] = 
  {
    0.890, 0.392, 0.250,	
    0.178, 4.25,  1.62,	
    0.940, 0.645, 0.480
  }; 

  //let's initialize of STL vectors with the help of private data tables above
 for(i=0; i < LOW_ENERGY_SIZE; i++){ 
  LowEnergy.push_back(LowEnergyTable[i]);
  ProtonProtonLowEnergyElasticXSec.push_back(PPLowEnergyElasticXSec[i]);
  NeutronProtonLowEnergyElasticXSec.push_back(NPLowEnergyElasticXSec[i]);
 }
 
  //let's count the total cross sections
  G4double sum;

  //proton-proton total cross section = single + double + elastic 
  for(i=0; i < NUCLEON_TOTAL_SIZE; i++){ 
    sum = ProtonProtonElasticXSec[i];
    if(i >= (NUCLEON_TOTAL_SIZE - NUCLEON_SINGLE_SIZE)) //the sizes of component vectors are not equal! 
      sum += ProtonProtonSingleProdXSec[i - (NUCLEON_TOTAL_SIZE - NUCLEON_SINGLE_SIZE)];
    if(i >= (NUCLEON_TOTAL_SIZE - DOUBLE_SIZE))
      sum += ProtonProtonDoubleProdXSec[i - (NUCLEON_TOTAL_SIZE - DOUBLE_SIZE)];
    ProtonProtonTotalXSec.push_back(sum);
  }

  
  //neutron-proton total cross section = single + double + elastic 
  for(i=0; i < NUCLEON_TOTAL_SIZE; i++){ 
    sum = NeutronProtonElasticXSec[i]; 
    if(i >= (NUCLEON_TOTAL_SIZE - NUCLEON_SINGLE_SIZE))
      sum += NeutronProtonSingleProdXSec[i - (NUCLEON_TOTAL_SIZE - NUCLEON_SINGLE_SIZE)];
    if(i >=(NUCLEON_TOTAL_SIZE - DOUBLE_SIZE))
      sum += NeutronProtonDoubleProdXSec[i - (NUCLEON_TOTAL_SIZE - DOUBLE_SIZE)];
    NeutronProtonTotalXSec.push_back(sum);
  }

  //pion+ - proton total cross section = single + elastic
  for(i=0; i < PION_TOTAL_SIZE; i++){ 
    sum = PionPlusProtonElasticXSec[i]; 
    if(i >= (PION_TOTAL_SIZE - PION_SINGLE_SIZE))
      sum += PionPlusProtonSingleProdXSec[i - (PION_TOTAL_SIZE - PION_SINGLE_SIZE)];
    PionPlusProtonTotalXSec.push_back(sum);
  }  
  //pion- -proton total cross section = single + exchange + elastic
  for(i=0; i < PION_TOTAL_SIZE; i++){  
    sum = PionMinusProtonElasticXSec[i] + PionMinusProtonExchangeXSec[i]; 
    if(i >= (PION_TOTAL_SIZE - PION_SINGLE_SIZE))
      sum += PionMinusProtonSingleProdXSec[i - (PION_TOTAL_SIZE - PION_SINGLE_SIZE)];
    PionMinusProtonTotalXSec.push_back(sum);
  }    
 //pion0 - proton total cross section = single + elastic 
  for(i=0; i < PION_TOTAL_SIZE; i++){ 
    sum = PionZeroProtonElasticXSec[i]; 
    if(i >= (PION_TOTAL_SIZE - PION_SINGLE_SIZE))
      sum += PionZeroProtonSingleProdXSec[i - (PION_TOTAL_SIZE - PION_SINGLE_SIZE)];
    PionZeroProtonTotalXSec.push_back(sum);
  }    
  //pion- neutron total cross section = single + elastic
  for(i=0; i < PION_TOTAL_SIZE; i++){ 
    sum = PionMinusNeutronElasticXSec[i]; 
    if(i >= (PION_TOTAL_SIZE - PION_SINGLE_SIZE))
      sum += PionMinusNeutronSingleProdXSec[i - (PION_TOTAL_SIZE - PION_SINGLE_SIZE)];
    PionMinusNeutronTotalXSec.push_back(sum);
  }  
}

G4BertiniData::~G4BertiniData(){
  delete theInstance;
}  

void G4BertiniData::SetVerboseLevel(G4int level){
  verboseLevel = level;
}
   
G4double G4BertiniData::GetCrossSection(G4Interaction interaction,
				      G4double particleEnergy){

 //index for upper tab.value corresponding to 'particleEnergy' 
  G4int i = VectorIndex(interaction, particleEnergy);
 
 //lower tabulated energy corresponding to 'particleEnergy'
  G4double prevTableEnergy = Energy(interaction, i);

  //a number that tells how the 'particleEnergy' is situated compared to lower and 
  //upper tabulated energy values
  G4double fraction = (particleEnergy - prevTableEnergy) / ENERGY_SPACING;

  switch(interaction){
  case G4ProtonProtonSingleProd:
      return ProtonProtonSingleProdXSec[i - 1] + fraction
	* (ProtonProtonSingleProdXSec[i] - ProtonProtonSingleProdXSec[i - 1]); //linear interpolation
      
   
  case G4ProtonProtonDoubleProd:
      return ProtonProtonDoubleProdXSec[i - 1] + fraction
	* (ProtonProtonDoubleProdXSec[i] - ProtonProtonDoubleProdXSec[i - 1]);      
      
  
  case G4ProtonProtonElastic:
    if(particleEnergy < ENERGY_SPACING){ //use low energy cross sections
      G4int ie;
        for (ie  = 1; ie < LOW_ENERGY_SIZE; ie++) {
	  if (particleEnergy <= LowEnergy[ie]) {
	    G4double temp = log(ProtonProtonLowEnergyElasticXSec[ie - 1]) +
	      (particleEnergy - LowEnergy[ie - 1]) / (LowEnergy[ie] - LowEnergy[ie - 1]) *
	      (log(ProtonProtonLowEnergyElasticXSec[ie]) - log(ProtonProtonLowEnergyElasticXSec[ie - 1])); 
	    return exp(temp) * 1.0e-24;
	  }
 	}
    }
    else{
      return ProtonProtonElasticXSec[i - 1] + fraction
	* (ProtonProtonElasticXSec[i] - ProtonProtonElasticXSec[i - 1]); 
    }      
  
  case G4ProtonProtonTotal:
      return ProtonProtonTotalXSec[i - 1] + fraction
	* (ProtonProtonTotalXSec[i] - ProtonProtonTotalXSec[i - 1]);
      


  case G4NeutronProtonSingleProd:
      return NeutronProtonSingleProdXSec[i - 1] + fraction
	* (NeutronProtonSingleProdXSec[i] - NeutronProtonSingleProdXSec[i - 1]);
      
  
  case G4NeutronProtonDoubleProd:
      return NeutronProtonDoubleProdXSec[i - 1] + fraction
	* (NeutronProtonDoubleProdXSec[i] - NeutronProtonDoubleProdXSec[i - 1]);      
      

  case G4NeutronProtonElastic:
    if(particleEnergy < ENERGY_SPACING){
      G4int ie;
        for (ie  = 1; ie < LOW_ENERGY_SIZE; ie++) {
	  if (particleEnergy <= LowEnergy[ie]) {
	    G4double temp = log(NeutronProtonLowEnergyElasticXSec[ie - 1]) +
	      (particleEnergy - LowEnergy[ie - 1]) / (LowEnergy[ie] - LowEnergy[ie - 1]) *
	      (log(NeutronProtonLowEnergyElasticXSec[ie]) - log(NeutronProtonLowEnergyElasticXSec[ie - 1])); 
	    return exp(temp) * 1.0e-24;
	  }
 	}
    }
    else{
      return NeutronProtonElasticXSec[i - 1] + fraction
	* (NeutronProtonElasticXSec[i] - NeutronProtonElasticXSec[i - 1]); 
    }
      
  case G4NeutronProtonTotal:
      return NeutronProtonTotalXSec[i - 1] + fraction
	* (NeutronProtonTotalXSec[i] - NeutronProtonTotalXSec[i - 1]);
    
      
      
  case G4PionPlusProtonSingleProd:
      return PionPlusProtonSingleProdXSec[i - 1] + fraction
	* (PionPlusProtonSingleProdXSec[i] - PionPlusProtonSingleProdXSec[i - 1]);
      
  
  
  case G4PionPlusProtonElastic: 
      return PionPlusProtonElasticXSec[i - 1] + fraction
	* (PionPlusProtonElasticXSec[i] - PionPlusProtonElasticXSec[i - 1]); 
      

  case G4PionPlusProtonTotal:
      return PionPlusProtonTotalXSec[i - 1] + fraction
	* (PionPlusProtonTotalXSec[i] - PionPlusProtonTotalXSec[i - 1]);
      


  case G4PionMinusProtonSingleProd:
      return PionMinusProtonSingleProdXSec[i - 1] + fraction
	* (PionMinusProtonSingleProdXSec[i] - PionMinusProtonSingleProdXSec[i - 1]);
      

  case G4PionMinusProtonExchange:
      return PionMinusProtonExchangeXSec[i - 1] + fraction
	* (PionMinusProtonExchangeXSec[i] - PionMinusProtonExchangeXSec[i - 1]);   
      

  case G4PionMinusProtonElastic:
      return PionMinusProtonElasticXSec[i - 1] + fraction
	* (PionMinusProtonElasticXSec[i] - PionMinusProtonElasticXSec[i - 1]); 
      

  case G4PionMinusProtonTotal:
      return PionMinusProtonTotalXSec[i - 1] + fraction
	* (PionMinusProtonTotalXSec[i] - PionMinusProtonTotalXSec[i - 1]);
      

 
  case G4PionZeroProtonSingleProd:
      return PionZeroProtonSingleProdXSec[i - 1] + fraction
	* (PionZeroProtonSingleProdXSec[i] - PionZeroProtonSingleProdXSec[i - 1]);
      

  case G4PionZeroProtonElastic:
      return PionZeroProtonElasticXSec[i - 1] + fraction
	* (PionZeroProtonElasticXSec[i] - PionZeroProtonElasticXSec[i - 1]); 
      

  case G4PionZeroProtonTotal:    
      return PionZeroProtonTotalXSec[i - 1] + fraction
	* (PionZeroProtonTotalXSec[i] - PionZeroProtonTotalXSec[i - 1]);
      
 

  case G4PionMinusNeutronSingleProd:
      return PionMinusNeutronSingleProdXSec[i - 1] + fraction
	* (PionMinusNeutronSingleProdXSec[i] - PionMinusNeutronSingleProdXSec[i - 1]);
      

  case G4PionMinusNeutronElastic:
      return PionMinusNeutronElasticXSec[i - 1] + fraction
	* (PionMinusNeutronElasticXSec[i] - PionMinusNeutronElasticXSec[i - 1]); 
      

  case G4PionMinusNeutronTotal:
      return PionMinusNeutronTotalXSec[i - 1] + fraction
	* (PionMinusNeutronTotalXSec[i] - PionMinusNeutronTotalXSec[i - 1]);
      
  }

}  
  
G4int G4BertiniData::VectorIndex(G4Interaction interaction,
				 G4double particleEnergy){

    if(interaction == G4ProtonProtonSingleProd || interaction == G4NeutronProtonSingleProd)
	return static_cast<G4int>((particleEnergy - NUCLEON_SINGLE_START_ENERGY) / ENERGY_SPACING + 1.0);

    if(interaction == G4ProtonProtonDoubleProd || interaction == G4NeutronProtonDoubleProd )
	return static_cast<G4int>((particleEnergy - NUCLEON_DOUBLE_START_ENERGY) / ENERGY_SPACING + 1.0);
    
    if(interaction==G4PionPlusProtonSingleProd || interaction==G4PionMinusProtonSingleProd ||
       interaction==G4PionZeroProtonSingleProd || interaction==G4PionMinusNeutronSingleProd )
      return static_cast<G4int>((particleEnergy - PION_SINGLE_START_ENERGY) / ENERGY_SPACING + 1.0);
    
    else
      return static_cast<G4int>(particleEnergy / ENERGY_SPACING + 1.0); 
}

G4double G4BertiniData::Energy(G4Interaction interaction, G4int index){
    if(interaction == G4ProtonProtonSingleProd || interaction == G4NeutronProtonSingleProd)
	return (index-1)*ENERGY_SPACING + NUCLEON_SINGLE_START_ENERGY;

    if(interaction == G4ProtonProtonDoubleProd || interaction == G4NeutronProtonDoubleProd )
	return (index-1)*ENERGY_SPACING + NUCLEON_DOUBLE_START_ENERGY;
    
    if(interaction==G4PionPlusProtonSingleProd || interaction==G4PionMinusProtonSingleProd ||
       interaction==G4PionZeroProtonSingleProd || interaction==G4PionMinusNeutronSingleProd )
      return (index-1)*ENERGY_SPACING + PION_SINGLE_START_ENERGY;
    
    else
      return (index-1)*ENERGY_SPACING;
}








