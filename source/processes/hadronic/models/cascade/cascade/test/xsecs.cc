#include <fstream>
#include "../cascade/src/G4BertiniData.cc" 

//test method for all the data vectors in G4BertiniData
//particle1 = 1,2,3,4,5; particle2 = 1,2; process = 4 (total) or x (all) 

int main(){

  G4BertiniData* dataBase;
  G4BertiniData* instance = dataBase->Instance();
  G4int particle1, particle2, process, j=0;
  G4double energy[176];

  for(j=0; j < 176; j++)
    energy[j] = 20.0*(j+1);

  ifstream in;

  in.open("xsecs.in", ios::in);

  in >> particle1 >> particle2 >> process;
  in.close();

  if(particle1 < 1 || particle1 > 5 || particle2 < 1 || particle2 > 2){
    cout << "not valid parameter values" << endl;
    return 1;
  }
  
  iterator i;

  if(particle2 == 1){
  switch(particle1){
  case 1:
      if(process == 4){
	j=0;
      for(i=instance->ProtonProtonTotalXSec.begin();
	  i < instance->ProtonProtonTotalXSec.end(); i++){
	cout << energy[j] <<"\t"<<  *i << endl;
	j++;
	  }
      break;
      }
      else{
        for(i = instance->ProtonProtonSingleProdXSec.begin();
        i < instance->ProtonProtonSingleProdXSec.end(); i++){
	cout << *i << endl;
	}
	

   
      for(i=instance->ProtonProtonDoubleProdXSec.begin();
	  i < instance->ProtonProtonDoubleProdXSec.end(); i++){
	cout << *i << endl;
       }

      for(i=instance->ProtonProtonElasticXSec.begin();
	  i < instance->ProtonProtonElasticXSec.end(); i++){
	cout << *i << endl;
      }

      for(i=instance->ProtonProtonTotalXSec.begin();
	  i < instance->ProtonProtonTotalXSec.end(); i++){
	cout << *i << endl;
	  }
      
      break;
      }
  case 2:
    
	  for(i=instance->NeutronProtonSingleProdXSec.begin();
	  i < instance->NeutronProtonSingleProdXSec.end(); i++){
	    cout << *i << endl;
	    //j++;
	  }
	  	
	  //j=0;
	  for(i=instance->NeutronProtonDoubleProdXSec.begin();
	      i < instance->NeutronProtonDoubleProdXSec.end(); i++){
	    cout << *i << endl;
	    //j++;
	      }
	  //j=0;

	  for(i=instance->NeutronProtonElasticXSec.begin();
	      i < instance->NeutronProtonElasticXSec.end(); i++){
	    cout << *i << endl;
	    //j++;
	      }
	  // j=0;

	  for(i=instance->NeutronProtonTotalXSec.begin();
	      i < instance->NeutronProtonTotalXSec.end(); i++){
	    cout << *i << endl;
	    //j++;
	      }
	  
  	break;
  case 3:
        
	  for(i=instance->PionPlusProtonSingleProdXSec.begin();
	  i < instance->PionPlusProtonSingleProdXSec.end(); i++){
	    cout << *i << endl;
	   }
	 	
	   for(i=instance->PionPlusProtonElasticXSec.begin();
	    i < instance->PionPlusProtonElasticXSec.end(); i++){
	    cout << *i << endl;
	    }
	 
	
	  for(i=instance->PionPlusProtonTotalXSec.begin();
	      i < instance->PionPlusProtonTotalXSec.end(); i++){
	    cout << *i << endl;
	      }
	
	 break;
  case 4:
       
       if(process == 4){
	j=0;
      for(i=instance->PionMinusProtonTotalXSec.begin();
	  i < instance->PionMinusProtonTotalXSec.end(); i++){
	cout << energy[j] <<"\t"<<  *i << endl;
	j++;
      }
      break;
       }

       else{
	  for(i=instance->PionMinusProtonSingleProdXSec.begin();
	  i < instance->PionMinusProtonSingleProdXSec.end(); i++){
	    cout << *i << endl;
	      }
	
	  for(i=instance->PionMinusProtonExchangeXSec.begin();
	      i < instance->PionMinusProtonExchangeXSec.end(); i++){
	    cout << *i << endl;
	      }
	
	  for(i=instance->PionMinusProtonElasticXSec.begin();
	      i < instance->PionMinusProtonElasticXSec.end(); i++){
	    cout << *i << endl;
	      }
	
	
	  for(i=instance->PionMinusProtonTotalXSec.begin();
	      i < instance->PionMinusProtonTotalXSec.end(); i++){
	    cout << *i << endl;
	      }
	
	break;
       }
  case 5:
        
	
	  for(i=instance->PionZeroProtonSingleProdXSec.begin();
	  i < instance->PionZeroProtonSingleProdXSec.end(); i++){
	    cout << *i << endl;
	      }
	
	
	  for(i=instance->PionZeroProtonElasticXSec.begin();
	      i < instance->PionZeroProtonElasticXSec.end(); i++){
	    cout << *i << endl;
	      }
	
	  for(i=instance->PionZeroProtonTotalXSec.begin();
	      i < instance->PionZeroProtonTotalXSec.end(); i++){
	    cout << *i << endl;

	      }
	  break;
  } 
      
  
 
  }
  else{
    
	  for(i=instance->PionMinusNeutronSingleProdXSec.begin();
	  i < instance->PionMinusNeutronSingleProdXSec.end(); i++){
	    cout << *i << endl;
	      }
	
	
	  for(i=instance->PionMinusNeutronElasticXSec.begin();
	      i < instance->PionMinusNeutronElasticXSec.end(); i++){
	    cout << *i << endl;
	      }
	
	  for(i=instance->PionMinusNeutronTotalXSec.begin();
	      i < instance->PionMinusNeutronTotalXSec.end(); i++){
	    cout << *i << endl;
	      }
	
  }

}





