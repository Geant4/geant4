#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <fstream>

int main(){
  char name[128], temp[128], temp2[128], temp3[128], *oldout, *newout, command1[128], *command2, *command3,  c;
  int i=0, j=0;

  cout << "Give the name of the test: " << endl;

  while((c=getchar())!='\n' && i<127){
    name[i]=c;
    i++;
  }
  name[i]='\0';

  ifstream file(name);
  if(!file){
    cout << "The file " << name << " does not exist! " << endl; 
    return 1;
  }
  
  for(j=0; j < i+1; j++) temp[j] = name[j];
  for(j=0; j < i+1; j++) temp2[j] = name[j];
  for(j=0; j < i+1; j++) temp3[j] = name[j];
  
  oldout = name;
  newout = temp;
  

  strcat(oldout, ".old");
  strcat(newout, ".out");
  ifstream oldFile(newout);

  //cout << oldout << endl;
  //cout << newout << endl;

  //command1 = new char(128);
  if(oldFile){
    strcpy(command1, "cp "); 
    strcat(command1, newout); 
    strcat(command1, " "); 
    strcat(command1, oldout); 
    //cout << command1 << endl;
  }
  
  //command2 = new char(128);
  //strcpy(command2, temp2);
  command2=temp2;
  strcat(command2, " > ");
  strcat(command2, newout);
  //cout << command2 << endl;

  command3 = new char(128);
  strcpy(command3, "diff ");
  strcat(command3, oldout);
  strcat(command3, " ");
  strcat(command3, newout); 
  //cout << command3 << endl;
  
  cout << "-------------------------------" << endl;
  if(oldFile){
    cout << command1 << endl;
    system(command1);
  }
  system(temp3);
  cout << "-------------------------------------------" << endl;
  cout << command2 << endl;
  system(command2);
  if(oldFile){
    cout << command3 << endl;
    system(command3);
  }

  char plot;
  cout << "-------------------------------------------" << endl;
  cout << "Do you want to visualize the results (y/n)?" << endl;
  cin >> c;
  if(c == 'y'){
    cout << "Histogram or x-y plot (h/x)?" << endl;
    cin >> plot;
    if(plot == 'x'){
      char* command4 = new char(128);
      strcpy(command4, "cp ");
      strcat(command4, newout); 
      strcat(command4, " draw2.out");  
      cout << command4 << endl;
      system(command4);
      cout << "----------------------------------------------------------------------------------" << endl;
      cout << "Start ROOT by command 'root' and visualize results by command '.x draw2.C' in ROOT" << endl; 
      delete(command4);
    }
  }

  delete(command3);
  return 0;
}


