#include <fstream>
#include <iostream>
#include <string>


extern "C" system(const char *);

main()
{
  int count = 0;
  for(;;)
  {
    count ++;
    cout << "iteration "<<count<<endl;
    string i1, i2, dummy;

    system("rm  ...i1; ls -lt > ...i1");
    ifstream it1("...i1");
    it1 >> dummy >> dummy>> dummy>> dummy>> dummy>> dummy;
    it1 >> i1;
    cout << i1 << endl;
    system("sleep 600");
    
    system("rm  ...i2; ls -lt > ...i2");
    ifstream it2("...i2");
    it2 >> dummy >> dummy>> dummy>> dummy>> dummy>> dummy;
    it2 >> i2;
    cout << i2 << endl;   
    if(i1 == i2) 
    {
      system("rm  ...i3; ps -uhpw | grep robustness > ...i3");
      string i3;
      ifstream it3("...i3");
      it3 >> i3;
      i3 = "kill "+i3;;
      system(i3.c_str());
    }
  }
}
