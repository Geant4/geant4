#include "TROOT.h"
#include "TApplication.h"

#include "mkplmainframe.h"

using namespace std;
 
int main(int argc, char ** argv)
{
  TApplication theApp("App", &argc, argv);

  mkplmainframe MainWindow(gClient->GetRoot(), 700, 500);

  theApp.Run();

  return 0;
}
