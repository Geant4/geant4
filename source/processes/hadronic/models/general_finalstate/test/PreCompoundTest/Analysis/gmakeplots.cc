#include "TROOT.h"
#include "TApplication.h"

#include "mkplmakeplots.h"

using namespace std;

int main(int argc, char ** argv)
{
  TApplication theApp("App", &argc, argv);

  //  mkplmakeplots MainWindow(gClient->GetRoot(), 700, 500);
  mkplmakeplots MainWindow(gClient->GetRoot(), 500, 500);

  theApp.Run();

  return 0;
}
