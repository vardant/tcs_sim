#include "deepgen.h"
#include <TSystem.h>

using namespace std;

//void run_deepgen(string fname) {
int main(int argc, char* argv[]) {

  string fname = argv[1];

  //gSystem->Load("/work/hallc/nps/vardan/G4.nps/g4.10.03/tcs/aux/tcsgen_data.test/deepgen.C");
  deepgen d(fname);
  d.Loop();

}
