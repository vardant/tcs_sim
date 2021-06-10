#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <iostream>

using namespace std;

//Generate a pair of random integer seeds to be fed to the G4 simulation.

int main ()
{
  /* initialize random seed: */
  srand (time(NULL));

  cout << rand() << " " << rand() << endl;

  return 0;
}
