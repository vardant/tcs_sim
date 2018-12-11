#include <TTree.h>
#include <TFile.h>
#include <TROOT.h>
#include <iostream>

using namespace std;

// Generate class files for the given root file.

void make_class(string root_file="RootFiles.tracked/DEEPGen_10000_tracked.root")
{
  //  TFile* r_file = new TFile(root_file.c_str());

  //  vector<string> tree_name {"beam", "target", "calo", "hodox", "hodoy",
  //      "tracker", "kin"};

  // for (vector<string>::iterator it = tree_name.begin(); it < tree_name.end();
  //       it++) {
  //    TTree *T = (TTree*)gROOT->FindObject((*it).c_str());
  //    T->MakeClass((*it).c_str());
  //  }

  TFile* tcs_file = new TFile("RootFiles/DEEPGen_10000.root");
  TTree *T = (TTree*)gROOT->FindObject("TCS_Tree");
  T->MakeClass("sim_tcs");

}
