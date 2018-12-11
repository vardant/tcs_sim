//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Dec  7 22:32:02 2018 by ROOT version 6.14/04
// from TTree hodox/TCS X hodoscopes' per event hit collections
// found on file: RootFiles.tracked/DEEPGen_10000_tracked.root
//////////////////////////////////////////////////////////

#ifndef hodox_h
#define hodox_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"

class hodox {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   vector<int>     *detcont;
   vector<unsigned int> *chancont;
   vector<double>  *edepcont;
   vector<int>     *pidcont;

   // List of branches
   TBranch        *b_detcont;   //!
   TBranch        *b_chancont;   //!
   TBranch        *b_edepcont;   //!
   TBranch        *b_pidcont;   //!

   hodox(TTree *tree=0);
   hodox(string file_name, TTree *tree=0);   ////
   virtual ~hodox();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   ////   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   bool CheckEdep(long entry, int pid, int quarter);   ////
};

#endif

#ifdef hodox_cxx
hodox::hodox(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("RootFiles.tracked/DEEPGen_10000_tracked.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("RootFiles.tracked/DEEPGen_10000_tracked.root");
      }
      f->GetObject("hodox",tree);

   }
   Init(tree);
}

hodox::~hodox()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t hodox::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t hodox::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void hodox::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   detcont = 0;
   chancont = 0;
   edepcont = 0;
   pidcont = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("detcont", &detcont, &b_detcont);
   fChain->SetBranchAddress("chancont", &chancont, &b_chancont);
   fChain->SetBranchAddress("edepcont", &edepcont, &b_edepcont);
   fChain->SetBranchAddress("pidcont", &pidcont, &b_pidcont);
   Notify();
}

Bool_t hodox::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void hodox::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t hodox::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

//==============================================================================

//Additional methods

hodox::hodox(string file_name, TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
     TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(file_name.c_str());
      if (!f || !f->IsOpen()) {
	f = new TFile(file_name.c_str());
      }
      f->GetObject("hodox",tree);

   }
   Init(tree);
}

bool hodox::CheckEdep(long entry, int pid, int quarter) {

  const double Ethr = 1.5; //[MeV]

  int tentry = GetEntry(entry);

  double edep = 0.;
  for (uint j = 0; j < detcont->size(); j++) {
    if (pidcont->at(j) == pid && detcont->at(j) == quarter)
      edep += edepcont->at(j);
  }

  return (edep > Ethr ? true : false);
}

#endif // #ifdef hodox_cxx
