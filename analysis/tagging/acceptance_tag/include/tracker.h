//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Dec  7 22:32:02 2018 by ROOT version 6.14/04
// from TTree tracker/TCS trackers' per event hit collections
// found on file: RootFiles.tracked/DEEPGen_10000_tracked.root
//////////////////////////////////////////////////////////

#ifndef tracker_h
#define tracker_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class tracker {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   vector<double>  *xcont;
   vector<double>  *ycont;
   vector<double>  *edepcont;
   vector<double>  *lengthcont;
   vector<int>     *detcont;
   vector<int>     *layercont;
   vector<int>     *pidcont;
   vector<int>     *pidorigcont;
   vector<int>     *trackidcont;
   vector<int>     *nstepcont;

   // List of branches
   TBranch        *b_xcont;   //!
   TBranch        *b_ycont;   //!
   TBranch        *b_edepcont;   //!
   TBranch        *b_lengthcont;   //!
   TBranch        *b_detcont;   //!
   TBranch        *b_layercont;   //!
   TBranch        *b_pidcont;   //!
   TBranch        *b_pidorigcont;   //!
   TBranch        *b_trackidcont;   //!
   TBranch        *b_nstepcont;   //!

   tracker(TTree *tree=0);
   tracker(string file_name, TTree *tree=0);   ////
   virtual ~tracker();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   ////   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   bool CheckTrack(long entry, int pid, int layer, int quarter);   ////
};

#endif

#ifdef tracker_cxx
tracker::tracker(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("RootFiles.tracked/DEEPGen_10000_tracked.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("RootFiles.tracked/DEEPGen_10000_tracked.root");
      }
      f->GetObject("tracker",tree);

   }
   Init(tree);

   cout << "A tracker class object created" << endl;
}

tracker::~tracker()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t tracker::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t tracker::LoadTree(Long64_t entry)
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

void tracker::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   xcont = 0;
   ycont = 0;
   edepcont = 0;
   lengthcont = 0;
   detcont = 0;
   layercont = 0;
   pidcont = 0;
   pidorigcont = 0;
   trackidcont = 0;
   nstepcont = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("xcont", &xcont, &b_xcont);
   fChain->SetBranchAddress("ycont", &ycont, &b_ycont);
   fChain->SetBranchAddress("edepcont", &edepcont, &b_edepcont);
   fChain->SetBranchAddress("lengthcont", &lengthcont, &b_lengthcont);
   fChain->SetBranchAddress("detcont", &detcont, &b_detcont);
   fChain->SetBranchAddress("layercont", &layercont, &b_layercont);
   fChain->SetBranchAddress("pidcont", &pidcont, &b_pidcont);
   fChain->SetBranchAddress("pidorigcont", &pidorigcont, &b_pidorigcont);
   fChain->SetBranchAddress("trackidcont", &trackidcont, &b_trackidcont);
   fChain->SetBranchAddress("nstepcont", &nstepcont, &b_nstepcont);
   Notify();
   cout << "A tracker class object initiated" << endl;
}

Bool_t tracker::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void tracker::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t tracker::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

//==============================================================================

//Additional methods

tracker::tracker(string file_name, TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
     TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(file_name.c_str());
      if (!f || !f->IsOpen()) {
	f = new TFile(file_name.c_str());
      }
      f->GetObject("tracker",tree);

   }
   Init(tree);

   cout << "A tracker class object created on file" << file_name << endl;
}

bool tracker::CheckTrack(long entry, int pid, int layer, int quarter) {

  int trackid = 0;
  switch (pid) {
  case   11: trackid = 1; break;
  case  -11: trackid = 2; break;
  case 2212: trackid = 3; break;;
  default: cout << "*** tracker::IsTrack: wrong pid = " << pid << " ***"
		<< endl;
  }

  Long64_t tentry = GetEntry(entry);

  bool flag = false;
  for (UInt_t j = 0; j < detcont->size(); ++j) {
    if (detcont->at(j) == quarter && layercont->at(j) == layer &&
	pidorigcont->at(j) == pid && trackidcont->at(j) == trackid) {
      flag = true;
      break;
    }
  }

  return flag;
}

#endif // #ifdef tracker_cxx
