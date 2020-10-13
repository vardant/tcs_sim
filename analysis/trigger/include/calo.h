
//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Dec  7 22:32:02 2018 by ROOT version 6.14/04
// from TTree calo/TCS Calorimeters' per event hit collections
// found on file: RootFiles.tracked/DEEPGen_10000_tracked.root
//////////////////////////////////////////////////////////

#ifndef calo_h
#define calo_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"

class calo {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   vector<int>     *detcont;
   vector<unsigned int> *colcont;
   vector<unsigned int> *rowcont;
   vector<double>  *edepcont;
   vector<int>     *pidcont;

   // List of branches
   TBranch        *b_detcont;   //!
   TBranch        *b_colcont;   //!
   TBranch        *b_rowcont;   //!
   TBranch        *b_edepcont;   //!
   TBranch        *b_pidcont;   //!

   calo(TTree *tree=0);
   calo(string file_name, TTree *tree=0);   ////
   virtual ~calo();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t GetEntriesFast();
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   ////   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   bool CheckEdep(long entry, int pid, int quarter);
   double Edep(long entry, int pid, int quarter);
   //   bool Accept(int iq, int ihit);
   bool CheckAccept(long entry, int pid, int quarter);
};

#endif

#ifdef calo_cxx
calo::calo(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("RootFiles.tracked/DEEPGen_10000_tracked.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("RootFiles.tracked/DEEPGen_10000_tracked.root");
      }
      f->GetObject("calo",tree);

   }
   Init(tree);
}

calo::~calo()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t calo::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t calo::LoadTree(Long64_t entry)
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

void calo::Init(TTree *tree)
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
   colcont = 0;
   rowcont = 0;
   edepcont = 0;
   pidcont = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("detcont", &detcont, &b_detcont);
   fChain->SetBranchAddress("colcont", &colcont, &b_colcont);
   fChain->SetBranchAddress("rowcont", &rowcont, &b_rowcont);
   fChain->SetBranchAddress("edepcont", &edepcont, &b_edepcont);
   fChain->SetBranchAddress("pidcont", &pidcont, &b_pidcont);
   Notify();
}

Bool_t calo::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void calo::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t calo::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

//==============================================================================

//Additional methods

calo::calo(string file_name, TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
     TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(file_name.c_str());
      if (!f || !f->IsOpen()) {
	f = new TFile(file_name.c_str());
      }
      f->GetObject("calo",tree);

   }
   Init(tree);
}

//------------------------------------------------------------------------------

bool calo::CheckEdep(long entry, int pid, int quarter) {

  //  const double Ethr = 1000.; //[MeV], for e+, e-.
  const double Ethr = 0.;

  double edep = Edep(entry, pid, quarter);

  return (edep > Ethr ? true : false);
}

//------------------------------------------------------------------------------

double calo::Edep(long entry, int pid, int quarter) {

  int tentry = GetEntry(entry);

  double edep = 0.;
  for (uint j = 0; j < detcont->size(); j++) {
    //if (pidcont->at(j) == pid && detcont->at(j) == quarter && Accept(quarter,j))
    if (pidcont->at(j) == pid && detcont->at(j) == quarter)
      edep += edepcont->at(j);
  }

  return edep;
}
//------------------------------------------------------------------------------

//Check calorimeter acceptance. Calculate weighted coordinates of
//shower (in terms of row and column numbers). Accept event if
//coordinates are within calorimeter (the outer rim of one module
//width excluded).

bool calo::CheckAccept(long entry, int pid, int quarter) {

  const int ncol=23;
  const int nrow=23;

  int tentry = GetEntry(entry);

  double x = 0.;
  double y = 0.;
  double edep = 0.;
  for (uint ihit = 0; ihit < detcont->size(); ihit++) {
    if (pidcont->at(ihit) == pid && detcont->at(ihit) == quarter) {
      double e = edepcont->at(ihit);
      x += (colcont->at(ihit)+1+0.5)*e;
      y += (rowcont->at(ihit)+1+0.5)*e;
      edep += e;
    }
  }

  if (edep>0.) {
    x /= edep;
    y /= edep;
  }

  return (x>1. && x<ncol && y>0. && y<nrow);
  //  return (x>2. && x<ncol-1 && y>2. && y<nrow-1);
}

//------------------------------------------------------------------------------

//Remove 1/8-th corner close to beam pipe.
//Should be applied to the primary tracks.
/*
bool calo::Accept(int iq, int ihit) {

  int j = rowcont->at(ihit) + 1;
  int k = colcont->at(ihit) + 1;

  bool acc = false;

  float x0, y0, x1, y1;
  if (iq==0) {
    x0=23-11.5;
    y0=1;
    x1=23;
    y1=11.5;
    acc = j > y0+(y1-y0)/(x1-x0)*(k-x0);
  }
  if (iq==1) {
    x0=1;
    y0=11.5;
    x1=11.5;
    y1=1;
    acc = j > y0+(y1-y0)/(x1-x0)*(k-x0);
  }
  if (iq==3) {
    x0=23-11.5;
    y0=23;
    x1=23;
    y1=23-11.5;
    acc = j < y0+(y1-y0)/(x1-x0)*(k-x0);
  }
  if (iq==2) {
    x0=1;
    y0=23-11.5;
    x1=11.5;
    y1=23;
    acc = j < y0+(y1-y0)/(x1-x0)*(k-x0);
  }

  return acc;
}
*/
//------------------------------------------------------------------------------

Long64_t calo::GetEntriesFast()
{
   if (!fChain) return 0;
   return fChain->GetEntriesFast();
}

#endif // #ifdef calo_cxx
