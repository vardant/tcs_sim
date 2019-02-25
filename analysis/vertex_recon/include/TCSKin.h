//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Feb  6 19:36:52 2019 by ROOT version 6.02/08
// from TChain kin/
//////////////////////////////////////////////////////////

#ifndef TCSKin_h
#define TCSKin_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

struct KinVar {
  Double_t        Q2;
  Double_t        t;
  Double_t        s;
  Double_t        xi;
  Double_t        tau;
  Double_t        eta;
  Double_t        phi_cm;
  Double_t        the_cm;
  Double_t        psf;
  Double_t        flux_factor;
  Double_t        crs_BH;
  Double_t        Eg;
  Double_t        pminus[4];
  Double_t        pplus[4];
  Double_t        precoil[4];
  Double_t        vertex[3];
};

// Header file for the classes stored in the TTree if any.

class TCSKin {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        Q2;
   Double_t        t;
   Double_t        s;
   Double_t        xi;
   Double_t        tau;
   Double_t        eta;
   Double_t        phi_cm;
   Double_t        the_cm;
   Double_t        psf;
   Double_t        flux_factor;
   Double_t        crs_BH;
   Double_t        Eg;
   Double_t        pminus[4];
   Double_t        pplus[4];
   Double_t        precoil[4];
   Double_t        vertex[3];

   // List of branches
   TBranch        *b_Q2;   //!
   TBranch        *b_t;   //!
   TBranch        *b_s;   //!
   TBranch        *b_xi;   //!
   TBranch        *b_tau;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_phi_cm;   //!
   TBranch        *b_the_cm;   //!
   TBranch        *b_psf;   //!
   TBranch        *b_flux_factor;   //!
   TBranch        *b_crs_BH;   //!
   TBranch        *b_Eg;   //!
   TBranch        *b_pminus;   //!
   TBranch        *b_pplus;   //!
   TBranch        *b_precoil;   //!
   TBranch        *b_vertex;   //!

   TCSKin(TTree *tree=0);
   virtual ~TCSKin();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   ///   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   TCSKin(int runlo, int runhi, TTree *tree=0);
   KinVar GetKinVar();
};

#endif

#ifdef TCSKin_cxx
TCSKin::TCSKin(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("kin",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("kin","");
      chain->Add("RootFiles.tracked/DEEPGen_1000_tracked.root/kin");
      chain->Add("RootFiles.tracked/DEEPGen_1001_tracked.root/kin");
      chain->Add("RootFiles.tracked/DEEPGen_1003_tracked.root/kin");
      chain->Add("RootFiles.tracked/DEEPGen_1005_tracked.root/kin");
      chain->Add("RootFiles.tracked/DEEPGen_1006_tracked.root/kin");
      chain->Add("RootFiles.tracked/DEEPGen_1007_tracked.root/kin");
      chain->Add("RootFiles.tracked/DEEPGen_1008_tracked.root/kin");
      chain->Add("RootFiles.tracked/DEEPGen_1009_tracked.root/kin");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

TCSKin::~TCSKin()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t TCSKin::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TCSKin::LoadTree(Long64_t entry)
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

void TCSKin::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Q2", &Q2, &b_Q2);
   fChain->SetBranchAddress("t", &t, &b_t);
   fChain->SetBranchAddress("s", &s, &b_s);
   fChain->SetBranchAddress("xi", &xi, &b_xi);
   fChain->SetBranchAddress("tau", &tau, &b_tau);
   fChain->SetBranchAddress("eta", &eta, &b_eta);
   fChain->SetBranchAddress("phi_cm", &phi_cm, &b_phi_cm);
   fChain->SetBranchAddress("the_cm", &the_cm, &b_the_cm);
   fChain->SetBranchAddress("psf", &psf, &b_psf);
   fChain->SetBranchAddress("flux_factor", &flux_factor, &b_flux_factor);
   fChain->SetBranchAddress("crs_BH", &crs_BH, &b_crs_BH);
   fChain->SetBranchAddress("Eg", &Eg, &b_Eg);
   fChain->SetBranchAddress("pminus", pminus, &b_pminus);
   fChain->SetBranchAddress("pplus", pplus, &b_pplus);
   fChain->SetBranchAddress("precoil", precoil, &b_precoil);
   fChain->SetBranchAddress("vertex", vertex, &b_vertex);
   Notify();
}

Bool_t TCSKin::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void TCSKin::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t TCSKin::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

//==============================================================================

KinVar TCSKin::GetKinVar() {

  KinVar k;

  k.Q2 = Q2;
  k.t = t;
  k.s = s;
  k.xi = xi;
  k.tau = tau;
  k.eta = eta;
  k.phi_cm = phi_cm;
  k.the_cm = the_cm;
  k.psf = psf;
  k.flux_factor = flux_factor;
  k.crs_BH = crs_BH;
  k.Eg = Eg;
  for (int i=0; i<4; i++) {
    k.pminus[i] = pminus[i];
    k.pplus[i] = pplus[i];
    k.precoil[i] = precoil[i];
    k.vertex[i] = vertex[i];
  }

  return k;
}

//-----------------------------------------------------------------------------

TCSKin::TCSKin(int runlo, int runhi, TTree *tree) : fChain(0) 
{

// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.

  if (tree == 0) {
    TChain * chain = new TChain("kin","");
    int nchain = 0;
    for (int runno=runlo; runno<runhi; runno++) {
      string rootfile=Form("RootFiles.tracked/DEEPGen_%d_tracked.root",runno);
      if (FILE *file = fopen(rootfile.c_str(), "r")) {
	fclose(file);
	chain->Add(rootfile.c_str(),0);
	//	cout << rootfile << " chained" << endl;
	nchain++;
      }
    }
    cout << "TCSKin::TCSKin: " << nchain << " root files chained." << endl;
    tree = chain;
  }

  Init(tree);
}

#endif // #ifdef TCSKin_cxx
