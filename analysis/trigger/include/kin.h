//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Oct  8 10:45:23 2020 by ROOT version 6.14/04
// from TTree kin/TCS kinematics
// found on file: RootFiles.tracked/DEEPGen_0-99_tracked.root
//////////////////////////////////////////////////////////

#ifndef kin_h
#define kin_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class kin {
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

   kin(TTree *tree=0);
   kin(string file_name, TTree *tree=0);
   virtual ~kin();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   //   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   double GetWeight(int ientry);
   double dsigdy_br(double y);
   double sig_br(double y);
};

#endif

#ifdef kin_cxx
kin::kin(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("RootFiles.tracked/DEEPGen_0-99_tracked.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("RootFiles.tracked/DEEPGen_0-99_tracked.root");
      }
      f->GetObject("kin",tree);

   }
   Init(tree);
}

kin::~kin()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t kin::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t kin::LoadTree(Long64_t entry)
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

void kin::Init(TTree *tree)
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

Bool_t kin::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void kin::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t kin::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
  const float Q2min = 4.;
  const float Q2max = 9.;
  const float tmin = -1.;
  const float tmax = 0.;
  const float tQ2max = 0.3;
  GetEntry(entry);
  if (Q2>Q2min && Q2<9. && t>tmin && t<tmax && -t/Q2<tQ2max)
    return 1;
  else
    return 0;
}

//==============================================================================

kin::kin(string file_name, TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
     TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(file_name.c_str());
      if (!f || !f->IsOpen()) {
	f = new TFile(file_name.c_str());
      }
      f->GetObject("kin",tree);

   }
   Init(tree);
}

//------------------------------------------------------------------------------

//TCS event weight according to Marie's prescriptions.
//The W_tot_unpol is replaced by crs_BH.
//Shall be devided by total number of sampled events to get rates.

double kin::GetWeight(int ientry) {

  const double luminosity_exp = 1.47e23*1.5e12; //cm^-2 s^-1, UVA tar.+ CPS beam

  const double phasespace_MC = 1.0*6.4*2.094*6.283;
  //  const int TrueEventNumber = 16700;               //approximate number
  const float TrueEventNumberFactor = 1.6700;

  const double Emax = 11.;  //GeV
  const double Emin = 5.5;  //GeV

  const double pb = 1.e-36;   //cm^2

  const double f_br = (Emax-Emin)/Emax/(sig_br(Emax/Emax)-sig_br(Emin/Emax));

  GetEntry(ientry);

  double BeamProfileRescale_bmr = dsigdy_br(Eg/Emax) * f_br;
  double weight = luminosity_exp * phasespace_MC * sin(the_cm) *
                  BeamProfileRescale_bmr * crs_BH * pb / TrueEventNumberFactor;

  //  cout << "f_br = " << f_br << endl;
  //  cout << "BeamProfileRescale_bmr = " << BeamProfileRescale_bmr << endl;
  //  cout << "the_cm = " << the_cm << endl;
  //  cout << "crs_BH = " << crs_BH << endl;
  //  cout << "weight = " << weight << endl;
  //  getchar();

  return weight;
}

//dsigdy reduced brem. cross section.

double kin::dsigdy_br(double y) {
  return 1./y*(4./3.-4./3*y+y*y);
}

double kin::sig_br(double y) {
  return 4./3.*log(y)-4./3.*y+y*y/2.;
}

#endif // #ifdef kin_cxx
