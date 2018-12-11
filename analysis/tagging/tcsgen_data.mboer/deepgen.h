//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Feb  6 14:53:01 2018 by ROOT version 6.10/02
// from TTree TCS_Tree/TCS generated events
// found on file: RootFiles/DEEPGen_10000.root
//////////////////////////////////////////////////////////

#ifndef deepgen_h
#define deepgen_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class deepgen {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        ALV_minus_lab[4];
   Double_t        ALV_plus_lab[4];
   Double_t        ALV_gamma_in[4];
   Double_t        ALV_el_in[4];
   Double_t        ALV_Recoil_lab[4];
   Double_t        ALV_el_out[4];
   Double_t        Q2;
   Double_t        theta_gamma;
   Double_t        theta_beam;
   Double_t        phi_beam;
   Double_t        Psi_s;
   Double_t        yy;
   Double_t        Qp2;
   Double_t        tt;
   Double_t        RotAxis;
   Double_t        ttmin;
   Double_t        Phi_CMV;
   Double_t        Theta_CMV;
   Double_t        Egamma;
   Float_t         cross_tot_unpol;
   Float_t         cross_BH;
   Float_t         cross_TCS;
   Float_t         W_tot_unpol;
   Float_t         W_BH;
   Float_t         W_TCS;
   Int_t           VirtualFlag;
   Double_t        Flux_qr;
   Double_t        Flux_bmr;
   Double_t        epsilon;
   Float_t         BSA;
   Float_t         TSA;
   Float_t         BTSA;
   Long64_t        EventNumber;
   Long64_t        TrueEventNumber;
   Float_t         cross_tot_pol;
   Float_t         W_tot_pol;
   Float_t         cross_tot_pol_beam;
   Float_t         cross_tot_pol_beam_circ;
   Float_t         cross_tot_pol_beam_lin;
   Float_t         W_tot_pol_beam;
   Float_t         W_tot_pol_beam_circ;
   Float_t         W_tot_pol_beam_lin;
   Float_t         cross_tot_pol_target;
   Float_t         W_tot_pol_target;
   Double_t        phi_s;
   Double_t        theta_s;
   Int_t           target_spindir;
   Int_t           beam_spindir;
   Double_t        poltargetdeg;
   Double_t        polbeamdeg_circ;
   Double_t        polbeamdeg_lin;
   Double_t        polbeamdeg;
   Float_t         thetamin_nocut;
   Int_t           FlagSing;

   // List of branches
   TBranch        *b_ALV_minus_lab;   //!
   TBranch        *b_ALV_plus_lab;   //!
   TBranch        *b_ALV_gamma_in;   //!
   TBranch        *b_ALV_el_in;   //!
   TBranch        *b_ALV_Recoil_lab;   //!
   TBranch        *b_ALV_el_out;   //!
   TBranch        *b_Q2;   //!
   TBranch        *b_theta_gamma;   //!
   TBranch        *b_theta_beam;   //!
   TBranch        *b_phi_beam;   //!
   TBranch        *b_Psi_s;   //!
   TBranch        *b_yy;   //!
   TBranch        *b_Qp2;   //!
   TBranch        *b_tt;   //!
   TBranch        *b_RotAxis;   //!
   TBranch        *b_ttmin;   //!
   TBranch        *b_Phi_CMV;   //!
   TBranch        *b_Theta_CMV;   //!
   TBranch        *b_Egamma;   //!
   TBranch        *b_cross_tot_unpol;   //!
   TBranch        *b_cross_BH;   //!
   TBranch        *b_cross_TCS;   //!
   TBranch        *b_W_tot_unpol;   //!
   TBranch        *b_W_BH;   //!
   TBranch        *b_W_TCS;   //!
   TBranch        *b_VirtualFlag;   //!
   TBranch        *b_Flux_qr;   //!
   TBranch        *b_Flux_bmr;   //!
   TBranch        *b_epsilon;   //!
   TBranch        *b_BSA;   //!
   TBranch        *b_TSA;   //!
   TBranch        *b_BTSA;   //!
   TBranch        *b_EventNumber;   //!
   TBranch        *b_TrueEventNumber;   //!
   TBranch        *b_cross_tot_pol;   //!
   TBranch        *b_W_tot_pol;   //!
   TBranch        *b_cross_tot_pol_beam;   //!
   TBranch        *b_cross_tot_pol_beam_circ;   //!
   TBranch        *b_cross_tot_pol_beam_lin;   //!
   TBranch        *b_W_tot_pol_beam;   //!
   TBranch        *b_W_tot_pol_beam_circ;   //!
   TBranch        *b_W_tot_pol_beam_lin;   //!
   TBranch        *b_cross_tot_pol_target;   //!
   TBranch        *b_W_tot_pol_target;   //!
   TBranch        *b_phi_s;   //!
   TBranch        *b_theta_s;   //!
   TBranch        *b_target_spindir;   //!
   TBranch        *b_beam_spindir;   //!
   TBranch        *b_poltargetdeg;   //!
   TBranch        *b_polbeamdeg_circ;   //!
   TBranch        *b_polbeamdeg_lin;   //!
   TBranch        *b_polbeamdeg;   //!
   TBranch        *b_thetamin_nocut;   //!
   TBranch        *b_FlagSing;   //!

   deepgen(TTree *tree=0);
   deepgen(string root_file, TTree *tree=0);
   virtual ~deepgen();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef deepgen_cxx
deepgen::deepgen(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("RootFiles/DEEPGen_10000.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("RootFiles/DEEPGen_10000.root");
      }
      f->GetObject("TCS_Tree",tree);

   }
   Init(tree);
}

deepgen::deepgen(string root_file, TTree *tree) : fChain(0) 
{
  std::cout << "Root file " << root_file << std::endl;
// if parameter tree is not specified (or zero), connect the file.
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(root_file.c_str());
    if (!f || !f->IsOpen()) {
      f = new TFile(root_file.c_str());
    }
    f->GetObject("TCS_Tree",tree);

  }
  Init(tree);
}

deepgen::~deepgen()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t deepgen::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t deepgen::LoadTree(Long64_t entry)
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

void deepgen::Init(TTree *tree)
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

   fChain->SetBranchAddress("ALV_minus_lab", ALV_minus_lab, &b_ALV_minus_lab);
   fChain->SetBranchAddress("ALV_plus_lab", ALV_plus_lab, &b_ALV_plus_lab);
   fChain->SetBranchAddress("ALV_gamma_in", ALV_gamma_in, &b_ALV_gamma_in);
   fChain->SetBranchAddress("ALV_el_in", ALV_el_in, &b_ALV_el_in);
   fChain->SetBranchAddress("ALV_Recoil_lab", ALV_Recoil_lab, &b_ALV_Recoil_lab);
   fChain->SetBranchAddress("ALV_el_out", ALV_el_out, &b_ALV_el_out);
   fChain->SetBranchAddress("Q2", &Q2, &b_Q2);
   fChain->SetBranchAddress("theta_gamma", &theta_gamma, &b_theta_gamma);
   fChain->SetBranchAddress("theta_beam", &theta_beam, &b_theta_beam);
   fChain->SetBranchAddress("phi_beam", &phi_beam, &b_phi_beam);
   fChain->SetBranchAddress("Psi_s", &Psi_s, &b_Psi_s);
   fChain->SetBranchAddress("yy", &yy, &b_yy);
   fChain->SetBranchAddress("Qp2", &Qp2, &b_Qp2);
   fChain->SetBranchAddress("tt", &tt, &b_tt);
   fChain->SetBranchAddress("RotAxis", &RotAxis, &b_RotAxis);
   fChain->SetBranchAddress("ttmin", &ttmin, &b_ttmin);
   fChain->SetBranchAddress("Phi_CMV", &Phi_CMV, &b_Phi_CMV);
   fChain->SetBranchAddress("Theta_CMV", &Theta_CMV, &b_Theta_CMV);
   fChain->SetBranchAddress("Egamma", &Egamma, &b_Egamma);
   fChain->SetBranchAddress("cross_tot_unpol", &cross_tot_unpol, &b_cross_tot_unpol);
   fChain->SetBranchAddress("cross_BH", &cross_BH, &b_cross_BH);
   fChain->SetBranchAddress("cross_TCS", &cross_TCS, &b_cross_TCS);
   fChain->SetBranchAddress("W_tot_unpol", &W_tot_unpol, &b_W_tot_unpol);
   fChain->SetBranchAddress("W_BH", &W_BH, &b_W_BH);
   fChain->SetBranchAddress("W_TCS", &W_TCS, &b_W_TCS);
   fChain->SetBranchAddress("VirtualFlag", &VirtualFlag, &b_VirtualFlag);
   fChain->SetBranchAddress("Flux_qr", &Flux_qr, &b_Flux_qr);
   fChain->SetBranchAddress("Flux_bmr", &Flux_bmr, &b_Flux_bmr);
   fChain->SetBranchAddress("epsilon", &epsilon, &b_epsilon);
   fChain->SetBranchAddress("BSA", &BSA, &b_BSA);
   fChain->SetBranchAddress("TSA", &TSA, &b_TSA);
   fChain->SetBranchAddress("BTSA", &BTSA, &b_BTSA);
   fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
   fChain->SetBranchAddress("TrueEventNumber", &TrueEventNumber, &b_TrueEventNumber);
   fChain->SetBranchAddress("cross_tot_pol", &cross_tot_pol, &b_cross_tot_pol);
   fChain->SetBranchAddress("W_tot_pol", &W_tot_pol, &b_W_tot_pol);
   fChain->SetBranchAddress("cross_tot_pol_beam", &cross_tot_pol_beam, &b_cross_tot_pol_beam);
   fChain->SetBranchAddress("cross_tot_pol_beam_circ", &cross_tot_pol_beam_circ, &b_cross_tot_pol_beam_circ);
   fChain->SetBranchAddress("cross_tot_pol_beam_lin", &cross_tot_pol_beam_lin, &b_cross_tot_pol_beam_lin);
   fChain->SetBranchAddress("W_tot_pol_beam", &W_tot_pol_beam, &b_W_tot_pol_beam);
   fChain->SetBranchAddress("W_tot_pol_beam_circ", &W_tot_pol_beam_circ, &b_W_tot_pol_beam_circ);
   fChain->SetBranchAddress("W_tot_pol_beam_lin", &W_tot_pol_beam_lin, &b_W_tot_pol_beam_lin);
   fChain->SetBranchAddress("cross_tot_pol_target", &cross_tot_pol_target, &b_cross_tot_pol_target);
   fChain->SetBranchAddress("W_tot_pol_target", &W_tot_pol_target, &b_W_tot_pol_target);
   fChain->SetBranchAddress("phi_s", &phi_s, &b_phi_s);
   fChain->SetBranchAddress("theta_s", &theta_s, &b_theta_s);
   fChain->SetBranchAddress("target_spindir", &target_spindir, &b_target_spindir);
   fChain->SetBranchAddress("beam_spindir", &beam_spindir, &b_beam_spindir);
   fChain->SetBranchAddress("poltargetdeg", &poltargetdeg, &b_poltargetdeg);
   fChain->SetBranchAddress("polbeamdeg_circ", &polbeamdeg_circ, &b_polbeamdeg_circ);
   fChain->SetBranchAddress("polbeamdeg_lin", &polbeamdeg_lin, &b_polbeamdeg_lin);
   fChain->SetBranchAddress("polbeamdeg", &polbeamdeg, &b_polbeamdeg);
   fChain->SetBranchAddress("thetamin_nocut", &thetamin_nocut, &b_thetamin_nocut);
   fChain->SetBranchAddress("FlagSing", &FlagSing, &b_FlagSing);
   Notify();
}

Bool_t deepgen::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void deepgen::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t deepgen::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef deepgen_cxx
