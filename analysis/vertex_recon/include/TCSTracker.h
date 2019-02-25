//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Feb  6 19:36:52 2019 by ROOT version 6.02/08
// from TChain tracker/
//////////////////////////////////////////////////////////

#ifndef TCSTracker_h
#define TCSTracker_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

#include <TVector3.h>
#include <TMath.h>
#include <TRandom.h>

#define TRACKER1_DIST 1200.   //mm, consistent with G4 coding.
#define TRACKER2_DIST 1300.
#define TRACKER3_DIST 1400.

#define TILT_ANGLE 13.835     //deg, consistent with G4 coding.
#define ROT_ANGLE  10.034

#define GEM_ACCURACY 100.E-3   //100um

struct hit {
  double x;
  double y;
  int det;
  int layer;
};

class TCSTracker {
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

   TCSTracker(TTree *tree=0);
   virtual ~TCSTracker();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   ///   virtual void     Loop();
   virtual Bool_t   Notify(string comment);
   virtual void     Show(Long64_t entry = -1);

   TCSTracker(int runlo, int runhi, TTree *tree=0);
   virtual Long64_t GetEntries();
   bool GoodTrack(int ientry, int particle_id, int prong_id, int quarter,
		  TVector3 &vdir);
   double FitStraightTrack(const vector<hit> &hitlist, const int axis,
			   double &offset, double &slope);
   void RotateToLab(TVector3 &v, const int ndet);
   double GEMSmear(const double xtrack);
   TRandom RandomGen;
};

#endif

#ifdef TCSTracker_cxx
TCSTracker::TCSTracker(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
     cout << "TCSTracker::TCSTracker: SINGLE_TREE defined" << endl;
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("tracker",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      cout << "TCSTracker::TCSTracker: SINGLE_TREE not defined" << endl;
      TChain * chain = new TChain("tracker","");
      chain->Add("RootFiles.tracked/DEEPGen_1000_tracked.root/tracker");
      chain->Add("RootFiles.tracked/DEEPGen_1001_tracked.root/tracker");
      chain->Add("RootFiles.tracked/DEEPGen_1003_tracked.root/tracker");
      chain->Add("RootFiles.tracked/DEEPGen_1005_tracked.root/tracker");
      chain->Add("RootFiles.tracked/DEEPGen_1006_tracked.root/tracker");
      chain->Add("RootFiles.tracked/DEEPGen_1007_tracked.root/tracker");
      chain->Add("RootFiles.tracked/DEEPGen_1008_tracked.root/tracker");
      chain->Add("RootFiles.tracked/DEEPGen_1009_tracked.root/tracker");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
  ////  RandomGen.SetSeed(0);   //seed from comp. clock.
}

TCSTracker::~TCSTracker()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t TCSTracker::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TCSTracker::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify("Tree loaded.");
   }
   return centry;
}

void TCSTracker::Init(TTree *tree)
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
   Notify("Chain initialized.");
}

Bool_t TCSTracker::Notify(string comment)
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

  cout << comment << endl;

   return kTRUE;
}

void TCSTracker::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t TCSTracker::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef TCSTracker_cxx

//==============================================================================

Long64_t TCSTracker::GetEntries() {
  return fChain->GetEntries();
}

//-----------------------------------------------------------------------------

TCSTracker::TCSTracker(int runlo, int runhi, TTree *tree) : fChain(0) 
{

// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.

  if (tree == 0) {
    TChain * chain = new TChain("tracker","");
    for (int runno=runlo; runno<runhi; runno++) {
      string rootfile=Form("RootFiles.tracked/DEEPGen_%d_tracked.root",runno);
      if (FILE *file = fopen(rootfile.c_str(), "r")) {
	fclose(file);
	chain->Add(rootfile.c_str());
      }
    }
    tree = chain;
  }

  Init(tree);
}

//------------------------------------------------------------------------------

bool TCSTracker::GoodTrack( int ientry, int particle_id, int prong_id,
			    int quarter, TVector3 &vdir) {

  vector<hit> hitlist;

  //    cout << "Filling hitlist:" << endl;
  //    cout << " detcont size = " << detcont->size() << endl;
  //    cout << " pidcont size = " << pidcont->size() << endl;
  //    cout << " pidorigcont size = " << pidorigcont->size() << endl;
  //    cout << " trackidcont size = " << trackidcont->size() << endl;

  GetEntry(ientry);

  for (UInt_t j = 0; j < detcont->size(); ++j) {
    if (pidorigcont->at(j) == particle_id &&
	pidcont->at(j) == pidorigcont->at(j) && trackidcont->at(j) == prong_id
	&& detcont->at(j) == quarter) {
       hit goodHit = {GEMSmear(xcont->at(j)), GEMSmear(ycont->at(j)),
                       detcont->at(j), layercont->at(j)};

      hitlist.push_back(goodHit);
    }
  }

  //    cout << " hitlist size = " << hitlist.size() << endl;

  bool good_hitlist = false;

  for (uint j=0; j<hitlist.size(); j++) {
    if (hitlist.at(j).layer == 0) {
      good_hitlist = true;
      break;
    }
  }

  if (good_hitlist) {

    for (uint j=0; j<hitlist.size(); j++) {
      if (!good_hitlist) break;
      for (uint k=j+1; k<hitlist.size(); k++)
	if (hitlist.at(k).layer == hitlist.at(j).layer) {
	  good_hitlist = false;
	  break;
	}
    }

  }

  good_hitlist = good_hitlist && hitlist.size() > 1;

  if (good_hitlist) {

    double xslope = 0.;
    double xoffset = 0.;
    double xChi2 = FitStraightTrack(hitlist, 0, xoffset, xslope);

    double yslope = 0.;
    double yoffset = 0.;
    double yChi2 = FitStraightTrack(hitlist, 1, yoffset, yslope);

    //      cout << "hitlist size = " << hitlist.size() << endl;
    //      cout << "xslope  = " << xslope << endl;
    //      cout << "xoffset = " << xoffset << endl;
    //      cout << "xChi2   = " << xChi2 << endl;
    //      cout << "yslope  = " << yslope << endl;
    //      cout << "yoffset = " << yoffset << endl;
    //      cout << "yChi2   = " << yChi2 << endl;
    //      getchar();

    //      TVector3 voff(xoffset, yoffset, 0.);
    //      RotateToLab(voff, ndet);

    ////      TVector3 vdir(0., 0., 1.);
    vdir.SetXYZ(xslope, yslope, 1.);

    RotateToLab(vdir, quarter);
  }

  return good_hitlist;
}

//------------------------------------------------------------------------------

double TCSTracker::FitStraightTrack(const vector<hit> &hitlist, const int axis,
				    double &offset, double &slope) {

  const double zpos[] {TRACKER1_DIST, TRACKER2_DIST, TRACKER3_DIST};

  double sum_zx = 0.;
  double sum_x = 0.;
  double sum_z = 0.;
  double sum_zz = 0.;

  uint n = hitlist.size();
  for (uint i=0; i<n; i++) {
    double z = zpos[hitlist.at(i).layer];
    double x = (axis==0 ? hitlist.at(i).x : hitlist.at(i).y);
    sum_zx += z*x;
    sum_x += x;
    sum_z += z;
    sum_zz += z*z;
  }

  slope = (n*sum_zx - sum_x*sum_z)/(n*sum_zz - sum_z*sum_z);
  offset = (sum_x - sum_z*slope)/n;

  double x2 = 0.;
  for (uint i=0; i<n; i++) {
    double z = zpos[hitlist.at(i).layer];
    double x = (axis==0 ? hitlist.at(i).x : hitlist.at(i).y);
    x2 += (x-slope*z-offset)*(x-slope*z-offset)/(GEM_ACCURACY*GEM_ACCURACY);
  };
  x2 /= 2.;

  return x2;
}

//------------------------------------------------------------------------------

void TCSTracker::RotateToLab(TVector3 &v, const int ndet) {

  double xangle = 0.;
  double yangle = 0.;
  //  int xflip = 0;
  //  int yflip = 0;

  switch (ndet) {
    case 0:
      xangle = -TILT_ANGLE;
      yangle =  ROT_ANGLE;
      //      xflip = +1;
      //      yflip = +1;
      break;
    case 1:
      xangle = -TILT_ANGLE;
      yangle = -ROT_ANGLE;
      //      xflip = -1;
      //      yflip = +1;
      break;
    case 2:
      xangle =  TILT_ANGLE;
      yangle = -ROT_ANGLE;
      //      xflip = -1;
      //      yflip = -1;
      break;
    case 3:
      xangle =  TILT_ANGLE;
      yangle =  ROT_ANGLE;
      //      xflip = +1;
      //      yflip = -1;
      break;
  default:
    cout << "*** RotateToLab: wrong ndet = " << ndet << " ! ***" << endl;
  }

  //  cout << "RotateToLab:" << endl;
  //  cout << " ndet   = " << ndet << endl;
  //  cout << " xangle = " << xangle << endl;
  //  cout << " yangle = " << yangle << endl;
  //  cout << " xflip  = " << xflip << endl;
  //  cout << " yflip  = " << yflip << endl;

  //  cout << " vector before flip:" << endl;
  //  v.Print();

  //  v.SetX(xflip*v.X());
  //  v.SetY(yflip*v.Y());

  //  cout << " vector after flip:" << endl;
  //  v.Print();

  v.RotateX(xangle*TMath::DegToRad());

  //  cout << " vector after X-rotataion:" << endl;
  //  v.Print();

  v.RotateY(yangle*TMath::DegToRad());

  //  cout << " vector after XY-rotataion:" << endl;
  //  v.Print();
  //  getchar();
}

//------------------------------------------------------------------------------

double TCSTracker::GEMSmear(const double xtrack) {
  return RandomGen.Gaus(xtrack, GEM_ACCURACY);
  //  return RandomGen.Gaus(xtrack, 0);
}
