#include "TSystem.h"
#include "TString.h"
#include "TEnv.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cstring>
#include <cstdlib>

using namespace std;

Double_t lumi_scale;
Int_t bgno;

  double BGnumberofevents;
  double BGcrosssection;
  double BGweight;

void errormsg(TString var) {
    cout << endl;
    cout << "  ERROR : Environment variable(s) not properly set or missing : " << var << endl;
    cout << "          Exiting." << endl;
    cout << endl;
    std::exit(0);
}

void get_environment(TString background) {

  TEnv *xs = new TEnv("INPUT");

  TString *base = new TString("/user/mukherjee/limits_LFV/scripts/xs.cfg");
  //TString *base = new TString("/home/home1/institut_3a/gueth/ANALYSIS_skimmer_v104/macro/SirPlotAlot_thomas/SirPlotAlot/xs.cfg");

  const char* ret;
  ret = gSystem->FindFile(".", *base);

  /*
  cout << *base << endl;
  cout << ret << endl;
  */

  /*
  if (!strstr(ret, *base)) {
    cout << endl;
    cout << "  ERROR : Configuration file '" << *bbase << "' not found in current path !" << endl;
    cout << "          Exiting." << endl;
    cout << endl;
    exit(0);
  }  
  */

  xs->ReadFile("/user/mukherjee/limits_LFV/scripts/xs.cfg", kEnvChange); 
  //xs->ReadFile("/home/home1/institut_3a/gueth/ANALYSIS_skimmer_v104/macro/SirPlotAlot_thomas/SirPlotAlot/xs.cfg", kEnvChange);

  lumi_scale = xs->GetValue("lumi", -1.);
  if (lumi_scale<0) errormsg("lumi"); 

  //  BGnumberofevents = xs->GetValue(background+".Nev", -1.);
  BGcrosssection   = xs->GetValue(background+".xs", -1.);
  BGweight         = xs->GetValue(background+".weight", -1.);    


  // if (BGnumberofevents<0 || BGcrosssection<0 || BGweight<0)
  if (BGcrosssection<0 || BGweight<0)                                                                                         
    {errormsg(background);}

  delete base;

}
