//#include "create_input_histos.h"

#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TString.h"
#include <sstream>
#include "TH1.h"
#include "TF1.h"
#include "TString.h"
#include "Environment.h"
#include <iostream>
#include <fstream>
#include <TMath.h>
#include "TKey.h"
#include "TLegend.h"
#include "TString.h"
#include "TGraphErrors.h"
#include "TIterator.h"

//using namespace RooFit ;
using namespace std;
bool debug = 0;

void create_input_histos(int width_signal_region){

  //TAG: Change to file name 
  TFile* infile; 

  //get the backgrounds
  //names of the samples

  //change to names of directories / bkgs in datacard

  TString rootfilenames[] = {
"/net/scratch_cms/institut_3a/mukherjee/LFV/ttbar_tot.root",
"/net/scratch_cms/institut_3a/mukherjee/LFV/WW_tot.root",
"/net/scratch_cms/institut_3a/mukherjee/LFV/SingleTop_tot.root",
"/net/scratch_cms/institut_3a/mukherjee/LFV/DY_tot.root",
"/net/scratch_cms/institut_3a/mukherjee/LFV/WZ_tot.root",
"/net/scratch_cms/institut_3a/mukherjee/LFV/ZZ_tot.root",
"/net/scratch_cms/institut_3a/mukherjee/LFV/Wjet_QCD_tot.root"
};

  TString sample_names[] = {"TT_tot","WW_tot","single_top_tot","DY_tot","WZ_tot","ZZ_tot","Wjet_tot"};
  const int arraySize = sizeof(sample_names)/sizeof(sample_names[0]);  

  //names of systematics for shape-based input
  TString syst_names[] ={"Ele_syst_Scale","Muon_syst_Scale","Muon_syst_Resolution"};
  const int arraySize_systs = sizeof(syst_names)/sizeof(syst_names[0]);    

  std::cout << "arraySize_systs " << arraySize_systs << " arraySize " << arraySize << std::endl;

  Char_t file_title[400];
  Char_t dir_title[400];  
  Char_t dir_title_syst[400];  

  //outfile for all histos
  TFile *outfile = new TFile("root_out/out.root","recreate");

  //outfile for all text to be used as input for the limit cfgs
  ofstream myfile;
  myfile.open ("txt_out/normalization.txt");

  double mass_min=200;
  double mass_max=6500;

  int step_size=100;

  int N_points=(int)(mass_max-mass_min)/(double)step_size;

  double masses[400];
  double num_total[400];

  double num_signal[400];
  double num_observed[400];

  if (debug) std::cout << "will start initializer loop" << std::endl;

  for(int n=0; n<400; n++)
  {
    masses[n]=0.;
    num_total[n]=0.;
    num_signal[n]=0.;
    num_observed[n]=0.;
  }  

  TF1* fit_resolution=new TF1("fit_resolution","[0]+[1]*x+[2]*x*x+[3]*x*x*x",mass_min,mass_max);
  fit_resolution->SetParameter(0,0.012830);
  fit_resolution->SetParameter(1,0.000015);  
  fit_resolution->SetParameter(2,0.0);
  fit_resolution->SetParameter(3,0.0);   
 
  //###############
  //loop over backgrounds
  //###############  
  if (debug) std::cout << "will start background loop"<< std::endl;
  
  for (int i = 0; i < arraySize; ++i){
    if (debug) std::cout << "will get the bkg root file"<< std::endl;
    std::cout << "rootfilename " << rootfilenames[i] << "  sample_name " << sample_names[i] << std::endl;
    TFile* infile = new TFile(rootfilenames[i]);
    
    //TAG change to one input file per bkg
    TH1F* hist_ori;
    //    cout << ((TString)(sample_names[i]) + "_7_0/"+(TString)(sample_names[i])+ "_ori").Data() << endl;
    if (debug) std::cout << "will get the hist h1_0_emu_Mass"<< std::endl;
    
    hist_ori=(TH1F*)infile->Get("emu/Stage_0/h1_0_emu_Mass");
    hist_ori->SetName(sample_names[i]);
    
    outfile->cd();
    hist_ori->Write(sample_names[i]);
    
    //TAG check which binning is used for the histograms (here: 1-10000, 1 GeV binning)
    myfile << "bkg " << sample_names[i] << " " << hist_ori->Integral(1,6000) << "\n";
    
    //###############
    //loop over systematics
    //###############  
    if (debug) std::cout << "will start the syst loop"<< std::endl;
    
    
    for(int k=0; k<arraySize_systs; k++)
      {
	TH1F* hist_syst_up;
	TH1F* hist_syst_down;	
	
	//hist_syst_up=(TH1F*)infile->Get(((TString)(sample_names[i]) + "_7_0/"+(TString)(sample_names[i])+ "_" + syst_names[k] + "Up").Data());
	if (debug) std::cout << "will get syst hist UP"<< std::endl;
	
	hist_syst_up   = (TH1F*)infile->Get("emu/Stage_0/sys/h1_0_emu_Mass_"+syst_names[k]+"Up");
	
	if (debug) std::cout << "will get syst hist DOWN" << std::endl;
	hist_syst_down = (TH1F*)infile->Get("emu/Stage_0/sys/h1_0_emu_Mass_"+syst_names[k]+"Down");
	if (debug) std::cout << "Successfully  got syst hist DOWN"<< std::endl;
	
	if (debug) std::cout << "Will set name up" << std::endl;
	if (debug) std::cout << "hist_syst_up " << hist_syst_up << std::endl;
	
	hist_syst_up->SetName(sample_names[i]+"_"+syst_names[k]+"Up");
	if (debug) std::cout << "Successfully  set name up, will write hist up" << std::endl;
	hist_syst_up->Write(); 
	
	if (debug) std::cout << "Successfully  written hist_syst_up"<< std::endl;
	
	myfile << "bkg " << sample_names[i]+"_"+syst_names[k]+"Up" << " " << hist_syst_up->Integral(1,6000) << "\n";
	hist_syst_down->SetName(sample_names[i]+"_"+syst_names[k]+"Down");
	hist_syst_down->Write(); 
	myfile << "bkg " << sample_names[i]+"_"+syst_names[k]+"Down" << " " << hist_syst_down->Integral(1,6000) << "\n";
	
	if (debug) std::cout << "Will delete hist_syst_up hist_syst_down "<< std::endl;
	
	delete hist_syst_up;
	delete hist_syst_down;	
	if (debug) std::cout << "Successfully deleted hist_syst_up hist_syst_down "<< std::endl;
	
      }
    
    if (debug) std::cout << "Will delete hist_ori "<< std::endl;
    delete hist_ori;
    if (debug) std::cout << "Successfully deleted hist_ori "<< std::endl;
    
    
  }
 

  //###############
  //end loop over backgrounds
  //###############    

  //get datahist

  //TAG get the file with the data histogram
  if (debug) std::cout << "will get data root file"<< std::endl;
  TFile* data_file = new TFile("/net/scratch_cms/institut_3a/mukherjee/LFV/allData.root");
  TH1F* data;
  //  data=(TH1F*)data_file->Get("h1_inv_mass_1mu_1tau_aligned_7_0");
  if (debug) std::cout << "will get data hist"<< std::endl;
  data=(TH1F*)data_file->Get("emu/Stage_0/h1_0_emu_Mass"); 

  outfile->cd();
  data->SetName("data_obs");
  data->Write();

  myfile << "data " << "data" << " " << data->Integral(1,6000) << "\n";

  //###############
  //prepare input for the individual mass points
  //###############

  //exit(0);

  const int num_samples=(int)((mass_max-mass_min)/(double)step_size+1);
  
  double mass_sig=(double)mass_min;
  double resolution=5.;
  //double width=1.;
  
  TF1* gauss = new TF1("gauss","TMath::Gaus(x,[0],[1])",mass_sig-8*resolution,mass_sig+8*resolution);

  TF1* fit_acceff=new TF1("fit_acceff","[0]+[1]/(x+[2])+[3]*x",300.,6000.);
 
  fit_acceff->SetParameter(0,0.8248);
  fit_acceff->SetParameter(1,-161.9);  
  fit_acceff->SetParameter(2,143.777824);
  fit_acceff->SetParameter(3,-0.000025); 

  double acceff=0.7;

  double xsec=1.0;

  double limit_lower=mass_min;
  double limit_upper=mass_max;  

  //TAG change according to the binning used in the histograms
  int bin_limit_lower=1;
  int bin_limit_upper=6000;  

  for(int k=0;k<num_samples;k++)
    {
      TH1F* signal_temp;
      
      sprintf(file_title,"root_out/out_mass_%d.root",(int)mass_sig);

      sprintf(dir_title,"txt_out/normalization_Mass_%d_input_histos.txt",(int)mass_sig);

      ofstream myfile_temp;
     
      //cout << "Mass point: " << mass_sig << endl;

      myfile_temp.open(dir_title);

      TFile* outfile_signal = new TFile(file_title,"recreate");
      outfile_signal->cd();

      resolution=mass_sig*(fit_resolution->Eval(mass_sig));
      acceff=fit_acceff->Eval(mass_sig);

      limit_lower=mass_sig-width_signal_region*resolution;
      limit_upper=mass_sig+width_signal_region*resolution;
      
      //cout << "limit low: " << limit_lower << endl;

      //if(mass_sig>=800.)limit_upper=9950.;
      bin_limit_lower=(int)limit_lower;
      bin_limit_upper=(int)limit_upper;	

      int counter_histos=0;

      TKey *key;
      TIter next(outfile->GetListOfKeys());

      while ((key = (TKey*)next())) {

	TH1 *h1 = (TH1*)key->ReadObj();
	TString *name = new TString(h1->GetName());
	
        //cout << "TString " << name->Data() << endl;

	/*
	for(int n=1;n<(bin_limit_lower);n++)
	  {
	    h1->SetBinContent(n,0.);
	    h1->SetBinError(n,0.);	    
	  }
	for(int m=(bin_limit_upper+1);m<10000;m++)
	  {
	    h1->SetBinContent(m,0.);
	    h1->SetBinError(m,0.);	    
	  }	
	*/

	//cout << h1->GetName() << endl;

	if(strcmp(h1->GetName(),"data_obs") ==0 )
	  {
	    myfile_temp << "data " << h1->GetName() << " " << h1->Integral(1,6000) << "\n";	    
	  }	
	else
	  {
	    myfile_temp << "bkg " << h1->GetName() << " " << h1->Integral(1,6000) << "\n"; 
	  }
	
	h1->Write();

	counter_histos++;

	delete name;
	delete h1;

      }


      //set the signal PDFs here

      //TAG get signal cross section
      sprintf(file_title,"LQD_001_LLE_001_MSnl_scale_down_%d",(int)mass_sig);
      get_environment(file_title);
      xsec=BGcrosssection;

      gauss->SetParameter(0,mass_sig);
      gauss->SetParameter(1,resolution);      
	  
      signal_temp=new TH1F("signal_temp","",6000,0.,6000.); 	  
      signal_temp->FillRandom("gauss",10000);
      signal_temp->Sumw2();
      signal_temp->Scale(xsec*acceff*lumi_scale/10000.);
	  
      cout << "signal mass: " << mass_sig << " N events: " << signal_temp->Integral(1,6000) << endl; 

      outfile_signal->cd();
      signal_temp->Write("signal");

      myfile << "signal " << mass_sig << " " << signal_temp->Integral(1,6000) << "\n";
      myfile_temp << "signal " << mass_sig << " " << signal_temp->Integral(1,6000) << "\n";

      masses[k]=mass_sig;

      mass_sig += step_size;   
      delete signal_temp;
  
      myfile_temp.close();

    }


      //close the output txt file:
      myfile.close();
    
      delete gauss;

}
