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
  std::cout << "Starting .... " << std::endl;
  TFile* infile; 

  TString rootfilenames[] = {
    "/net/scratch_cms/institut_3a/erdweg/public/13TeV_rpv/Jan_04/ttbar_tot.root",
    "/net/scratch_cms/institut_3a/erdweg/public/13TeV_rpv/Jan_04/WW_tot.root",
    "/net/scratch_cms/institut_3a/erdweg/public/13TeV_rpv/Jan_04/SingleTop_tot.root",
    "/net/scratch_cms/institut_3a/erdweg/public/13TeV_rpv/Jan_04/DY_tot.root",
    "/net/scratch_cms/institut_3a/erdweg/public/13TeV_rpv/Jan_04/WZ_tot.root",
    "/net/scratch_cms/institut_3a/erdweg/public/13TeV_rpv/Jan_04/ZZ_tot.root",
    "/net/scratch_cms/institut_3a/erdweg/public/13TeV_rpv/Jan_04/Wgamma.root",
    "/net/scratch_cms/institut_3a/erdweg/public/13TeV_rpv/Jan_04/QCD.root"
};

  //TString sample_names[] = {"TT_tot","WW_tot","single_top_tot","WZ_tot","ZZ_tot","datadriven"};
  TString sample_names[] = {"TT_tot","WW_tot","single_top_tot","DY_tot","WZ_tot","ZZ_tot","Wgamma","datadriven"};
  const int arraySize = sizeof(sample_names)/sizeof(sample_names[0]);  

  //names of systematics for shape-based input
  TString syst_names[] ={"pileup_syst_","eleID_syst_","Ele_syst_Scale","Muon_syst_Scale","Muon_syst_Resolution","muoID_syst_"};
  const int arraySize_systs = sizeof(syst_names)/sizeof(syst_names[0]);    

  std::cout << "No- of systematics " << arraySize_systs << "  No. of bkg MC samples " << arraySize << std::endl;

  Char_t file_title[400];
  Char_t dir_title[400];  
  Char_t dir_title_syst[400];  

  //outfile for all histos
  TFile *outfile = new TFile("root_out/out.root","recreate");

  //outfile for all text to be used as input for the limit cfgs
  ofstream myfile;
  myfile.open ("txt_out/normalization.txt");

  double mass_min=200;
  double mass_max=1000;

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
  fit_resolution->SetParameter(0,0.013);
  fit_resolution->SetParameter(1,1.5e-05);  
  fit_resolution->SetParameter(2,-3.8e-09);
  fit_resolution->SetParameter(3,4.4e-13);   

  
  //
  //  PDF Systematics
  //
  TString pdffilename = "/net/scratch_cms/institut_3a/erdweg/public/13TeV_rpv/Nov_24/for_limit.root";
  TFile* pdffile = new TFile(pdffilename);
  std::cout << "File for PDF uncert : " << pdffilename << std::endl;
  TH1D* hist_pdf_mean=(TH1D*)pdffile->Get("PDF4LHC15_nnlo_mc_mean");
  TH1D* hist_pdf_up=(TH1D*)pdffile->Get("PDF4LHC15_nnlo_mc_up");
  TH1D* hist_pdf_down=(TH1D*)pdffile->Get("PDF4LHC15_nnlo_mc_down");

  TH1D* hist_pdf_rel_up = new TH1D("reluppdf", "pdfuprel", 6200, 0, 6200);   //  = ((*hist_pdf_up)/(*hist_pdf_mean));
  TH1D* hist_pdf_rel_down = new TH1D("reldownpdf", "pdfdownrel", 6200, 0, 6200);//= ((*hist_pdf_down)/(*hist_pdf_mean));

  for (int i=1; i<6201; i++) {
    double up=hist_pdf_up->GetBinContent(i);
    double down=hist_pdf_down->GetBinContent(i);
    double mean=hist_pdf_mean->GetBinContent(i);
    double rel_up = 0;   
    double rel_down = 0; 
    if (mean!=0.0) rel_up =   up/mean;
    if (mean!=0.0) rel_down = down/mean;

    //std::cout << "bin " << i << "rel_up " << rel_up << "  rel_down " << rel_down << std::endl;
    hist_pdf_rel_up->SetBinContent(i,rel_up);
    hist_pdf_rel_down->SetBinContent(i,rel_down);
  }

  //###############
  //loop over backgrounds
  //###############  
  if (debug) std::cout << "will start background loop"<< std::endl;

  //////////////////////////////////////
  //              UPDATE              //
  double Lumi_bkg = (2464.0/1000.0);
  std::cout << "Lumi scale factor for bkg : " << Lumi_bkg << std::endl;
  /////////////////////////////////////
  

  double ttbar_kfact = 1.138;
  for (int i = 0; i < arraySize; ++i){
    if (debug) std::cout << "will get the bkg root file"<< std::endl;
    std::cout << "\n\nrootfilename " << rootfilenames[i] << "  sample_name " << sample_names[i] << std::endl;
    TFile* infile = new TFile(rootfilenames[i]);
    
    //TAG change to one input file per bkg
    TH1D* hist_ori;
    //    cout << ((TString)(sample_names[i]) + "_7_0/"+(TString)(sample_names[i])+ "_ori").Data() << endl;
    if (debug) std::cout << "will get the hist h1_0_emu_Mass"<< std::endl;
    
    hist_ori=(TH1D*)infile->Get("emu/Stage_0/h1_0_emu_Mass");
    //if (sample_names[i] != "datadriven") 
    hist_ori->Scale(Lumi_bkg);
    if (sample_names[i]=="TT_tot") {
      std::cout << "Extra k-factor scaling will be done for TTbar background" << std::endl;
      hist_ori->Scale(ttbar_kfact);
    }
      
    hist_ori->SetName(sample_names[i]);
    
    outfile->cd();
    hist_ori->Write(sample_names[i]);
    
    //TAG check which binning is used for the histograms (here: 1-10000, 1 GeV binning)
    std::cout  << "bkg " << sample_names[i] << " " << hist_ori->Integral(1,6000) << "\n";
    myfile << "bkg " << sample_names[i] << " " << hist_ori->Integral(1,6000) << "\n";
    
    
    if (debug) std::cout << "Now multiply mass_mean_hist_bkg with pdf_rel (up and down) hist" << std::endl;

    //if (sample_names[i] != "datadriven") {
      TH1D* hist_pdf_thisBkg_Up = new TH1D("bkgUPpdf", "BkgUpPDF", 6200, 0, 6200);   //  = hist_pdf_rel_up *(*hist_ori);
      TH1D* hist_pdf_thisBkg_Down = new TH1D("bkgDOWNpdf", "BkgDOWNPDF", 6200, 0, 6200); ; //= hist_pdf_rel_down *(*hist_ori);

      for (int ii=1; ii<6201; ii++) {
	double my_up=hist_pdf_rel_up->GetBinContent(ii);
	double my_down=hist_pdf_rel_down->GetBinContent(ii);
	double my_mean=hist_ori->GetBinContent(ii);
	double mul_my_up = my_up*my_mean;
	double mul_my_down = my_down*my_mean;

	hist_pdf_thisBkg_Up->SetBinContent(ii,mul_my_up);
	hist_pdf_thisBkg_Down->SetBinContent(ii,mul_my_down);
      }
      
      //  hist_pdf_thisBkg_Up.Scale( (hist_ori->Integral()) / (hist_pdf_thisBkg_Up.Integral()) );
      hist_pdf_thisBkg_Up->SetName(sample_names[i]+"_"+"pdf_syst"+"Up");
      hist_pdf_thisBkg_Up->Write();
      std::cout << "MEAN hist classname " << hist_ori->ClassName() << std::endl;
      std::cout << "UP hist classname " << hist_pdf_thisBkg_Up->ClassName() << std::endl;
      myfile << "bkg " << sample_names[i]+"_"+"pdf_syst"+"Up" << " " << hist_pdf_thisBkg_Up->Integral(1,6000) << "\n";
      
      // hist_pdf_thisBkg_Down.Scale( (hist_ori->Integral()) / (hist_pdf_thisBkg_Down.Integral()) );
      hist_pdf_thisBkg_Down->SetName(sample_names[i]+"_"+"pdf_syst"+"Down");
      hist_pdf_thisBkg_Down->Write(); 
      myfile << "bkg " << sample_names[i]+"_"+"pdf_syst"+"Down" << " " << hist_pdf_thisBkg_Down->Integral(1,6000) << "\n";
      //}
    //###############
    //loop over other systematics
    //###############  
    if (debug) std::cout << "will start the syst loop"<< std::endl;
    
    //if (sample_names[i] != "datadriven") {
      for(int k=0; k<arraySize_systs; k++)
	{
	  TH1D* hist_syst_up;
	  TH1D* hist_syst_down;	
	  
	  //hist_syst_up=(TH1F*)infile->Get(((TString)(sample_names[i]) + "_7_0/"+(TString)(sample_names[i])+ "_" + syst_names[k] + "Up").Data());
	  if (debug) std::cout << "will get syst hist UP"<< std::endl;
	  
	  hist_syst_up   = (TH1D*)infile->Get("emu/Stage_0/sys/h1_0_emu_Mass_"+syst_names[k]+"Up");
	  std::cout << "Lumi scaling for " << sample_names[i] << " and " << syst_names[k] << " UP"  << std::endl;
	  hist_syst_up->Scale(Lumi_bkg);
	  hist_syst_down = (TH1D*)infile->Get("emu/Stage_0/sys/h1_0_emu_Mass_"+syst_names[k]+"Down");
	  std::cout << "Lumi scaling for " << sample_names[i]  << " and " << syst_names[k] << " DOWN"  << std::endl;
	  hist_syst_down->Scale(Lumi_bkg);
	  hist_syst_up->SetName(sample_names[i]+"_"+syst_names[k]+"Up");
	  hist_syst_up->Write(); 
	  myfile << "bkg " << sample_names[i]+"_"+syst_names[k]+"Up" << " " << hist_syst_up->Integral(1,6000) << "\n";
	  
	  hist_syst_down->SetName(sample_names[i]+"_"+syst_names[k]+"Down");
	  hist_syst_down->Write(); 
	  myfile << "bkg " << sample_names[i]+"_"+syst_names[k]+"Down" << " " << hist_syst_down->Integral(1,6000) << "\n";
	  
	  if (debug) std::cout << "Will delete hist_syst_up hist_syst_down "<< std::endl;
	  
	  delete hist_syst_up;
	  delete hist_syst_down;	
	if (debug) std::cout << "Successfully deleted hist_syst_up hist_syst_down "<< std::endl;
	
      }
      // }
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
  TFile* data_file = new TFile("/net/scratch_cms/institut_3a/erdweg/public/13TeV_rpv/Nov_24/allData.root");
  TH1D* data;
  //  data=(TH1F*)data_file->Get("h1_inv_mass_1mu_1tau_aligned_7_0");
  if (debug) std::cout << "will get data hist"<< std::endl;
  data=(TH1D*)data_file->Get("emu/Stage_0/h1_0_emu_Mass"); 

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

  TF1* fit_acceff=new TF1("fit_acceff","[0]+[1]/(x+[2])+[3]*x",200.,6000.);
 
  fit_acceff->SetParameter(0,0.8184);
  fit_acceff->SetParameter(1,-124.6);  
  fit_acceff->SetParameter(2,100.090490);
  fit_acceff->SetParameter(3,-0.000024); 

  double acceff=0.7;

  double xsec=1.0;

  double limit_lower=mass_min;
  double limit_upper=mass_max;  

  //TAG change according to the binning used in the histograms
  int bin_limit_lower=1;
  int bin_limit_upper=6200;  

  for(int k=0;k<num_samples;k++)
    {
      TH1D*  signal_temp;
      TH1D*  signal_temp_massScale_up;  
      TH1D*  signal_temp_massScale_down;
      TH1D*  signal_temp_massRes_up;
      TH1D*  signal_temp_massRes_down;

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
	for(int m=(bin_limit_upper+1);m<6000;m++)
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
      sprintf(file_title,"LQD_01_LLE_01_MSnl_scale_down_%d",(int)mass_sig);
      get_environment(file_title);
      xsec=BGcrosssection;

      gauss->SetParameter(0,mass_sig);
      gauss->SetParameter(1,resolution);      
      signal_temp=new TH1D("signal_temp","",6200,0.,6200.); 	  
      signal_temp->FillRandom("gauss",10000);
      signal_temp->Sumw2();
      signal_temp->Scale((xsec/1000.0)*acceff*lumi_scale/10000.);
      myfile << "signal " << mass_sig << " " << signal_temp->Integral(1,6000) << "\n";
      /*
      for(int n=1;n<(bin_limit_lower);n++)
	{
	  signal_temp->SetBinContent(n,0.);
	  signal_temp->SetBinError(n,0.);
	}
      for(int m=(bin_limit_upper+1);m<6000;m++)
	{
	  signal_temp->SetBinContent(m,0.);
	  signal_temp->SetBinError(m,0.);
	}
      */
      cout << "signal mass: " << mass_sig << " N events: " << signal_temp->Integral(1,6000) << endl; 

      // Signal mass scale systematic //
      gauss->SetParameter(0,mass_sig*(1.0+0.02));
      gauss->SetParameter(1,resolution);
      signal_temp_massScale_up=new TH1D("signal_massScale_systUp","",6200,0.,6200.);
      signal_temp_massScale_up->FillRandom("gauss",10000);
      signal_temp_massScale_up->Sumw2();
      signal_temp_massScale_up->Scale((xsec/1000.0)*acceff*lumi_scale/10000.);
      myfile << "signal syst " << mass_sig << " signal_massScale_systUp " << signal_temp_massScale_up->Integral(1,6000) << "\n";

      gauss->SetParameter(0,mass_sig*(1.0-0.02));
      gauss->SetParameter(1,resolution);
      signal_temp_massScale_down=new TH1D("signal_massScale_systDown","",6200,0.,6200.);
      signal_temp_massScale_down->FillRandom("gauss",10000);
      signal_temp_massScale_down->Sumw2();
      signal_temp_massScale_down->Scale((xsec/1000.0)*acceff*lumi_scale/10000.);
      myfile << "signal syst " << mass_sig << " signal_massScale_systDown " << signal_temp_massScale_down->Integral(1,6000) << "\n";

      // Signal mass resolution systematic //
      gauss->SetParameter(0,mass_sig);
      gauss->SetParameter(1,resolution*(1.0+0.05));
      signal_temp_massRes_up=new TH1D("signal_massRes_systUp","",6200,0.,6200.);
      signal_temp_massRes_up->FillRandom("gauss",10000);
      signal_temp_massRes_up->Sumw2();
      signal_temp_massRes_up->Scale((xsec/1000.0)*acceff*lumi_scale/10000.);
      myfile << "signal syst " << mass_sig << " signal_massRes_systUp " << signal_temp_massRes_up->Integral(1,6000) << "\n";

      gauss->SetParameter(0,mass_sig);
      gauss->SetParameter(1,resolution*(1.0-0.05));
      signal_temp_massRes_down=new TH1D("signal_massRes_systDown","",6200,0.,6200.);
      signal_temp_massRes_down->FillRandom("gauss",10000);
      signal_temp_massRes_down->Sumw2();
      signal_temp_massRes_down->Scale((xsec/1000.0)*acceff*lumi_scale/10000.);
      myfile << "signal syst " << mass_sig << " signal_massRes_systDown " << signal_temp_massRes_down->Integral(1,6000) << "\n";
      
      // Signal pdf systematic //
      TH1D* hist_pdf_thisSig_Up = new TH1D("SigPDFUP", "SigPDFUP", 6200, 0, 6200); //  = hist_pdf_rel_up *(*signal_temp);
      TH1D* hist_pdf_thisSig_Down  = new TH1D("SigPDFDOWN", "SigPDFDOWN", 6200, 0, 6200); //= hist_pdf_rel_down *(*signal_temp);
      for (int ii=1; ii<6201; ii++) {
        double mySig_up=hist_pdf_rel_up->GetBinContent(ii);
        double mySig_down=hist_pdf_rel_down->GetBinContent(ii);
        double mySig_mean=signal_temp->GetBinContent(ii);
	double sig_mul_up = mySig_up*mySig_mean;
	double sig_mul_down = mySig_down*mySig_mean;
        hist_pdf_thisSig_Up->SetBinContent(ii,sig_mul_up);
        hist_pdf_thisSig_Down->SetBinContent(ii,sig_mul_down);
      }
      // std::cout << "Will set name now (UP)" << std::endl;
      TString mass_sig_name = "";
      mass_sig_name += int(mass_sig); 
      std::cout << "mass_sig_name " << mass_sig_name << std::endl;
      hist_pdf_thisSig_Up->SetName("signal_pdf_systUp");
      myfile << "signal syst " << mass_sig  << " signal_pdf_systUp " <<   hist_pdf_thisSig_Up->Integral(1,6000) << "\n";
      // std::cout << "Will set name now (DOWN)" << std::endl;
      hist_pdf_thisSig_Down->SetName("signal_pdf_systDown");
      myfile << "signal syst " << mass_sig << " signal_pdf_systDown " <<  hist_pdf_thisSig_Down->Integral(1,6000) << "\n";

      // Write in output file
      outfile_signal->cd();
      signal_temp->Write("signal");
      signal_temp_massScale_up->Write("signal_massScale_systUp");
      signal_temp_massScale_down->Write("signal_massScale_systDown");
      signal_temp_massRes_up->Write("signal_massRes_systUp");
      signal_temp_massRes_down->Write("signal_massRes_systDown");
      
      hist_pdf_thisSig_Up->Write();
      hist_pdf_thisSig_Down->Write();
      myfile_temp << "signal " << mass_sig << " " << signal_temp->Integral(1,6000) << "\n";
      myfile_temp << "signal syst " << mass_sig << " signal_massScale_systUp " << signal_temp_massScale_up->Integral(1,6000) << "\n";
      myfile_temp << "signal syst " << mass_sig << " signal_massScale_systDown " << signal_temp_massScale_down->Integral(1,6000) << "\n";
      myfile_temp << "signal syst " << mass_sig << " signal_massRes_systUp " << signal_temp_massRes_up->Integral(1,6000) << "\n";
      myfile_temp << "signal syst " << mass_sig << " signal_massRes_systDown " << signal_temp_massRes_down->Integral(1,6000) << "\n";
      myfile_temp << "signal syst " << mass_sig << " signal_pdf_systUp " << hist_pdf_thisSig_Up->Integral(1,6000) << "\n";
      myfile_temp << "signal syst " << mass_sig << " signal_pdf_systDown " << hist_pdf_thisSig_Down->Integral(1,6000) << "\n";

      masses[k]=mass_sig;

      mass_sig += step_size;   
      delete signal_temp;
  
      myfile_temp.close();

    }


      //close the output txt file:
      myfile.close();
    
      delete gauss;

}
