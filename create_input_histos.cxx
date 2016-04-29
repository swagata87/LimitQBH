//#include "create_input_histos.h"

#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TString.h"
#include <sstream>
#include "TH1.h"
#include "TF1.h"
#include "TString.h"
//#include "Environment.h"
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

  int bin_limit_lower=1;
  int bin_limit_upper=10000;

  TString rootfilenames[] = {
    "/net/scratch_cms/institut_3a/mukherjee/March24_Envelope/ttbar_tot.root",
    "/net/scratch_cms/institut_3a/mukherjee/March24_Envelope/WW_tot.root",
    "/net/scratch_cms/institut_3a/mukherjee/March24_Envelope/SingleTop_tot.root",
    "/net/scratch_cms/institut_3a/mukherjee/March24_Envelope/DY_tot.root",
    "/net/scratch_cms/institut_3a/mukherjee/March24_Envelope/WZ_tot.root",
    "/net/scratch_cms/institut_3a/mukherjee/March24_Envelope/ZZ_tot.root",
    "/net/scratch_cms/institut_3a/mukherjee/March24_Envelope/Wgamma.root",
    "/net/scratch_cms/institut_3a/mukherjee/NewJson/Bkg_DataDriven_MuJet.root"
};

  //TString sample_names[] = {"TT_tot","WW_tot","single_top_tot","WZ_tot","ZZ_tot","datadriven"};
  TString sample_names[] = {"TT_tot","WW_tot","single_top_tot","DY_tot","WZ_tot","ZZ_tot","Wgamma","datadriven"};
  const int arraySize = sizeof(sample_names)/sizeof(sample_names[0]);  

  //names of systematics for shape-based input
  TString syst_names[] ={"pileup_syst_","eleID_syst_","Ele_syst_Scale","Muon_syst_Scale","Muon_syst_Resolution","muoID_syst_","TopEnvelope_syst_"};
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

  //  double mass_min=500;
  //  double mass_max=4500;
  //  int step_size=1000;
  //  int N_points=(int)(mass_max-mass_min)/(double)step_size;

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

  //
  //  PDF Systematics
  //
  TString pdffilename = "/net/scratch_cms/institut_3a/erdweg/public/13TeV_rpv/76_test/for_limit.root";
  TFile* pdffile = new TFile(pdffilename);
  std::cout << "File for PDF uncert : " << pdffilename << std::endl;
  TH1D* hist_pdf_mean=(TH1D*)pdffile->Get("PDF4LHC15_nnlo_mc_mean");
  TH1D* hist_pdf_up=(TH1D*)pdffile->Get("PDF4LHC15_nnlo_mc_up");
  TH1D* hist_pdf_down=(TH1D*)pdffile->Get("PDF4LHC15_nnlo_mc_down");

  std::cout << "Nbin pdf bkg = " << hist_pdf_up->GetXaxis()->GetNbins() << std::endl;


  TH1D* hist_pdf_rel_up = new TH1D("reluppdf", "pdfuprel", 6000, 0, 6000);   //  = ((*hist_pdf_up)/(*hist_pdf_mean));
  TH1D* hist_pdf_rel_down = new TH1D("reldownpdf", "pdfdownrel", 6000, 0, 6000);//= ((*hist_pdf_down)/(*hist_pdf_mean));

  for (int i=1; i<6000; i++) {
    double up=hist_pdf_up->GetBinContent(i);
    double down=hist_pdf_down->GetBinContent(i);
    double mean=hist_pdf_mean->GetBinContent(i);
    double rel_up = 0;   
    double rel_down = 0; 
    if (mean!=0.0) rel_up =   up/mean;
    if (mean!=0.0) rel_down = down/mean;

    //    std::cout << "bin " << i << "rel_up " << rel_up << "  rel_down " << rel_down << std::endl;
    hist_pdf_rel_up->SetBinContent(i,rel_up);
    hist_pdf_rel_down->SetBinContent(i,rel_down);
  }

  //###############
  //loop over backgrounds
  //###############  
  if (debug) std::cout << "will start background loop"<< std::endl;

  //////////////////////////////////////
  //              UPDATE              //
  double Lumi_bkg = (2671.0/1000.0);
  std::cout << "Lumi scale factor for bkg : " << Lumi_bkg << std::endl;
  /////////////////////////////////////
  

  double ttbar_kfact = 1.116;
  double WW_kfact = 1.1619;

  for (int i = 0; i < arraySize; ++i){
    if (debug) std::cout << "will get the bkg root file"<< std::endl;
    std::cout << "\n\nrootfilename " << rootfilenames[i] << "  sample_name " << sample_names[i] << std::endl;
    TFile* infile = new TFile(rootfilenames[i]);
    
    //TAG change to one input file per bkg
    TH1D* hist_ori;
    //    cout << ((TString)(sample_names[i]) + "_7_0/"+(TString)(sample_names[i])+ "_ori").Data() << endl;
    if (debug) std::cout << "will get the hist h1_0_emu_Mass"<< std::endl;
    
    hist_ori=(TH1D*)infile->Get("emu/Stage_0/h1_0_emu_Mass");
    if (sample_names[i] != "datadriven") hist_ori->Scale(Lumi_bkg);
    if (sample_names[i]=="TT_tot") {
      std::cout << "Extra k-factor scaling will be done for TTbar background" << std::endl;
      hist_ori->Scale(ttbar_kfact);
    }

    if (sample_names[i]=="WW_tot") {
      std::cout << "Extra k-factor scaling will be done for WW background" << std::endl;
      hist_ori->Scale(WW_kfact);
    }

    hist_ori->SetName(sample_names[i]);
    

    /*    for(int n=1;n<(bin_limit_lower);n++)
      {
	hist_ori->SetBinContent(n,0.);
	hist_ori->SetBinError(n,0.);    
      }
    for(int m=(bin_limit_upper+1);m<10000;m++)
      {
	hist_ori->SetBinContent(m,0.);
    	hist_ori->SetBinError(m,0.);    
      }
    */
    outfile->cd();
    hist_ori->Write(sample_names[i]);
    
    //TAG check which binning is used for the histograms (here: 1-10000, 1 GeV binning)
    std::cout  << "bkg " << sample_names[i] << " " << hist_ori->Integral(bin_limit_lower,bin_limit_upper) << "\n";
    myfile << "bkg " << sample_names[i] << " " << hist_ori->Integral(bin_limit_lower,bin_limit_upper) << "\n";
    
    
    if (debug) std::cout << "Now multiply mass_mean_hist_bkg with pdf_rel (up and down) hist" << std::endl;

    if (sample_names[i] != "datadriven") {
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

      /*
      for(int n=1;n<(bin_limit_lower);n++)
	{
	  hist_pdf_thisBkg_Up->SetBinContent(n,0.);
	  hist_pdf_thisBkg_Up->SetBinError(n,0.);    
	}
      for(int m=(bin_limit_upper+1);m<10000;m++)
	{
	  hist_pdf_thisBkg_Up->SetBinContent(m,0.);
	  hist_pdf_thisBkg_Up->SetBinError(m,0.);    
	}
      */

      hist_pdf_thisBkg_Up->Write();
      //std::cout << "MEAN hist classname " << hist_ori->ClassName() << std::endl;
      // std::cout << "UP hist classname " << hist_pdf_thisBkg_Up->ClassName() << std::endl;
      myfile << "bkg " << sample_names[i]+"_"+"pdf_syst"+"Up" << " " << hist_pdf_thisBkg_Up->Integral(bin_limit_lower,bin_limit_upper) << "\n";
      
      // hist_pdf_thisBkg_Down.Scale( (hist_ori->Integral()) / (hist_pdf_thisBkg_Down.Integral()) );
      hist_pdf_thisBkg_Down->SetName(sample_names[i]+"_"+"pdf_syst"+"Down");

      /*
      for(int n=1;n<(bin_limit_lower);n++)
	{
          hist_pdf_thisBkg_Down->SetBinContent(n,0.);
          hist_pdf_thisBkg_Down->SetBinError(n,0.);
        }
      for(int m=(bin_limit_upper+1);m<10000;m++)
        {
          hist_pdf_thisBkg_Down->SetBinContent(m,0.);
          hist_pdf_thisBkg_Down->SetBinError(m,0.);
        }
      */

      hist_pdf_thisBkg_Down->Write(); 
      myfile << "bkg " << sample_names[i]+"_"+"pdf_syst"+"Down" << " " << hist_pdf_thisBkg_Down->Integral(bin_limit_lower,bin_limit_upper) << "\n";
    }
    //###############
    //loop over other systematics
    //###############  
    if (debug) std::cout << "will start the syst loop"<< std::endl;
    
    if (sample_names[i] != "datadriven") {
      for(int k=0; k<arraySize_systs; k++)
	{
	  TH1D* hist_syst_up;
	  TH1D* hist_syst_down;	
	  
	  //hist_syst_up=(TH1F*)infile->Get(((TString)(sample_names[i]) + "_7_0/"+(TString)(sample_names[i])+ "_" + syst_names[k] + "Up").Data());
	  if (debug) std::cout << "will get syst hist UP"<< std::endl;
	  
	  hist_syst_up   = (TH1D*)infile->Get("emu/Stage_0/sys/h1_0_emu_Mass_"+syst_names[k]+"Up");
	  //std::cout << "Lumi scaling for " << sample_names[i] << " and " << syst_names[k] << " UP"  << std::endl;
	  hist_syst_up->Scale(Lumi_bkg);
	  hist_syst_down = (TH1D*)infile->Get("emu/Stage_0/sys/h1_0_emu_Mass_"+syst_names[k]+"Down");
	  //std::cout << "Lumi scaling for " << sample_names[i]  << " and " << syst_names[k] << " DOWN"  << std::endl;
	  hist_syst_down->Scale(Lumi_bkg);
	  hist_syst_up->SetName(sample_names[i]+"_"+syst_names[k]+"Up");
	  /*
	  for(int n=1;n<(bin_limit_lower);n++)
	    {
	      hist_syst_up->SetBinContent(n,0.);
	      hist_syst_up->SetBinError(n,0.);
	    }
	  for(int m=(bin_limit_upper+1);m<10000;m++)
	    {
	      hist_syst_up->SetBinContent(m,0.);
	      hist_syst_up->SetBinError(m,0.);
	    }
	  */

	  hist_syst_up->Write(); 
	  myfile << "bkg " << sample_names[i]+"_"+syst_names[k]+"Up" << " " << hist_syst_up->Integral(bin_limit_lower,bin_limit_upper) << "\n";
	  
	  hist_syst_down->SetName(sample_names[i]+"_"+syst_names[k]+"Down");

	  /*
	  for(int n=1;n<(bin_limit_lower);n++)
            {
              hist_syst_down->SetBinContent(n,0.);
              hist_syst_down->SetBinError(n,0.);
            }
          for(int m=(bin_limit_upper+1);m<10000;m++)
            {
              hist_syst_down->SetBinContent(m,0.);
              hist_syst_down->SetBinError(m,0.);
            }
	  */


	  hist_syst_down->Write(); 
	  myfile << "bkg " << sample_names[i]+"_"+syst_names[k]+"Down" << " " << hist_syst_down->Integral(bin_limit_lower,bin_limit_upper) << "\n";
	  
	  if (debug) std::cout << "Will delete hist_syst_up hist_syst_down "<< std::endl;
	  
	  delete hist_syst_up;
	  delete hist_syst_down;	
	if (debug) std::cout << "Successfully deleted hist_syst_up hist_syst_down "<< std::endl;
	
      }
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
  TFile* data_file = new TFile("/net/scratch_cms/institut_3a/mukherjee/NewJson/allData.root");
  TH1D* data;
  //  data=(TH1F*)data_file->Get("h1_inv_mass_1mu_1tau_aligned_7_0");
  if (debug) std::cout << "will get data hist"<< std::endl;
  data=(TH1D*)data_file->Get("emu/Stage_0/h1_0_emu_Mass"); 

  outfile->cd();
  data->SetName("data_obs");
  /*
  for(int n=1;n<(bin_limit_lower);n++)
    {
      data->SetBinContent(n,0.);
      data->SetBinError(n,0.);
    }
  for(int m=(bin_limit_upper+1);m<10000;m++)
    {
      data->SetBinContent(m,0.);
      data->SetBinError(m,0.);
    }
  */


  data->Write();

  myfile << "data " << "data" << " " << data->Integral(bin_limit_lower,bin_limit_upper) << "\n";

  //###############
  //prepare input for the individual mass points
  //###############
  //exit(0);
  
  TFile* infile_sig;
  TString rootfilenames_sig[] = {
    "/net/scratch_cms/institut_3a/erdweg/public/13TeV_rpv/Feb_02/QBHToEMu_M-500_n1_RS-QBH_13TeV_P8-skimid3601.root",
    "/net/scratch_cms/institut_3a/erdweg/public/13TeV_rpv/Feb_02/QBHToEMu_M-1000_n1_RS-QBH_13TeV_P8-skimid3609.root",
    "/net/scratch_cms/institut_3a/erdweg/public/13TeV_rpv/Feb_02/QBHToEMu_M-1500_n1_RS-QBH_13TeV_P8-skimid3607.root",
    "/net/scratch_cms/institut_3a/erdweg/public/13TeV_rpv/Feb_02/QBHToEMu_M-2500_n1_RS-QBH_13TeV_P8-skimid3613.root",
    "/net/scratch_cms/institut_3a/erdweg/public/13TeV_rpv/Feb_02/QBHToEMu_M-3000_n1_RS-QBH_13TeV_P8-skimid3600.root",
    "/net/scratch_cms/institut_3a/erdweg/public/13TeV_rpv/Feb_02/QBHToEMu_M-3500_n1_RS-QBH_13TeV_P8-skimid3646.root",
    "/net/scratch_cms/institut_3a/erdweg/public/13TeV_rpv/Feb_02/QBHToEMu_M-4000_n1_RS-QBH_13TeV_P8-skimid3597.root",
    "/net/scratch_cms/institut_3a/erdweg/public/13TeV_rpv/Feb_02/QBHToEMu_M-4500_n1_RS-QBH_13TeV_P8-skimid3643.root",
    "/net/scratch_cms/institut_3a/erdweg/public/13TeV_rpv/Feb_02/QBHToEMu_M-5000_n1_RS-QBH_13TeV_P8-skimid3594.root",
  };
  
  TFile* infile_sig_all;
  TString rootfilenames_sig_all[] = {
    "/net/scratch_cms/institut_3a/mukherjee/March17_pT_smear/QBHToEMu_M-500_n1_RS-QBH_13TeV_P8-skimid3601.root",
    "/net/scratch_cms/institut_3a/mukherjee/March17_pT_smear/QBHToEMu_M-1000_n1_RS-QBH_13TeV_P8-skimid3609.root",
    "/net/scratch_cms/institut_3a/mukherjee/March17_pT_smear/QBHToEMu_M-1500_n1_RS-QBH_13TeV_P8-skimid3607.root",
    "/net/scratch_cms/institut_3a/mukherjee/March17_pT_smear/QBHToEMu_M-2500_n1_RS-QBH_13TeV_P8-skimid3613.root",
    "/net/scratch_cms/institut_3a/mukherjee/March17_pT_smear/QBHToEMu_M-3000_n1_RS-QBH_13TeV_P8-skimid3600.root",
    "/net/scratch_cms/institut_3a/mukherjee/March17_pT_smear/QBHToEMu_M-3500_n1_RS-QBH_13TeV_P8-skimid3646.root",
    "/net/scratch_cms/institut_3a/mukherjee/March17_pT_smear/QBHToEMu_M-4000_n1_RS-QBH_13TeV_P8-skimid3597.root",
    "/net/scratch_cms/institut_3a/mukherjee/March17_pT_smear/QBHToEMu_M-4500_n1_RS-QBH_13TeV_P8-skimid3643.root",
    "/net/scratch_cms/institut_3a/mukherjee/March17_pT_smear/QBHToEMu_M-5000_n1_RS-QBH_13TeV_P8-skimid3594.root",
  };
  

  std::string sample_names_sig[] = {"500" , "1000", "1500", "2500", "3000", "3500", "4000", "4500", "5000"};
  //std::string sample_names_sig[] = {"5000"};
  const int arraySize_sig = sizeof(sample_names_sig)/sizeof(sample_names_sig[0]);

  TString syst_names_sig[] ={"pileup_syst_","eleID_syst_","Ele_syst_Scale","Muon_syst_Scale","Muon_syst_Resolution","muoID_syst_"};
  const int arraySize_systs_sig = sizeof(syst_names_sig)/sizeof(syst_names_sig[0]);


  TF1* fit_acceff=new TF1("fit_acceff","[0]+[1]/(x+[2])+[3]*x",200.,6000.);
  fit_acceff->SetParameter(0,0.7961);
  fit_acceff->SetParameter(1,-130.8);
  fit_acceff->SetParameter(2,212.174487);
  fit_acceff->SetParameter(3,-0.000022);
  double acceff=0.7;
  /*
  for(int mass_sig=500;mass_sig<5000;mass_sig+=500) {
    acceff=fit_acceff->Eval(mass_sig);
    std::cout <<  mass_sig << "&   " <<     acceff << " \\\\" << std::endl;


  }
  */


  std::cout << "No- of systematics " << arraySize_systs_sig << "  No. of signal MC samples " << arraySize_sig << std::endl;

  for (int k = 0; k < arraySize_sig; ++k){
    TFile* infile_sig_all = new TFile(rootfilenames_sig_all[k]);
    TFile* infile_sig = new TFile(rootfilenames_sig[k]);

    TH1D* signal_temp;
    sprintf(file_title,"root_out/out_mass_%s.root",sample_names_sig[k].c_str());
    sprintf(dir_title,"txt_out/normalization_Mass_%s_input_histos.txt",sample_names_sig[k].c_str());
    ofstream myfile_temp;
    myfile_temp.open(dir_title);
    TFile* outfile_signal = new TFile(file_title,"recreate");
    outfile_signal->cd();
    
    int counter_histos=0;
    
    TKey *key;
    TIter next(outfile->GetListOfKeys());
    
  while ((key = (TKey*)next())) {
    
    TH1 *h1 = (TH1*)key->ReadObj();
    TString *name = new TString(h1->GetName());

    
    if(strcmp(h1->GetName(),"data_obs") ==0 )
      {
	myfile_temp << "data " << h1->GetName() << " " << h1->Integral(bin_limit_lower,bin_limit_upper) << "\n";	    
      }	
    else
      {
	myfile_temp << "bkg " << h1->GetName() << " " << h1->Integral(bin_limit_lower,bin_limit_upper) << "\n"; 
      }
    
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

    h1->Write();
    
    counter_histos++;
    
    delete name;
    delete h1;
    
  }
  
  signal_temp=(TH1D*)infile_sig_all->Get("emu/Stage_0/h1_0_emu_Mass");
  //std::cout << "Nbin = " << signal_temp->GetXaxis()->GetNbins() << std::endl;
  std::cout << "signal mass: " << sample_names_sig[k] << "  Mass resolution = " << signal_temp->GetRMS() << std::endl;
  signal_temp->Scale(Lumi_bkg);
  myfile << "signal " << sample_names_sig[k] << " " << signal_temp->Integral(bin_limit_lower,bin_limit_upper) << "\n";
  //cout << "Nbin = " << signal_temp->GetNbinsX() << endl;
  //cout << "Simple (1-6000):: signal mass: " << sample_names_sig[k] << " N events: " << signal_temp->Integral(1,6000) << endl;
  //cout << "Overflow (1-7000):: signal mass: " << sample_names_sig[k] << " N events: " << signal_temp->Integral(1,7000) << endl;
  outfile_signal->cd();
  myfile_temp << "signal " << sample_names_sig[k] << " " << signal_temp->Integral(bin_limit_lower,bin_limit_upper) << "\n";

  TH1D* hist_pdf_mean_sig=(TH1D*)infile_sig->Get("PDF4LHC15_nnlo_mc_mean");
  TH1D* hist_pdf_up_sig=(TH1D*)infile_sig->Get("PDF4LHC15_nnlo_mc_up");
  TH1D* hist_pdf_down_sig=(TH1D*)infile_sig->Get("PDF4LHC15_nnlo_mc_down");

  //  std::cout << "Nbin1 = " << hist_pdf_up_sig->GetXaxis()->GetNbins() << std::endl;
  // std::cout << "Nbin2 = " << hist_pdf_down_sig->GetXaxis()->GetNbins() << std::endl;


  TH1D* hist_pdf_rel_up_sig = new TH1D("reluppdf", "pdfuprel", 6000, 0, 6000);
  TH1D* hist_pdf_rel_down_sig = new TH1D("reldownpdf", "pdfdownrel", 6000, 0, 6000);

  
  for (int i=1; i<6000; i++) {
    double up=hist_pdf_up_sig->GetBinContent(i);
    double down=hist_pdf_down_sig->GetBinContent(i);
    double mean=hist_pdf_mean_sig->GetBinContent(i);
    double rel_up = 0;
    double rel_down = 0;
    if (mean!=0.0) rel_up =   up/mean;
    if (mean!=0.0) rel_down = down/mean;
    //if (i<100) std::cout << "SIGNAL bin " << i  << "  mean " << mean << " up " << up << " down " << down << std::endl;
    hist_pdf_rel_up_sig->SetBinContent(i,rel_up);
    hist_pdf_rel_down_sig->SetBinContent(i,rel_down);
  
  }

  TH1D* hist_pdf_thisSig_Up = new TH1D("SigPDFUP", "SigPDFUP", 6000, 0, 6000); //  = hist_pdf_rel_up *(*signal_temp);                                                                          
  TH1D* hist_pdf_thisSig_Down  = new TH1D("SigPDFDOWN", "SigPDFDOWN", 6000, 0, 6000); //= hist_pdf_rel_down *(*signal_temp);                                                                   
 
  for (int ii=1; ii<6000; ii++) {
    double mySig_up=hist_pdf_rel_up_sig->GetBinContent(ii);
    double mySig_down=hist_pdf_rel_down_sig->GetBinContent(ii);
    double mySig_mean=signal_temp->GetBinContent(ii);
    double sig_mul_up = mySig_up*mySig_mean;
    double sig_mul_down = mySig_down*mySig_mean;
    //if (ii<100) 
    //std::cout << "bin " << ii << " :  " << mySig_up << "*" << mySig_mean << "=" << sig_mul_up << std::endl; 
    hist_pdf_thisSig_Up->SetBinContent(ii,sig_mul_up);
    hist_pdf_thisSig_Down->SetBinContent(ii,sig_mul_down);
  }

  hist_pdf_thisSig_Up->SetName("signal_pdf_systUp");
  myfile << "signal syst " << sample_names_sig[k]  << " signal_pdf_systUp " <<   hist_pdf_thisSig_Up->Integral(bin_limit_lower,bin_limit_upper) << "\n";
  myfile_temp << "signal syst " << sample_names_sig[k]  << " signal_pdf_systUp " <<   hist_pdf_thisSig_Up->Integral(bin_limit_lower,bin_limit_upper) << "\n";

  hist_pdf_thisSig_Down->SetName("signal_pdf_systDown");
  myfile << "signal syst " << sample_names_sig[k] << " signal_pdf_systDown " <<  hist_pdf_thisSig_Down->Integral(bin_limit_lower,bin_limit_upper) << "\n";
  myfile_temp << "signal syst " << sample_names_sig[k] << " signal_pdf_systDown " <<  hist_pdf_thisSig_Down->Integral(bin_limit_lower,bin_limit_upper) << "\n";

  /*
  for(int n=1;n<(bin_limit_lower);n++)
    {
      signal_temp->SetBinContent(n,0.);
      signal_temp->SetBinError(n,0.);
    }
  for(int m=(bin_limit_upper+1);m<10000;m++)
    {
      signal_temp->SetBinContent(m,0.);
      signal_temp->SetBinError(m,0.);
    }
  */


  signal_temp->Write("signal");

  /*

  for(int n=1;n<(bin_limit_lower);n++)
    {
      hist_pdf_thisSig_Up->SetBinContent(n,0.);
      hist_pdf_thisSig_Up->SetBinError(n,0.);
    }
  for(int m=(bin_limit_upper+1);m<10000;m++)
    {
      hist_pdf_thisSig_Up->SetBinContent(m,0.);
      hist_pdf_thisSig_Up->SetBinError(m,0.);
    }
  */

  hist_pdf_thisSig_Up->Write();

  /*
  for(int n=1;n<(bin_limit_lower);n++)
    {
      hist_pdf_thisSig_Down->SetBinContent(n,0.);
      hist_pdf_thisSig_Down->SetBinError(n,0.);
    }
  for(int m=(bin_limit_upper+1);m<10000;m++)
    {
      hist_pdf_thisSig_Down->SetBinContent(m,0.);
      hist_pdf_thisSig_Down->SetBinError(m,0.);
    }
  */

  hist_pdf_thisSig_Down->Write();

  

  // Systematics
  //if (debug) 
  // std::cout << "will start the syst_sig loop"<< std::endl;
  for(int kk=0; kk<arraySize_systs_sig; kk++)
    {
      TH1D* sig_hist_syst_up;
      TH1D* sig_hist_syst_down;
      sig_hist_syst_up   = (TH1D*)infile_sig_all->Get("emu/Stage_0/sys/h1_0_emu_Mass_"+syst_names_sig[kk]+"Up");
      //std::cout << "Lumi scaling for " << sample_names_sig[k] << " and " << syst_names_sig[kk] << " UP"  << std::endl;
      sig_hist_syst_up->Scale(Lumi_bkg);
      sig_hist_syst_down = (TH1D*)infile_sig_all->Get("emu/Stage_0/sys/h1_0_emu_Mass_"+syst_names_sig[kk]+"Down");
      // std::cout << "Lumi scaling for " << sample_names_sig[k]  << " and " << syst_names_sig[kk] << " DOWN"  << std::endl;
      sig_hist_syst_down->Scale(Lumi_bkg);
      sig_hist_syst_up->SetName("signal_"+syst_names_sig[kk]+"Up");


      /*
      for(int n=1;n<(bin_limit_lower);n++)
	{
	  sig_hist_syst_up->SetBinContent(n,0.);
	  sig_hist_syst_up->SetBinError(n,0.);
	}
      for(int m=(bin_limit_upper+1);m<10000;m++)
	{
	  sig_hist_syst_up->SetBinContent(m,0.);
	  sig_hist_syst_up->SetBinError(m,0.);
	}
      */

      sig_hist_syst_up->Write();
      myfile << "signal syst " << sample_names_sig[k] << " " << "signal_"+syst_names_sig[kk]+"Up" << " " << sig_hist_syst_up->Integral(bin_limit_lower,bin_limit_upper) << "\n";
      myfile_temp << "signal syst " << sample_names_sig[k] << " " << "signal_"+syst_names_sig[kk]+"Up" << " " << sig_hist_syst_up->Integral(bin_limit_lower,bin_limit_upper) << "\n";

      sig_hist_syst_down->SetName("signal_"+syst_names_sig[kk]+"Down");


      /*
      for(int n=1;n<(bin_limit_lower);n++)
	{
          sig_hist_syst_down->SetBinContent(n,0.);
          sig_hist_syst_down->SetBinError(n,0.);
        }
      for(int m=(bin_limit_upper+1);m<10000;m++)
        {
          sig_hist_syst_down->SetBinContent(m,0.);
          sig_hist_syst_down->SetBinError(m,0.);
        }
      */

      sig_hist_syst_down->Write();
      myfile << "signal syst " << sample_names_sig[k] << " " << "signal_"+syst_names[kk]+"Down" << " " << sig_hist_syst_down->Integral(bin_limit_lower,bin_limit_upper) << "\n";
      myfile_temp << "signal syst " << sample_names_sig[k] << " " << "signal_"+syst_names[kk]+"Down" << " " << sig_hist_syst_down->Integral(bin_limit_lower,bin_limit_upper) << "\n";

      if (debug) std::cout << "Will delete hist_syst_up hist_syst_down "<< std::endl;

      delete sig_hist_syst_up;
      delete sig_hist_syst_down;
      if (debug) std::cout << "Successfully deleted hist_syst_up hist_syst_down "<< std::endl;

    }

  delete signal_temp;
  myfile_temp.close();
  }
  //close the output txt file:                                                                                                                                                                     
  myfile.close();
  //  delete gauss;                                                                                                                                                                                

}
