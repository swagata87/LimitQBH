#include "TH1.h"
#include "TTree.h"
#include "TBranch.h"
#include <algorithm>
#include "TGraphAsymmErrors.h"
#include "/user/mukherjee/limits_LFV/scripts/Environment.h"
#include "TSystem.h"
#include "TString.h"
#include "TEnv.h"

#include "TH1F.h"
#include "TFile.h"
#include "TString.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "TVectorD.h"
#include "TLine.h"

#include "TPad.h"
#include "TLatex.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include "TSpline.h"
#include <cstdlib>
//#include "/user/gueth/CMSSW_5_3_11_patch1/src/limits/Environment.h"

using std::cout;
using std::endl;
using namespace std;
using namespace TMath;

void limit_scan()
{

  gROOT->Reset();
  gROOT->SetStyle("Plain");

  /*  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.15);  
  gStyle->SetPadTopMargin(0.05);   
  gStyle->SetTitleXSize(0.04);
  gStyle->SetTitleXOffset(1.05);
  gStyle->SetTitleYSize(0.04);
  gStyle->SetTitleYOffset(1.05);   
  */
  //expected including MC stat. uncertainty
  double observed_xsec[200];
  double mass[200];
  double kM[200];
  double kM_pb[200];

  int counter_array=0;

  for(int k=0; k<200; k++)
    {
      observed_xsec[k]=0.;
      mass[k]=0.;
      kM_pb[k]=0.;
      kM[k]=0.;      
    }

  string line;
  ifstream observed_in ("expected_limit_mass_kM.txt");
  //ifstream observed_in ("observed_limit_mass_kM.txt");
  //ifstream observed_in ("observed_limit_mass_kM_add_cuts.txt");
  if (observed_in.is_open())
  {
    while(!observed_in.eof())
    {

	  observed_in >> mass[counter_array];
	  observed_in >> kM[counter_array];      
	  observed_in >> observed_xsec[counter_array];      
     
	  std::cout << "\n" << mass[counter_array] << endl;
	  std::cout << kM[counter_array] << endl;      
	  std::cout << observed_xsec[counter_array] << endl;      

      counter_array++;
      std::cout << counter_array << endl;

      if(mass[counter_array-1]==4000)break;

    }
    observed_in.close();
  }

  std::cout << counter_array << endl;

  double lambda_fix[4];
  lambda_fix[0]=0.07;  
  lambda_fix[1]=0.05;
  lambda_fix[2]=0.01;
  lambda_fix[3]=0.007;  

  //remove after testing

  double xsec_max[200];

  int counter_points_plot[4];

  double lambdaP_solution[4][200];

  for(int k=0; k<200; k++)
    {
      xsec_max[k]=0.;
    }

  for(int j=0;j<4;j++)
    {
      counter_points_plot[j]=0.;
      for(int f=0;f<counter_array;f++)
	{
	  lambdaP_solution[j][f]=0.;
	}
    }
  
  for(int j=0;j<4;j++)
    {
      for(int f=0;f<counter_array;f++)
	{
	  
	  //kM[f]=kM_pb[f]*1000.;
	  xsec_max[f]=kM[f]*2./3.*pow(lambda_fix[j],2);

	  //check the solvability condition
	  if(xsec_max[f]<observed_xsec[f])break;
	  
	  //get the solution and keep counting ...
	  lambdaP_solution[j][f]=sqrt((2*pow(lambda_fix[j],2))/((2*pow(lambda_fix[j],2)*kM[f]/observed_xsec[f])-3));
	  
	  /*
	    if(j==0)
	    {
	      std::cout << f << endl;
	      std::cout << mass[f] << endl;
	      std::cout << observed_xsec[f] << " " << xsec_max[f] << endl;	    
	      std::cout << lambdaP_solution[j][f] << endl;
	      std::cout << endl;
	    }
	  */

	  /*
	    std::cout << mass[f] << endl;
	    std::cout << kM[f] << endl;
	    std::cout << observed_xsec[f] << endl;      
	    std::cout << lambdaP_solution[j][f] << endl;
	    std::cout << endl;
	  */

	  counter_points_plot[j]++;	
	}
      
      if(counter_points_plot[j]<counter_array)
	{
	  lambdaP_solution[j][counter_points_plot[j]]=10.;
	  counter_points_plot[j]++;
	}
    }

  /*
  for(int fg=0;fg<counter_points_plot[0];fg++)
    {
	    std::cout << mass[fg] << endl;
	    std::cout << observed_xsec[fg] << endl;      
	    std::cout << lambdaP_solution[3][fg] << endl;
	    std::cout << endl;      
    }
  */

  TCanvas *c_lambdaP=new TCanvas("c_lambdaP","",800,800);
  c_lambdaP->SetGridx();
  c_lambdaP->SetGridy();  

  //lambdaP=0.07
  double lambdaP_solution_0[200];
  for(int f=0;f<counter_points_plot[0];f++)
    {
      lambdaP_solution_0[f]=lambdaP_solution[0][f];
    }

  TGraph *graph_lambdaP_0 = new TGraph(counter_points_plot[0],mass,lambdaP_solution_0);  

 
  c_lambdaP->cd();
  gPad->SetLogy();

  graph_lambdaP_0->SetLineWidth(3);
  graph_lambdaP_0->SetLineStyle(10);  
  graph_lambdaP_0->SetLineColor(kRed);  
  graph_lambdaP_0->SetMarkerStyle(22);
  graph_lambdaP_0->SetMarkerSize(1.2); 
  graph_lambdaP_0->SetMarkerColor(kRed);  
  graph_lambdaP_0->GetYaxis()->SetRangeUser(0.0005,0.2);
  graph_lambdaP_0->GetXaxis()->SetRangeUser(200.,5000.);

  graph_lambdaP_0->GetXaxis()->SetTitleFont(42);
  graph_lambdaP_0->GetYaxis()->SetTitleFont(42);
  graph_lambdaP_0->GetXaxis()->SetLabelFont(42);
  graph_lambdaP_0->GetYaxis()->SetLabelFont(42);  

  graph_lambdaP_0->GetYaxis()->SetTitle("#lambda'_{311}");
  graph_lambdaP_0->GetXaxis()->SetTitle("M_{#tilde{#nu}_{#tau}} (GeV)");  
  graph_lambdaP_0->GetXaxis()->SetLabelSize(0.03);

  //lambdaP=0.05
  double lambdaP_solution_1[200];
  for(int f=0;f<counter_points_plot[1];f++)
    {
      lambdaP_solution_1[f]=lambdaP_solution[1][f];
    }

  TGraph *graph_lambdaP_1 = new TGraph(counter_points_plot[1],mass,lambdaP_solution_1);  

  graph_lambdaP_1->SetLineWidth(3);
  graph_lambdaP_1->SetLineStyle(7); 
  graph_lambdaP_1->SetLineColor(kBlue);  
  graph_lambdaP_1->SetMarkerStyle(22);
  graph_lambdaP_1->SetMarkerSize(1.2); 
  graph_lambdaP_1->SetMarkerColor(kBlue);   
  graph_lambdaP_1->GetYaxis()->SetRangeUser(0.001,0.4);  
  graph_lambdaP_1->GetXaxis()->SetRangeUser(200.,5000.);

  //lambdaP=0.01
  double lambdaP_solution_2[200];
  for(int f=0;f<counter_points_plot[2];f++)
    {
      lambdaP_solution_2[f]=lambdaP_solution[2][f];
    }

  TGraph *graph_lambdaP_2 = new TGraph(counter_points_plot[2],mass,lambdaP_solution_2);  

  graph_lambdaP_2->SetLineWidth(3);
  graph_lambdaP_2->SetLineStyle(9); 
  graph_lambdaP_2->SetLineColor(418);  
  graph_lambdaP_2->SetMarkerStyle(22);
  graph_lambdaP_2->SetMarkerSize(1.2);
  graph_lambdaP_2->SetMarkerColor(kGreen);    
  graph_lambdaP_2->GetYaxis()->SetRangeUser(0.001,0.4); 
  graph_lambdaP_2->GetXaxis()->SetRangeUser(200.,5000.);

  //lambdaP=0.007
  double lambdaP_solution_3[200];
  for(int f=0;f<counter_points_plot[3];f++)
    {
      lambdaP_solution_3[f]=lambdaP_solution[3][f];
    }

  TGraph *graph_lambdaP_3 = new TGraph(counter_points_plot[3],mass,lambdaP_solution_3);  

  graph_lambdaP_3->SetLineWidth(3);
  graph_lambdaP_3->SetLineStyle(1); 
  graph_lambdaP_3->SetLineColor(kBlack);  
  graph_lambdaP_3->SetMarkerStyle(22);
  graph_lambdaP_3->SetMarkerSize(1.2);
  graph_lambdaP_3->SetMarkerColor(kBlack);    
  graph_lambdaP_3->GetYaxis()->SetRangeUser(0.001,0.4); 
  graph_lambdaP_3->GetXaxis()->SetRangeUser(200.,5000.);

  //Let's draw it

  graph_lambdaP_0->Draw("Al");  
  graph_lambdaP_1->Draw("l,same");
  graph_lambdaP_2->Draw("l,same");
  graph_lambdaP_3->Draw("l,same");  

  TLegend *leg_example = new TLegend(0.45,0.15,0.95,0.40);
  leg_example->SetFillColor(0);
  leg_example->SetTextFont(42);
  
  leg_example->AddEntry(graph_lambdaP_0, "95% CL limit #lambda_{132}=#lambda_{231}=0.07","l");
  leg_example->AddEntry(graph_lambdaP_1, "95% CL limit #lambda_{132}=#lambda_{231}=0.05","l"); 
  leg_example->AddEntry(graph_lambdaP_2, "95% CL limit #lambda_{132}=#lambda_{231} =0.01","l"); 
  leg_example->AddEntry(graph_lambdaP_3, "95% CL limit #lambda_{132}=#lambda_{231}=0.007","l"); 

  leg_example->Draw("same");  


    TLatex* CMS_text = new TLatex(0.14,0.85,"CMS");
    CMS_text->SetNDC();
    CMS_text->SetTextSize(0.04);
    CMS_text->SetTextAngle(0);
    CMS_text->Draw("same");
    
    TLatex* CMS_text_2 = new TLatex(0.14,0.81,"Preliminary");
    //TLatex* CMS_text_2 = new TLatex(0.20,0.83,"own work in");
    CMS_text_2->SetNDC();
    CMS_text_2->SetTextFont(52);
    CMS_text_2->SetTextSize(0.04);
    CMS_text_2->SetTextAngle(0);
    CMS_text_2->Draw("same");    
    
    /*
    TLatex* CMS_text_3 = new TLatex(0.20,0.78,"progress");
    CMS_text_3->SetNDC();
    CMS_text_3->SetTextFont(52);
    CMS_text_3->SetTextSize(0.05);
    CMS_text_3->SetTextAngle(0);
    CMS_text_3->Draw("same");        
    */

    TLatex* lumiText = new TLatex(0.93,0.92,"2.5 fb^{-1} (13 TeV)");
    lumiText->SetNDC();
    lumiText->SetTextFont(42);
    lumiText->SetTextSize(0.03);
    lumiText->SetTextAlign(32);
    lumiText->Draw("same");     

    /*
    TLatex* expectedText = new TLatex(0.62,0.43,"expected limits");
    expectedText->SetNDC();
    expectedText->SetTextFont(42);
    expectedText->SetTextSize(0.04);
    expectedText->SetTextAngle(0);
    expectedText->Draw("same");       
    */
}
