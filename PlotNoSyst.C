#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TF1.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TH3F.h"

#include <iostream>
#include <fstream>
using namespace std;


int PlotNoSyst(TString file, TString output, string filename, int bin, TString yaxis, TString xaxis,TString yaxis1, TString xaxis1, double min, double max, TString title, string title1) {
	gROOT->SetBatch(kTRUE);
	
	
	
	
	TImage *img = TImage::Open("Claslogo.pdf");
	
        gStyle->SetPaintTextFormat("4.1f");
        gStyle->SetOptStat(1);
        gStyle->SetPalette(55);
        gStyle->SetLabelSize(.04, "xyz");
        gStyle->SetTitleSize(.04, "xyz");
        gStyle->SetTitleSize(.07,"t");
        gStyle->SetFrameLineWidth(1);
        gStyle->SetLineWidth(1);
        gStyle->SetHistLineWidth(1);
        gStyle->SetMarkerStyle(13);
        gStyle->SetTitleW(0.8);  //per cent of the pad width
        gStyle->SetTitleH(0.1); //per cent of the pad height

	int a = 9;
	int b = 0;

	//1.5GeV
	TGraph* hist1 = new TGraph("TheoryMarc/tcs_fb_diff_z_7p0_qps2p25_theta65_phi 0_rximrst02_1_bv1p0_bs5p0_dt_e1_t.dat","%lg %*lg %*lg %*lg %lg ");
	hist1->SetLineColor(kGreen);
	
	//1.8GeV
	TGraph* hist2 = new TGraph("TheoryMarc/tcs_fb_diff_z_7p0_qps3p24_theta50_phi 0_rximrst02_1_bv1p0_bs5p0_dt_e1_t.dat","%lg %*lg %*lg %*lg %lg ");
	hist2->SetLineColor(kRed);
	hist2->SetLineStyle(7);
	TGraph* hist3 = new TGraph("TheoryMarcNew/tcs_fb_diff_z_7p2_qps3p24_theta65_phi0_rximrst02_1_bv1p0_bs5p0_dt_e1_t.dat","%lg %*lg %*lg %lg ");
	hist3->SetLineColor(kRed+1);
	
	TGraph* hist4 = new TGraph("TheoryMarcNew/tcs_fb_diff_z_7p2_qps3p24_theta65_phi0_rximrst02_1_bv1p0_bs5p0_t.dat","%lg %*lg %*lg %lg ");
	hist4->SetLineColor(kViolet+b);
	
	TGraph* hist3p = new TGraph("TheoryMarc/tcs_fb_diff_z_7p0_qps3p24_theta65_phi0_rximrst02_1_bv1p0_bs5p0_dt_e1_t.dat","%lg %*lg %*lg %*lg %lg ");
	hist3->SetLineColor(kRed+1);
	TGraph* hist4p = new TGraph("TheoryMarc/tcs_fb_diff_z_7p0_qps3p24_theta65_phi0_rximrst02_1_bv1p0_bs5p0_t.dat","%lg %*lg %*lg %*lg %lg ");
	hist4->SetLineColor(kViolet+b);
	TGraph* hist5 = new TGraph("TheoryMarc/tcs_fb_diff_z_7p0_qps3p24_theta80_phi 0_rximrst02_1_bv1p0_bs5p0_dt_e1_t.dat","%lg %*lg %*lg %*lg %lg ");
	hist5->SetLineColor(kRed+2);
	hist5->SetLineStyle(8);
	
	//2GeV
	TGraph* hist6 = new TGraph("TheoryMarc/tcs_fb_diff_z_7p0_qps4p0_theta65_phi 0_rximrst02_1_bv1p0_bs5p0_dt_e1_t.dat","%lg %*lg %*lg %*lg %lg ");
	hist6->SetLineColor(kRed);
	hist6->SetLineStyle(8);
	//2.2GeV
	TGraph* hist7 = new TGraph("TheoryMarcNew/tcs_fb_diff_z_8p1_qps5p6_theta65_phi0_rximrst02_1_bv1p0_bs5p0_dt_e1_t.dat","%lg %*lg %*lg %lg ");
	hist7->SetLineColor(kRed+1);
	
	TGraph* hist8 = new TGraph("TheoryMarcNew/tcs_fb_diff_z_8p1_qps5p6_theta65_phi0_rximrst02_1_bv1p0_bs5p0_t.dat","%lg %*lg %*lg %lg ");
	hist8->SetLineColor(kViolet+b);
	
	  
	  
	//BSA
	TGraph* hist9 = new TGraph("TheoryMarcNew/tcs_diff_z_7p2_qps3p24_theta90_phi90_rximrst02_1_bv1p0_bs5p0_dt_e1_t.dat","%lg %*lg %lg");
	hist9->SetLineColor(kCyan);
	TGraph* hist10 = new TGraph("TheoryMarcNew/tcs_diff_z_7p2_qps3p24_theta90_phi90_rximrst02_1_bv1p0_bs5p0_dt_e1_t.dat","%lg %*lg %lg");
	hist10->SetLineColor(kViolet+b);
	
	
	
	
	
	////Results Pawel
	// BSA t
	TGraph* hist11 = new TGraph("TheoryPawelNew/ACU_E_7.29_M_1.8_phi_90_theta_92.dat","%lg %lg");
	hist11->SetLineColor(kTeal+a);
	// BSA M
	TGraph* hist15 = new TGraph("TheoryPawel/ACU_dep_M_E_6.8_t_-0.33","%lg %lg");
	hist15->SetLineColor(kTeal+a);
	
	// AFB t
	TGraph* hist12 = new TGraph("TheoryPawel/AFB_dep_t_E_6.5_M_1.69","%lg %lg");
	hist12->SetLineColor(kTeal+a);
	
	TGraph* hist13 = new TGraph("TheoryPawelNew/AFB_E_7.24_M_1.8_phi_m40_40_theta_50_80.dat","%lg %lg");
	hist13->SetLineColor(kTeal+a);
	
	TGraph* hist14 = new TGraph("TheoryPawelNew/AFB_E_8.13_M_2.25_phi_m40_40_theta_50_80.dat","%lg %lg");
	hist14->SetLineColor(kTeal+a);
	
	//AFB_M
	
	TGraph* hist16 = new TGraph("TheoryPawel/AFB_dep_M_E_6.8_t_-0.37","%lg %lg");
	hist16->SetLineColor(kTeal+a);
	
	
	TGraph* hist17 = new TGraph("TheoryPawel/AFB_dep_M_E_6.8_t_-0.33","%lg %lg");
	hist17->SetLineColor(kTeal+5);
	
	//AFB_E
	
	TGraph* hist18 = new TGraph("TheoryPawel/AFB_dep_E_M_1.79_t_-0.37","%lg %lg");
	hist18->SetLineColor(kOrange);
	
	int tailleLigne = 4;
	hist2->SetLineWidth(tailleLigne);
  	hist3->SetLineWidth(tailleLigne);
   	hist4->SetLineWidth(tailleLigne);
   	hist5->SetLineWidth(tailleLigne);
   	hist12->SetLineWidth(tailleLigne);
   	hist13->SetLineWidth(tailleLigne);
   	hist14->SetLineWidth(tailleLigne);
	hist7->SetLineWidth(tailleLigne);
   	hist8->SetLineWidth(tailleLigne);
   	hist14->SetLineWidth(tailleLigne);
	hist10->SetLineWidth(tailleLigne);
   	hist11->SetLineWidth(tailleLigne);
	   
	TFile *FileResultat = new TFile("../Base/"+file+".root");//new TFile("../results/"+file+".root");
	TGraphAsymmErrors *Resultat =(TGraphAsymmErrors*)FileResultat->Get("Graph");
	Resultat->SetLineColor(kBlue);
	Resultat->SetMarkerColor(kBlue);
	Resultat->SetMarkerSize(2);
	Resultat->SetLineWidth(3);
	Resultat->GetHistogram()->SetMaximum(max);   // along          
   	Resultat->GetHistogram()->SetMinimum(min);
   	Resultat->GetXaxis()->SetTitle(xaxis);
   	Resultat->SetTitle(title);
  	
  	double xres,yres,xres1,yres1;
  	Resultat->GetPoint(1,xres,yres); 
  	Resultat->GetPoint(bin,xres1,yres1); 
  	double minA=xres-(Resultat->GetErrorXlow(0));
  	cout<<xres1<<endl;
  	double maxA=xres1+(Resultat->GetErrorXhigh(0));
  	cout<<"min A "<<minA<<" "<<maxA<<endl;
        
        TFile *systFileGrape = new TFile("../flat/"+file+"Systflat.root");
        //TFile *systFileGrape = new TFile("../resultsGRAPEFullEff/"+file+"SystGRAPE.root");
	TGraphAsymmErrors *systGrape =(TGraphAsymmErrors*)systFileGrape->Get("Graph");
	systGrape->SetLineColor(kRed);
	systGrape->SetMarkerColor(kRed);
	systGrape->SetMarkerSize(1);
	
	TFile *systFilePositron = new TFile("../08/"+file+"SystPositron.root");
	TGraphAsymmErrors *systPositron =(TGraphAsymmErrors*)systFilePositron->Get("Graph");
	systPositron->SetLineColor(kRed);
	systPositron->SetMarkerColor(kRed);
	systPositron->SetMarkerSize(1);
	
	TFile *systFileEff = new TFile("../NoEff/"+file+"SystEff.root");
	TGraphAsymmErrors *systEff =(TGraphAsymmErrors*)systFileEff->Get("Graph");
	systEff->SetLineColor(kRed);
	systEff->SetMarkerColor(kRed);
	systEff->SetMarkerSize(1);
	
	TFile *systFileMethod = new TFile("/vol0/pierre/Bureau/Hipo4Ana/Systematics/TCSGen0/"+file+"SystMethod.root");
	TGraphAsymmErrors *systMethod =(TGraphAsymmErrors*)systFileMethod->Get("Graph");
	systMethod->SetLineColor(kRed);
	systMethod->SetMarkerColor(kRed);
	systMethod->SetMarkerSize(1);
	
	TFile *systFileBin = new TFile("../ExcluCuts/"+file+"SystExclu.root");
	TGraphAsymmErrors *systBin =(TGraphAsymmErrors*)systFileBin->Get("Graph");
	systBin->SetLineColor(kRed);
	systBin->SetMarkerColor(kRed);
	systBin->SetMarkerSize(1);
	
	TFile *systFileChi2 = new TFile("../chi2/"+file+"Systchi2.root");
	TGraphAsymmErrors *systChi2 =(TGraphAsymmErrors*)systFileChi2->Get("Graph");
	systChi2->SetLineColor(kRed);
	systChi2->SetMarkerColor(kRed);
	systChi2->SetMarkerSize(1);
	
	
	TFile *BHFile; 
	TGraphAsymmErrors *BHR;
	
	if(filename.rfind("R", 0) == 0 || filename.rfind("A", 0) == 0 || filename.rfind("B", 0) == 0){
	BHFile= new TFile("/vol0/pierre/Bureau/Hipo4Ana/SimuAnalysis/bigStat/"+file+".root");
	BHR =(TGraphAsymmErrors*)BHFile->Get("Graph");
	BHR->SetLineColor(kRed);
	BHR->SetMarkerColor(kRed);
	BHR->SetMarkerSize(6);
	BHR->SetMarkerStyle(22);
	BHR->SetLineWidth(3);
	}
	
	
  	
  	TGraphAsymmErrors *Systematic =new TGraphAsymmErrors(bin);
  	for(int i=0; i<bin;i++){
  	
  	double x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7;
  	systEff->GetPoint(i,x1,y1); 
  	systGrape->GetPoint(i,x2,y2); 
  	systPositron->GetPoint(i,x3,y3); 
  	systMethod->GetPoint(i,x4,y4); 
  	systBin->GetPoint(i,x6,y6); 
  	systChi2->GetPoint(i,x7,y7); 
  	Resultat->GetPoint(i,x5,y5); 
  	double errorStat = Resultat->GetErrorY(i); 
  	//double value1 = y1;
  	

	double e1Err2=0.0;if(i!=0)e1Err2=systGrape->GetErrorYlow(i-1);
  	double e1Err3=0.0;if(i!=0)e1Err3=systPositron->GetErrorYlow(i-1);
  	double e1Err4=0.0;if(i!=0)e1Err4=systMethod->GetErrorYlow(i-1);
  	//double e1Err5=0.0;if(i!=0)e1Err5=systBin->GetErrorYlow(i-1);
  	double e1Err6=0.0;if(i!=0)e1Err6=systChi2->GetErrorYlow(i-1);
  	double e1Err7=0.0;if(i!=0)e1Err7=systEff->GetErrorYlow(i-1);


  	double e1Err2H=0.0;if(i!=0)e1Err2H=systGrape->GetErrorYhigh(i-1);
  	double e1Err3H=0.0;if(i!=0)e1Err3H=systPositron->GetErrorYhigh(i-1);
  	double e1Err4H=0.0;if(i!=0)e1Err4H=systMethod->GetErrorYhigh(i-1);
  	//double e1Err5H=0.0;if(i!=0)e1Err5H=systBin->GetErrorYhigh(i-1);
  	double e1Err6H=0.0;if(i!=0)e1Err6H=systChi2->GetErrorYhigh(i-1);
  	double e1Err7H=0.0;if(i!=0)e1Err7H=systEff->GetErrorYhigh(i-1);

	double e2Err2=0.0;if(i!=(bin-1))e2Err2=systGrape->GetErrorYlow(i+1);
  	double e2Err3=0.0;if(i!=(bin-1))e2Err3=systPositron->GetErrorYlow(i+1);
  	double e2Err4=0.0;if(i!=(bin-1))e2Err4=systMethod->GetErrorYlow(i+1);
  	//double e2Err5=0.0;if(i!=(bin-1))e2Err5=systBin->GetErrorYlow(i+1);
  	double e2Err6=0.0;if(i!=(bin-1))e2Err6=systChi2->GetErrorYlow(i+1);
  	double e2Err7=0.0;if(i!=(bin-1))e2Err7=systEff->GetErrorYlow(i+1);
  	

  	double e2Err2H=0.0;if(i!=(bin-1))e2Err2H=systGrape->GetErrorYhigh(i+1);
  	double e2Err3H=0.0;if(i!=(bin-1))e2Err3H=systPositron->GetErrorYhigh(i+1);
  	double e2Err4H=0.0;if(i!=(bin-1))e2Err4H=systMethod->GetErrorYhigh(i+1);
  	//double e2Err5H=0.0;if(i!=(bin-1))e2Err5H=systBin->GetErrorYhigh(i+1);
  	double e2Err6H=0.0;if(i!=(bin-1))e2Err6H=systChi2->GetErrorYhigh(i+1);
  	double e2Err7H=0.0;if(i!=(bin-1))e2Err7H=systEff->GetErrorYhigh(i+1);
	 
	double AverageNB = 3.; if(i==(bin-1) || i==0)AverageNB=2.;
  	
  	double CErr2=systGrape->GetErrorYlow(i);
  	double CErr3=systPositron->GetErrorYlow(i);
  	double CErr4=systMethod->GetErrorYlow(i);
  	//double CErr5=systBin->GetErrorYlow(i);
  	double CErr6=systChi2->GetErrorYlow(i);
  	double CErr7=systEff->GetErrorYlow(i);
  	

  	double CErr2H=systGrape->GetErrorYhigh(i);
  	double CErr3H=systPositron->GetErrorYhigh(i);
  	double CErr4H=systMethod->GetErrorYhigh(i);
  	//double CErr5H=systBin->GetErrorYhigh(i);
  	double CErr6H=systChi2->GetErrorYhigh(i);
  	double CErr7H=systEff->GetErrorYhigh(i);
  	
  	double Err2=(e1Err2+e2Err2+CErr2)/AverageNB;
  	double Err3=(e1Err3+e2Err3+CErr3)/AverageNB;
  	double Err4=(e1Err4+e2Err4+CErr4)/AverageNB;
  	//double Err5=(e1Err5+e2Err5+CErr5)/AverageNB;
  	double Err6=(e1Err6+e2Err6+CErr6)/AverageNB;
  	double Err7=(e1Err7+e2Err7+CErr7)/AverageNB;
  	

  	double Err2H=(e1Err2H+e2Err2H+CErr2H)/AverageNB;
  	double Err3H=(e1Err3H+e2Err3H+CErr3H)/AverageNB;
  	double Err4H=(e1Err4H+e2Err4H+CErr4H)/AverageNB;
  	//double Err5H=(e1Err5H+e2Err5H+CErr5H)/AverageNB;
  	double Err6H=(e1Err6H+e2Err6H+CErr6H)/AverageNB;
  	double Err7H=(e1Err7H+e2Err7H+CErr7H)/AverageNB;
  	
  	//Special treatment for exclu cut syst
  	cout<<" y 5 "<<y5<<" "<<(systBin->GetErrorYlow(0))<<endl;
  	double Err5=(systBin->GetErrorYlow(0))*y5;
  	double Err5H=(systBin->GetErrorYhigh(0))*y5;
  	
  	/*
  	double Err2=systGrape->GetErrorYlow(i);
  	double Err3=systPositron->GetErrorYlow(i);
  	double Err4=systMethod->GetErrorYlow(i);
  	double Err5=systBin->GetErrorYlow(i);
  	double Err6=systChi2->GetErrorYlow(i);
  	double Err7=systEff->GetErrorYlow(i);
  	

  	double Err2H=systGrape->GetErrorYhigh(i);
  	double Err3H=systPositron->GetErrorYhigh(i);
  	double Err4H=systMethod->GetErrorYhigh(i);
  	double Err5H=systBin->GetErrorYhigh(i);
  	double Err6H=systChi2->GetErrorYhigh(i);
  	double Err7H=systEff->GetErrorYhigh(i);
  	*/
  	
  	
  	cout<<"decomp. low "<<" Grape "<<Err2<<" Posi "<<Err3<<" Method "<<Err4<<endl;
  	cout<<"decomp. "<< "Grape "<<Err2H<<" Posi "<<Err3H<<" Method "<<Err4H<<endl;
  	
  	double width1=systEff->GetErrorXlow(i);
  	double width2=systEff->GetErrorXhigh(i);
  	Systematic->SetPoint(i,x5,y5);
  	double errorYlow=TMath::Sqrt(Err2*Err2+Err3*Err3+Err4*Err4+Err5*Err5+Err6*Err6+Err7*Err7);
  	double errorYhigh=TMath::Sqrt(Err2H*Err2H+Err3H*Err3H+Err4H*Err4H+Err5H*Err5H+Err6H*Err6H+Err7H*Err7H);
  	cout<<y1<<" "<<errorYlow<<" "<<errorYhigh<<endl;
  	Systematic->SetPointError(i,0.02,0.02,errorYlow,errorYhigh);
  	//Systematic->SetPointError(i,width1,width2,errorYlow,errorYhigh);
  	
  	  	
  	
  	}
  	Systematic->SetFillColorAlpha(kBlue, 0.35);
  	//Systematic->SetLineColor(kBlue);
  	
  	TCanvas *cancG0  = new TCanvas("canG0","Assym",3000,1800);//1500,2000);//1500,2000);    
  	TPad *pad1 = new TPad("pad1", "pad1", 0, 0., 1, 1.0);
  	pad1->SetBottomMargin(0.11); // Upper and lower plot are joined
   	//pad1->SetGridx();         // Vertical grid
  	//pad1->SetGridy();
  	pad1->Draw();             // Draw the upper pad: pad1
 	pad1->cd();               // pad1 becomes the current pad
 	
 	Resultat->SetMarkerSize(5);
 	Resultat->GetHistogram()->SetMaximum(max);
 	Resultat->GetHistogram()->SetMinimum(min);
        //Resultat->GetHistogram()->SetMinimum(-0.13);
 	Resultat->GetXaxis()->SetRangeUser(minA,maxA);
 	Resultat->SetTitle("");
  	Resultat->Draw("A P ");
  	
  	if(filename.rfind("R", 0) == 0 || filename.rfind("A", 0) == 0 || filename.rfind("B", 0) == 0)BHR->Draw(" P ");
  	
  	auto legend = new TLegend();
  	if(title1.rfind("AFB_M", 0) == 0 ){ legend = new TLegend(0.10,0.3,0.6,0.11); }
  	else legend = new TLegend(0.42,0.30,0.9,0.11); 
  	legend->SetNColumns(2);
  	legend->SetBorderSize(1);
  	legend->SetTextFont(42);
  	//legend->SetBorderSize(0);
  	//legend->SetFillStyle(0);
         legend->AddEntry(Resultat,"DATA","lep");  
          legend->AddEntry(Systematic,"Tot. Syst.","f1");  
      if(filename.rfind("R", 0) == 0 || filename.rfind("A", 0) == 0 || filename.rfind("B", 0) == 0) legend->AddEntry(BHR,"BH","lp");  
  	
  	if(title1.rfind("1.69", 0) == 0 || title1.rfind("1.79", 0) == 0 ){
  	//hist1->Draw(" L ");
  	//hist2->Draw(" L ");
  	hist3->Draw(" L ");
  	//hist3p->Draw(" L ");
  	hist4->Draw(" L ");
  	hist12->SetLineStyle(4);
  	hist13->SetLineStyle(4);
  	hist4->SetLineStyle(7);
  	//hist5->Draw(" L ");
  	if(title1.rfind("1.69", 0) == 0)hist12->Draw(" L ");
  	if(title1.rfind("1.79", 0) == 0)hist13->Draw(" L ");
  	if(title1.rfind("1.69", 0) == 0)legend->AddEntry(hist12,"GK, no D-term","l");  
  	if(title1.rfind("1.79", 0) == 0)legend->AddEntry(hist13,"GK, no D-term","l");  
  	//legend->AddEntry(hist1,"VGG, 1.5 GeV/65#circ","l");  
  	//legend->AddEntry(hist2,"VGG, 1.8 GeV/50#circ","l");  
  	//legend->AddEntry(hist3,"VGG, 1.8 GeV/65#circ","l");  
  	//legend->AddEntry(hist4,"VGG, 1.8 GeV/65#circ/No D.","l");  
  	legend->AddEntry(hist3,"VGG","l");  
  	legend->AddEntry(hist4,"VGG, no D-term","l");  
  	//legend->AddEntry(hist5,"VGG, 1.8 GeV/80#circ","l");  
  	}
  	
  	if(title1.rfind("2.21", 0) == 0 ){
  	//hist6->Draw(" L ");
  	hist7->Draw(" L ");
  	hist8->Draw(" L ");
  	hist14->Draw(" L ");
  	//legend->AddEntry(hist6,"VGG, 2 GeV/65#circ","l");  
  	hist14->SetLineStyle(4);
  	hist8->SetLineStyle(7);
  	legend->AddEntry(hist14,"GK, no D-term","l");
  	legend->AddEntry(hist7,"VGG","l");    
  	legend->AddEntry(hist8,"VGG, no D-term","l"); 
  	
  	 
  	
  	}
  	
  	if(title1.rfind("BSA_t", 0) == 0 ){
  	legend->SetNColumns(2);
  	//hist9->Draw(" L ");
  	hist10->SetLineStyle(7);
  	hist11->SetLineStyle(4);
  	hist10->Draw(" L ");
  	hist11->Draw(" L ");
  	//legend->AddEntry(hist9,"VGG, b_{SEA}=1","l");  
  	legend->AddEntry(hist11,"GK","l");  
  	legend->AddEntry(hist10,"VGG","l"); 
  	Resultat->GetYaxis()->SetTitle("A_{#odot U}");
  	
  	}
  	
  	if(title1.rfind("BSA_M", 0) == 0 ){
  	hist15->Draw(" L ");
  	legend->AddEntry(hist15,"GK","l");  
  	}
  	
  	
  	if(title1.rfind("AFB_M", 0) == 0 ){
  	legend->SetNColumns(2);
  	hist16->Draw(" L ");
  	legend->AddEntry(hist16,"GK, No D.","l");  
  	hist17->Draw(" L ");
  	legend->AddEntry(hist17,"GK, No D., -t=0.33 GeV^{2}","l");  
  	}
  	
  	if(title1.rfind("AFB_Eg", 0) == 0 ){
  	hist18->Draw(" L ");
  	legend->AddEntry(hist18,"GK, No D.","l");  
  	}
  	
  	Systematic->Draw("E2 same ");
  	//img->Draw("same");
  	
  	   TText *t = new TText();
 	t->SetTextAlign(22);
   t->SetTextColorAlpha(kGray, 0.50);
  // t->SetTextColor(kRed+2);
   t->SetTextFont(40);
   t->SetTextSize(0.25);
   t->SetTextAngle(25);
  	
 	// t->DrawTextNDC(.5,.53, "Preliminary");
  
	/*
	TText *t1 = new TText();
 	t1->SetTextAlign(22);
  	t1->SetTextColor(kBlack);
  	t1->SetTextFont(40);
  	t1->SetTextSize(0.06);
  	if(title1.rfind("AFB", 0) == 0)t1->DrawTextNDC(.8,.8, "60^#circ<#theta<80^#circ");
  	*/
  	
  	TLatex T1, T2;
  	 if(title1.rfind("2.21", 0) == 0 ){T1.DrawLatexNDC(.14,.85, "Forward angular bin:");T2.DrawLatexNDC(.14,.8, "#theta #in [50#circ,80#circ], #phi #in [-40#circ, 40#circ]");}
  	 if(title1.rfind("1.79", 0) == 0 || title1.rfind("1.69", 0) == 0  ){T1.DrawLatexNDC(.5,.85, "Forward angular bin:");T2.DrawLatexNDC(.5,.8, "#theta #in [50#circ,80#circ], #phi #in [-40#circ, 40#circ]");}
  	
	//legend->AddEntry(Resultat3,"No Acc","l");     	    	
   	legend->Draw("");
   	 
   	  //Resultat->Draw("P");
   	
   

  	cancG0->SaveAs(output+"0_NoSyst.pdf");
  	
  	
  	gApplication->Terminate();

return 0;
}
