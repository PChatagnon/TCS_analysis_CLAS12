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


int Plot(TString file, TString output, string filename, int bin, TString yaxis, TString xaxis,TString yaxis1, TString xaxis1, double min, double max, TString title) {
	gROOT->SetBatch(kTRUE);
	
	
	ofstream myfile;
  	myfile.open (filename+".txt");
  	myfile << "\\begin{center}\n";
  	myfile << "\\begin{tabular}{ |c|c|c|c|c| } \n";
  	myfile << "\\hline\n";
  	myfile << "$"<<xaxis1<<"$"<<" & "<<"$"<<yaxis1<<"$"<<" & Stat. error & Low Syst. Error & High Syst. Error \\\\ \n";
  	myfile << "\\hline\n";
	//Int_t argc=gApplication->Argc();
	//char** argv=gApplication->Argv();
	
	//TString file = TString(argv[3]);//"BSAversust";//
	cout<<file<<endl;
	//TString output = TString(argv[5]);//"BSA";//
	
	
	TImage *img = TImage::Open("Claslogo.pdf");
	
        gStyle->SetPaintTextFormat("4.1f");
        gStyle->SetOptStat(1);
        gStyle->SetPalette(55);
        gStyle->SetLabelSize(.04, "xyz");
        gStyle->SetTitleSize(.04, "xyz");
        gStyle->SetTitleSize(.07,"t");
        gStyle->SetFrameLineWidth(2);
        gStyle->SetLineWidth(2);
        gStyle->SetHistLineWidth(2);
        gStyle->SetMarkerStyle(13);
        gStyle->SetTitleW(0.8);  //per cent of the pad width
        gStyle->SetTitleH(0.1); //per cent of the pad height

	
	   
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
	BHR->SetMarkerSize(4);
	BHR->SetLineWidth(3);
	}
	
	
  	
  	TGraphAsymmErrors *Systematic =new TGraphAsymmErrors(bin);
  	for(int i=0; i<bin;i++){
  	
  	double x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7,x8,y8;
  	systEff->GetPoint(i,x1,y1); 
  	systGrape->GetPoint(i,x2,y2); 
  	systPositron->GetPoint(i,x3,y3); 
  	systMethod->GetPoint(i,x4,y4); 
  	systBin->GetPoint(i,x6,y6); 
  	systChi2->GetPoint(i,x6,y6); 
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
  	double CErr6=systChi2->GetErrorYlow(i);
  	double CErr7=systEff->GetErrorYlow(i);
  	
	//cout<<file<<" "<<i<<" "<<e1Err5<<" "<<e2Err5<<" "<<CErr5<<" "<<endl;

  	double CErr2H=systGrape->GetErrorYhigh(i);
  	double CErr3H=systPositron->GetErrorYhigh(i);
  	double CErr4H=systMethod->GetErrorYhigh(i);
  	double CErr6H=systChi2->GetErrorYhigh(i);
  	double CErr7H=systEff->GetErrorYhigh(i);
  	
  	double Err2=(e1Err2+e2Err2+CErr2)/AverageNB;
  	double Err3=(e1Err3+e2Err3+CErr3)/AverageNB;
  	double Err4=(e1Err4+e2Err4+CErr4)/AverageNB;
  	double Err6=(e1Err6+e2Err6+CErr6)/AverageNB;
  	double Err7=(e1Err7+e2Err7+CErr7)/AverageNB;
  	

  	double Err2H=(e1Err2H+e2Err2H+CErr2H)/AverageNB;
  	double Err3H=(e1Err3H+e2Err3H+CErr3H)/AverageNB;
  	double Err4H=(e1Err4H+e2Err4H+CErr4H)/AverageNB;
  	double Err6H=(e1Err6H+e2Err6H+CErr6H)/AverageNB;
  	double Err7H=(e1Err7H+e2Err7H+CErr7H)/AverageNB;
  	
  	//Special treatment for exclu cut syst
  	cout<<" y 5 "<<y5<<" "<<(systBin->GetErrorYlow(0))<<endl;
  	double Err5=(systBin->GetErrorYlow(0))*y5;
  	double Err5H=(systBin->GetErrorYhigh(0))*y5;
  	
  	
  	double width1=systEff->GetErrorXlow(i);
  	double width2=systEff->GetErrorXhigh(i);
  	Systematic->SetPoint(i,x1,y5);
  	double errorYlow=TMath::Sqrt(Err2*Err2+Err3*Err3+Err4*Err4+Err5*Err5+Err6*Err6+Err7*Err7);
  	double errorYhigh=TMath::Sqrt(Err2H*Err2H+Err3H*Err3H+Err4H*Err4H+Err5H*Err5H+Err6H*Err6H+Err7H*Err7H);
  	cout<<y1<<" "<<errorYlow<<" "<<errorYhigh<<endl;
  	Systematic->SetPointError(i,width1,width2,errorYlow,errorYhigh);
  	
  	  	myfile << setprecision (3) << x5<<" & "<<y5<<" & "<<errorStat<<" & "<<errorYlow <<" & "<< errorYhigh <<" \\\\"<<"\n";
  	
  	}
  	Systematic->SetFillColorAlpha(kBlack, 0.35);
  	
  	TCanvas *cancG0  = new TCanvas("canG0","Assym",6000,8000);//1500,2000);//1500,2000);   
  	TPad *pad0 = new TPad("pad0", "pad0", 0, 0.0, 1, 1.0); 
  	TPad *pad1 = new TPad("pad1", "pad1", 0, 0.5, 1, 1.0);
  	pad1->SetBottomMargin(0.11); // Upper and lower plot are joined
   	//pad1->SetGridx();         // Vertical grid
  	//pad1->SetGridy();
  	pad1->Draw();             // Draw the upper pad: pad1
 	pad1->cd();               // pad1 becomes the current pad
 	
 	//Resultat->UseCurrentStyle();
 		Resultat->SetMarkerSize(5);
 	Resultat->GetHistogram()->SetMaximum(0.7);
        //Resultat->GetHistogram()->SetMinimum(-0.13);
 	Resultat->GetXaxis()->SetRangeUser(minA,maxA);
  	Resultat->Draw("A P ");
  	Systematic->Draw("E2 same ");
  	//img->Draw("same");
  	if(filename.rfind("R", 0) == 0 || filename.rfind("A", 0) == 0 || filename.rfind("B", 0) == 0)BHR->Draw(" P ");
  	   TText *t = new TText();
 	t->SetTextAlign(22);
   t->SetTextColorAlpha(kGray, 0.50);
  // t->SetTextColor(kRed+2);
   t->SetTextFont(40);
   t->SetTextSize(0.25);
   t->SetTextAngle(25);
  	
 	 t->DrawTextNDC(.5,.53, "Preliminary");
  

  	
  	auto legend = new TLegend(0.8,0.25,0.9,0.11); 
  	//legend->SetBorderSize(0);
  	//legend->SetFillStyle(0);
         legend->AddEntry(Resultat,"DATA","lp");  
          legend->AddEntry(Systematic,"Tot. Syst.","f1");  
      if(filename.rfind("R", 0) == 0 || filename.rfind("A", 0) == 0 || filename.rfind("B", 0) == 0) legend->AddEntry(BHR,"BH","lp");  
	//legend->AddEntry(Resultat3,"No Acc","l");     	    	
   	legend->Draw("");
   	
   	pad0->cd(); 
   	//pad0->SetGridx();         // Vertical grid
  	//pad0->SetGridy(); 
   	Resultat->Draw("A P ");
   	//Systematic->Smooth(1);
  	Systematic->Draw("E2 same ");
  	BHR->Draw(" P ");
   	t->DrawTextNDC(.5,.53, "Preliminary");
   	legend->Draw("");
   	
   	
   	//gStyle->SetTitleH(0.3);
   	
  	cancG0->cd();
  	TPad *pad2 = new TPad("pad2", "pad2", 0, 0.00, 1, 0.1);
   	pad2->SetTopMargin(0);
   	pad2->SetGridy();
   	pad2->SetBottomMargin(0.1);
   	//pad2->SetGridx(); // vertical grid
   	pad2->Draw();
   	pad2->cd();       // pad2 becomes the current pad
   	
 	systGrape->SetTitle("");// 	systGrape->SetTitle("Acceptance");
   	//systGrape->SetTitle("Acceptance");
   	//systGrape->GetYaxis()->SetTitle("Acceptance");
   	systGrape->GetHistogram()->SetTitleSize(1.);
   	//systGrape->SetTitleSize(25.);
   	systGrape->GetHistogram()->SetMaximum(0.13);
        systGrape->GetHistogram()->SetMinimum(-0.13);
   	systGrape->SetFillColorAlpha(kRed, 0.5);
   	systGrape->GetXaxis()->SetRangeUser(minA,maxA);
  	systGrape->Draw(" A E2");
  	systGrape->GetXaxis()->SetLabelSize(0.08);
	systGrape->GetYaxis()->SetLabelSize(0.08);
  	
  	 TText *t1 = new TText();
 	t1->SetTextAlign(22);
   t1->SetTextColorAlpha(kBlack, 1);
  // t->SetTextColor(kRed+2);
   t1->SetTextFont(40);
   t1->SetTextSize(0.18);
  	
 	 t1->DrawTextNDC(.5,.9, "Acceptance");
 	 t1->SetTextSize(0.2);
  	cancG0->cd();
  	TPad *pad6 = new TPad("pad6", "pad6", 0, 0.1, 1, 0.18);
   	pad6->SetTopMargin(0.);
   	pad6->SetGridy();
   	pad6->SetBottomMargin(0.0);
   	//pad6->SetGridx(); // vertical grid
   	pad6->Draw();
   	pad6->cd();       // pad2 becomes the current pad
	systBin->SetTitle("");//	systBin->SetTitle("Exclusivity");
   	systBin->GetHistogram()->SetMaximum(0.35);
        systBin->GetHistogram()->SetMinimum(-0.35);
   	systBin->SetFillColorAlpha(kViolet, 0.5);
   	systBin->GetXaxis()->SetRangeUser(minA,maxA);
  	systBin->Draw(" A E2");
  	systBin->GetXaxis()->SetLabelSize(0.08);
	systBin->GetYaxis()->SetLabelSize(0.08);
  	 t1->DrawTextNDC(.5,.9, "Exclusivity (relative uncertainty)");
  	
  	cancG0->cd();
  	TPad *pad4 = new TPad("pad4", "pad4", 0, 0.18, 1, 0.26);
   	pad4->SetTopMargin(0);
   	pad4->SetGridy();
   	pad4->SetBottomMargin(0.0);
   	//pad4->SetGridx(); // vertical grid
   	pad4->Draw();
   	pad4->cd();       // pad2 becomes the current pad
   	systEff->SetTitle("");//systEff->SetTitle("Efficiency");
   	systEff->GetHistogram()->SetMaximum(0.13);
        systEff->GetHistogram()->SetMinimum(-0.13);
   	systEff->SetFillColorAlpha(kOrange, 0.5);
   	systEff->GetXaxis()->SetRangeUser(minA,maxA);
  	systEff->Draw(" A E2");
  	systEff->GetXaxis()->SetLabelSize(0.08);
	systEff->GetYaxis()->SetLabelSize(0.08);
  	t1->DrawTextNDC(.5,.9, "Efficiency");
	cancG0->cd();
  	TPad *pad3 = new TPad("pad3", "pad3", 0, 0.26, 1, 0.34);
   	pad3->SetTopMargin(0);
   	pad3->SetGridy();
   	pad3->SetBottomMargin(0.0);
   	//pad3->SetGridx(); // vertical grid
   	pad3->Draw();
   	pad3->cd();       // pad2 becomes the current pad
   	systPositron->SetTitle("");//systPositron->SetTitle("Positron ID");
   	systPositron->GetHistogram()->SetMaximum(0.13);
        systPositron->GetHistogram()->SetMinimum(-0.13);
   	systPositron->SetFillColorAlpha(kBlue, 0.5);
   	systPositron->GetXaxis()->SetRangeUser(minA,maxA);
  	systPositron->Draw(" A E2");
  	systPositron->GetXaxis()->SetLabelSize(0.08);
	systPositron->GetYaxis()->SetLabelSize(0.08);
  	  	t1->DrawTextNDC(.5,.9, "Positron ID");
  
	cancG0->cd();
  	TPad *pad7 = new TPad("pad7", "pad7", 0, 0.34, 1, 0.42);
   	pad7->SetTopMargin(0.);
   	pad7->SetGridy();
   	pad7->SetBottomMargin(0.0);
   	//pad7->SetGridx(); // vertical grid
   	pad7->Draw();
   	pad7->cd();       // pad2 becomes the current pad
   	systChi2->SetTitle("");//systChi2->SetTitle("#Chi^{2} proton");
   	systChi2->GetHistogram()->SetMaximum(0.13);
        systChi2->GetHistogram()->SetMinimum(-0.13);
   	systChi2->SetFillColorAlpha(kGreen, 0.5);
   	systChi2->GetXaxis()->SetRangeUser(minA,maxA);
  	systChi2->Draw(" A E2");
  	systChi2->GetXaxis()->SetLabelSize(0.08);
	systChi2->GetYaxis()->SetLabelSize(0.08);
	
	
	t1->DrawTextNDC(.5,.9, "Proton identification");
	cancG0->cd();
  	TPad *pad5 = new TPad("pad5", "pad5", 0, 0.42, 1, 0.5);
   	pad5->SetTopMargin(0.);
   	pad5->SetGridy();
   	pad5->SetBottomMargin(0.0);
   	//pad5->SetGridx(); // vertical grid
   	pad5->Draw();
   	pad5->cd();       // pad2 becomes the current pad
   	systMethod->SetTitle("");//systMethod->SetTitle("Method");
   	systMethod->GetHistogram()->SetMaximum(0.13);
        systMethod->GetHistogram()->SetMinimum(-0.13);
   	systMethod->SetFillColorAlpha(kYellow, 0.5);
   	systMethod->GetXaxis()->SetRangeUser(minA,maxA);
  	systMethod->Draw(" A E2");
  	systMethod->GetXaxis()->SetLabelSize(0.08);
	systMethod->GetYaxis()->SetLabelSize(0.08);
	t1->DrawTextNDC(.5,.9, "Method");
	
	
	
	
	/*cancG0->cd();
	TPad *l = new TPad("l","l",0.4,0.4,0.8,0.8);
   	//gPad->cd(0);
   	l->SetFillStyle(4000);
   	l->SetFrameFillColor(0);
   	l->SetFrameLineColor(0);
   l->SetFrameFillStyle(0);
   l->SetFrameBorderMode(0);
  	l->Draw();
   	l->cd();
  	img->Draw();*/
	

  	cancG0->SaveAs(output+"0.pdf");
  	
  	TCanvas *cancG01  = new TCanvas("canG01","Assym1",3000,1800);//1500,2000);//1500,2000);  
  	pad0->SetBottomMargin(0.11);  
  	pad0->Draw();             // Draw the upper pad: pad1
  	cancG01->SaveAs(output+"0_sansSyst.pdf");
  	
  	myfile << "\\hline\n";
  	myfile << "\\end{tabular} \n";
  	myfile << "\\end{center}\n";
  	myfile.close();
  	
  	Resultat->SaveAs(output+"File.root");
  	
  	gApplication->Terminate();

return 0;
}
