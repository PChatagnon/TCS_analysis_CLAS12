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
#include "bib/TCSclass.h"
#include "bib/TCSfunc.h"
using namespace std;


Particle RadiativeCorr(Particle vParticle, Particle Photons[], double thetaWin,double thetaWin1, int np){

	for(int i=0;i<np;i++){
		if(Photons[i].Vector.Angle(vParticle.Vector.Vect())*TMath::RadToDeg()<thetaWin && abs(Photons[i].Vector.Theta()-vParticle.Vector.Theta())*TMath::RadToDeg()<thetaWin1 && Photons[i].Vector.P()>0.){
		vParticle.Vector=(vParticle.Vector+Photons[i].Vector);
		//cout<<"in rad corr "<<endl;
		}
	}
 return vParticle;
}

int pippim() {
	gStyle->SetOptFit(1111);
	//gStyle->SetOptStat(10);
	gStyle->SetPalette(55);	
	gStyle->SetLabelSize(.05, "xyz");
	gStyle->SetTitleSize(.05, "xyz");
	gStyle->SetTitleSize(.07,"t");
	gStyle->SetFrameLineWidth(1); 
	gStyle->SetLineWidth(1);
	gStyle->SetHistLineWidth(1);
	gStyle->SetMarkerStyle(13);                              
	//gStyle->SetMarkerStyle(13);  
	//TFile *file = new TFile("FilterPP5532.root");

	gROOT->SetBatch(kTRUE);
	//tc->Add("out_out_Proton_67.root");
	//tc->Add("SimuPPFinal.root");
	//tc->Add("SimuFinalFinal.root");
	//TTree *Tree =(TTree*)file->Get("hipo2root");
	//tc->GetEntry(0);                                                                                                                                                       
	// TTree *Tree = (TTree*)tc->GetTree();          


	TH2F *plotQ2W = new TH2F("plotQ2","plotQ2",50,0.8,4.2,50,1.8,4.2);
	plotQ2W->GetXaxis()->SetTitle("Q2");

	TH2F *MCplotQ2W = new TH2F("MCplotQ2","MCplotQ2",100,0,10,100,0,10);
	MCplotQ2W->GetXaxis()->SetTitle("Q2");


	TH1F *InvariantMassPP = new TH1F("InvariantMassPP","InvariantMassPP",100,0,2);
	InvariantMassPP->GetXaxis()->SetTitle("M (GeV)");
	TH1F *MMass = new TH1F("MMass","Missing mass proton",100,-5,5);
	MMass->GetXaxis()->SetTitle("M^{2} (GeV^{2})");
	TH2F *MMassTheta = new TH2F("MMassTheta","Theta vs Missing mass proton",100,0,2,100,5,70);
	MMassTheta->GetXaxis()->SetTitle("M^{2} (GeV^{2})");
	TH1F *MMassCentral = new TH1F("MMassC","Missing mass proton in central",100,0.6,1.1);
	MMassCentral->GetXaxis()->SetTitle("M (GeV)");
	TH1F *MMassCentralRaw = new TH1F("MMassCR","Missing mass proton in central Raw",100,0,2);
	MMassCentralRaw->GetXaxis()->SetTitle("M (GeV)");
	TH2F *MMassW = new TH2F("MMassW","W vs Missing mass proton",100,0,2,100,2,5);
	MMassW->GetXaxis()->SetTitle("M^{2} (GeV^{2})");
	TH2F *MMassThetaE = new TH2F("MMassThetaE","ThetaE vs Missing mass proton",100,0,2,100,5,15);
	MMassThetaE->GetXaxis()->SetTitle("M^{2} (GeV^{2})");
	TH2F *MMassThetaPip = new TH2F("MMassThetaPip","ThetaPip vs Missing mass proton",100,0,2,100,0,50);
	MMassThetaPip->GetXaxis()->SetTitle("M^{2} (GeV^{2})");
	TH2F *MMassThetaPim = new TH2F("MMassThetaPim","ThetaPim vs Missing mass proton",100,0,2,100,0,50);
	MMassThetaPim->GetXaxis()->SetTitle("M^{2} (GeV^{2})");
	TH2F *MMassPE = new TH2F("MMassPE","PE vs Missing mass proton",100,0,2,100,0,10);
	MMassPE->GetXaxis()->SetTitle("M^{2} (GeV^{2})");
	TH2F *MMassPPip = new TH2F("MMassPPip","PPip vs Missing mass proton",100,0,2,100,0,10);
	MMassPPip->GetXaxis()->SetTitle("M^{2} (GeV^{2})");
	TH2F *MMassPPim = new TH2F("MMassPPim","PPim vs Missing mass proton",100,0,2,100,0,10);
	MMassPPim->GetXaxis()->SetTitle("M^{2} (GeV^{2})");
	TH2F *MMassMmP = new TH2F("MMassMmP","Missing P vs Missing mass proton",100,0,2,100,0,5);
	MMassMmP->GetXaxis()->SetTitle("M^{2} (GeV^{2})");

	TH2F *ThetaPW = new TH2F("ThetaPW","ThetaPW",100,0,90,100,0,5);
	ThetaPW->GetXaxis()->SetTitle("Theta P");
	TH2F *ThetaPQ2 = new TH2F("ThetaPQ2","ThetaPQ2",100,0,90,100,0,5);
	ThetaPQ2->GetXaxis()->SetTitle("Theta P");

	TH1F *Chi2Pp = new TH1F("Chi2Pp","Chi2Pp",100,-10,5);
	Chi2Pp->GetXaxis()->SetTitle("Chi2Pp");
	TH1F *Chi2Pm = new TH1F("Chi2Pm","Chi2Pm",100,-10,5);
	Chi2Pm->GetXaxis()->SetTitle("Chi2Pm");
	TH1F *Chi2E = new TH1F("Chi2E","Chi2E",100,-10,5);
	Chi2E->GetXaxis()->SetTitle("Chi2E");

	TH2F *MMassChi2E = new TH2F("MMassChi2E","Chi2E vs MMass",100,0,2,100,-10,8);
	MMassChi2E->GetXaxis()->SetTitle("MMass");
	TH2F *MMassChi2Pp = new TH2F("MMassChi2Pp","Chi2Pp vs MMass",100,0,2,100,-10,8);
	MMassChi2Pp->GetXaxis()->SetTitle("MMass");
	TH2F *MMassChi2Pm = new TH2F("MMassChi2Pm","Chi2Pm vs MMass",100,0,2,100,-10,8);
	MMassChi2Pm->GetXaxis()->SetTitle("MMass");



	TH1F *MMassTot = new TH1F("MMassTot","Missing mass total",100,-0.2,0.2);
	MMassTot->GetXaxis()->SetTitle("M^{2} (GeV^{2})");
	TH1F *MMassPTot = new TH1F("MMassPTot","Missing mass total proton",100,0,2);
	MMassPTot->GetXaxis()->SetTitle("M^{2} (GeV^{2})");
	TH1F *MMassPpTot = new TH1F("MMassPpTot","Missing mass pi+ total",100,-1,1);
	MMassPpTot->GetXaxis()->SetTitle("M^{2} (GeV^{2})");
	TH1F *MMassPmTot = new TH1F("MMassPmTot","Missing mass pi- total",100,-1,1);
	MMassPmTot->GetXaxis()->SetTitle("M^{2} (GeV^{2})");

	TH1F *bMMassTot = new TH1F("bMMassTot","Missing mass total",100,-0.1,0.1);
	bMMassTot->GetXaxis()->SetTitle("M^{2} (GeV^{2})");
	bMMassTot->SetLineColor(3);
	TH1F *bMMassTot1 = new TH1F("bMMassTot1","Missing mass total",100,-0.1,0.1);
	bMMassTot1->GetXaxis()->SetTitle("M^{2} (GeV^{2})");
	bMMassTot1->SetLineColor(3);
	
	
	TH1F *bPitheta = new TH1F("bPitheta","",100,-5.,5.);
	bPitheta->GetXaxis()->SetTitle("#Delta #theta");
	bPitheta->SetLineColor(3);
	TH1F *bPitheta1 = new TH1F("bPitheta1","",100,-5.,5.);
	bPitheta1->GetXaxis()->SetTitle("#Delta #theta");
	bPitheta1->SetLineColor(3);
	
	TH1F *bMMassPTot = new TH1F("bMMassPTot","Missing mass total proton",100,0,2);
	bMMassPTot->GetXaxis()->SetTitle("M^{2} (GeV^{2})");
	bMMassPTot->SetLineColor(3);
	TH1F *bMMassPpTot = new TH1F("bMMassPpTot","Missing mass pi+ total",100,-1,1);
	bMMassPpTot->GetXaxis()->SetTitle("M^{2} (GeV^{2})");
	bMMassPpTot->SetLineColor(3);
	TH1F *bMMassPmTot = new TH1F("bMMassPmTot","Missing mass pi- total",100,-1,1);
	bMMassPmTot->GetXaxis()->SetTitle("M^{2} (GeV^{2})");
	bMMassPmTot->SetLineColor(3);


	TH1F *Pbins = new TH1F("Pbins","Pbins",4,0.4,1.5);
	Double_t ThetaBinstab [3] = {37.,45.,65.};
	TH1F *Thetabins = new TH1F("Thetabins","Thetabins",2,ThetaBinstab);
	TH1F styleEff("styleEff","styleEff",18,-180,180);
	TH1F *Eff[4][2];
	for(int Pbin=0;Pbin<4;Pbin++){
		for(int thbin=0;thbin<2;thbin++){
			Eff[Pbin][thbin]= new TH1F(styleEff);
			Eff[Pbin][thbin]->SetName(Form("%d%d",Pbin,thbin));
			Eff[Pbin][thbin]->SetTitle(Form("P bin %d / th bin %d ; #phi ",Pbin,thbin));
			Eff[Pbin][thbin]->Sumw2();
			
		}
	}
	
	TH1F *mmEff[4][2];
	for(int Pbin=0;Pbin<4;Pbin++){
		for(int thbin=0;thbin<2;thbin++){
			mmEff[Pbin][thbin]= new TH1F(styleEff);
			mmEff[Pbin][thbin]->SetName(Form("mm%d%d",Pbin,thbin));
			mmEff[Pbin][thbin]->SetTitle(Form("MM: P bin %d / th bin %d ; #phi ",Pbin,thbin));
			mmEff[Pbin][thbin]->Sumw2();
			
		}
	}

	TH1F *MomP = new TH1F("MomP","Momentum Proton",30,0,2);
	MomP->GetXaxis()->SetTitle(" P (GeV)");
	MomP->Sumw2();
	TH1F *PhiP = new TH1F("PhiP","Phi Proton",30,-180,180);
	PhiP->GetXaxis()->SetTitle(" Phi ");
	PhiP->Sumw2();
	TH1F *ThetaP = new TH1F("ThetaP","Theta Proton",40,30,70);
	ThetaP->GetXaxis()->SetTitle(" Theta ");
	ThetaP->Sumw2();
	
	
	
	TH2F *MomPhiP = new TH2F("MomPhiP",";#phi;P",50,-180,180,50,0,3);
	TH2F *ThetaPhiP = new TH2F("ThetaPhiP",";#phi;Theta",20,-180,180,50,30,70);
	TH2F *ThetaPhiP1 = new TH2F("ThetaPhiP1",";#phi;Theta",20,-180,180,50,30,70);
	TH2F *ThetaPhiP2 = new TH2F("ThetaPhiP2",";#phi;Theta",20,-180,180,50,30,70);
	TH2F *ThetaPP = new TH2F("ThetaPP",";#P;Theta",20,0,3,20,30,70);

	TH1F *mmMomP = new TH1F("mmMomP","Momentum mProton",30,0,2);
	mmMomP->GetXaxis()->SetTitle(" P (GeV)");
	mmMomP->Sumw2();
	TH1F *mmPhiP = new TH1F("mmPhiP","Phi mProton",30,-180,180);
	mmPhiP->GetXaxis()->SetTitle(" Phi ");
	mmPhiP->Sumw2();
	TH1F *mmThetaP = new TH1F("mmThetaP","Theta mProton",40,30,70);
	mmThetaP->GetXaxis()->SetTitle(" Theta ");
	mmThetaP->Sumw2();
	
	
	TH2F *mmMomPhiP = new TH2F("mmMomPhiP",";#phi;P",50,-180,180,50,0,3);
	TH2F *mmThetaPhiP = new TH2F("mmThetaPhiP",";#phi;Theta",20,-180,180,50,30,70);
	TH2F *mmThetaPhiP1 = new TH2F("mmThetaPhiP1",";#phi;Theta",20,-180,180,50,30,70);
	TH2F *mmThetaPhiP2 = new TH2F("mmThetaPhiP2",";#phi;Theta",20,-180,180,50,30,70);
	TH2F *mmThetaPP = new TH2F("ThetaPP",";#P;Theta",20,0,3,20,30,70);

	TH1F *ResoP = new TH1F("ResoP","Resolution P",100,-1,1);
	ResoP->GetXaxis()->SetTitle("Reso P (GeV)");
	TH1F *ResoTheta = new TH1F("ResoTheta","Resolution theta",100,-30,30);
	ResoTheta->GetXaxis()->SetTitle("Reso theta ");
	TH1F *ResoPhi = new TH1F("ResoPhi","Resolution Phi",100,-20,20);
	ResoPhi->GetXaxis()->SetTitle("Reso Phi");

	TH2F *ResoThetaVSThetaP = new TH2F("ResoThetaVSThetaP","Resolution theta VS theta Proton",80,40,60,80,-20,30);
	ResoThetaVSThetaP->GetXaxis()->SetTitle("Theta");
	TH2F *ResoThetaVSPhiP = new TH2F("ResoThetaVSPhiP","Resolution theta VS phi Proton",80,-180,180,80,-20,30);
	ResoThetaVSPhiP->GetXaxis()->SetTitle("Phi");
	TH2F *ResoThetaVSPP = new TH2F("ResoThetaVSPP","Resolution theta VS p Proton",80,0,3,80,-20,30);
	ResoThetaVSPP->GetXaxis()->SetTitle("P");

	TH2F *ResoThetaVSThetaE = new TH2F("ResoThetaVSThetaE","Resolution theta VS theta E",80,0,15,80,-20,30);
	ResoThetaVSThetaE->GetXaxis()->SetTitle("Theta");
	TH2F *ResoThetaVSPhiE = new TH2F("ResoThetaVSPhiE","Resolution theta VS phi E",80,-180,180,80,-20,30);
	ResoThetaVSPhiE->GetXaxis()->SetTitle("Phi");
	TH2F *ResoThetaVSPE = new TH2F("ResoThetaVSPE","Resolution theta VS p E",80,0,9,80,-20,30);
	ResoThetaVSPE->GetXaxis()->SetTitle("P");

	TH2F *ResoThetaVSThetaPip = new TH2F("ResoThetaVSThetaPip","Resolution theta VS theta pip",80,0,50,80,-20,30);
	ResoThetaVSThetaPip->GetXaxis()->SetTitle("Theta");
	TH2F *ResoThetaVSPhiPip = new TH2F("ResoThetaVSPhiPip","Resolution theta VS phi pip",80,-180,180,80,-20,30);
	ResoThetaVSPhiPip->GetXaxis()->SetTitle("Phi");
	TH2F *ResoThetaVSPPip = new TH2F("ResoThetaVSPPip","Resolution theta VS p pip",80,0,9,80,-20,30);
	ResoThetaVSPPip->GetXaxis()->SetTitle("P");

	TH2F *ResoThetaVSThetaPim = new TH2F("ResoThetaVSThetaPim","Resolution theta VS theta pim",80,0,50,80,-20,30);
	ResoThetaVSThetaPim->GetXaxis()->SetTitle("Theta");
	TH2F *ResoThetaVSPhiPim = new TH2F("ResoThetaVSPhiPim","Resolution theta VS phi pim",80,-180,180,80,-20,30);
	ResoThetaVSPhiPim->GetXaxis()->SetTitle("Phi");
	TH2F *ResoThetaVSPPim = new TH2F("ResoThetaVSPPim","Resolution theta VS p pim",80,0,9,80,-20,30);
	ResoThetaVSPPim->GetXaxis()->SetTitle("P");

	//Check Fiducial cuts
	TH2F *UPionP = new TH2F("UPionP","UPionP",100,0,420,100,0,10);
	TH2F *VPionP = new TH2F("VPionP","VPionP",100,0,420,100,0,10);
	TH2F *WPionP = new TH2F("WPionP","WPionP",100,0,420,100,0,10);
	TH2F *UPionM = new TH2F("UPionM","UPionM",100,0,420,100,0,10);
	TH2F *VPionM = new TH2F("VPionM","VPionM",100,0,420,100,0,10);
	TH2F *WPionM = new TH2F("WPionM","WPionM",100,0,420,100,0,10);
	TH2F *UElec = new TH2F("UElec","UElec",100,0,420,100,0,10);
	TH2F *VElec = new TH2F("VElec","VElec",100,0,420,100,0,10);
	TH2F *WElec = new TH2F("WElec","WElec",100,0,420,100,0,10);

	TH2F *ThetaPproton = new TH2F("ThetaPproton","ThetaPproton",100,30,70,100,0,3);


	//Momentum correction
	TH2F *PhiDeltaPproton = new TH2F("PhiDeltaPproton","PhiDeltaPproton",50,-180,180,50,-0.5,0.5);
	TH2F *ThetaDeltaPproton = new TH2F("ThetaDeltaPproton","ThetaDeltaPproton",50,40,60,50,-0.5,0.5);


	TH2F *ResoPVSThetaE = new TH2F("ResoPVSThetaE","Resolution P VS theta E",80,0,15,50,-0.5,0.5);
	TH2F *ResoPVSPhiE = new TH2F("ResoPVSPhiE","Resolution P VS phi E",80,-180,180,50,-0.5,0.5);
	TH2F *ResoPVSPE = new TH2F("ResoPVSPE","Resolution P VS p E",80,0,9,50,-0.5,0.5);
	TH2F *ResoPVSThetaPip = new TH2F("ResoPVSThetaPip","Resolution P VS thetapip",80,0,50,50,-0.5,0.5);
	TH2F *ResoPVSPhiPip = new TH2F("ResoPVSPhiPip","Resolution P VS phi pip",80,-180,180,50,-0.5,0.5);
	TH2F *ResoPVSPPip = new TH2F("ResoPVSPPip","Resolution P VS p pip",80,0,9,50,-0.5,0.5);
	TH2F *ResoPVSThetaPim = new TH2F("ResoPVSThetaPim","Resolution P VS thetapim",80,0,50,50,-0.5,0.5);
	TH2F *ResoPVSPhiPim = new TH2F("ResoPVSPhiPim","Resolution P VS phi pim",80,-180,180,50,-0.5,0.5);
	TH2F *ResoPVSPPim = new TH2F("ResoPVSPPim","Resolution P VS p pim",80,0,9,50,-0.5,0.5);
	
	
	TH2F *ResoPVSTheta = new TH2F("ResoPVSTheta","; #theta;#Delta P /P ",100,35,60,100,-0.5,0.5);
	TH2F *ResoPVSPhi = new TH2F("ResoPVSPhi","; #phi_{CVT}(#circ) + 90#circ ;#Delta P /P ",50,-180,180,100,-0.5,0.5);
	TH2F *ResoPVSPhi1 = new TH2F("ResoPVSPhi1","; #phi_{CVT}(#circ);#Delta P /P ",50,-180,180,100,-0.5,0.5);
	TH2F *ResoPVSPhi2 = new TH2F("ResoPVSPhi2","; #phi_{CVT}(#circ);#Delta P /P ",50,-180,180,100,-0.5,0.5);
	TH2F *ResoPVSP = new TH2F("ResoPVSP","; P;#Delta P /P ",100,0,1.5,100,-0.5,0.5);
	TH2F *PVSPhi = new TH2F("PVSPhi","PVSPhi ",100,-180,180,100,0.2,1.5);
	


	//check proton momentum dependance
	TH2F *PVSThetaE = new TH2F("PVSThetaE"," P VS theta E",80,0,15,50,0,3);
	TH2F *PVSPhiE = new TH2F("PVSPhiE"," P VS phi E",80,-180,180,50,0,3);
	TH2F *PVSPE = new TH2F("PVSPE"," P VS p E",80,0,9,50,0,3);
	TH2F *PVSThetaPip = new TH2F("PVSThetaPip"," P VS thetapip",80,0,50,50,0,3);
	TH2F *PVSPhiPip = new TH2F("PVSPhiPip"," P VS phi pip",80,-180,180,50,0,3);
	TH2F *PVSPPip = new TH2F("PVSPPip"," P VS p pip",80,0,9,50,0,3);
	TH2F *PVSThetaPim = new TH2F("PVSThetaPim"," P VS thetapim",80,0,50,50,0,3);
	TH2F *PVSPhiPim = new TH2F("PVSPhiPim"," P VS phi pim",80,-180,180,50,0,3);
	TH2F *PVSPPim = new TH2F("PVSPPim"," P VS p pim",80,0,9,50,0,3);



	//fidtrack
	TH1F* thetatrack=new TH1F("thetatrack","thetatrack",100,30,60);
	TH1F* Chi2thetatrack=new TH1F("Chi2thetatrack","Chi2thetatrack",100,30,60);
	
	TH1F* phitrack=new TH1F("phitrack","phitrack",100,-180,180);
	TH1F* Chi2phitrack=new TH1F("Chi2phitrack","Chi2phitrack",100,-180,180);
	
	TH1F* Resophitrack=new TH1F("Resophitrack","Resophitrack",100,-180,180);


	TH1F *t = new TH1F("t","t",100,0,2);
	TH1F *essai = new TH1F("essai","essai",100,0,40);
	TH2F *tmass = new TH2F("tmass","tmass",100,0,2,100,0,2);

	//Check electron kinematics
	TH2F *PThetaE = new TH2F("PThetaE"," P VS theta E",80,0,10,80,0,30);
	
	
	TF1 *f1P = new TF1("f1P","[0]+[1]*x+[2]*x*x",0,1);
	f1P->SetParameters(0.153319,-0.298968,0.1607);
	TF1 *f2P = new TF1("f1P","[0]+[1]*x+[2]*x*x",0,1);
	f2P->SetParameters(0.0398946,-0.0748125,0.0395764);
	TF1 *f3P = new TF1("f1P","[0]+[1]*x",0,1);
	f3P->SetParameters(0.0292947,-0.0577956);
	
	
	/*
	TF1 *f1Pproton = new TF1("f1Pproton","[0]+[1]*x",-180.,180.);
	f1Pproton->SetParameters(0.0111159,-6.83753e-05);
	TF1 *f2Pproton = new TF1("f2Pproton","[0]+[1]*x",-180.,180.);
	f2Pproton->SetParameters(-0.047977,0.00084448);
	TF1 *f3Pproton = new TF1("f3Pproton","[0]+[1]*x",-180.,180.);
	f3Pproton->SetParameters(-0.0543518,0.00045238);
*/


	TF1 *f1Pproton = new TF1("f1Pproton","[0]+[1]*x",-130.,10.);
	f1Pproton->SetParameters(-0.00146068,-2.13735e-05);
	TF1 *f2Pproton = new TF1("f2Pproton","[0]+[1]*x",-100.,40.);
	f2Pproton->SetParameters(-0.0608671,0.000849025);
	TF1 *f3Pproton = new TF1("f3Pproton","[0]+[1]*x",25.,160.);
	f3Pproton->SetParameters(-0.0670748,0.000419003);
	/*
	TF1 *f1Pproton = new TF1("f1Pproton","[0]+[1]*x",-130.,10.);
	f1Pproton->SetParameters(-0.0314093,-5.85409e-05);
	TF1 *f2Pproton = new TF1("f2Pproton","[0]+[1]*x",-100.,40.);
	f2Pproton->SetParameters(-0.0924035,0.000822015);
	TF1 *f3Pproton = new TF1("f3Pproton","[0]+[1]*x",25.,160.);
	f3Pproton->SetParameters(-0.0973604,0.00046165);
	*/
	
	TString nameFiles="/vol0/pierre/Bureau/DATA/FilterppessaiApril2020new.hipo";//BGeffCentral_new_3.hipo";//centralP.hipo";//////simuRho.hipo";//SimuRhoOk.hipo";////
	hipo::reader  reader;
	reader.open(nameFiles);

	hipo::dictionary  factory;
	reader.readDictionary(factory);
	factory.show();
	hipo::event      event;

	hipo::bank EVENT(factory.getSchema("REC::Event"));
	hipo::bank PART(factory.getSchema("REC::Particle"));
	hipo::bank TRACK(factory.getSchema("REC::Track"));
	hipo::bank TRAJ(factory.getSchema("REC::Traj"));
	/*hipo::bank SCIN(factory.getSchema("REC::Scintillator"));
	hipo::bank CHE(factory.getSchema("REC::Cherenkov"));
	hipo::bank CALO(factory.getSchema("REC::Calorimeter"));*/

	hipo::bank RUN(factory.getSchema("RUN::config"));
	int run=0;

	int n=0;
	int n1=0;
	while (reader.next() /*&& n<100000*/){

		reader.read(event);
		event.getStructure(PART);
		/*event.getStructure(SCIN);
		event.getStructure(CHE);
		event.getStructure(CALO);*/
		event.getStructure(EVENT);
		event.getStructure(TRAJ);
		event.getStructure(TRACK);
		//event.getStructure(RUN);

		int IDrun=0;//RUN.getInt("run",0);
		
		
		
		if(IDrun!=run){cout<<IDrun<<endl;
		run=IDrun;}

		n++;
		int pid;
		float massN=0.938;
		float massPI=0.139;
		float me =0.000511;
		float ebeam=10.6;

		Particle PionP;
		Particle PionM;
		Particle Electron;
		Particle Proton;
		TLorentzVector RestProton;
		RestProton.SetPxPyPzE(0,0,0,massN);
		TLorentzVector Beam;
		TLorentzVector MissingP;
		TLorentzVector MissingTot;
		TLorentzVector MissingPTot;
		TLorentzVector MissingPpTot;
		TLorentzVector MissingPmTot;
		Beam.SetPxPyPzE(0,0,ebeam,ebeam);

		int recE=0;
		int recPim=0;
		int recPip=0;
		int recP=0;
		int np=PART.getRows();
		Particle Photons[np];

		Track Region1Track;

		for( int i=0; i < np; i++ ){
			//cout<<"particle"<<endl;
			int pid=PART.getInt("pid",i);
			float  px = PART.getFloat("px",i);
			float  py = PART.getFloat("py",i);
			float  pz = PART.getFloat("pz",i);
			float  beta = PART.getFloat("beta",i);
			int status = abs(PART.getInt("status",i));
			int charge = PART.getInt("charge",i);
			float  chi2 = PART.getFloat("chi2pid",i);
			float  vx = PART.getFloat("vx",i);
			float  vy = PART.getFloat("vy",i);
			float  vz = PART.getFloat("vz",i);

			if(pid==11 && status>2000){
				Electron.Vector.SetXYZM(px,py,pz,me);
				Electron.index=i;
				Electron.pid=11;
				Electron.status=status;
				Electron.chi2=chi2;
				Electron.vertex.z=vz;
				recE++;
			}

			if(pid==211 && status<4000){
				PionP.Vector.SetXYZM(px,py,pz,massPI);
				PionP.index=i;
				PionP.pid=211;
				PionP.status=status;
				PionP.chi2=chi2;
				recPip++;
			}

			if(pid==-211 && status<4000){
				PionM.Vector.SetXYZM(px,py,pz,massPI);
				PionM.index=i;
				PionM.pid=-211;
				PionM.status=status;
				PionM.chi2=chi2;
				recPim++;
			}

			if(pid==2212){
				Proton.Vector.SetXYZM(px,py,pz,massN);
				Proton.index=i;
				Proton.pid=2212;
				Proton.status=status;
				Proton.chi2=chi2;
				recP++;
			}
			
			if(pid==22 ){
				Photons[i].Vector.SetXYZM(px,py,pz,0.);
				Photons[i].index=i;
				Photons[i].pid=22;
				Photons[i].status=status;
				Photons[i].beta=beta;
			}


		}

		if(recE==1)essai->Fill(Electron.Vector.Theta()*TMath::RadToDeg());

		vector<Particle> Particles={Electron,Proton,PionM,PionP};

		bool PIDflag=false;
		bool FiducialCut=true;



		if(recP==1){
		
		for(int i=0;i<TRACK.getRows();i++){
		
			


			int sector = 0;
			float chi2 = -1000.0;
			int NDF = 1;
			int indexTrack=-5000;
		
			if(TRACK.getInt("status",i)>0 && TRACK.getInt("detector",i)==5 ){
				sector=TRACK.getInt("sector",i);
				chi2=TRACK.getFloat("chi2",i);
				NDF=TRACK.getInt("NDF",i);
				indexTrack=TRACK.getInt("index",i);
			}




			float Nchi2=chi2/NDF;


			for(int c=0;c<TRAJ.getRows();c++){

				int detector=0;
				int pindex = TRAJ.getInt("pindex",c);
				int index = TRAJ.getInt("index",c);
				int layer = TRAJ.getInt("layer",c);
				detector = TRAJ.getInt("detector",c);
				float x = TRAJ.getFloat("x",c);
				float y = TRAJ.getFloat("y",c);
				float z = TRAJ.getFloat("z",c);

				//Region1
				if(index==(indexTrack) && detector==5 && layer==12 && pindex==Proton.index ){
				//cout<<layer<<endl;
					Region1Track.detector=detector;
					Region1Track.pindex=index;
					Region1Track.sector=sector;
					Region1Track.layer=layer;
					Region1Track.x=x;
					Region1Track.y=y;
					Region1Track.z=z;
					Region1Track.chi2=Nchi2;

				}


			}
			}
		
		
		}

		double x=Region1Track.x;
			double y=Region1Track.y;
			double z=Region1Track.z;
//cout<<x<<" "<<y<<" "<<z<<endl;
			double Nchi2=Region1Track.chi2;

			float d = sqrt(x*x+y*y);
			//cout<<d<<endl;
			float phi = (TMath::ATan2(y,x)*TMath::RadToDeg());
				//if(sector>3 || sector<7)phi = (TMath::ATan2(y,x)*TMath::RadToDeg())+(sector*(-60))+420;
			float theta = TMath::ATan(d/z)*TMath::RadToDeg();
		//if((phi<160. && phi>140.) || (phi<40. && phi>20.) || (phi<-80. && phi>-100.))continue;
		
					if(recP==1 && Proton.status>4000){
					phitrack->Fill(phi);Chi2phitrack->Fill(phi,Region1Track.chi2);
					thetatrack->Fill(theta);Chi2thetatrack->Fill(theta,Region1Track.chi2);
					}

			
/*
		if(PIDflag){
			CalorimeterResp Calo;
			CheResp Che;
			ScinResp Scin;
			for (int i=0;i<4 ; i++){

				for(int c=0;c<SCIN.getRows();c++){
					int Calopindex = CALO.getInt("pindex",c);                        
					int Calosector = CALO.getInt("sector",c);                        
					int Calolayer = CALO.getInt("layer",c);                          
					int Calodetector = CALO.getInt("detector",c);                    
					float Caloenergy = CALO.getFloat("energy",c);                    
					float Calox = CALO.getFloat("x",c);                              
					float Caloy = CALO.getFloat("y",c);                              
					float Caloz = CALO.getFloat("z",c);                              
					float Calou = CALO.getFloat("lu",c);                             
					float Calov = CALO.getFloat("lv",c);                             
					float Calow = CALO.getFloat("lw",c);                             
					float Calodu = CALO.getFloat("du",c);                            
					float Calodv = CALO.getFloat("dv",c);                            
					float Calodw = CALO.getFloat("dw",c);                            
					float Calom2u = CALO.getFloat("m2u",c);                          
					float Calom2v = CALO.getFloat("m2v",c);                          
					float Calom2w = CALO.getFloat("m2w",c);                          
					float Calom3u = CALO.getFloat("m3u",c);                          
					float Calom3v = CALO.getFloat("m3v",c);                          
					float Calom3w = CALO.getFloat("m3w",c);
					if(Calopindex==(Particles[i].index)){
						Calo.detector=Calodetector;
						Calo.pindex=Calopindex;
						Calo.sector=Calosector;
						Calo.layer=Calolayer;
						Calo.energy=Caloenergy;
						Calo.x=Calox;
						Calo.y=Caloy;
						Calo.z=Caloz;
						Calo.u=Calou;
						Calo.v=Calov;
						Calo.w=Calow;
						Calo.du=Calodu;
						Calo.dv=Calodv;
						Calo.dw=Calodw;
						Calo.m2u=Calom2u;
						Calo.m2v=Calom2v;
						Calo.m2w=Calom2w;
						Calo.m3u=Calom3u;
						Calo.m3v=Calom3v;
						Calo.m3w=Calom3w;
						Particles[i].Calorimeter.push_back(Calo);
					}
				}

			}
			Electron=Particles[0];
			Proton=Particles[1];
			PionM=Particles[2];
			PionP=Particles[3];

		}	
*/

		double avant=Electron.Vector.P();
		bool RadCorr = true;//true;//
		if(RadCorr){
		//cout<<"here "<<endl;
		Electron = RadiativeCorr(Electron,Photons,10.,1.5,np);
		}
		
		if(avant!=Electron.Vector.P()){
		//cout<<"avant "<<avant<<"après "<<Electron.Vector.P()<<" "<<(avant-Electron.Vector.P())/avant<<" status "<<Electron.status<<endl;
		}
		
		
		/*if(recP==1){
		cout<<"event"<<endl;
		cout<<Proton.Vector.P()<<endl;
		}*/
		//MC correction
		double momprotonavant= Proton.Vector.P();
		if(recP==1){
		double PP=Proton.Vector.P();
		double newPP = 0.0;
		if(Proton.status<4000 && Proton.Vector.Theta()*TMath::RadToDeg()>27.)newPP=PP+(f1P->Eval(PP));
		if(Proton.status<4000 && Proton.Vector.Theta()*TMath::RadToDeg()<27.)newPP=PP+(f2P->Eval(PP));
		if(Proton.status>4000)newPP=PP+(f3P->Eval(PP));
		
		
		//cout<<"avant "<<momprotonavant<<endl;
		Proton.Vector.SetRho(newPP);
		}
		//if(recP==1)cout<<(Proton.Vector.P()-momprotonavant)/momprotonavant<<endl;

		double W = (RestProton+Beam-Electron.Vector).M(); 
		double Q2=(Beam-Electron.Vector).M2();
		if(recPip==1 && recPim==1 && recE==1 && Electron.status>2000 /*&& (PionP.Vector+PionM.Vector).M()>0.6 && (PionP.Vector+PionM.Vector).M()<1.*/){
		
		
			
		
		
			n1++;
			MissingP=Beam+RestProton-Electron.Vector-PionP.Vector-PionM.Vector;
			PThetaE->Fill(Electron.Vector.P(),Electron.Vector.Theta()*TMath::RadToDeg());
			//if(MissingP.Phi()*TMath::RadToDeg()<25. && MissingP.Phi()*TMath::RadToDeg()>0.)continue;
			//MMass->Fill(MissingP.M2());
			MMass->Fill(Electron.vertex.z);

			if(MissingP.Theta()*TMath::RadToDeg()>37.0 && MissingP.Theta()*TMath::RadToDeg()<65.0){
				MMassCentralRaw->Fill(MissingP.M());
				if( PionP.Vector.P()>1.0 && PionM.Vector.P()>1.0 && PionP.chi2<4. && PionM.chi2<4. && PionP.chi2>-4. && PionM.chi2>-4. && PionP.status<3000 && PionM.status<3000)MMassW->Fill(MissingP.M2(),W);
				if( PionP.Vector.P()>1.0 && PionM.Vector.P()>1.0 && PionP.chi2<4. && PionM.chi2<4. && PionP.chi2>-4. && PionM.chi2>-4. && PionP.status<3000 && PionM.status<3000)MMassThetaE->Fill(MissingP.M2(),Electron.Vector.Theta()*TMath::RadToDeg());
				if( PionP.Vector.P()>1.0 && PionM.Vector.P()>1.0 && PionP.chi2<4. && PionM.chi2<4. && PionP.chi2>-4. && PionM.chi2>-4. && PionM.status<3000)MMassThetaPip->Fill(MissingP.M2(),PionP.Vector.Theta()*TMath::RadToDeg());
				if( PionP.Vector.P()>1.0 && PionM.Vector.P()>1.0 && PionP.chi2<4. && PionM.chi2<4. && PionP.chi2>-4. && PionM.chi2>-4. && PionP.status<3000)MMassThetaPim->Fill(MissingP.M2(),PionM.Vector.Theta()*TMath::RadToDeg());
				if( PionP.Vector.P()>1.0 && PionM.Vector.P()>1.0 && PionP.chi2<4. && PionM.chi2<4. && PionP.chi2>-4. && PionM.chi2>-4. && PionP.status<3000 && PionM.status<3000)MMassPE->Fill(MissingP.M2(),Electron.Vector.P());
				if(PionM.Vector.P()>1.0 && PionP.chi2<4. && PionM.chi2<4. && PionP.chi2>-4. && PionM.chi2>-4. && PionP.status<3000 && PionM.status<3000)MMassPPip->Fill(MissingP.M2(),PionP.Vector.P());
				if( PionP.Vector.P()>1.0 && PionP.chi2<4. && PionM.chi2<4. && PionP.chi2>-4. && PionM.chi2>-4. && PionP.status<3000 && PionM.status<3000)MMassPPim->Fill(MissingP.M2(),PionM.Vector.P());
				if( PionP.Vector.P()>1.0 && PionM.Vector.P()>1.0 && PionP.chi2<4. && PionM.chi2<4. && PionP.chi2>-4. && PionM.chi2>-4. && PionP.status<3000 && PionM.status<3000)MMassChi2E->Fill(MissingP.M2(),Electron.chi2);
				if( PionP.Vector.P()>1.0 && PionM.Vector.P()>1.0 && PionM.chi2<4. && PionM.chi2>-4. && PionP.status<3000 && PionM.status<3000)MMassChi2Pp->Fill(MissingP.M2(),PionP.chi2);
				if( PionP.Vector.P()>1.0 && PionM.Vector.P()>1.0 && PionP.chi2<4. && PionP.chi2>-4. && PionP.status<3000 && PionM.status<3000)MMassChi2Pm->Fill(MissingP.M2(),PionM.chi2);

			}

			if(/*MissingP.M2()>0.6 && MissingP.M2()<1.2 &&*/
			
			 MissingP.Theta()*TMath::RadToDeg()>37.0 && MissingP.Theta()*TMath::RadToDeg()<65. //65.0 
			&& MissingP.P()>0.3 && MissingP.P()<1.5
			
			&& PionP.Vector.P()>1.0 && PionM.Vector.P()>1.0 && PionP.chi2<4. && PionM.chi2<4. && PionP.chi2>-4. && PionM.chi2>-4. && PionP.status<3000 && PionM.status<3000){
				MMassTheta->Fill(MissingP.M2(),MissingP.Theta()*TMath::RadToDeg());
				plotQ2W->Fill(-Q2,W);
				ThetaPproton->Fill(MissingP.Theta()*TMath::RadToDeg(),MissingP.P());
				MMassCentral->Fill(MissingP.M());
				mmMomP->Fill(MissingP.P(),1.);
				mmMomPhiP->Fill(MissingP.Phi()*TMath::RadToDeg(),MissingP.P(),1.);
				
				
				if(MissingP.P()<0.66)mmThetaPhiP->Fill(MissingP.Phi()*TMath::RadToDeg(),MissingP.Theta()*TMath::RadToDeg(),1.);
				if(MissingP.P()>0.66 && MissingP.P()<1.)mmThetaPhiP1->Fill(MissingP.Phi()*TMath::RadToDeg(),MissingP.Theta()*TMath::RadToDeg(),1.);
				if(MissingP.P()>1.)mmThetaPhiP2->Fill(MissingP.Phi()*TMath::RadToDeg(),MissingP.Theta()*TMath::RadToDeg(),1.);
				
				mmThetaPP->Fill(MissingP.P(),MissingP.Theta()*TMath::RadToDeg(),1.);
				mmPhiP->Fill(MissingP.Phi()*TMath::RadToDeg(),1.);
				mmThetaP->Fill(MissingP.Theta()*TMath::RadToDeg(),1.);
				MMassMmP->Fill(MissingP.M2(),MissingP.P());

				int Pbinning = Pbins->GetXaxis()->FindBin(MissingP.P())-1;
				int Thetabinning = Thetabins->GetXaxis()->FindBin(MissingP.Theta()*TMath::RadToDeg())-1;
				if(Pbinning>-1 && Pbinning<4 && Thetabinning>-1 && Thetabinning<2)mmEff[Pbinning][Thetabinning]->Fill(MissingP.Phi()*TMath::RadToDeg());


				PVSThetaE->Fill(Electron.Vector.Theta()*TMath::RadToDeg(),MissingP.P());
				PVSPhiE->Fill(Electron.Vector.Phi()*TMath::RadToDeg(),MissingP.P());
				PVSPE->Fill(Electron.Vector.P(),MissingP.P());

				PVSThetaPip->Fill(PionP.Vector.Theta()*TMath::RadToDeg(),MissingP.P());
				PVSPhiPip->Fill(PionP.Vector.Phi()*TMath::RadToDeg(),MissingP.P());
				PVSPPip->Fill(PionP.Vector.P(),MissingP.P());

				PVSThetaPim->Fill(PionM.Vector.Theta()*TMath::RadToDeg(),MissingP.P());
				PVSPhiPim->Fill(PionM.Vector.Phi()*TMath::RadToDeg(),MissingP.P());
				PVSPPim->Fill(PionM.Vector.P(),MissingP.P());

				t->Fill(-(MissingP-RestProton).M2());
				tmass->Fill(-(MissingP-RestProton).M2(),(PionM.Vector+PionP.Vector).M());


				if(/*recPip==1 && recPim==1 && recE==1 && */recP==1 /*&& Electron.status>2000*/){
//if((recP==1 && Proton.status>4000) && ((phi<170. && phi>130.) || (phi<50. && phi>10.) || (phi<-70. && phi>-110.)))continue;
					MissingTot=Beam+RestProton-Electron.Vector-PionP.Vector-PionM.Vector-Proton.Vector;
					MissingPTot=Beam+RestProton-Electron.Vector-PionP.Vector-PionM.Vector;
					MissingPpTot=Beam+RestProton-Electron.Vector-PionM.Vector-Proton.Vector;
					MissingPmTot=Beam+RestProton-Electron.Vector-PionP.Vector-Proton.Vector;

//if(Proton.status<4000)cout<<Proton.status<<endl;

					double MMtot = MissingTot.M2();
					double MMPtot = MissingPTot.M2();
					double MMPptot = MissingPpTot.M2();
					double MMPmtot = MissingPmTot.M2();

					if(Proton.Vector.Theta()*TMath::RadToDeg()>37.0 && Proton.Vector.Theta()*TMath::RadToDeg()<65.0){
						MMassTot->Fill(MMtot);
						MMassPTot->Fill(MMPtot);
						MMassPpTot->Fill(MMPptot);
						MMassPmTot->Fill(MMPmtot);
					}

					if(PionP.Vector.P()>1.0 && PionM.Vector.P()>1.0 && PionP.chi2<4. && PionM.chi2<4. && PionP.chi2>-4. && PionM.chi2>-4. && PionP.status<3000 && PionM.status<3000){

						if(Proton.status>4000 /*&& (phi>-90. && phi<30.)*/ && Proton.Vector.Theta()*TMath::RadToDeg()>37.0 && Proton.Vector.Theta()*TMath::RadToDeg()<65.0 /*&& MMPtot>0. && MMPtot<1.5 && MMPptot>-0.2 && MMPptot<0.2 && MMPmtot>-0.2 && MMPmtot<0.2 */)bMMassTot->Fill(MMtot);
						
						if(Proton.status>4000 && (phi>-90. && phi<30.) && Proton.Vector.Theta()*TMath::RadToDeg()>37.0 && Proton.Vector.Theta()*TMath::RadToDeg()<65.0 )bPitheta->Fill(MissingPmTot.Theta()*TMath::RadToDeg()-PionM.Vector.Theta()*TMath::RadToDeg());
						
						if(Proton.Vector.Theta()*TMath::RadToDeg()>37.0 && Proton.Vector.Theta()*TMath::RadToDeg()<65.0 /*&& MMtot<0.2 && MMtot>-0.2 && MMPptot>-0.2 && MMPptot<0.2 && MMPmtot>-0.2 && MMPmtot<0.2*/ )bMMassPTot->Fill(MMPtot);
						if(MissingPTot.Theta()*TMath::RadToDeg()>37.0 && Proton.Vector.Theta()*TMath::RadToDeg()<65.0/*&& MMtot<0.2 && MMtot>-0.2 && MMPtot>0. && MMPtot<1.5 && MMPmtot>-0.2 && MMPmtot<0.2 */)bMMassPpTot->Fill(MMPptot);
						if(MissingPTot.Theta()*TMath::RadToDeg()>37.0 && Proton.Vector.Theta()*TMath::RadToDeg()<65.0/*&& MMtot<0.2 && MMtot>-0.2 && MMPtot>0. && MMPtot<1.5 && MMPptot>-0.2 && MMPptot<0.2 */)bMMassPmTot->Fill(MMPmtot);

						if( Proton.status>4000 && Proton.Vector.Theta()*TMath::RadToDeg()>40.0 && Proton.Vector.Theta()*TMath::RadToDeg()<65.0 /*&& MissingPTot.M2()>0.6 && MissingPTot.M2()<1.2*/ /* && MMtot<0.2 && MMtot>-0.2 && MMPtot>0. && MMPtot<1.5 && MMPptot>-0.2 && MMPptot<0.2 && MMPmtot>-0.2 && MMPmtot<0.2 */){



							
							double resoP=(Proton.Vector.P()-MissingPTot.P())/MissingPTot.P();
							double resoP1=(Proton.Vector.P()-MissingPTot.P())/Proton.Vector.P();
							double diffP=(Proton.Vector.P()-MissingPTot.P());
							
							Resophitrack->Fill(phi,resoP);
							//momentum correction
							PhiDeltaPproton->Fill(Proton.Vector.Phi()*TMath::RadToDeg(),resoP);
							ThetaDeltaPproton->Fill(Proton.Vector.Theta()*TMath::RadToDeg(),resoP);

							ResoPVSThetaE->Fill(Electron.Vector.Theta()*TMath::RadToDeg(),resoP);
							ResoPVSPhiE->Fill(Electron.Vector.Phi()*TMath::RadToDeg(),resoP);
							ResoPVSPE->Fill(Electron.Vector.P(),resoP);

							ResoPVSThetaPip->Fill(PionP.Vector.Theta()*TMath::RadToDeg(),resoP);
							ResoPVSPhiPip->Fill(PionP.Vector.Phi()*TMath::RadToDeg(),resoP);
							ResoPVSPPip->Fill(PionP.Vector.P(),resoP);

							ResoPVSThetaPim->Fill(PionM.Vector.Theta()*TMath::RadToDeg(),resoP);
							ResoPVSPhiPim->Fill(PionM.Vector.Phi()*TMath::RadToDeg(),resoP);
							ResoPVSPPim->Fill(PionM.Vector.P(),resoP);

							
							bool correctMom=false;
							
							double PP=Proton.Vector.P();
							double newPP = 0.0;
							if(phi>150. || phi<-90.){
							double phi1=phi;
										if(phi>0.){phi1=phi-270.;}
										else {phi1=phi+90.;}
										newPP=PP*(1.-(f1Pproton->Eval(phi1)));
										}
							if(phi>-90. && phi<30.){newPP=PP*(1.-(f2Pproton->Eval(phi)));}
							if(phi<150. && phi>30.){newPP=PP*(1.-(f3Pproton->Eval(phi)));}
		
							if(correctMom)Proton.Vector.SetRho(newPP);
							
							MissingTot=Beam+RestProton-Electron.Vector-PionP.Vector-PionM.Vector-Proton.Vector;
							MMtot = MissingTot.M2();
							if(Proton.Vector.Theta()*TMath::RadToDeg()>37.0 && Proton.Vector.Theta()*TMath::RadToDeg()<65.0 && (phi>-90. && phi<30.))bMMassTot1->Fill(MMtot);
							
							MissingPmTot=Beam+RestProton-Electron.Vector-PionP.Vector-Proton.Vector;
							if(Proton.status>4000 && (phi>-90. && phi<30.) && Proton.Vector.Theta()*TMath::RadToDeg()>37.0 && Proton.Vector.Theta()*TMath::RadToDeg()<65.0 )bPitheta1->Fill(MissingPmTot.Theta()*TMath::RadToDeg()-PionM.Vector.Theta()*TMath::RadToDeg());
							
							resoP1=(Proton.Vector.P()-MissingPTot.P())/Proton.Vector.P();
							diffP=(Proton.Vector.P()-MissingPTot.P());
							
							if(phi>150. || phi<-90.){
							double phi1=phi;
										if(phi>0.){phi1=phi-270.;}
										else {phi1=phi+90.;}
										ResoPVSPhi->Fill(phi1,resoP1);
										}
							if(phi>-90. && phi<30.){ResoPVSPhi1->Fill(phi,resoP1);}
							if(phi<150. && phi>30.){ResoPVSPhi2->Fill(phi,resoP1);}
							//cout<<Region1Track.sector<<endl;
							
							ResoPVSTheta->Fill(Proton.Vector.Theta()*TMath::RadToDeg(),diffP);
							ResoPVSP->Fill(Proton.Vector.P(),diffP);
						
							PVSPhi->Fill(Proton.Vector.Phi()*TMath::RadToDeg(),Proton.Vector.P());
							

							//end momentum correction

							ThetaPW->Fill(MissingP.Theta()*TMath::RadToDeg(),W);
							ThetaPQ2->Fill(MissingP.Theta()*TMath::RadToDeg(),-Q2);

							InvariantMassPP->Fill((PionP.Vector+PionM.Vector).M());

							/*MomP->Fill(Proton.Vector.P(),1.);
							MomPhiP->Fill(Proton.Vector.Phi()*TMath::RadToDeg(),Proton.Vector.P(),1.);
							ThetaPhiP->Fill(Proton.Vector.Phi()*TMath::RadToDeg(),Proton.Vector.Theta()*TMath::RadToDeg(),1.);
							ThetaPP->Fill(Proton.Vector.P(),Proton.Vector.Theta()*TMath::RadToDeg(),1.);
							PhiP->Fill(Proton.Vector.Phi()*TMath::RadToDeg(),1.);
							ThetaP->Fill(Proton.Vector.Theta()*TMath::RadToDeg(),1.);*/
							
							MomP->Fill(MissingP.P(),1.);
							MomPhiP->Fill(MissingP.Phi()*TMath::RadToDeg(),MissingP.P(),1.);
							
							if(MissingP.P()<0.66)ThetaPhiP->Fill(MissingP.Phi()*TMath::RadToDeg(),MissingP.Theta()*TMath::RadToDeg(),1.);
							if(MissingP.P()>0.66  && MissingP.P()<1.)ThetaPhiP1->Fill(MissingP.Phi()*TMath::RadToDeg(),MissingP.Theta()*TMath::RadToDeg(),1.);
							if(MissingP.P()>1.)ThetaPhiP2->Fill(MissingP.Phi()*TMath::RadToDeg(),MissingP.Theta()*TMath::RadToDeg(),1.);
							
							ThetaPP->Fill(MissingP.P(),MissingP.Theta()*TMath::RadToDeg(),1.);
							PhiP->Fill(MissingP.Phi()*TMath::RadToDeg(),1.);
							ThetaP->Fill(MissingP.Theta()*TMath::RadToDeg(),1.);
							
							if(Pbinning>-1 && Pbinning<4 && Thetabinning>-1 && Thetabinning<2)Eff[Pbinning][Thetabinning]->Fill(MissingP.Phi()*TMath::RadToDeg());
							
							ResoP->Fill(resoP);
							ResoTheta->Fill(Proton.Vector.Theta()*TMath::RadToDeg()-MissingPTot.Theta()*TMath::RadToDeg());
							ResoPhi->Fill(Proton.Vector.Phi()*TMath::RadToDeg()-MissingPTot.Phi()*TMath::RadToDeg());

							ResoThetaVSThetaP->Fill(Proton.Vector.Theta()*TMath::RadToDeg(),Proton.Vector.Theta()*TMath::RadToDeg()-MissingPTot.Theta()*TMath::RadToDeg());
							ResoThetaVSPhiP->Fill(Proton.Vector.Phi()*TMath::RadToDeg(),Proton.Vector.Theta()*TMath::RadToDeg()-MissingPTot.Theta()*TMath::RadToDeg());
							ResoThetaVSPP->Fill(Proton.Vector.P(),Proton.Vector.Theta()*TMath::RadToDeg()-MissingPTot.Theta()*TMath::RadToDeg());

							ResoThetaVSThetaE->Fill(Electron.Vector.Theta()*TMath::RadToDeg(),Proton.Vector.Theta()*TMath::RadToDeg()-MissingPTot.Theta()*TMath::RadToDeg());
							ResoThetaVSPhiE->Fill(Electron.Vector.Phi()*TMath::RadToDeg(),Proton.Vector.Theta()*TMath::RadToDeg()-MissingPTot.Theta()*TMath::RadToDeg());
							ResoThetaVSPE->Fill(Electron.Vector.P(),Proton.Vector.Theta()*TMath::RadToDeg()-MissingPTot.Theta()*TMath::RadToDeg());

							ResoThetaVSThetaPip->Fill(PionP.Vector.Theta()*TMath::RadToDeg(),Proton.Vector.Theta()*TMath::RadToDeg()-MissingPTot.Theta()*TMath::RadToDeg());
							ResoThetaVSPhiPip->Fill(PionP.Vector.Phi()*TMath::RadToDeg(),Proton.Vector.Theta()*TMath::RadToDeg()-MissingPTot.Theta()*TMath::RadToDeg());
							ResoThetaVSPPip->Fill(PionP.Vector.P(),Proton.Vector.Theta()*TMath::RadToDeg()-MissingPTot.Theta()*TMath::RadToDeg());

							ResoThetaVSThetaPim->Fill(PionM.Vector.Theta()*TMath::RadToDeg(),Proton.Vector.Theta()*TMath::RadToDeg()-MissingPTot.Theta()*TMath::RadToDeg());
							ResoThetaVSPhiPim->Fill(PionM.Vector.Phi()*TMath::RadToDeg(),Proton.Vector.Theta()*TMath::RadToDeg()-MissingPTot.Theta()*TMath::RadToDeg());
							ResoThetaVSPPim->Fill(PionM.Vector.P(),Proton.Vector.Theta()*TMath::RadToDeg()-MissingPTot.Theta()*TMath::RadToDeg());
						}

					}
				}
			}
		}

	}


	TFile* EffNew=new TFile("EffNew.root","recreate");
	EffNew->cd();

	TCanvas *canEfficiency = new TCanvas("canEfficiency","",12000,6000);
	canEfficiency->Divide(4,2);
	for(int Pbin=0;Pbin<4;Pbin++){
		for(int thbin=0;thbin<2;thbin++){
			TH1F *effMomPbis= (TH1F*)Eff[Pbin][thbin]->Clone(Form("mm%d%d",Pbin,thbin));
			effMomPbis->Divide(mmEff[Pbin][thbin]);
			canEfficiency->cd(1+Pbin+thbin*4);effMomPbis->Draw();
			effMomPbis->Write();
			
		}
	}
	canEfficiency->SaveAs("canEfficiency.pdf");
	EffNew->Close();

	// Fit slices projected along Y fron bins in X [7,32] with more than 20 bins  in Y filled
	PhiDeltaPproton->FitSlicesY();
	ThetaDeltaPproton->FitSlicesY();

	ResoPVSThetaE->FitSlicesY(0,0,-1,20);
	ResoPVSPhiE->FitSlicesY(0,0,-1,20);
	ResoPVSPE->FitSlicesY(0,0,-1,20);
	ResoPVSThetaPip->FitSlicesY(0,0,-1,20);
	ResoPVSPhiPip->FitSlicesY(0,0,-1,20);
	ResoPVSPPip->FitSlicesY(0,0,-1,20);
	ResoPVSThetaPim->FitSlicesY(0,0,-1,20);
	ResoPVSPhiPim->FitSlicesY(0,0,-1,20);
	ResoPVSPPim->FitSlicesY(0,0,-1,20);

	TH1D *PhiDeltaPproton_0 = (TH1D*)gDirectory->Get("PhiDeltaPproton_0");
	TH1D *PhiDeltaPproton_1 = (TH1D*)gDirectory->Get("PhiDeltaPproton_1");

	TH1D *ThetaDeltaPproton_0 = (TH1D*)gDirectory->Get("ThetaDeltaPproton_0");
	TH1D *ThetaDeltaPproton_1 = (TH1D*)gDirectory->Get("ThetaDeltaPproton_1");

	TH1D *ResoPVSThetaE_1 = (TH1D*)gDirectory->Get("ResoPVSThetaE_1");
	TH1D *ResoPVSPhiE_1 = (TH1D*)gDirectory->Get("ResoPVSPhiE_1");
	TH1D *ResoPVSPE_1 = (TH1D*)gDirectory->Get("ResoPVSPE_1");
	TH1D *ResoPVSThetaPip_1 = (TH1D*)gDirectory->Get("ResoPVSThetaPip_1");
	TH1D *ResoPVSPhiPip_1 = (TH1D*)gDirectory->Get("ResoPVSPhiPip_1");
	TH1D *ResoPVSPPip_1 = (TH1D*)gDirectory->Get("ResoPVSPPip_1");
	TH1D *ResoPVSThetaPim_1 = (TH1D*)gDirectory->Get("ResoPVSThetaPim_1");
	TH1D *ResoPVSPhiPim_1 = (TH1D*)gDirectory->Get("ResoPVSPhiPim_1");
	TH1D *ResoPVSPPim_1 = (TH1D*)gDirectory->Get("ResoPVSPPim_1");


	TH1F *effMomP= (TH1F*)MomP->Clone("MomP");
	TH2F *effMomPhiP= (TH2F*)MomPhiP->Clone("MomPhiP");
	TH2F *effThetaPhiP= (TH2F*)ThetaPhiP->Clone("ThetaPhiP");
	TH2F *effThetaPhiP1= (TH2F*)ThetaPhiP1->Clone("ThetaPhiP1");
	TH2F *effThetaPhiP2= (TH2F*)ThetaPhiP2->Clone("ThetaPhiP2");
	TH2F *effThetaPP= (TH2F*)ThetaPP->Clone("ThetaPP");
	TH1F *effThetaP= (TH1F*)ThetaP->Clone("ThetaP");
	TH1F *effPhiP= (TH1F*)PhiP->Clone("PhiP");
	effMomP->Divide(mmMomP); 
	effMomPhiP->Divide(mmMomPhiP); 
	effThetaPhiP->Divide(mmThetaPhiP); 
	effThetaPhiP1->Divide(mmThetaPhiP1); 
	effThetaPhiP2->Divide(mmThetaPhiP2); 
	effThetaPP->Divide(mmThetaPP); 
	cout<<effMomP->GetBinContent(1)<<endl;
	effThetaP->Divide(mmThetaP); 
	effPhiP->Divide(mmPhiP); 
	/*ResoP->Fit("gaus");
	  ResoPhi->Fit("gaus");
	  ResoTheta->Fit("gaus");*/ 
	//MMassCentralRaw->Scale((MMassCentral->GetMaximum())/(MMassCentralRaw->GetMaximum()));
	//MMassCentral->Fit("gaus"); 
	bMMassPTot->Fit("gaus"); 
	bMMassPpTot->Fit("gaus"); 
	bMMassPmTot->Fit("gaus"); 

	/*new TCanvas;
	  tmass->Draw();

	  new TCanvas;
	  t->Draw();


	  new TCanvas;
	  essai->Draw();*/

	TCanvas *cantrack = new TCanvas("cantrack","",6000,6000);
	cantrack->Divide(3,3);
	cantrack->cd(1);phitrack->Draw();
	cantrack->cd(2);Chi2phitrack->Draw();
	TH1F * copy = (TH1F*)Chi2phitrack->Clone("Chi2phitrack");
	cantrack->cd(3);copy->Divide(phitrack);copy->Draw("hist");
	cantrack->cd(4);Resophitrack->Draw();
	TH1F * copy2 = (TH1F*)Resophitrack->Clone("Resophitrack");
	cantrack->cd(5);copy2->Divide(phitrack);copy2->Draw("hist");
	cantrack->cd(7);thetatrack->Draw();
	cantrack->cd(8);Chi2thetatrack->Draw();
	TH1F * copytheta = (TH1F*)Chi2thetatrack->Clone("Chi2thetatrack");
	cantrack->cd(9);copytheta->Divide(thetatrack);copytheta->Draw("hist");
	cantrack->SaveAs("cantrack.pdf");

	TCanvas *canCC  = new TCanvas("canCC","",20000,10000);
	canCC->Divide(3,3);
	canCC->cd(1);PhiDeltaPproton->Draw("col");
	canCC->cd(3);PhiDeltaPproton_0->Draw();
	canCC->cd(2);PhiDeltaPproton_1->Draw();
	canCC->cd(4);ThetaDeltaPproton->Draw("col");
	canCC->cd(6);ThetaDeltaPproton_0->Draw();
	canCC->cd(5);ThetaDeltaPproton_1->Draw();


	TCanvas *canCD  = new TCanvas("canCD","",20000,10000);
	canCD->Divide(3,3);
	canCD->cd(1);ResoPVSThetaE->Draw("col");
	canCD->cd(2);ResoPVSPhiE->Draw("col");
	canCD->cd(3);ResoPVSPE->Draw("col");
	canCD->cd(4);ResoPVSThetaPip->Draw("col");
	canCD->cd(5);ResoPVSPhiPip->Draw("col");
	canCD->cd(6);ResoPVSPPip->Draw("col");
	canCD->cd(7);ResoPVSThetaPim->Draw("col");
	canCD->cd(8);ResoPVSPhiPim->Draw("col");
	canCD->cd(9);ResoPVSPPim->Draw("col");
	canCD->SaveAs("MomentumCorrectionCheck.pdf");
	
	
	TCanvas *canCDProton  = new TCanvas("canCDProton","",5000,3000);
	canCDProton->Divide(3,3);
	canCDProton->cd(1);ResoPVSTheta->Draw("col");
	
	canCDProton->cd(2);ResoPVSP->Draw("col");
	canCDProton->cd(3);PVSPhi->Draw("col");
	
	ResoPVSPhi->FitSlicesY(0,0,-1,20);
	ResoPVSPhi1->FitSlicesY(0,0,-1,20);
	ResoPVSPhi2->FitSlicesY(0,0,-1,20);
	f1Pproton->SetLineWidth(2);
	f1Pproton->SetLineColor(kBlue);
	f2Pproton->SetLineWidth(2);
	f2Pproton->SetLineColor(kBlue);
	f3Pproton->SetLineWidth(2);
	f3Pproton->SetLineColor(kBlue);
	TH1D *ResoPVSPhi_1 = (TH1D*)gDirectory->Get("ResoPVSPhi_1");
	TH1D *ResoPVSPhi1_1 = (TH1D*)gDirectory->Get("ResoPVSPhi1_1");
	TH1D *ResoPVSPhi2_1 = (TH1D*)gDirectory->Get("ResoPVSPhi2_1");
	canCDProton->cd(7);ResoPVSPhi_1->SetMaximum(0.5);ResoPVSPhi_1->SetMinimum(-0.5);ResoPVSPhi_1->Fit("pol1","","",-90.,-20.);f1Pproton->Draw("same");
	canCDProton->cd(8);ResoPVSPhi1_1->SetMaximum(0.5);ResoPVSPhi1_1->SetMinimum(-0.5);ResoPVSPhi1_1->Fit("pol1","","",-90.,30.);f2Pproton->Draw("same");
	canCDProton->cd(9);ResoPVSPhi2_1->SetMaximum(0.5);ResoPVSPhi2_1->SetMinimum(-0.5);ResoPVSPhi2_1->Fit("pol1","","",50.,130.);f3Pproton->Draw("same");
	canCDProton->cd(4);ResoPVSPhi->Draw("col");
	canCDProton->cd(5);ResoPVSPhi1->Draw("col");
	canCDProton->cd(6);ResoPVSPhi2->Draw("col");
	
	canCDProton->SaveAs("MomentumCorrectionProton.pdf");
	
	TCanvas *canCDProton1  = new TCanvas("canCDProton1","",5000,4000);
	canCDProton1->Divide(2,2);

	
	canCDProton1->cd(1);ResoPVSPhi->Draw("col");f1Pproton->Draw("same");
	canCDProton1->cd(2);ResoPVSPhi1->Draw("col");f2Pproton->Draw("same");
	canCDProton1->cd(3);ResoPVSPhi2->Draw("col");f3Pproton->Draw("same");
	
	canCDProton1->SaveAs("MomentumCorrectionProton1.pdf");
	canCDProton1->SaveAs("MomentumCorrectionProton1.root");
	
	TCanvas *canCDa  = new TCanvas("canCDa","",20000,10000);
	canCDa->Divide(3,3);
	canCDa->cd(1);ResoPVSThetaE_1->Draw();
	canCDa->cd(2);ResoPVSPhiE_1->Draw();
	canCDa->cd(3);ResoPVSPE_1->Draw();
	canCDa->cd(4);ResoPVSThetaPip_1->Draw();
	canCDa->cd(5);ResoPVSPhiPip_1->Draw();
	canCDa->cd(6);ResoPVSPPip_1->Draw();
	canCDa->cd(7);ResoPVSThetaPim_1->Draw();
	canCDa->cd(8);ResoPVSPhiPim_1->Draw();
	canCDa->cd(9);ResoPVSPPim_1->Draw();
	canCDa->SaveAs("MomentumCorrectionCheck1.pdf");
		

	TCanvas *canCDp  = new TCanvas("canCDp","",20000,10000);
	canCDp->Divide(3,3);
	canCDp->cd(1);PVSThetaE->Draw("col");
	canCDp->cd(2);PVSPhiE->Draw("col");
	canCDp->cd(3);PVSPE->Draw("col");
	canCDp->cd(4);PVSThetaPip->Draw("col");
	canCDp->cd(5);PVSPhiPip->Draw("col");
	canCDp->cd(6);PVSPPip->Draw("col");
	canCDp->cd(7);PVSThetaPim->Draw("col");
	canCDp->cd(8);PVSPhiPim->Draw("col");
	canCDp->cd(9);PVSPPim->Draw("col");
	canCDp->SaveAs("CheckMomentumPDpendance.pdf");

	/*	

		TCanvas *canF  = new TCanvas("canF","",10000,12000);
		canF->Divide(3,3);
		canF->cd(1);UElec->Draw("col");
		canF->cd(2);VElec->Draw("col");
		canF->cd(3);WElec->Draw("col");
		canF->cd(4);UPionP->Draw("col");
		canF->cd(5);VPionP->Draw("col");
		canF->cd(6);WPionP->Draw("col");
		canF->cd(7);UPionM->Draw("col");
		canF->cd(8);VPionM->Draw("col");
		canF->cd(9);WPionM->Draw("col");

	 */
	TCanvas *cant  = new TCanvas("cant","",14000,12000);
	cant->Divide(2,1);
	cant->cd(1);MMass->Draw();
	MMassCentral->SetLineColor(3);
	cant->cd(2);MMassCentralRaw->Draw("hist");MMassCentral->Draw("same");
	
	TCanvas *canMass  = new TCanvas("canMass","");
	MMassCentral->Fit("gaus","","",0.89,0.99);MMassCentral->SaveAs("MassProton.root");
	canMass->SaveAs("MassProton.pdf");

	TCanvas *can  = new TCanvas("can","",14000,12000);
	can->Divide(4,5);
	can->cd(1);MMass->Draw();
	MMassCentral->SetLineColor(3);
	can->cd(2);MMassCentralRaw->Draw("hist");MMassCentral->Draw("same");
	can->cd(3);MMassTheta->Draw("col");
	can->cd(4);MMassW->Draw("col");
	can->cd(5);MMassThetaE->Draw("col");
	can->cd(6);MMassThetaPip->Draw("col");
	can->cd(7);MMassThetaPim->Draw("col");
	can->cd(9);MMassPE->Draw("col");
	can->cd(10);MMassPPip->Draw("col");
	can->cd(11);MMassPPim->Draw("col");
	can->cd(14);MMassChi2Pp->Draw("col");
	can->cd(15);MMassChi2Pm->Draw("col");
	can->cd(13);MMassChi2E->Draw("col");
	can->cd(17);MMassMmP->Draw("col");
	can->cd(18);ThetaPW->Draw("col");
	can->cd(19);ThetaPQ2->Draw("col");
	can->cd(20);ThetaPproton->Draw("col");
	can->SaveAs("pippim.pdf");


	TCanvas *can0  = new TCanvas("can0","",20000,10000);
	can0->Divide(3,3);
	can0->cd(1);MMassTot->Draw();bMMassTot->Draw("same");
	can0->cd(2);MMassPTot->Draw();bMMassPTot->Draw("same");
	can0->cd(5);MMassPpTot->Draw();bMMassPpTot->Draw("same");
	can0->cd(4);MMassPmTot->Draw();bMMassPmTot->Draw("same");
	can0->cd(3);InvariantMassPP->Draw();
	can0->cd(6);plotQ2W->Draw("col");
	can0->SaveAs("pippim0.pdf");
bMMassTot->SaveAs("MMass.root");
bMMassTot1->SaveAs("MMass1.root");

bPitheta->SaveAs("ThetaPi.root");
bPitheta1->SaveAs("ThetaPi1.root");

	TCanvas *can1  = new TCanvas("can1","",20000,10000);
	can1->Divide(4,3);
	can1->cd(1);MomP->Draw();
	can1->cd(2);ThetaP->Draw();
	can1->cd(3);PhiP->Draw();
	can1->cd(4);MomPhiP->Draw("colz");
	can1->cd(5);mmMomP->Draw();
	can1->cd(6);mmThetaP->Draw();
	can1->cd(7);mmPhiP->Draw();
	can1->cd(8);mmMomPhiP->Draw("colz");
	/*effMomP->GetYaxis()->SetRangeUser(0., 1.);
	effThetaP->GetYaxis()->SetRangeUser(0., 1.);
	effPhiP->GetYaxis()->SetRangeUser(0., 1.);*/
	can1->cd(9);effMomP->Draw("E1");
	can1->cd(10);effThetaP->Draw("E1");
	can1->cd(11);effPhiP->Draw("E1");
	can1->cd(12);effMomPhiP->Draw("colz");
	can1->SaveAs("pippimRatio1D.pdf");
	effThetaP->SaveAs("effThetaP.root");
	effPhiP->SaveAs("effPhiP.root");
	effMomP->SaveAs("effMomP.root");

	TCanvas *candd  = new TCanvas("candd","",20000,10000);
	candd->Divide(3,3);
	/*effMomPhiP->GetZaxis()->SetRangeUser(0., 1.);
	effThetaPhiP->GetZaxis()->SetRangeUser(0., 1.);
	effThetaPP->GetZaxis()->SetRangeUser(0., 1.);*/
	candd->cd(1);MomPhiP->Draw("colz");
	candd->cd(2);ThetaPhiP->Draw("colz");
	candd->cd(3);ThetaPP->Draw("colz");
	candd->cd(4);mmMomPhiP->Draw("colz");
	candd->cd(5);mmThetaPhiP->Draw("colz");
	candd->cd(6);mmThetaPP->Draw("colz");
	candd->cd(7);effMomPhiP->Draw("colz");
	candd->cd(8);effThetaPhiP->Draw("colz");
	candd->cd(9);effThetaPP->Draw("colz");
	candd->SaveAs("pippimRatio2D.pdf");
	
	TCanvas *candd1  = new TCanvas("candd1","",20000,10000);
	candd1->Divide(3,3);
	/*effMomPhiP->GetZaxis()->SetRangeUser(0., 1.);
	effThetaPhiP->GetZaxis()->SetRangeUser(0., 1.);
	effThetaPP->GetZaxis()->SetRangeUser(0., 1.);*/
	candd1->cd(1);ThetaPhiP->Draw("colz");
	candd1->cd(2);ThetaPhiP1->Draw("colz");
	candd1->cd(3);ThetaPhiP2->Draw("colz");
	candd1->cd(4);mmThetaPhiP->Draw("colz");
	candd1->cd(5);mmThetaPhiP1->Draw("colz");
	candd1->cd(6);mmThetaPhiP2->Draw("colz");
	candd1->cd(7);effThetaPhiP->Draw("colz");
	candd1->cd(8);effThetaPhiP1->Draw("colz");
	candd1->cd(9);effThetaPhiP2->Draw("colz");
	candd1->SaveAs("pippimRatio2DbinInMom.pdf");

	TCanvas *canElec  = new TCanvas("canElec","",2000,1500);
	PThetaE->Draw("colz");
	canElec->SaveAs("ElecKine.pdf");
	
	TCanvas *can3  = new TCanvas("can3","",2000,1000);
	can3->Divide(3,1);
	can3->cd(1);ResoP->Fit("gaus");
	can3->cd(2);ResoTheta->Fit("gaus");
	can3->cd(3);ResoPhi->Fit("gaus");
	can3->SaveAs("resolutions.pdf");
	/*TCanvas *can2  = new TCanvas("can2","",10000,12000);
	  can2->Divide(3,5);
	  can2->cd(1);ResoP->Draw();
	  can2->cd(2);ResoTheta->Draw();
	  can2->cd(3);ResoPhi->Draw();
	  can2->cd(4);ResoThetaVSThetaP->Draw("col");
	  can2->cd(5);ResoThetaVSPhiP->Draw("col");
	  can2->cd(6);ResoThetaVSPP->Draw("col");
	  can2->cd(7);ResoThetaVSThetaE->Draw("col");
	  can2->cd(8);ResoThetaVSPhiE->Draw("col");
	  can2->cd(9);ResoThetaVSPE->Draw("col");
	  can2->cd(10);ResoThetaVSThetaPip->Draw("col");
	  can2->cd(11);ResoThetaVSPhiPip->Draw("col");
	  can2->cd(12);ResoThetaVSPPip->Draw("col");
	  can2->cd(13);ResoThetaVSThetaPim->Draw("col");
	  can2->cd(14);ResoThetaVSPhiPim->Draw("col");
	  can2->cd(15);ResoThetaVSPPim->Draw("col");
	  can2->SaveAs("pippimCheck.png");
	 */
	cout<<n<<endl;
	cout<<n1<<endl;
	return 1;
}
