/*Author J. Van de Walle CERN July 2007
 edited by Liam Gaffney
 This program has NO guarantee of being BUG FREE !!!*/
#include <iostream>
#include "stdio.h"
#include "TMath.h"
#include "TRandom.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TVector3.h"
#include "TTree.h"
#include "TStyle.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"
using namespace std;

void plotkinematics(double e_proj, double Eex_p, double Aproj, double Atarget)
{
	
	gStyle->SetOptStat(0);
	TH2F* EvsThCM=new TH2F("EvsThCM","Energy [MeV] vs CM Scattering Angle [deg]",180,0,180,1000,0,1000);
	TH2F* EvsThLAB=new TH2F("EvsThLAB","Energy [MeV] vs LAB Scattering Angle [deg]",180,0,180,1000,0,1000);
	TH2F* ThLABvsThCM=new TH2F("ThLABvsThCM","Theta LAB vs Theta CM",180,0,180,180,0,180);
	TH2F* ProjLABTargLAB=new TH2F("ProjLABTargLAB","Projectile LAB vs Target LAB;Projectile angle;Target angle",180,0,180,90,0,90);
	
	double pi=TMath::Pi();
	double r2d=180./pi;
	double Etilde=e_proj-(Eex_p*(1+(Aproj/Atarget)));
	double tau=(Aproj/Atarget)*sqrt(e_proj/Etilde);
	double tau_t=sqrt(e_proj/Etilde);
	double e_ejec=0.;
	
	TGraph* gEvsThCMProj=new TGraph();
	TGraph* gEvsThLABProj=new TGraph();
	TGraph* gEvsThCMEjec=new TGraph();
	TGraph* gEvsThLABEjec=new TGraph();
	TGraph* gThLABProjvsThCM=new TGraph();
	TGraph* gThLABEjecvsThCM=new TGraph();
	TGraph* gProjLABTargLAB=new TGraph();
	gEvsThCMEjec->SetLineColor(2);
	gEvsThLABEjec->SetLineColor(2);
	gEvsThCMEjec->SetLineWidth(2);
	gEvsThLABEjec->SetLineWidth(2);
	gEvsThCMProj->SetLineWidth(2);
	gEvsThLABProj->SetLineWidth(2);
	gThLABEjecvsThCM->SetLineColor(2);
	gThLABEjecvsThCM->SetLineWidth(2);
	gThLABProjvsThCM->SetLineWidth(2);
	gProjLABTargLAB->SetLineWidth(2);
	
	for( int i=0;i<=180;i++ ){
		double theta_cm_proj=i*1.*pi/180.;
		double theta_cm_ejec=pi-theta_cm_proj;
		
		double theta_lab_proj=TMath::ATan(TMath::Sin(theta_cm_proj)/(TMath::Cos(theta_cm_proj)+tau));
		double theta_lab_ejec=TMath::ATan(TMath::Sin(theta_cm_ejec)/(TMath::Cos(theta_cm_ejec)+tau_t));
		
		if( theta_lab_proj<0. ){ theta_lab_proj=theta_lab_proj+pi; };
		
		e_proj=pow(Atarget/(Atarget+Aproj),2)*(1+pow(tau,2)+(2*tau*TMath::Cos(theta_cm_proj)))*Etilde;
		e_ejec=(Atarget*Aproj/pow(Atarget+Aproj,2))*(1+pow(tau_t,2)+(2*tau_t*TMath::Cos(theta_cm_ejec)))*Etilde;
		
		gThLABEjecvsThCM->SetPoint(i,theta_cm_proj*r2d,theta_lab_ejec*r2d);
		gThLABProjvsThCM->SetPoint(i,theta_cm_proj*r2d,theta_lab_proj*r2d);
		
		gEvsThCMProj->SetPoint(i,theta_cm_proj*r2d,e_proj);
		gEvsThLABProj->SetPoint(i,theta_lab_proj*r2d,e_proj);
		gEvsThCMEjec->SetPoint(i,theta_cm_proj*r2d,e_ejec);
		gEvsThLABEjec->SetPoint(i,theta_lab_ejec*r2d,e_ejec);
		gProjLABTargLAB->SetPoint(i,theta_lab_proj*r2d,theta_lab_ejec*r2d);
		
	};
	
	TLine* l1=new TLine(0,16.4,180,16.4);
	TLine* l2=new TLine(0,53,180,53);
	TLine* l3=new TLine(16.4,0,16.4,180);
	TLine* l4=new TLine(53,0,53,180);
	
	TLegend* leg1= new TLegend(0.7,0.7,0.9,0.9);
	TLegend* leg2= new TLegend(0.7,0.7,0.9,0.9);
	TLegend* leg3= new TLegend(0.7,0.7,0.9,0.9);
	
	leg1->SetFillColor(0);
	leg2->SetFillColor(0);
	leg3->SetFillColor(0);
	
	leg1->SetTextSize(0.04);
	leg2->SetTextSize(0.04);
	leg3->SetTextSize(0.04);
	
	leg1->AddEntry(gEvsThCMProj,"Projectile","l");
	leg1->AddEntry(gEvsThCMEjec,"Ejectile","l");
	
	leg2->AddEntry(gEvsThLABProj,"Projectile","l");
	leg2->AddEntry(gEvsThLABEjec,"Ejectile","l");
	
	leg3->AddEntry(gThLABProjvsThCM,"Projectile","l");
	leg3->AddEntry(gThLABEjecvsThCM,"Ejectile","l");
	leg3->AddEntry(l1,"CD inner strip","l");
	leg3->AddEntry(l2,"CD outer strip","l");
	
	l1->SetLineWidth(2);
	l2->SetLineWidth(2);
	l1->SetLineStyle(2);
	l2->SetLineStyle(2);
	
	TCanvas* CM=new TCanvas("CM","CM",700,500);
	CM->cd();
	EvsThCM->Draw();
	gEvsThCMProj->Draw("L");
	gEvsThCMEjec->Draw("L");
	leg1->Draw();
	
	TCanvas* LAB=new TCanvas("LAB","LAB",700,500);
	LAB->cd();
	EvsThLAB->Draw();
	gEvsThLABProj->Draw("L");
	gEvsThLABEjec->Draw("L");
	leg2->Draw();
	
	TCanvas* ANG=new TCanvas("ANGLES","ANGLES",700,500);
	ANG->cd();
	ThLABvsThCM->Draw();
	gThLABProjvsThCM->Draw("L");
	gThLABEjecvsThCM->Draw("L");
	l1->Draw("L");
	l2->Draw("L");
	leg3->Draw();
	
	TCanvas* LABvsLAB=new TCanvas("LABvsLAB","LABvsLAB",700,500);
	LABvsLAB->cd();
	ProjLABTargLAB->Draw();
	gProjLABTargLAB->Draw("L");
	l1->Draw("L");
	l2->Draw("L");
	l3->Draw("L");
	l4->Draw("L");
	
	
}
void AngleKinemtaics(double e_proj, double Eex_p, double Aproj, double Atarget, double thlp, double thle)
{
	
	
	double pi=TMath::Pi();
	double r2d=180./pi;
	double d2r=pi/180.;
	double Etilde=e_proj-(Eex_p*(1+(Aproj/Atarget)));
	double tau=(Aproj/Atarget)*sqrt(e_proj/Etilde);
	double tau_t=sqrt(e_proj/Etilde);
	double e_ejec=0.;
	
	double thetacm=0.;
	
	if( thlp!=0 ){
		double det=pow(1/TMath::Tan(thlp*d2r),2);+1-(tau*tau);
		double costhetacm=-tau+(1/TMath::Tan(thlp*d2r))*sqrt(det);
		costhetacm=costhetacm/(1+pow(1/TMath::Tan(thlp*d2r),2));
		thetacm=TMath::ACos(costhetacm);
	};
	
	cout<<thetacm*r2d<<endl;
	
};
