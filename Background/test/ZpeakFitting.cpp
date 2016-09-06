#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <utility>
#include <math.h>

#include "TFile.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TStopwatch.h"
#include "TLatex.h"
#include "TLine.h"
#include "TF1.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooAddPdf.h"
#include "RooGaussian.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooMsgService.h"
#include "RooMinuit.h"
#include "RooPolynomial.h"
#include "RooRandom.h"

#include "boost/program_options.hpp"
#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string/classification.hpp"
#include "boost/algorithm/string/predicate.hpp"
#include "../interface/RooPowerLaw.h"
#include "../interface/RooPowerLawSum.h"
#include "RooGenericPdf.h"
#include "../interface/PdfModelBuilder.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooMultiPdf.h"
#include "HiggsAnalysis/GBRLikelihood/interface/RooDoubleCBFast.h"   

using namespace std;
using namespace RooFit;
using namespace boost;
namespace po = boost::program_options;


void runFit(RooAbsPdf *pdf, RooDataSet *data, double *NLL, int *stat_t, int MaxTries, int mhLow, int mhHigh){

	int ntries=0;
	int stat=1;
	double minnll=10e8;
	while (stat!=0){
	  if (ntries>=MaxTries) break;
	  RooFitResult *fitTest = pdf->fitTo(*data,RooFit::Save(1),Range(mhLow,mhHigh));
          stat = fitTest->status();
	  minnll = fitTest->minNll();
	  ntries++; 
	}
	*stat_t = stat;
	*NLL = minnll;
}


void GoodnessOfFit (RooRealVar *mass, RooAbsPdf *pdf, RooDataSet *data, int nFitParam, double &chi2, int &ndf, double &pval, int bins){

RooPlot *plot = mass->frame();//to convert pdf to curve
data->plotOn(plot,Name("data"),Binning(bins));
pdf->plotOn(plot,Name("pdf"));

  // Calculate the chi^2/NDOF of this curve with respect to the histogram
  // 'hist' accounting nFitParam floating parameters in case the curve
  // was the result of a fit
  // Find curve object
  RooCurve* curve = (RooCurve*) plot->findObject("pdf", RooCurve::Class());
  // Find histogram object
  RooHist* hist = (RooHist*) plot->findObject("data", RooHist::Class()) ;
  Int_t i,np = hist->GetN();
  Double_t x,y,xl,xh ;
  // Find starting and ending bin of histogram based on range of RooCurve
  Double_t xstart,xstop ;
#if ROOT_VERSION_CODE >= ROOT_VERSION(4,0,1)
  curve->GetPoint(0,xstart,y) ;
  curve->GetPoint(curve->GetN()-1,xstop,y) ;
#else
  const_cast<RooCurve*>(curve)->GetPoint(0,xstart,y) ;
  const_cast<RooCurve*>(curve)->GetPoint(curve->GetN() - 1,xstop,y) ;
#endif
  Int_t nbin(0) ;
  Int_t non0bin(0) ;
  chi2=0;
  for (i=0 ; i<np ; i++) {   
    // Retrieve histogram contents
    hist->GetPoint(i,x,y) ;
    xl = x - hist->GetEXlow()[i] ;
    xh = x + hist->GetEXhigh()[i] ;
    // Check if the whole bin is in range of curve
    if (xl < xstart || xstop < xh) continue ;
    nbin++ ;
    // Integrate function over this bin.
    // Start a hack to work around a bug in RooCurve::interpolate
    // that sometimes gives a wrong result.
    Double_t avg = curve->average(xl, xh);
    Double_t avg2 = 0.5 * (curve->average(xl, x) + curve->average(x, xh));
    if (avg + avg2 > 0 &&
	(avg2 - avg) / (avg2 + avg) > 0.1) {
      avg = curve->interpolate(x);
    }
    // End of hack around the bug in RooCurve::interpolate
    // JV: Adjust observed and expected number of events for bin width to represent
    // number of events.
    Double_t norm = (xh - xl) / plot->getFitRangeBinW();
    y *= norm;
    avg *= norm;
    double dy=hist->GetErrorY(i);
    // JV: Use the expected number of events for the y uncertainty,
    // See (33.34) of http://pdg.lbl.gov/2011/reviews/rpp2011-rev-statistics.pdf
    // Add pull^2 to chisq
    if (y != 0) {      
      Double_t resid = y - avg;
      chi2 += pow(resid/dy,2) ;
      non0bin++;
  //cout<<"bin "<<i<<" x "<<x<<" y "<<y<<" pdf "<<avg<<" dy "<<dy<<" (y-pdf)Â²/dyÂ² "<<pow(resid/dy,2)<<" chi2 "<<chi2<<endl;
    }
  }
  ndf=non0bin -nFitParam;
  pval=TMath::Prob(chi2,ndf);
   
}



void plot(RooRealVar *mass, RooAbsPdf *pdf, RooDataSet *data, string catname, string outdir, int bins, int mhLow, int mhHigh, bool runGOF=false, int nFitParams=1, bool isData=false){
 
 
  double chi2=0;
  double pval=0;
  int ndf=0;

  if(runGOF)   GoodnessOfFit(mass,pdf,data,nFitParams,chi2,ndf,pval,bins);    
 
  RooPlot *plot = mass->frame();
  mass->setRange(mhLow,mhHigh);

  data->plotOn(plot,Binning(bins));

  TCanvas *canv = new TCanvas();
  TPad *pad =new TPad("haut","haut",0,0.25,1,1);
  pad->SetNumber(1);
  pad->SetBottomMargin(0);
  pad->Draw();
  TPad *pad2 =new TPad("bas","bas",0,0,1,0.23);
  pad2->SetNumber(2);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.35);
  pad2->Draw();

  pdf->plotOn(plot);//,RooFit::NormRange("fitdata_1,fitdata_2"));
  pdf->paramOn(plot,RooFit::Layout(0.35,0.89,0.89),RooFit::Format("NEA",AutoPrecision(1)));  
  canv->cd(1);
  plot->SetMaximum(plot->GetMaximum()*1.5);    
  //plot->SetMinimum(0.0001);
  plot->SetTitle("");
  plot->SetXTitle("");
  plot->SetYTitle(Form("Events / %.2fGeV",float(mhHigh-mhLow)/float(bins)));
  //plot->GetYaxis()->SetLabelSize(0.05);
  plot->SetTitleSize(0.045, "Y");
  plot->SetTitleOffset(1.,"Y");
  plot->SetLabelColor(0,"X");
  plot->Draw();

  TLatex *lat = new TLatex();
  lat->SetNDC();
  lat->SetTextFont(42);
  if(runGOF)   {
		lat->DrawLatex(0.65,0.43,Form("chi2 / ndf = %.1f / %d",chi2,ndf));
		lat->DrawLatex(0.65,0.33,Form("pval = %.3f",pval));
	        }
  lat->DrawLatex(0.75,0.93,"2.7fb^{-1} (13TeV)"); 
  lat->DrawLatex(0.12,0.85,"CMS Preliminary");   
  lat->DrawLatex(0.1,0.93, Form("%s", catname.c_str()));  

  TLegend *leg = new TLegend(0.12,0.7,0.32,0.8);
  leg->SetFillColor(0);
  leg->SetLineColor(1);
  if(!isData)   leg->AddEntry(data,"Simulation","lep");
  else    leg->AddEntry(data,"Data","lep");
  leg->Draw("same");

  canv->cd(2);
  RooHist* hPull = plot->pullHist();
  TLine *line = new TLine(mhLow, 0., mhHigh, 0.);  
  line->SetLineWidth(1.5);
  RooPlot *plotPull = mass->frame();
  plotPull->addPlotable(hPull,"P");
  plotPull->addObject(line);
  plotPull->GetXaxis()->SetLabelSize(0.15);
  plotPull->GetYaxis()->SetLabelSize(0.15);
  plotPull->SetXTitle("m_{#gamma#gamma}(GeV)");
  plotPull->SetTitleSize(0.15, "XY");
  plotPull->SetTitleOffset(0.3,"Y");
  plotPull->SetYTitle("Pull");
  plotPull->SetTitle("");
  plotPull->GetYaxis()->SetNdivisions(6);
  //plotPull->SetMarkerSize(0.05);
  plotPull->Draw();


  //pull plots
  TCanvas *canv_Pull = new TCanvas();
  canv_Pull->cd(1);

  const int histoBins = bins;
  int yMaxIndex[histoBins];
  TMath::Sort(histoBins, hPull->GetY(), yMaxIndex, true);
  int MaxIndexNum = 0;
  if(TMath::Abs(hPull->GetY()[yMaxIndex[0]]) < TMath::Abs(hPull->GetY()[yMaxIndex[bins-1]]))  MaxIndexNum = bins-1;
  TH1F *h_pull = new TH1F("h_pull", "h_pull", int(2*TMath::Abs(hPull->GetY()[yMaxIndex[MaxIndexNum]])/0.5), -TMath::Abs(hPull->GetY()[yMaxIndex[MaxIndexNum]]), TMath::Abs(hPull->GetY()[yMaxIndex[MaxIndexNum]]));
  for(int pd=0; pd<=bins-1; pd++){
       double xp_d=0; double yp_d=0;
       hPull->GetPoint(pd, xp_d, yp_d);
       h_pull->Fill(yp_d);
  }
  h_pull->SetXTitle("pull");
  //h_pull->Rebin(2);
  //h_pull->Fit("gaus");
  TF1 *f1_gau = new TF1("gaus","gaus",-TMath::Abs(hPull->GetY()[yMaxIndex[MaxIndexNum]]), TMath::Abs(hPull->GetY()[yMaxIndex[MaxIndexNum]]));
  h_pull->Fit(f1_gau);
  h_pull->SetMaximum(h_pull->GetMaximum()*1.5);     
  f1_gau->SetLineColor(4);
  h_pull->Draw();
  f1_gau->Draw("same");

  TLatex *lat_pull = new TLatex();
  lat_pull->SetNDC();
  lat_pull->SetTextFont(42);
  lat_pull->DrawLatex(0.32,0.82,""); 
  
  string name=Form("%s/DYee_%s_DCB",outdir.c_str(),catname.c_str());

  canv->SaveAs(Form("%s.png",name.c_str()));
  canv->SaveAs(Form("%s.pdf",name.c_str()));
  canv_Pull->SaveAs(Form("%s_pull.png", name.c_str())); 
  canv_Pull->SaveAs(Form("%s_pull.pdf", name.c_str())); 
 
  delete canv;
  delete lat;

  cout<<"chi2 / ndf = "<<chi2/ndf<<" pval = "<<pval<<endl;  

}



int main(int argc, char* argv[]){

  string fileName;
  string flashggCatsStr;
  vector<string> flashggCats;
  string logfile;
  string outDir;
  string outfilename;
  int bins; 
  int mhLow; 
  int mhHigh; 
  bool runGOF=false;  
  bool isData=false;  
  bool addExp=false;  
  bool verbose=true;


po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",                                                                                  "Show help")
    ("infilename,i", po::value<string>(&fileName)->default_value("output_Zee.root"),                                              "In file name")
    ("flashggCats,f", po::value<string>(&flashggCatsStr)->default_value("UntaggedTag_0,UntaggedTag_1,UntaggedTag_2,VBFTag_0"),       "Flashgg category names to consider")
    ("logfile,l", po::value<string>(&logfile)->default_value("Zee_Yield_Zfit.log"),                  "log file of fit results")
    ("outDir,D", po::value<string>(&outDir)->default_value("plotsZpeak"),                      "Out directory for plots")
    ("outfile,o", po::value<string>(&outfilename)->default_value("output_Zee_fits.root"),         		"output file with workspace")
    ("mhLow,L", po::value<int>(&mhLow)->default_value(60),                                                                                                                 "Starting point for scan") 
    ("mhHigh,H", po::value<int>(&mhHigh)->default_value(120),                                                                                                               "End point for scan") 
    ("bins,B", po::value<int>(&bins)->default_value(25),                                                                                                                 "Bins for the plot") 
    ("runGOF",                                                                                  "Run GOF")   
    ("isData",                                                                                  "Run zee data")   
    ("addExp",                                                                                  "add Exponential")   
    ("verbose,v",                                                                               "Run with more output")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc,argv,desc),vm);
  po::notify(vm);
  if (vm.count("help")) { cout << desc << endl; exit(1); }
  if (vm.count("runGOF")) runGOF=true;  
  if (vm.count("isData")) isData=true;   
  if (vm.count("addExp")) addExp=true;   

  if (vm.count("verbose")) verbose=true;

 if (!verbose) {
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    RooMsgService::instance().setSilentMode(true);
    gErrorIgnoreLevel=kWarning;
  }

  split(flashggCats,flashggCatsStr,boost::is_any_of(","));

  string ext = "13TeV";

  TFile *outputfile;
  outputfile = new TFile(outfilename.c_str(),"RECREATE");
  RooWorkspace *outputws;
  outputws = new RooWorkspace(); 
  outputws->SetName("Zpeak");
  //FIXME : Import Lumi variable?

  ofstream logfile_stream(logfile.c_str());  

  system(Form("mkdir -p %s",outDir.c_str()));

  TFile *inFile = TFile::Open(fileName.c_str());
  TDirectoryFile *dir=(TDirectoryFile*) inFile->Get("tagsDumper");
  RooWorkspace *inWS = (RooWorkspace*)dir->Get(Form("cms_hgg_%s",ext.c_str()));

  //FIXME : To re-optimize ! Or let alpha float ...
  vector<pair<float,float> > alphaCB12;   
  alphaCB12.push_back(pair<float,float>(1.5, 1.5)); //cat 0
  alphaCB12.push_back(pair<float,float>(1.5, 1.5));   //cat 1
  alphaCB12.push_back(pair<float,float>(1.5, 1.8));   //cat 2
  alphaCB12.push_back(pair<float,float>(2., 2.)); //vbf 0

  PdfModelBuilder pdfsModel;
  RooRealVar *mass = (RooRealVar*)inWS->var("CMS_hgg_mass");
  mass->setRange(mhLow,mhHigh); 
  pdfsModel.setObsVar(mass);
	 
  unsigned int ncats=flashggCats.size();

for (unsigned int ifc =0; ifc<ncats ; ifc++){
  
    cout << "processing cat " << ifc <<" "<<flashggCats[ifc]<< endl; 

    //Get data
    RooDataSet *data = (RooDataSet*)inWS->data(Form("DY_%s_%s",ext.c_str(),flashggCats[ifc].c_str()));
    data->Print();
    //RooDataSet *data=(RooDataSet*)data0->reduce("leadPtMgg>0. && subleadPtMgg>0.");
    //data->Print();
    //FIXME : ok for Zee MC double-fake, check name in case of single-fake

    //Define pdf
    RooAbsPdf  *bkgPdf;

    float alphaCB1 = alphaCB12[ifc].first;
    float alphaCB2 = alphaCB12[ifc].second;

    if(!addExp)
	bkgPdf = pdfsModel.getDoubleCB(Form("DYee_%s_%s",ext.c_str(),flashggCats[ifc].c_str()), alphaCB1, alphaCB2);
//	bkgPdf = pdfsModel.getDoubleCB(Form("DYee_%s_%s",ext.c_str(),flashggCats[ifc].c_str()));//float alpha12
    //else 
	//bkgPdf = pdfsModel.getDoubleCBplusContinuum("Exponential",Form("DYee_%s_%s",ext.c_str(),flashggCats[ifc].c_str()),1,alphaCB1,alphaCB2,false).first;
	
    //Fit
    int fitStatus = 0; 
    double thisNll=0.;
    runFit(bkgPdf,data,&thisNll,&fitStatus,3,mhLow,mhHigh);

    //Plot
    int nFitParams;
    if(!addExp) nFitParams=4;
    //if(!addExp) nFitParams=6;//float alpha12
    else nFitParams=6;
    plot(mass,bkgPdf,data,flashggCats[ifc],outDir,bins,mhLow,mhHigh,runGOF,nFitParams,isData);

    //Write parameters walue in log file (FIXME : faster way?)

          RooArgSet *params = bkgPdf->getParameters(*data);
          float Mean = ((RooRealVar*)params->find(Form("DYee_%s_%s_DCB_mean",ext.c_str(),flashggCats[ifc].c_str())))->getValV();
          float MeanErrorL = ((RooRealVar*)params->find(Form("DYee_%s_%s_DCB_mean",ext.c_str(),flashggCats[ifc].c_str())))->getErrorLo();
          float MeanErrorH = ((RooRealVar*)params->find(Form("DYee_%s_%s_DCB_mean",ext.c_str(),flashggCats[ifc].c_str())))->getErrorHi(); 
          float Sigma = ((RooRealVar*)params->find(Form("DYee_%s_%s_DCB_sigma",ext.c_str(),flashggCats[ifc].c_str())))->getValV();
          float SigmaErrorL = ((RooRealVar*)params->find(Form("DYee_%s_%s_DCB_sigma",ext.c_str(),flashggCats[ifc].c_str())))->getErrorLo();
          float SigmaErrorH = ((RooRealVar*)params->find(Form("DYee_%s_%s_DCB_sigma",ext.c_str(),flashggCats[ifc].c_str())))->getErrorHi();
          float nCB1 = ((RooRealVar*)params->find(Form("DYee_%s_%s_DCB_nCB1",ext.c_str(),flashggCats[ifc].c_str())))->getValV();
          float nCB1ErrorL = ((RooRealVar*)params->find(Form("DYee_%s_%s_DCB_nCB1",ext.c_str(),flashggCats[ifc].c_str())))->getErrorLo();
          float nCB1ErrorH = ((RooRealVar*)params->find(Form("DYee_%s_%s_DCB_nCB1",ext.c_str(),flashggCats[ifc].c_str())))->getErrorHi();
          float nCB2 = ((RooRealVar*)params->find(Form("DYee_%s_%s_DCB_nCB2",ext.c_str(),flashggCats[ifc].c_str())))->getValV();
          float nCB2ErrorL = ((RooRealVar*)params->find(Form("DYee_%s_%s_DCB_nCB2",ext.c_str(),flashggCats[ifc].c_str())))->getErrorLo();
          float nCB2ErrorH = ((RooRealVar*)params->find(Form("DYee_%s_%s_DCB_nCB2",ext.c_str(),flashggCats[ifc].c_str())))->getErrorHi();
          //float alphaCB1 = ((RooRealVar*)params->find(Form("DYee_%s_%s_DCB_alphaCB1",ext.c_str(),flashggCats[ifc].c_str())))->getValV();
          //float alphaCB2 = ((RooRealVar*)params->find(Form("DYee_%s_%s_DCB_alphaCB2",ext.c_str(),flashggCats[ifc].c_str())))->getValV();

          logfile_stream << Form("cat: %s    ",flashggCats[ifc].c_str())<<"    nEntries     "<<data->sumEntries() << "    Mean:  " << Mean << "   MeanErrorL:  " << MeanErrorL << " MeanErrorH:  "<< MeanErrorH << "  Sigma:  " << Sigma << "  SigmaErrorL:  " << SigmaErrorL << " SigmaErrorH:  " << SigmaErrorH << "    nCB1:  " << nCB1 << "    nCB1ErrorL:   " << nCB1ErrorL << "   nCB1ErrorH:  " << nCB1ErrorH << "   nCB2:  " << nCB2 << "   nCB2ErrorL:  " << nCB2ErrorL << "    nCB2ErrorH:   " << nCB2ErrorH << "    alphaCB1:   " << alphaCB1 << "    alphaCB2:   " << alphaCB2 << endl;

	  if(addExp){
              float expParam = ((RooRealVar*)params->find(Form("DYee_%s_%s_DCB_exp_p1",ext.c_str(),flashggCats[ifc].c_str())))->getValV(); 
              float expParamErrorL = ((RooRealVar*)params->find(Form("DYee_%s_%s_DCB_exp_p1",ext.c_str(),flashggCats[ifc].c_str())))->getErrorLo();
              float expParamErrorH = ((RooRealVar*)params->find(Form("DYee_%s_%s_DCB_exp_p1",ext.c_str(),flashggCats[ifc].c_str())))->getErrorHi();
              float expCBFrac = ((RooRealVar*)params->find(Form("DYee_%s_%s_DCB_frac_sum1",ext.c_str(),flashggCats[ifc].c_str())))->getValV();
              float expCBFracErrorL = ((RooRealVar*)params->find(Form("DYee_%s_%s_DCB_frac_sum1",ext.c_str(),flashggCats[ifc].c_str())))->getErrorLo();
              float expCBFracErrorH = ((RooRealVar*)params->find(Form("DYee_%s_%s_DCB_frac_sum1",ext.c_str(),flashggCats[ifc].c_str())))->getErrorHi();
 
              logfile_stream << Form("cat: %s    ",flashggCats[ifc].c_str()) <<"   expParam:  " << expParam << "   expParamErrorL:  " << expParamErrorL << "   expParamErrorH: " << expParamErrorH << "   expCBFrac:  " << expCBFrac << "   expCBFracErrorL:  " << expCBFracErrorL << "   expCBFracErrorH:  " << expCBFracErrorH << endl;
	   }


    logfile_stream << endl;
    params->Print("v");

    //Save data and pdf in workspace
    outputws->import(*bkgPdf);
    outputws->import(*data);

  }

    outputfile->cd();
    outputws->Write();

    cout<<endl<<"+++++++++++++++++++++++++++++++++++++++++"<<endl;
    outputws->Print();
    cout<<"+++++++++++++++++++++++++++++++++++++++++"<<endl<<endl;
    
    outputfile->Close();

}

