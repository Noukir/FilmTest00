#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TNtuple.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>

#define THRESHOLD 20  //CAEN
#define HISTWIDTH 330 
#define HISTMINOS 30

static const int _DT5751DataSize = 2090; // byte 20170613_1415~~
static const int _DT5751Length = 1029;   // 20170613_1415~~

typedef struct{
	unsigned int header[8];
	unsigned short waveform[_DT5751Length];
} DT5751WFdata;

int main(int argc, char **argv){
	if (argc != 3)	{
		std::cerr << "Useage: " << argv[0] << " wfb-filename root-filename" << std::endl;
		return 1;
	}

	TFile *rootf = new TFile(argv[2], "RECREATE");
	std::ifstream *ifs = new std::ifstream(argv[1]);
	if (!ifs->good()) return 1;
	DT5751WFdata *wf = NULL;
	char *buf = new char[_DT5751DataSize];

	Float_t *x = new Float_t[_DT5751Length];
	Float_t *y = new Float_t[_DT5751Length];
	Float_t *z = new Float_t[_DT5751Length];
	Float_t Sareat[HISTWIDTH]={0};
	Float_t SSarea[HISTWIDTH]={0};
	Float_t dSSarea[HISTWIDTH]={0};
	Float_t t_std;
	
	int count=0;
	int ndate=0;
	
	for (int k = 0; k < _DT5751Length; k++){
		x[k] = static_cast<Float_t>(k) * 1.00; // [nsec]
	}
	TH1F *h = 0;
	TH1F *hped = 0;
	TGraph *g0 = 0;
	TGraph *g1 = 0;
	TGraph *s1 = 0;
	TGraph *s2 = 0;
	TGraph *s3 = 0;
	TGraph *std = 0;
	TGraph *std_nonfit = 0;
	TGraph *Mode = 0;
	TGraph *Landau = 0;

	std::ostringstream os;

	TNtuple *nt = new TNtuple("nt", "nt", "event:peak0:peakV0:rise0:fall0:area0");
	TNtuple *fitnt = new TNtuple("fitnt", "fitnt", "t:mu_gaus:sigma_gaus");

	os.str("");
	os << "hped";
	hped = new TH1F(os.str().c_str(), os.str().c_str(), 4096, 0, 4096);

	while (1){
		if (ifs->eof()) break;
		ifs->read(buf, _DT5751DataSize);
		wf = (DT5751WFdata *)(buf);

	//Pedestal
		if (h) delete h;
		os.str("");
		os << "" << std::setw(8) << std::setfill('0') << wf->header[4];
		h = new TH1F(os.str().c_str(), os.str().c_str(), 4096, 0, 4096);
		for(int k=0;k<_DT5751Length;k++){
			y[k]= wf-> waveform[k]; // tatejiku
			if((k<50||k>400)&&y[k]>970) h->Fill(y[k]);
		}
		h->Fit("gaus", "Q", "goff");
		Float_t pedestalVal   = h->GetFunction("gaus")->GetParameter(1);
		Float_t pedestalSigma = h->GetFunction("gaus")->GetParameter(2);
		hped->Fill( pedestalVal );

	// original wave
		if (g0)	delete g0;
		g0 = new TGraph(_DT5751Length, x, y);
		os.str("");
		os << "wf_00_" << std::setw(8) << std::setfill('0') << wf->header[4];
		g0->SetName(os.str().c_str());
		g0->SetTitle(os.str().c_str());

		const int N = 4;
		int peakPos[N] = {0};
		Float_t peakVal[N] = {0};
		int risePos[N] = {0}, fallPos[N] = {0};
		Float_t area[N] = {0};
		Float_t Sarea[HISTWIDTH]= {0};

	//Convert 凸 wave 
		for (int k = 0; k < _DT5751Length; k++){
				z[k] = pedestalVal - y[k];
		}

	// 凸 wave (pedestal - original)
		if (g1) delete g1;
		g1 = new TGraph(_DT5751Length, x, z);
		os.str("");
		os << "subtractedwf_00_" << std::setw(8) << std::setfill('0') << wf->header[4];
		g1->SetName(os.str().c_str());
		g1->SetTitle(os.str().c_str());

		for(int i=0; i<N; i++){ //multi-peak search
		//Find peak
			for (int k = 0; k < _DT5751Length; k++){		
				for( int j=0; j<i; j++ ){
					if( risePos[j] <= k && k <= fallPos[j] ){//ignore pre-peak
						k = fallPos[j];
					}
				}
				if ( z[k] > peakVal[i] ){
					peakVal[i] = z[k];
					peakPos[i] = k;
				}
			}
		//Find rise
			for (int k = peakPos[i]; k >= 0; k--){
				if (z[k] <= 2 * pedestalSigma){
					risePos[i] = k;
					break;
				}
			}
		//Find fall
			for (int k = peakPos[i]; k <= _DT5751Length; k++){
				if (z[k] <= 2 * pedestalSigma){
						fallPos[i] = k;
						break;
				}
			}
/* hanchihaba
			float xPlus[N] = {0}, xMinus[N] = {0};
			for (int jj = peakPos[i]; jj < peakPos[i] + 100; jj++){
				if (z[jj] > 0.5 * peakVal[i] && z[jj + 1] < 0.5 * peakVal[i]){
					//keisan
					xPlus[i] = (0.5 * peakVal[i] - (z[jj] - (z[jj + 1] - z[jj]) * jj)) / (z[jj + 1] - z[jj]);
					// xPlus ga deru
					break;
				}
			}
			for (int jj = peakPos[i]; jj > 0; jj--){
				if (z[jj] > 0.5 * peakVal[i] && z[jj - 1] < 0.5 * peakVal[i]){
					//keisan
					xMinus[i] = (0.5 * peakVal[i] - (z[jj] - (z[jj - 1] - z[jj]) * jj)) / (z[jj - 1] - z[jj]);
					// xMinus ga deru
					break;
				}
			}
*/
		//noise cut
			if( peakVal[i] <= THRESHOLD )		  						  peakPos[i] = 100000; //too small peak cut
			if( fallPos[i] - risePos[i] <= 20 )	  						  peakPos[i] = 100000; //too sharp peak cut
			if( risePos[i] - HISTMINOS < 0 )			  				  peakPos[i] = 100000; //left-hamidashi-wave cut (for sekibunhakei)
			if( risePos[i] - HISTMINOS + HISTWIDTH - 1 >= _DT5751Length ) peakPos[i] = 100000; //right-hamidashi-wave cut (for sekibunhakei)
			for(int j=0; j<N; j++){
				if( j!=i && (risePos[i] - HISTMINOS) <= peakPos[j] && peakPos[j] <= (risePos[i] - HISTMINOS + HISTWIDTH - 1) )
					  peakPos[i] = 100000; //wave kasanari cut
			}

		//Double_t area = 0; //0927kara
			/*for (int j = 80; j < 350; j++){
				//    for(int j=(peakPos[i]-startTotalAreaOffset[i]); j<(zeroPos);j++){}
				area[i] += z[j];
			}*/
		}//multi peak cut end
	
	//t-sort (peakV_jun) -> (jikan_jun)
		int p=0, r=0 ,f=0;
		Float_t v=0;
		for(int n=N-2; 0<=n; n--){
			for(int i=0; i<=n; i++){
				if( peakPos[i] >= peakPos[i+1] ){
					p=peakPos[i]; v=peakVal[i]; r = risePos[i]; f = fallPos[i];
					peakPos[i] = peakPos[i+1];
					peakVal[i] = peakVal[i+1];
					risePos[i] = risePos[i+1];
					fallPos[i] = fallPos[i+1];
					peakPos[i+1] = p; peakVal[i+1] = v; risePos[i+1] = r; fallPos[i+1] = f; 
				}
			}
		}

	//noise -> 0
		for(int i=0; i<N; i++){
			if( peakPos[i] == 100000 ) peakPos[i] = 0;
		}

	//sekibunhakei
		if( peakPos[0] != 0 ){
			ndate++;
			count = 0;
			for(int k=risePos[0] - HISTMINOS; count < HISTWIDTH; k++){
				area[0] += z[k];
				Sarea[count] = area[0];
				Sareat[count] = count;
				count++;
			}
			
		//sekibunhakei-normlize
			for(int k=0; k < HISTWIDTH; k++){
				Sarea[k] = Sarea[k]/Sarea[HISTWIDTH-1];
			}

		//sekibunhakei-goukei
			for(int i=0; i<HISTWIDTH; i++){
				SSarea[i] += Sarea[i];
			}

		//sekibunhakei->s1
			if (s1) delete s1;
			s1 = new TGraph(HISTWIDTH, Sareat, Sarea);
			os.str("");
			os << "Sarea_subtractedwf_00_" << std::setw(8) << std::setfill('0') << wf->header[4];
			s1->SetName(os.str().c_str());
			s1->SetTitle(os.str().c_str());
			s1->Write();
		}

		nt->Fill(wf->header[4], peakPos[0], peakVal[0], risePos[0], fallPos[0], area[0]);
		
		h->Write();
		g0->Write();
		g1->Write();
	}//while-end

//pedestal gaus-mean
	hped->Fit("gaus", "Q", "goff");
	Float_t pedmean = hped->GetFunction("gaus")->GetParameter(1);
	Float_t pedSigma = hped->GetFunction("gaus")->GetParameter(2);
	hped->Write();
	std::ofstream ofs("pedtStd.dat");
	ofs<< pedmean << std::endl;
	ofs<< pedSigma << std::endl;

//sekibunhakei-goukei mean
	for(int i=0; i< HISTWIDTH; i++){
		SSarea[i] = SSarea[i]/(Float_t)ndate;	
	}
	Float_t max=0;
	int maxt=0;
//sekibunhakei-goukei mean bibun
	for(int i=0; i<HISTWIDTH-1; i++){
		dSSarea[i] = ( SSarea[i+1] - SSarea[i] ) / ( Sareat[i+1] - Sareat[i] );
		if(dSSarea[i] >= max ) {
			maxt = Sareat[i];
			max = dSSarea[i];
		}
	}

	t_std =( (SSarea[maxt] + SSarea[maxt + 1])*0.5 ) / SSarea[HISTWIDTH-1];
	std::cout<<"t_std = "<< t_std << std::endl;
	std::cout<<"ndate = "<< ndate << std::endl;
	
	if (s2) delete s2;
	s2 = new TGraph(HISTWIDTH, Sareat, SSarea);
	os.str("");
	os << "Sarea_mean";
	s2->SetName(os.str().c_str());
	s2->SetTitle(os.str().c_str());
	s2->Write();

	for(int i=0; i<HISTWIDTH; i++){
		Sareat[i] = i + 0.5;
	}
	if (s3) delete s3;
	s3 = new TGraph(HISTWIDTH, Sareat, dSSarea);
	os.str("");
	os << "dSarea_mean";
	s3->SetName(os.str().c_str());
	s3->SetTitle(os.str().c_str());
	s3->Write();

	ifs->close();
	delete ifs;


//////////t_std wo koteishite StdWave wo tukuruyo//////////


	std::ifstream *iifs = new std::ifstream(argv[1]);
	if (!iifs->good()) return 1;

	Float_t StdWaveFit[HISTWIDTH]={0};
	Float_t StdWaveMean[HISTWIDTH]={0};
	Float_t StdWaveLandau[HISTWIDTH]={0};

//動的二次元配列
	Float_t **WaveHist;
	WaveHist = new Float_t*[ndate];
	for(int i=0; i < ndate; i++){
		WaveHist[i] = new Float_t[HISTWIDTH];
	}
	for(int i=0; i < ndate; i++){
		for(int j=0; j < HISTWIDTH; j++){
			WaveHist[i][j] = 0;
		}
	}
	
	ndate = 0;
	while (1){
		if (iifs->eof()) break;
		iifs->read(buf, _DT5751DataSize);
		wf = (DT5751WFdata *)(buf);
		const int N = 4;
		int peakPos[N] = {0};
		Float_t peakVal[N] = {0};
		int risePos[N] = {0}, fallPos[N] = {0};
		Float_t area[N] = {0};
		Float_t Sarea[HISTWIDTH]= {0};
		
		int t_std_Pos[N] = {0};
	
	//Pedestal
		if (h) delete h;
		os.str("");
		os << "" << std::setw(8) << std::setfill('0') << wf->header[4];
		h = new TH1F(os.str().c_str(), os.str().c_str(), 4096, 0, 4096);
		for(int k=0;k<_DT5751Length;k++){
			y[k]= wf-> waveform[k]; // tatejiku
			if((k<50||k>400)&&y[k]>970) h->Fill(y[k]);
		}
		h->Fit("gaus", "Q", "goff");
		Float_t pedestalVal   = h->GetFunction("gaus")->GetParameter(1);
		Float_t pedestalSigma = h->GetFunction("gaus")->GetParameter(2);

	//Convert 凸 wave 
		for (int k = 0; k < _DT5751Length; k++){
				z[k] = pedestalVal - y[k];
		}

		for(int i=0; i<N; i++){ //multi-peak search
		//Find peak
			for (int k = 0; k < _DT5751Length; k++){		
				for( int j=0; j<i; j++ ){
					if( risePos[j] <= k && k <= fallPos[j] ){//ignore pre-peak
						k = fallPos[j];
					}
				}
				if ( z[k] > peakVal[i] ){
					peakVal[i] = z[k];
					peakPos[i] = k;
				}
			}
		//Find rise
			for (int k = peakPos[i]; k >= 0; k--){
				if (z[k] <= 2 * pedestalSigma){
					risePos[i] = k;
					break;
				}
			}
		//Find fall
			for (int k = peakPos[i]; k <= _DT5751Length; k++){
				if (z[k] <= 2 * pedestalSigma){
						fallPos[i] = k;
						break;
				}
			}

		//noise cut
			if( peakVal[i] <= THRESHOLD )								  peakPos[i] = 100000; //small peak cut
			if( fallPos[i] - risePos[i] <= 20 )	  						  peakPos[i] = 100000; //too sharp peak cut
			if( risePos[i] - HISTMINOS < 0 )			  				  peakPos[i] = 100000; //left-hamidashi-wave cut (for sekibunhakei)
			if( risePos[i] - HISTMINOS + HISTWIDTH - 1 >= _DT5751Length ) peakPos[i] = 100000; //right-hamidashi-wave cut (for sekibunhakei)
			for(int j=0; j<N; j++){
				if( j!=i && (risePos[i] - HISTMINOS) <= peakPos[j] && peakPos[j] <= (risePos[i] - HISTMINOS + HISTWIDTH - 1) )
					  peakPos[i] = 100000; //wave kasanari cut
			}
		}//multi peak cut end
	
	//t-sort (peakV_jun) to (jikan_jun)
		int p=0, r=0 ,f=0;
		Float_t v=0;
		for(int n=N-2; 0<=n; n--){
			for(int i=0; i<=n; i++){
				if( peakPos[i] >= peakPos[i+1] ){
					p=peakPos[i]; v=peakVal[i]; r = risePos[i]; f = fallPos[i];
					peakPos[i] = peakPos[i+1];
					peakVal[i] = peakVal[i+1];
					risePos[i] = risePos[i+1];
					fallPos[i] = fallPos[i+1];
					peakPos[i+1] = p; peakVal[i+1] = v; risePos[i+1] = r; fallPos[i+1] = f; 
				}
			}
		}

	//noise -> 0
		for(int i=0; i<N; i++){
			if( peakPos[i] == 100000 ) peakPos[i] = 0;
		}

	//sekibunhakei
		if( peakPos[0] != 0 ){
			count = 0;
			for(int k=risePos[0] - HISTMINOS; count < HISTWIDTH; k++){
				area[0] += z[k];
				Sarea[count] = area[0];
				count++;
			}
		//sekibunhakei-normlize
			for(int k=0; k < HISTWIDTH; k++){
				Sarea[k] = Sarea[k]/Sarea[HISTWIDTH-1];
			}

		//find t_std_Pos[i]
			for(int k=0; k < HISTWIDTH; k++){
				if( Sarea[k] >= t_std ){
					t_std_Pos[0] = k + risePos[0] - HISTMINOS;
					break;
				}
			}

			if( t_std_Pos[0] - HISTMINOS < 0 ) peakPos[0] = 0; //left-hamidashi-wave cut
			if( t_std_Pos[0] - HISTMINOS + HISTWIDTH - 1 >= _DT5751Length ) peakPos[0] = 0; //right-hamidashi-wave cut
		//hakei wo kugitte tasuyo
			count = 0;
			if( peakPos[0] != 0){
				for(int k=t_std_Pos[0] - HISTMINOS; count < HISTWIDTH; k++){		
					WaveHist[ndate][count] = z[k];	
					count++;
				}
				ndate++;
			}
		}
	}//while end

	std::cout<<"Making StdWave"<<std::endl;
	std::cout<< ndate << std::endl;
	TH1F *wh=0;
	Float_t StdWaveSigma[HISTWIDTH] = {0};
	Float_t maxwave[HISTWIDTH]={0};
	Float_t bunsan[HISTWIDTH]={0};
	Float_t StdWaveMode[HISTWIDTH]={0};

	for(int i=0; i < HISTWIDTH; i++){
		if (wh) delete wh;
		os.str("");
		os << "WaveHist_"<< std::setw(3) << std::setfill('0') << i ;
		wh = new TH1F(os.str().c_str(), os.str().c_str(), 211, -10, 200);

		for(int k=0; k < ndate; k++){
			StdWaveMean[i] += WaveHist[k][i];
			wh->Fill( WaveHist[k][i] );
			if( WaveHist[k][i] >= maxwave[i] ) maxwave[i] = WaveHist[k][i];
		}
	//std wave Sample mean
		StdWaveMean[i] = StdWaveMean[i] / (Float_t)ndate;

	//Sample bunsan
		for(int k=0; k<ndate; k++){
			bunsan[i] += ( WaveHist[k][i] - StdWaveMean[i] ) * ( WaveHist[k][i] - StdWaveMean[i] );
		}
		bunsan[i] = bunsan[i] / ndate;

		wh->SetName(os.str().c_str());
		wh->SetTitle(os.str().c_str());
		wh->Write();

	//Fitting with Gaus
		TF1 *mygaus = new TF1("mygaus","[0]*exp( -0.5*( (x - [1]) / [2] )*( (x - [1]) / [2] ) ) ",-10, maxwave[i] );
				
		mygaus->SetParameter(1, wh->GetMaximumBin() - 10 );
		mygaus->SetParameter(2, sqrt(bunsan[i]) );
		mygaus->SetParameter(0, 5000);

		wh->Fit("mygaus", "Q", "goff", wh->FindFirstBinAbove( (wh->GetMaximum())*0.7 ) - 10, wh->FindLastBinAbove( (wh->GetMaximum())*0.5 ) -10 );
	//	WaveHist_027->Fit("mygaus", "Q", "goff", WaveHist_027->FindFirstBinAbove( (WaveHist_027->GetMaximum())*0.5 )/2.0 - 11, WaveHist_027->FindLastBinAbove( (WaveHist_027->GetMaximum())*0.5 )/2.0 -11 );
		
		StdWaveSigma[i] = wh->GetFunction("mygaus")->GetParameter(2);
		StdWaveFit[i] = wh->GetFunction("mygaus")->GetParameter(1);		

		wh->Fit("landau", "Q", "goff", -10, maxwave[i] );
		StdWaveLandau[i] = wh->GetFunction("landau")->GetParameter(1);

		StdWaveMode[i] = wh->GetMaximumBin()-10 ;

		
		//std::cout<< i <<" "<< "GausFit" <<" "<< "Mean" <<" "<< "Mode" <<" "<< "landauFit" <<std::endl;
		std::cout<< i <<" "<< StdWaveFit[i] <<" "<< StdWaveMean[i] <<" "<< StdWaveMode[i] <<" "<< StdWaveLandau[i] <<std::endl;

		fitnt->Fill( i, StdWaveFit[i], StdWaveSigma[i]);
	}

//Set yokojiku value
	count = 0;
	Float_t t[HISTWIDTH]={0};
	for(int i = -HISTMINOS; count < HISTWIDTH; i++){
		t[count] = i;
		count++;
	}

//StdWave--mean of fitting parameter
	std = new TGraph(HISTWIDTH, t, StdWaveFit);
	os.str("");
	os << "StdWaveFit";
	std->SetName(os.str().c_str());
	std->SetTitle(os.str().c_str());
	std->Write();

//StdWave--Mode
	Mode = new TGraph(HISTWIDTH, t, StdWaveMode);
	os.str("");
	os << "StdWaveMode";
	std->SetName(os.str().c_str());
	std->SetTitle(os.str().c_str());
	std->Write();

//StdWave--Landau
	Landau = new TGraph(HISTWIDTH, t, StdWaveLandau);
	os.str("");
	os << "StdWaveLandau";
	std->SetName(os.str().c_str());
	std->SetTitle(os.str().c_str());
	std->Write();

//StdWave--Sample mean 
	std_nonfit = new TGraph(HISTWIDTH, t, StdWaveMean);
	os.str("");
	os << "StdWaveMean";
	std_nonfit->SetName(os.str().c_str());
	std_nonfit->SetTitle(os.str().c_str());
	std_nonfit->Write();

	std::ofstream oofs("StdWaveMean.dat");
	for(int i=0; i<HISTWIDTH; i++){
		oofs<< t[i] <<" "<< StdWaveMean[i] <<" "<< StdWaveSigma[i] << std::endl;
	}oofs.close();

	ofs<< t_std ;
	ofs.close();

//2D map of StdWave
	TH2F *wh2d = new TH2F("wh2d", "wh2d", 330, -30, 300, 211, -10, 200);
	for(int i=0; i < HISTWIDTH; i++){
		for(int k=0; k < ndate; k++){
			wh2d->Fill( t[i], WaveHist[k][i] );
		}
	}
	wh2d->Write();

	iifs->close();
	delete iifs;
	
	nt->Write();
	fitnt->Write();
	rootf->Close();

	for(int i=0; i < ndate; i++){
		delete[] WaveHist[i];
	}
	delete[] WaveHist;

	return 0;
}