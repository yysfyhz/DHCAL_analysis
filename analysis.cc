
#include <fstream>
#include <iostream>
#include <string>
#include "TH1.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TH3.h"
#include "TH3I.h"
#include "TH3D.h"
#include "TF1.h"
#include "TROOT.h"
#include "TClonesArray.h"
#include "TProfile.h"
#include "TMath.h"
#include <TApplication.h>
#include <TDirectory.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TText.h>
#include <TTree.h>
#include <TLegend.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TError.h>
#include <vector>
#include <list>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sstream>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <math.h>
#include <ctime>
#include <stdlib.h>

using namespace std;


//The class Event contains the time and x, y, z positions vectors
//of an event as well as the statistics of the event.
class Event {
	//A vector of integers which contains the time at which each hit was registered
	std::vector<int> t;
	//A vector of integers which contains the time difference to the trigger of each hit that was registered
	vector<int> deltat;
	//Vectors of char which contain the (x, y, z) position at which hits were registered
	std::vector<char> x;
	std::vector<char> y;
	std::vector<char> z;
	std::vector<float> CalibH;
	//Number of hits in the event
	int numPoints;
	//Maximum Z layer index of the detector
	int maxChannelZ;
	//Maximum Z layer reached by the event
	int maxZ;

	//Timebin with the most hits for an event
	int t0;
	//Contains the timing information for an event

	vector<int> deltat0;
	//Vector containing the Density information
	vector<int> density;
	vector<int> density3D;


	//Contains the Cerenkov information
	int ckov;
	//Associated an ID to the Event
	int eventno;

	//A vector of integers which contains the number of hits registered by each Z layer and for the bottom middle and top RPCs
	std::vector<int> noHitsZ;
	std::vector<int> noHitsZb;
	std::vector<int> noHitsZm;
	std::vector<int> noHitsZt;
	//No Hits for an Events
	int noHits;
	//A vector of integers which contains the mean X (or Y) for each Z layer 
	std::vector<double> meanX;
	std::vector<double> meanY;

	//Mean x and y of the event
	float meanx;
	float meany;
	float meanl;
	//Number of active layers
	int nactivelayer;
	//calibrated Hitssum
	double hitssum;
	//Calibrated hits using weights
	double weightcalhits;
	vector<double> sum;
	//Boolean saying if Event showered and the layer it started
	bool didshower;
	int showerstart;
	//Defines type
	string type;


	//Z layer at which a shower started
        int depth;

	public:
	Event();
	Event(vector<int>, vector<int>, int);
	~Event();
	void deleteTime();
	void deletePoints();
	void removeHitsXChannel();

	//Accessors
	vector<int> getT();
	vector<int> getDeltaT();
	vector<int> getDensity();
	vector<int> getDensity3D();
	vector<char> getX();
	vector<char> getY();
	vector<char> getZ();
	vector<float> getcalibH();
	int getNumPoints();
	int getNoHits();
	int getNActiveLayer();
	vector<int> getNoHitsZ();
	vector<int> getNoHitsZb();
	vector<int> getNoHitsZm();
	vector<int> getNoHitsZt();
	double getHitssum();
	double getWeightHits();
	vector<double> getSum();
	double getLeakage();
	bool getdidShower();
	int getshowerStart();
	string getPID();
	vector<double> getMeanX();
	vector<double> getMeanY();
	int gett0();
	vector<int> getdeltat0();
	int getCkov();
	int getEventNo();
	int getMaxZ();
	//Mutators
	void setNumPoints(int);
	void setMaxChannelZ(int);
	int checkDoubleHits(vector<char>&,vector<char>&,vector<char>&,vector<int>&,vector<float>&, vector<int>&, vector<int>&,vector<int>&, int&);//Check if the same Position fired twice during an Event and erase the information of the event
	int checkTimingHits(vector<char>&,vector<char>&,vector<char>&,vector<int>&,vector<float>&, vector<int>&,vector<int>&, vector<int>&, int&);//Remove hits with wrong Timing 
	void removeDeadAsics(vector<char>&, vector<char>&, vector<char>&,vector<int>&,vector<int>&,vector<int>&, int&, vector<int>&, vector<int>&, vector<int>&, vector<double>&);
	bool checkMultiple(vector<char>& X, vector<char>& Y, vector<char>& Z);

	//Event discarding

	int isDoubleHits(int, int, int, int, int);

	void buildParameters();

	//Event processing
	int processTime(int, int, int*, int*);
	void computeMean(vector<char>&,vector<char>&,vector<char>&);
	void computeCalibHits(vector<float>, vector<float>, vector<float>, vector<float>&, vector<int>, vector<int>, vector<int>, int , bool , int,vector<char>&,vector<char>&);
	void computeWeightCalibHits(vector<float>&, vector<float>& ,vector<int>&);
	void computeCalibInt(vector<float>&);
	void computeMeanLayer(vector<char>&);
	void computedidShower(bool);
	void computeNActiveLayer(vector<char>&);
	void computeDensity(vector<int>&,vector<int>&,vector<char>,vector<char>,vector<char>);
	void computeNoHits(vector<char>);
	void computeNoHitsZ(vector<char>&,vector<char>&);
	void computeNoHitsX();
	void computet0(vector<int>);
	void computedeltat0(vector<int>&);
	void computecalibH(vector<float>&);
	void computeCkov(int);
	void computeTransverse(vector<double>,vector<char>&,vector<char>&,vector<char>&);
	void computeEventNo(int);
	void computeMaxZ(int);
	void computePID(bool, int, int,int&, int&, int&, double, double, double);
	void computePIDMC(bool, int, int,int&, int&, int&, bool);
};


//The run class contains a collection of events held in listOfEvents, 
//as well as various statistics of the Run, computed in the computeStat()
//method. Event rejection, plotting and so on is done in the
//process() method of this class.
class Run{

	//A list containing pointers to all the Events in a run
	std::list<Event*> listOfEvents;
	//A vector of integers containing the number of hits in each Z layer
	int noHits;
	std::vector<int> noHitsZ;
	std::vector<int> noHitsZb;
	std::vector<int> noHitsZm;
	std::vector<int> noHitsZt;
	int energy;
	std::string period;
	std::string runName;
	string runNo;
	int nomu,nopi,nopo;

	//A time vector where each bin is 0.1 second long
	std::vector<int> timeVector;
	//A vector of time vectors...
	std::vector<std::vector<int> > timeZVector;

	//Maximum possible channel numbers
	int maxChannelX;
	int maxChannelY;
	int maxChannelZ;

	int noEvents;


	public :
	Run();
	Run(int, std::string, std::string, std::string);
	~Run();

	//Accessors
	std::list<Event*> getListOfEvents();
	int getEnergy();
	int getNoEvents();
	std::string getPeriod();

	//Run processing
	void process(std::string , bool , bool , bool , bool , bool , bool , std::string , bool );
	void setNoEvents(int);
	void plotTimeAndMaxHits();

};


/*
	Functions Declarations
*/

void makeDirectory(std::string dirName);

void initializeTimeZVector(std::vector<std::vector<int > > &timeZ, int maxChannelZ);

int find_maximumChannels(std::string filename, int *maxX, int *maxY, int *maxZ, int searchRange);


void readcalibration(string currentLine, int energy, vector<float> &calibbottom, vector<float> &calibmiddle,  vector<float> &calibtop, double eff, double mult);


void readdeadasics(string currentLine, int energy, vector<int> &deadX, vector<int> &deadY,  vector<int> &deadZ);

/*
	Class Run Methods
*/

//Constructors
Run::Run(){
}

Run::Run(int Energy, std::string month, std::string run, string runno){
	noEvents = 0;
	energy = Energy;
	period = month;
	runName = run;
	runNo=runno;

}

//Destructor
Run::~Run(){
	for (std::list<Event*>::const_iterator it = listOfEvents.begin(), end=listOfEvents.end(); it!=end; ++it){
		delete (*it);
	}
	std::list<Event*>().swap(listOfEvents);
}

void Run::setNoEvents(int nE){
	noEvents = nE;
}


//Accessors
std::list<Event*> Run::getListOfEvents(){	return listOfEvents;}
std::string Run::getPeriod(){		return period;}
int Run::getEnergy(){			return energy;}
int Run::getNoEvents(){			return noEvents;}



void initializeTimeZVector(std::vector<std::vector<int > > &timeZ, int maxChannelZ){
	for (int j=0; j<maxChannelZ+1; j++){
		std::vector<int> newTimeVector;
		newTimeVector.push_back(0);
		timeZ.push_back(newTimeVector);
	}
}


void plotDensity(std::list<Event*> &listOfEvents, std::string period, std::string runName, int energy);


//This helper function finds maxChannelX, Y and Z by looking at the first searchRange events
int find_maximumChannels(std::string filename, int *maxX, int *maxY, int *maxZ, int searchRange){

	//Open the file
	std::fstream file;
	file.open(filename.c_str(), std::ios::in);
	if(!file.is_open()){
		cout<<"Error opening the file:"<<filename.c_str() << endl;
		return 0;
	}

	//Create a new current event
	//Event * curEvent = new Event;
	int tim, xpos, ypos, zpos;
	int i=0;
	int NoEvents=0;

	while (file >> tim && i<searchRange)
	{
		//file >> tim;
		file >> xpos;
		file >> ypos;
		file >> zpos;
		if(xpos<0) {
			NoEvents++;
		}
		if (xpos>*maxX && xpos < 105){
			*maxX=xpos;
		}
		if (ypos>*maxY && ypos < 105){
			*maxY=ypos;
		}
		if (zpos>*maxZ && zpos < 100){
			*maxZ=zpos;
		}
	}
	return NoEvents;
	file.close();
}


//This is the most important part of the program, where events are rejected,
// This is done by the read_options() function 
//which stores the data of each event in listOfEvents.
void Run::process(std::string filename, bool discardDouble, bool discardTiming, bool discardActive, bool discardMultiple, bool calibrate, bool PID, std::string currentDirectoryStr, bool ispion) 
{
	//Mean efficiency and multiplicity
	double eff=0.917;
	double mult=1.573;
	//Old calibration
	//double eff=0.953;
	//double mult=1.569;
	//Create folders
	std::stringstream dirNStream; 	
	dirNStream << energy << "Gev_"<< runName;	
	std::string dirName = dirNStream.str();

	std::string energyFolder;
	if (energy != 0){
		std::stringstream energyStream;
		energyStream << energy << "Gev";
		energyFolder=energyStream.str();
	}
	else{
		energyFolder="Muons";
	}

	chdir(currentDirectoryStr.c_str());
	cout << "~Making directory " << currentDirectoryStr<<"/"<<endl;


	//Find the dimensions of the detector for this run	
	maxChannelX=0;
	maxChannelY=0;
	maxChannelZ=0;
	int NoEvents=find_maximumChannels(filename, &maxChannelX, &maxChannelY, &maxChannelZ, 1000);
	cout<<"Number of Events "<<NoEvents<<endl;
	maxChannelZ=49;

	//Initialize timeZVector	
	initializeTimeZVector(timeZVector, maxChannelZ);

	//Read the calibration file

	vector<float> calibbottom;
	vector<float> calibmiddle;
	vector<float> calibtop;

	if (calibrate){
		chdir("/users/detector/bfreund/private/Analysis/New/Calibration/");
		readcalibration(runNo, energy, calibbottom, calibmiddle, calibtop, eff, mult);
		calibbottom.push_back(1);
		calibmiddle.push_back(1);
		calibtop.push_back(1);

	}
	
/*
	Read the file, discard and identify events
*/

	//Open the file
	std::fstream file;
	file.open(filename.c_str(), std::ios::in);
	if(!file.is_open()){
		cout << "Error opening file:"<< filename.c_str()<< endl;
		return;
	}


 	nomu=0;
	nopi=0;
	nopo=0;
	int tim, xpos, ypos, zpos;
	int EventNo=0;
	vector<char> X, Y, Z;
	vector<float> calibH;
	vector<int> T, DeltaT, Density, Density3D;
	vector<float> gausmean, gausmeancal;
	
	int d=0;
	int noPoints=0;
	int noEvents=0;
	int totalNoHits=0;
	int okayNoEvents=0;
	int connectorhit=0;
	int noTiming=0;
	int noActive=0;

	int initTime=0;
	int ttrigger=0;
	int Ckovtemp,Ckov;
	int counter=0;
	int last=0;
	int nodoublehits=0;
	int maxT=0;

	int lessThanDiscarded=0;
	int moreThanDiscarded=0;
	int noDoubleHits=0;
	int noTimingHits=0;

	int noEventsTiming=0;
	int noAvtrans=0;
	int multiple=0;
	double cutpi=0;
	double cutp=0;

	if (energy==1){
		cutpi=2.2;
		cutp=2.5;
	}
	else if (energy==2){
		cutp=4;
		cutpi=2.5;
	}
	else if (energy==3){
		cutp=4;
		cutpi=2.2;
	}
	else if (energy==4){
		cutp=6;
		cutpi=3;
	}
	else if (energy==6){
		cutp=7;
		cutpi=2.5;
	}
	else if (energy==8){
		cutp=10;
		cutpi=2.5;
	}
	else if (energy==10){
		cutp=11;
		cutpi=2.8;
	}
	else if (energy==32){
		cutp=15;
		cutpi=3;
	}

	cout << "~Processing all events... \n";

	time_t beforeP=time(0);


	//Read each row in the file
	while (file >> tim)
	{
		file >> xpos;
		file >> ypos;
		file >> zpos;
		totalNoHits++;

		//If this is the first time stamp...
		if (noEvents==0 && totalNoHits<2){
			initTime=tim;
		}

		//If we have a new event starting here, it is time to save the points
		//in X, Y, Z to an event and process it
		if (xpos == -1){
			ttrigger=tim;
			Ckovtemp=ypos;

			//Only keep events with more than 5 hits
			if (noPoints>5){
				noEvents++;
				listOfEvents.push_back(new Event(T,DeltaT, noPoints));
				bool discarded=false;
				//Process time stamps
				counter = listOfEvents.back()->processTime(initTime, counter, &maxT, &last);
				//Store the time information in the timeVector
				for (int i=0; i<T.size(); i++){
					int sizeBefore = timeVector.size();
					//Add more bins to the timeVector when needed
					if (T[i]/100000+1 >= sizeBefore){
						for (int k=0; k<=T[i]/100000-sizeBefore+1; k++){
							timeVector.push_back(0);
							for (int j=0; j<maxChannelZ+1; j++){
								(timeZVector[j]).push_back(0);
							}
						}
					}
					timeVector[T[i]/100000]++; //?
					(timeZVector[(int) Z[i]])[T[i]/100000]++;
				}
				if (discardDouble){
					listOfEvents.back()->checkDoubleHits(X, Y, Z, T, calibH, DeltaT, Density, Density3D, noDoubleHits);
				}
				if (discardTiming){
					listOfEvents.back()->checkTimingHits(X, Y, Z, T, calibH, DeltaT, Density, Density3D, noTimingHits);
				}

				//Initialize the event

				listOfEvents.back()->computeNoHitsZ(Z,Y);
				listOfEvents.back()->computeNoHits(Z);
				listOfEvents.back()->computeMean(X,Y,Z);
				listOfEvents.back()->computeMeanLayer(Z);
				listOfEvents.back()->setMaxChannelZ(maxChannelZ);
				listOfEvents.back()->computeMaxZ(noPoints/200);
				listOfEvents.back()->computet0(DeltaT);
				listOfEvents.back()->computeCkov(Ckov);
				listOfEvents.back()->computeEventNo(EventNo);
				listOfEvents.back()->computeDensity(Density,Density3D,X,Y,Z);
						
				//Compute calibrated hits
				if(calibrate){
					vector<int> nohitsb=listOfEvents.back()->getNoHitsZb();
					vector<int> nohitsm=listOfEvents.back()->getNoHitsZm();
					vector<int> nohitst=listOfEvents.back()->getNoHitsZt();
					bool didshower=listOfEvents.back()->getdidShower();
					if(didshower){d++;}
					listOfEvents.back()->computeCalibHits(calibbottom,calibmiddle,calibtop, calibH, nohitsb, nohitsm, nohitst,Ckov,didshower, d,Y,Z);
				}

				listOfEvents.back()->computedidShower(calibrate);

				bool diditshower=listOfEvents.back()->getdidShower();
				if(diditshower) {
					int showerstart=listOfEvents.back()->getshowerStart();
				}	



				if(PID && period!="MC"){
					bool didshower=listOfEvents.back()->getdidShower();
					int showerstart=listOfEvents.back()->getshowerStart();
					double div= listOfEvents.back()->getHitssum()/listOfEvents.back()->getNActiveLayer();
					listOfEvents.back()->computePID(didshower, showerstart, Ckov, nomu, nopi, nopo, div, cutpi, cutp);

				}

				//Remove events with wrong timing
				if(discardTiming){
					if ((listOfEvents.back()->gett0() <19 || listOfEvents.back()->gett0() >20) && discardTiming) {
						delete listOfEvents.back();
						listOfEvents.pop_back();
						noTiming++;
						discarded=true;
					}
				}


					
				//Remove events with more than four hits in first layer
				if(discardMultiple &&!discarded){
					if(listOfEvents.back()->checkMultiple(X, Y, Z)){
						delete listOfEvents.back();
						listOfEvents.pop_back();
						discarded=true;
						multiple++;
					}			
				}

				//Remove events with less than 5 active layers
				if(discardActive &&!discarded){
					if (listOfEvents.back()->getNActiveLayer()<5) {
						delete listOfEvents.back();
						listOfEvents.pop_back(); //erase the last element of the vector
						noActive++;
						discarded=true;
					}
				}
					
 
				//Delete points if the event wasn't already discarded
				if (!discarded){
					listOfEvents.back()->deletePoints();
					listOfEvents.back()->deleteTime();
				}

 
			}
 
			EventNo++;
			if(EventNo%100==0){
				cout<<"Progress "<<int(((double)EventNo/(double)NoEvents)*100)<<"%\r"<<flush;
			}
			//Prepare to read another event
			noPoints=0;
			std::vector<char>().swap(X);
			std::vector<int>().swap(T);
			std::vector<int>().swap(DeltaT);
			std::vector<float>().swap(calibH);
			std::vector<int>().swap(Density);
			vector<int>().swap(Density3D);
			std::vector<char>().swap(Y);
			std::vector<char>().swap(Z);
		}
		//If we don't have a new event starting here, add the points to the T, X, Y, Z vectors
		else if (!((xpos<2) && ((ypos>19 && ypos<25 )||(ypos>51 && ypos<57)||(ypos>83 && ypos<89))) && zpos<50) {

			noPoints++;
			DeltaT.push_back(ttrigger-tim);
			Density.push_back(0);
			Density3D.push_back(0);
			T.push_back(tim);
			X.push_back((char)xpos);
			Y.push_back((char)ypos);
			Z.push_back((char)zpos);
			calibH.push_back(0);
			Ckov=Ckovtemp;

		}
		else {
			connectorhit++;	
		}


	}
	cout<<"Progress "<<100<<"%\r"<<flush;
	cout<<endl;
	time_t afterP = time(0);
	cout << "(Took " << difftime(afterP, beforeP) << " seconds.)\n";
	//cout<<listOfEvents.back()->getEventNo()<<endl;
	file.close();


/*
	Output findings
*/
	okayNoEvents=noEvents-lessThanDiscarded-moreThanDiscarded-noTiming-noActive-multiple;
	cout << "Number of events :                               " << noEvents<<"\n";
	cout << "Total number of hits :                           " << totalNoHits<<"\n";
	if(PID){
	cout << "Total number of muons :                          " << nomu<<"\n";
	cout << "Total number of pions :                          " << nopi<<"\n";
	cout << "Total number of positrons :                      " << nopo<<"\n";
	}

	cout << "Maximum X channel :                              " << maxChannelX<<"\n";
	cout << "Maximum Y channel :                              " << maxChannelY<<"\n";
	cout << "Maximum Z channel :                              " << maxChannelZ<<"\n";
	cout << "Number of double hits :                          " << noDoubleHits<<"\n";
	cout << "Number of connector hits :                       " << connectorhit<<"\n";
	cout << "Number of bad Timing hits :                      " << noTimingHits<<" ("<<(100*(float)noTimingHits)/(float)totalNoHits<<" %)\n";
	cout << "Events with bad Timing :   	                 " << noTiming<<" ("<<(100*(float)noTiming)/(float)noEvents<<" %)\n";
	cout << "Multiple Events :   	                         " << multiple<<" ("<<(100*(float)multiple)/(float)noEvents<<" %)\n";
	cout << "Events with less than 5 active hits :            " << noActive<<" ("<<(100*(float)noActive)/(float)noEvents<<" %)\n";	
	cout << "Events with more than 4000 hits discarded :      " << moreThanDiscarded<<" ("<<(100*moreThanDiscarded)/noEvents<<" %)\n";		
	cout << "Remaining number of events:                      " << okayNoEvents<<" ("<<(100*(float)okayNoEvents)/(float)noEvents<<" %)\n";


	//Plot density

	makeDirectory("Density");
	chdir("Density");
	makeDirectory(energyFolder);
	chdir(energyFolder.c_str());
	cout << "~Making density plots.. \n";
	plotDensity(listOfEvents, period, runName, energy);

	stringstream titleStream;
	titleStream <<"Densitytest"<<runName<<".root";
	string name = titleStream.str();
	TFile* g=TFile::Open(name.c_str(),"RECREATE");


	TH1F *hDensityp = new TH1F("", "Density of positrons ",9, -0.5, 8.5);
	hDensityp -> GetXaxis()->SetTitle("Density");
	hDensityp -> GetYaxis()->SetTitle("Normalized by the # of entries");
	hDensityp -> SetLineColor(kViolet+4);

	TH1D *hDensitycalibp = new TH1D("", "Density of positrons ",9, -0.5, 8.5);
	hDensitycalibp -> GetXaxis()->SetTitle("Density");
	hDensitycalibp -> GetYaxis()->SetTitle("Normalized by the # of entries");
	hDensitycalibp -> SetLineColor(kViolet+4);

	TH1D *hDensitycalibp3D = new TH1D("", "Density of positrons 3D ",27, -0.5, 26.5);
	hDensitycalibp3D -> GetXaxis()->SetTitle("Density");
	hDensitycalibp3D -> GetYaxis()->SetTitle("Normalized by the # of entries");
	hDensitycalibp3D -> SetLineColor(kViolet+4);

	TH1F *hDensitypi = new TH1F("", "Density of pions ",9, -0.5, 8.5);
	hDensitypi -> GetXaxis()->SetTitle("Density");
	hDensitypi -> GetYaxis()->SetTitle("Normalized by the # of entries");
	hDensitypi -> SetLineColor(kViolet+4);

	TH1F *hDensitym = new TH1F("", "Density of muons ",9, -0.5, 8.5);
	hDensitym -> GetXaxis()->SetTitle("Density");
	hDensitym -> GetYaxis()->SetTitle("Normalized by the # of entries");
	hDensitym -> SetLineColor(kViolet+4);
	

	int l=0;
	nomu=0;
	nopo=0;
	nopi=0;
	vector<float> calibHtemp;

	//Fill the histograms
	for (std::list<Event*>::const_iterator it = listOfEvents.begin(), end=listOfEvents.end(); it!=end; ++it){
		if ((*it)->getPID()=="positron") {
			vector<int> density= (*it)->getDensity();
			vector<int> density3D= (*it)->getDensity3D();
			if(period!="MC"){
				calibHtemp= (*it)->getcalibH();
			}
			for (int i=0; i<density.size();i++){				
				hDensityp->Fill((Int_t) density[i]);
				if(period!="MC"){
					hDensitycalibp->Fill((Int_t) density[i],1/(calibHtemp[i]));
					hDensitycalibp3D->Fill((Int_t) density3D[i],1/(calibHtemp[i]));
				}
				else {
					hDensitycalibp->Fill((Int_t) density[i]);
					hDensitycalibp3D->Fill((Int_t) density3D[i]);
				}
				nopo++;
			}

			stringstream titleStream;
			titleStream <<"DensityDistribution_"<<l<<".root";
			stringstream titleStream3D;
			titleStream3D <<"DensityDistribution3D_"<<l<<".root";
			string name = titleStream.str();
			string name3D = titleStream3D.str();	
			//hDensitycalibp->Write(name.c_str(),TObject::kOverwrite);
			hDensitycalibp->Reset();
			hDensitycalibp3D->Write(name3D.c_str(),TObject::kOverwrite);
			hDensitycalibp3D->Reset();
			l++;
		}
		if ((*it)->getPID()=="muon") {
			vector<int> density= (*it)->getDensity();
			for (int i=0; i<density.size();i++){				
				hDensitym->Fill((Int_t) density[i]);
				nomu++;
			}
		}
		if ((*it)->getPID()=="pion") {
			vector<int> density= (*it)->getDensity();
			for (int i=0; i<density.size();i++){				
				hDensitypi->Fill((Int_t) density[i]);
				nopi++;
			}
		}
	}
	hDensityp->Sumw2();
	hDensityp->Scale(1/(double)nopo);
	hDensitym->Sumw2();
	hDensitym->Scale(1/(double)nomu);
	hDensitypi->Sumw2();
	hDensitypi->Scale(1/(double)nopi);
	//hDensityp->Write("hDensityp",TObject::kOverwrite);
	//hDensitym->Write("hDensitym",TObject::kOverwrite);
	//hDensitypi->Write("hDensitypi",TObject::kOverwrite);
	delete hDensityp;
	delete hDensitym;
	delete hDensitypi;
	g->Close();



	//All Events are deleted from listOfEvents and the subsequent memory is cleared.
	for (std::list<Event*>::const_iterator it = listOfEvents.begin(), end=listOfEvents.end(); it!=end; ++it){
		delete (*it);
	}
	listOfEvents.clear();
	list<Event*>().swap(listOfEvents);

	return;

}


/*
	Class Event Methods
*/

//Constructor
Event::Event(){


}

//Constructor
Event::Event(vector<int> T, vector<int> DeltaT, int noPoints){

	//Initialize noHitsZ and noHitsX
	for (int j=0; j<60; j++){
		noHitsZ.push_back(0);
		noHitsZb.push_back(0);
		noHitsZm.push_back(0);
		noHitsZt.push_back(0);
	}
	noHits=0;
	//for (int j=0; j<100; j++){
	//	noHitsX.push_back(0);
	//}
	//Initialize meanX and meanY
	for (int i=0; i<60; i++){
		meanX.push_back(0);
		meanY.push_back(0);
	}
	t=T;
	deltat=DeltaT;
	//x=X;
	//y=Y;
	//z=Z;
	numPoints=noPoints;
}


//Destructor
Event::~Event(){

	std::vector<int>().swap(t);
	std::vector<int>().swap(deltat);
	std::vector<char>().swap(x);
	std::vector<char>().swap(y);
	std::vector<char>().swap(z);
	std::vector<double>().swap(meanX);
	std::vector<double>().swap(meanY);
	std::vector<int>().swap(noHitsZ);
	std::vector<int>().swap(noHitsZb);
	std::vector<int>().swap(noHitsZm);
	std::vector<int>().swap(noHitsZt);
}

//Delete the time vector to save memory
void Event::deleteTime(){
	std::vector<int>().swap(t);
	std::vector<int>().swap(deltat);
}
//Delete the position vectors to save memory
void Event::deletePoints(){
	std::vector<char>().swap(x);
	std::vector<char>().swap(y);
	std::vector<char>().swap(z);
}

//Accessors
vector<int> Event::getT(){			return t;}
vector<int> Event::getDeltaT(){			return deltat;}
vector<int> Event::getDensity(){		return density;}
vector<int> Event::getDensity3D(){		return density3D;}
vector<char> Event::getX(){			return x;}
vector<char> Event::getY(){			return y;}
vector<char> Event::getZ(){			return z;}
vector<float> Event::getcalibH(){		return CalibH;}
int Event::getNumPoints(){			return numPoints;}
int Event::getNoHits(){				return noHits;}
vector<int> Event::getNoHitsZ(){		return noHitsZ;}
vector<int> Event::getNoHitsZb(){		return noHitsZb;}
vector<int> Event::getNoHitsZm(){		return noHitsZm;}
vector<int> Event::getNoHitsZt(){		return noHitsZt;}
int Event::getNActiveLayer(){			return nactivelayer;}
double Event::getHitssum(){			return hitssum;}
double Event::getWeightHits(){			return weightcalhits;}
vector<double> Event::getSum(){			return sum;}
bool Event::getdidShower(){			return didshower;}
int Event::getshowerStart(){			return showerstart;}
vector<double> Event::getMeanX(){		return meanX;}
vector<double> Event::getMeanY(){		return meanY;}
int Event::gett0(){				return t0;}
vector<int> Event::getdeltat0(){		return deltat0;}
int Event::getCkov(){				return ckov;}
int Event::getEventNo(){			return eventno;}
int Event::getMaxZ(){				return maxZ;}
string Event::getPID(){				return type;}

//Mutators
void Event::setNumPoints(int noPoints){		numPoints=noPoints;}
void Event::setMaxChannelZ(int noChannels){	maxChannelZ=noChannels;}


/*
	Event: Computations
*/

//The time stamp of hits recorded by the detector goes from 0 to 10 million (in units of 100 ns). 
//When 10 million is reached, the time index goes back to 0. This method, called successively on 
//each (ordered) event of a run, sets the time of the first hit of the first event to 0 and makes 
//sure time stamps continues increasing once the 10 million mark is it. In order to no exceed the 
//integer limit (2.147 billion) in the case of long runs, the time unit is brought from 100 ns 
//to 1 us. This is useful for plotting the number of hits as a function of time histograms.

int Event::processTime(int initTime, int counter, int *max, int *last){

	for(int i=0; i<t.size(); i++){
		//The +2 is there because it can happen that
		//a run has times smaller than init
		if (t[i]+1000000<*last){
				counter++;
		}
		*last=t[i];
		//The time is brought back to us; the +2 is there in case there
		//are times smaller than initTime in the first event (usually 1 unit
		//less). That way we can avoid negative times.
		t[i]=(t[i]-initTime+2)/10 + (counter)*1000000;
		if (t[i]>*max){
			*max=t[i];
		}
	}

	return counter;
	
}

//Counts how many hits per layer (x, y or z value) are in an event
void Event::computeNoHitsZ(vector<char>& Z,vector<char>& Y){

	//Count how many hits were registered in each z level for
	//this event
	if(Z.size()>0){
		for (int j=0; j<Z.size();j++){
			if((int)Z[j]<50){
				noHitsZ[(int) Z[j]]++;
				if(Y[j]<32){		
					noHitsZb[(int) Z[j]]++;
				}
				if(Y[j]>31 && Y[j]<64){
					noHitsZm[(int) Z[j]]++;
				}
				if(Y[j]>63 && Y[j]<96){
					noHitsZt[(int) Z[j]]++;
				}
			}
		}
	}
	return;

}

void Event::computeNoHits(vector<char> Z){

	//Count how many hits were registered for
	//this event
	for (int j=0; j<Z.size();j++){
		if((int)Z[j]<50){
		noHits++;}
	}
	return;

}



void Event::computeMean(vector<char>& X,vector<char>& Y,vector<char>& Z){

	meanx=0;
	meany=0;
	int mean0X=0;
	int mean0Y=0;
	int nohits[5]={0};
	//Find the meanX and meanY of each z layer


	for (int i=0; i<X.size(); i++){
		meanx=meanx+(int) X[i];
		meany=meany+(int) Y[i];

		meanX[(int)Z[i]]=meanX[(int)Z[i]]+ (int) X[i];
		meanY[(int)Z[i]]=meanY[(int)Z[i]]+(int) Y[i];

	}

	meanx=meanx/numPoints;
	meany=meany/numPoints;

	for (int i=0; i<60; i++){
		if (noHitsZ[i]>0){
			meanX[i]= (meanX[i]/((double)noHitsZ[i]));
			meanY[i]= (meanY[i]/((double)noHitsZ[i]));
		}
		else {
			meanX[i]=-1;
			meanY[i]=-1;
		}
	}

	return;
}


void Event::computeNActiveLayer(vector<char>& Z){

	nactivelayer=0;
	//Find the number of active layers
	for (int i=0; i<50; i++){
		bool active=false;
		for (int j=0; j<Z.size(); j++){
			if ((int) Z[j]==i){
				active=true;
			}
		}	
		if (active) {
			nactivelayer++;
		}
	}
	return;
}


void Event::computeCalibHits(vector<float> calibbottom,vector<float> calibmiddle, vector<float> calibtop, vector<float>& calibH, vector<int>  nohitsb, vector<int>   nohitsm, vector<int> nohitst, int Ckov, bool didShower, int d,vector<char>& Y,vector<char>& Z){

	//Initiate the sum vectors
	vector<float> sumb,summ,sumt;
	for (int i=0; i<calibbottom.size()+1; i++){
		sumb.push_back(0);
		summ.push_back(0);
		sumt.push_back(0);
		sum.push_back(0);
	}

	//Put the calibration Factors in a Vector
	if(Z.size()>0){
		for (int j=0; j<Z.size() && (int)Z[j]<50;j++){
			if(Y[j]<32){		
				calibH[j]=calibbottom[(int)Z[j]];

			}
			if(Y[j]>31 && Y[j]<64){
				calibH[j]=calibmiddle[(int)Z[j]];
			}
			if(Y[j]>63 && Y[j]<96){
				calibH[j]=calibtop[(int)Z[j]];
			}
		}
	}


	//Set the vectors to 0
	for(int i=1; i<calibbottom.size();i++) {
		sumb[i]=0;
		summ[i]=0;
		sumt[i]=0;
		sum[i]=0;	
	}

	for(int i=0; i<calibbottom.size();i++) {
		if (calibbottom[i]>0 && calibmiddle[i]>0 && calibtop[i]>0 && nohitsb[i]>=0&& nohitsm[i]>=0&& nohitst[i]>=0){
			sumb[i+1]+=((float)(nohitsb[i]))/calibbottom[i];
			summ[i+1]+=((float)(nohitsm[i]))/calibmiddle[i];
			sumt[i+1]+=((float)(nohitst[i]))/calibtop[i];
			sum[i+1]+=sumb[i+1]+summ[i+1]+sumt[i+1];
		}
	}
	hitssum=0;
	for(int i=1; i<51;i++) {
		if(sumb[i]>=0 &&summ[i]>=0 &&sumt[i]>=0){
			hitssum=hitssum+sumb[i]+summ[i]+sumt[i];
		}
	}
	CalibH=calibH;

	return;
}


void Event::computeMeanLayer(vector<char>& Z){

	meanl=0;
	int sum=0;
	int suml=0;

	//Find the meanX and meanY of each z layer
	for (int i=0; i<noHitsZ.size(); i++){
		if( noHitsZ[i]>0) {
			sum=sum+noHitsZ[i];
			suml++;
		}
	}
	meanl=float(sum)/suml;
	
	return;
}

void Event::computedidShower(bool calibrate){

	didshower=false;
	showerstart=-1;
		//Looks for two consecutive layers with more than four hits
	for (int i=1; i<maxChannelZ+1; i++){
		if(calibrate){
			if (sum[i]>4 && sum[i-1]>4) {
				didshower=true;
				if (showerstart==-1){
					showerstart=i-1;
				}
			}
		}
		else {
			if (noHitsZ[i]>4 && noHitsZ[i-1]>4) {
				didshower=true;
				if (showerstart==-1){
					showerstart=i-1;
				}
			}
		}
	}
	return;
}


void Event::computePID(bool didshower, int showerstart, int ckov, int& nomu, int& nopi, int& nopo, double div, double cutpi, double cutp){
		if (/*ckov==0 &&*/ didshower==false && nactivelayer>40) {		
			type="muon";
			nomu++;
		}
		else if (didshower==true && showerstart>1 && showerstart<11 && ckov==0 && div>cutpi/* &&div<6*/) {		
			type="pion";
			nopi++;
		}
		else if (ckov!=0 && didshower==true && div>cutp/* && div<8&& noHitsZ[maxChannelZ]<3*/) {
			type="positron";
			nopo++;
		}	

	return;
}

//Compute the maximum Z plane where hits are registered
void Event::computeMaxZ(int capHits){

      maxZ=0;
      for (int i=maxChannelZ; i>1; i--){

		if (noHitsZ[i]>capHits && noHitsZ[i-1]>capHits && noHitsZ[i-2]>capHits){
			maxZ=i;
			i=0;
		}
      }
      return;

}


void Event::computet0(vector<int> DeltaT){

	int t0temp[23]={0};
	int temp=0;
	int temp2=0;
	t0=0;
	for (unsigned i=0; i<DeltaT.size() ; i++) {
		deltat[i]=DeltaT[i];
		for(int j=0; j<23; j++){
			if (DeltaT[i]==j) {
				t0temp[j]++;
			}
		}
	}
	int l=0;
	while ( l<23) {
		temp=t0temp[l];
		if (temp>temp2) {
			temp2=temp;
			t0=l;
		}
		l++;
	}
	return;
}

void Event::computedeltat0(vector<int>& DeltaT){
	deltat0=DeltaT;
	return;
}

void Event::computecalibH(vector<float>& calibH){
	CalibH=calibH;
	return;
}

void Event::computeDensity(vector<int>& Density, vector<int>& Density3D,vector<char> X,vector<char> Y, vector<char> Z){
	for (unsigned i=0; i<Density.size() ; i++) {
		for (int j=0; j<Density.size(); j++) {
			if(j!=i){	
				if (Z[j]==Z[i]) {				
					if(abs(X[j]-X[i])<=1&&abs(Y[j]-Y[i])<=1){
						Density[i]++;
						Density3D[i]++;
					}					
				}
				if (Z[j]==Z[i]-1) {				
					if(abs(X[j]-X[i])<=1&&abs(Y[j]-Y[i])<=1){
						Density3D[i]++;
					}					
				}
				if (Z[j]==Z[i]+1) {				
					if(abs(X[j]-X[i])<=1&&abs(Y[j]-Y[i])<=1){
						Density3D[i]++;
					}					
				}
			}
		}

	}
	density=Density;
	density3D=Density3D;
	return;
}


void Event::computeCkov(int Ckov){
	ckov=Ckov;
	return;
}

void Event::computeEventNo(int EventNo){
	eventno=EventNo;
	return;
}


//Check if the same Position fired twice during an Event; 
int Event::checkDoubleHits(vector<char>& X, vector<char>& Y, vector<char>& Z, vector<int>& T,vector<float>& calibH,vector<int>& DeltaT, vector<int>& Density, vector<int>& Density3D, int& noDoubleHits){
 	for (unsigned i=0; i<X.size() ; i++) {
		for (unsigned j=i+1; j<X.size(); j++) {	
			if (X[i]==X[j] && Y[i]==Y[j] && Z[i]==Z[j]) {				
				X.erase(X.begin() + j);
				Y.erase(Y.begin() + j);
				Z.erase(Z.begin() + j);
				T.erase(T.begin() + j);
				calibH.erase(calibH.begin() + j);
				DeltaT.erase(DeltaT.begin() + j);
				Density.erase(Density.begin()+j);
				Density3D.erase(Density3D.begin()+j);
				noDoubleHits++;					
			}
		}
	}

	return 0;

}
//Check for multiple particles
bool Event::checkMultiple(vector<char>& X, vector<char>& Y, vector<char>& Z){
	int xtemp[5]={0};
	int ytemp[5]={0};
	int l=0;
	if(noHitsZ[0]>4){
		return true;
	}
	else {
 		for (unsigned i=0; i<Z.size() ; i++) {
			if(Z[i]==0){
				xtemp[l]=X[i];
				ytemp[l]=Y[i];
				l++;
			}	
		}
		for(int j=0;j<noHitsZ[0];j++){
			for(int k=j+1;k<noHitsZ[0];k++){
				if((abs(xtemp[j]-xtemp[k])+abs(ytemp[j]-ytemp[k]))>3){
					return true;
				}
			}
		}
		return false;
	}


}

//Remove hits with wrong Timing 
int Event::checkTimingHits(vector<char>& X, vector<char>& Y, vector<char>& Z, vector<int>& T,vector<float>& calibH, vector<int>& DeltaT,vector<int>& Density, vector<int>& Density3D, int& noTimingHits){

 	for (unsigned i=0; i<DeltaT.size() ; i++) {
		if (DeltaT[i]<19 ||DeltaT[i]>20) {
			X.erase(X.begin() + i);
			Y.erase(Y.begin() + i);
			Z.erase(Z.begin() + i);
			T.erase(T.begin() + i);
			calibH.erase(calibH.begin() + i);
			DeltaT.erase(DeltaT.begin() + i);
			Density.erase(Density.begin()+i);
			Density3D.erase(Density3D.begin()+i);

			noTimingHits++;	
			i--;					
		}
	}

	return 0;

}

/*
	File processing
*/


//Shortcut function to make a directory
void makeDirectory(std::string dirName){
	int temp;
	temp = umask(0);

	return;
}

void readcalibration(string runNo, int energy, vector<float> &calibbottom, vector<float> &calibmiddle,  vector<float> &calibtop, double eff, double mult){
	string fileName = "RunCalibrationConstants_n2_v2_Oct14_trfit.txt";
	cout << "~Reading calibration File "<< fileName<<endl;
	//Open the file
	fstream file;
	file.open(fileName.c_str(), std::ios::in);
	if(!file.is_open()){
		cout<<"Error opening the Calibration file"<<endl;
		return;
	}

	std::string currentLine;
	float calibtemp;
	float temp;
	float efftemp;
	float multtemp;
	int bmt;
	int layer;
	while (file >> currentLine)
	{

		if (currentLine == runNo){
			file >> layer;
			file >> bmt;
		//cout <<"bmt"<<bmt<<endl;
			if (bmt==0) {
				file >> efftemp;
		//cout <<"efftemp"<<efftemp<<endl;
				file >> multtemp;
		//cout <<"multtemp"<<multtemp<<endl;
				calibtemp=(efftemp*multtemp)/(eff*mult);
				calibbottom.push_back(calibtemp);
				file>>temp;
				//cout<<layer<<endl;
				//cout<<calibtemp<<endl;
				//cout<<temp<<endl;
			}
			if (bmt==1) {
				file >> efftemp;
				file >> multtemp;
				calibtemp=(efftemp*multtemp)/(eff*mult);
				calibmiddle.push_back(calibtemp);
				file>>temp;
			}
			if (bmt==2) {
				file >> efftemp;
				file >> multtemp;
				calibtemp=(efftemp*multtemp)/(eff*mult);
				calibtop.push_back(calibtemp);
				file>>temp;
			}
		}
	}

	file.close();

	//Plot calibration
	std::string energyFolder;
	if (energy != 0){
		std::stringstream energyStream;
		energyStream << energy << "Gev";
		energyFolder=energyStream.str();
	}
	std::stringstream titleStream;
	titleStream << "calibrationfactor"<< runNo<< energy;
	std::string titleStr = titleStream.str();

	//File save name
	std::stringstream fileNStream;
	fileNStream << runNo << "calibrationfactor" <<"_"<<energy<<".png";
	std::string fileName2 = fileNStream.str();

	//Declare objects 
	TCanvas* canv;
	TPad* distrPad;
	TPad* titlePad;

	TText *title=new TText(0.5, 0.5, titleStr.c_str());
	title->SetTextAlign(22);
	title->SetTextSize(0.4);


	TH1D *cal = new TH1D("cal","Distribution of the calibration factor",100, 0, 3);
	cal -> GetXaxis()->SetTitle("Calibration factor");
	cal -> GetYaxis()->SetTitle("Number calibration factors");

	//Fill the histograms
	for (int i=0; i<50; i++){
		cal->Fill(calibbottom[i]);
		cal->Fill(calibmiddle[i]);
		cal->Fill(calibtop[i]);
	}

	canv = new TCanvas("canv", "Timing ", 800, 800);


	distrPad = new TPad ("distrPad", "distrPad", 0.005, 0.005, 0.995, 0.945);
	titlePad = new TPad ("titlePad", "titlePad", 0.005, 0.950, 0.995, 0.995);

	titlePad->Draw();
	distrPad->Draw();

	titlePad->cd();
	title->Draw();

	//Draw the 1D distribution histograms
	distrPad ->cd();
	cal->Draw();
	chdir("..");
	chdir("batch");
	chdir("2Calibration");
	makeDirectory(energyFolder);
	chdir(energyFolder.c_str());
	//Save the file
	canv->SaveAs(fileName2.c_str());
	chdir("..");
	chdir("..");

	delete title;
	delete distrPad;
	delete titlePad;
	delete canv;
	delete cal;

	return;
}

/*
	Plotting methods
*/

void plotDensity(list<Event*> &listOfEvents, string period, string runName, int energy) {
	//Title
	std::stringstream titleStream;
	std::stringstream fileNStream;

	fileNStream << period << "_" << runName << "_density";

	titleStream <<"The Density for every hit ";

	titleStream << energy << " GeV ("<< period << ", "<<runName<< ") " ;
	fileNStream<< "_"<<energy;

	std::string titleStr = titleStream.str();
	fileNStream<<".png";

	std::string fileName = fileNStream.str();

	//Declare objects 
	TCanvas* canv;
	TPad* distrPad;
	TPad* titlePad;

	TText *title=new TText(0.5, 0.5, titleStr.c_str());
	title->SetTextAlign(22);
	title->SetTextSize(0.4);

	//Create the histogram

	TH1F *hDensity = new TH1F("", "Density ",9, -0.5, 8.5);
	hDensity -> GetXaxis()->SetTitle("Density");
	hDensity -> GetYaxis()->SetTitle("Number of Hits");
	hDensity -> SetLineColor(kViolet+4);
	hDensity -> SetLineWidth(3);
	hDensity -> GetYaxis()->SetTitleOffset(1.4);

	gStyle->SetOptStat(1111);

	TH1F *hDensitym = new TH1F("", "Density of muons ",9, -0.5, 8.5);
	hDensitym -> GetXaxis()->SetTitle("Density");
	hDensitym -> GetYaxis()->SetTitle("Number of Hits");
	hDensitym -> SetLineColor(kViolet+4);
	hDensitym -> SetLineWidth(3);
	hDensitym -> GetYaxis()->SetTitleOffset(1.4);

	TH1F *hDensitypi = new TH1F("", "Density of pions ",9, -0.5, 8.5);
	hDensitypi -> GetXaxis()->SetTitle("Density");
	hDensitypi -> GetYaxis()->SetTitle("Number of Hits");
	hDensitypi -> SetLineColor(kViolet+4);
	hDensitypi -> SetLineWidth(3);
	hDensitypi -> GetYaxis()->SetTitleOffset(1.4);

	TH1F *hDensitypi3D = new TH1F("", "Density of pions ",27, -0.5, 26.5);
	hDensitypi3D -> GetXaxis()->SetTitle("Density");
	hDensitypi3D -> GetYaxis()->SetTitle("Number of Hits");
	hDensitypi3D -> SetLineColor(kViolet+4);
	hDensitypi3D -> SetLineWidth(3);
	hDensitypi3D -> GetYaxis()->SetTitleOffset(1.4);

	TH1F *hDensityp = new TH1F("", "Density of positrons ",9, -0.5, 8.5);
	hDensityp -> GetXaxis()->SetTitle("Density");
	hDensityp -> GetYaxis()->SetTitle("Number of Hits");
	hDensityp -> SetLineColor(kViolet+4);
	hDensityp -> SetLineWidth(3);
	hDensityp -> GetYaxis()->SetTitleOffset(1.4);

	TH1F *hDensityp3D = new TH1F("", "Density of positrons ",27, -0.5, 26.5);
	hDensityp3D -> GetXaxis()->SetTitle("Density");
	hDensityp3D -> GetYaxis()->SetTitle("Number of Hits");
	hDensityp3D -> SetLineColor(kViolet+4);
	hDensityp3D -> SetLineWidth(3);
	hDensityp3D -> GetYaxis()->SetTitleOffset(1.4);

	//Fill the histograms
	for (std::list<Event*>::const_iterator it = listOfEvents.begin(), end=listOfEvents.end(); it!=end; ++it){
		vector<int> density= (*it)->getDensity();
		vector<int> density3D= (*it)->getDensity3D();
		for (int i=0; i<density.size();i++){				
			hDensity->Fill((Int_t) density[i]);
			if ((*it)->getPID()=="muon") {		
				hDensitym->Fill((Int_t) density[i]);
			}
			if ((*it)->getPID()=="pion") {		
				hDensitypi->Fill((Int_t) density[i]);
				hDensitypi3D->Fill((Int_t) density3D[i]);
			}
			if ((*it)->getPID()=="positron") {
				hDensityp->Fill((Int_t) density[i]);
				hDensityp3D->Fill((Int_t) density3D[i]);	
			}
		}
	}

	canv = new TCanvas("", "Density ", 1600, 1600);


	distrPad = new TPad ("distrPad", "distrPad", 0.005, 0.005, 0.995, 0.945);
	titlePad = new TPad ("titlePad", "titlePad", 0.005, 0.950, 0.995, 0.995);

	titlePad->Draw();
	distrPad->Draw();
	distrPad->Divide(2,2);


	titlePad->cd();
	title->Draw();

	//Draw the 1D distribution histograms
	distrPad ->cd(1);
	hDensity->Draw();
	distrPad ->cd(2);
	hDensitym->Draw();
	distrPad ->cd(3);
	hDensitypi->Draw();
	distrPad ->cd(4);
	hDensityp->Draw();

	//Save the file
             chdir("/home/inpactest/Yusen/dhcal_TB");
	canv->SaveAs(fileName.c_str());
       
	delete hDensity;
	delete hDensitym;
	delete hDensitypi;
	delete hDensityp;
	delete title;
	delete distrPad;
	delete titlePad;
	delete canv;

	//Title
	std::stringstream titleStream3D;
	std::stringstream fileNStream3D;

	fileNStream3D << period << "_" << runName << "_density3D";

	titleStream3D <<"The 3D Density for positrons for every hit ";

	titleStream3D << energy << " GeV ("<< period << ", "<<runName<< ") " ;
	fileNStream3D<< "_"<<energy;

	std::string titleStr3D = titleStream3D.str();
	fileNStream3D<<".png";

	std::string fileName3D = fileNStream3D.str();

	//Declare objects 
	TCanvas* canv3D;
	TPad* distrPad3D;
	TPad* titlePad3D;

	TText *title3D=new TText(0.5, 0.5, titleStr3D.c_str());
	title3D->SetTextAlign(22);
	title3D->SetTextSize(0.4);
	canv3D = new TCanvas("", "Density ", 1600, 1600);


	distrPad3D = new TPad ("distrPad", "distrPad", 0.005, 0.005, 0.995, 0.945);
	titlePad3D = new TPad ("titlePad", "titlePad", 0.005, 0.950, 0.995, 0.995);

	titlePad3D->Draw();
	distrPad3D->Draw();
	distrPad3D->Divide(1,2);

	titlePad3D->cd();
	title3D->Draw();

	//Draw the 1D distribution histograms
	distrPad3D ->cd(1);
	hDensityp3D->Draw();
	distrPad3D ->cd(2);
	hDensitypi3D->Draw();

	//Save the file
	chdir("/home/inpactest/Yusen/dhcal_TB");
	canv3D->SaveAs(fileName3D.c_str());


	delete hDensityp3D;
	delete hDensitypi3D;
	delete title3D;
	delete distrPad3D;
	delete titlePad3D;
	delete canv3D;
}

int main()
{

		//Get the path of the current directory
	char currentDirectory[200];
	getcwd(currentDirectory, sizeof(currentDirectory));
	std::string currentDirectoryStr(currentDirectory);

	std::string month="Nov2011";
	bool ispion=false;
	bool doMuons;
	bool discardDouble;
	bool discardTiming;
	bool discardActive;
	bool discardMultiple;
	bool calibrate;
	bool PID;


		doMuons=true;
		discardDouble=true;
		discardTiming=true;
		discardActive=true;
		discardMultiple=true;
		calibrate=true;
		PID=true;



	//Open the period file

	std::stringstream searchFileStream;
	searchFileStream << month;
	searchFileStream << ".txt";
	std::string searchFile = searchFileStream.str();

	//Open the file

	std::fstream periodFile;
	periodFile.open(searchFile.c_str(), std::ios::in);
	if(!periodFile.is_open()){
		cout<<"Error opening the period file "<<searchFile<<endl;
		return -1;
	}

	std::string fileName;
	int energy;
	std::string run;

	//If we did select to process all runs contained in the period, the program will look for 
	//the period file which contains every run name we want to be processed in a period, 
	//listed under their energy. It will read it and make a new Run object for each run. 
	//Pointers to these objects will be stored in the listOfRuns vector. Each run is then 
	//processed by its member function process(). Once each run has been taken care of, 
	//the plotRunTiming() method is called.

	//cout << "The 'All' option was selected." << endl;


	std::vector<Run*> listOfRuns;


	chdir("batch");

	std::string currentLine;
	
	int noProcessed=0;
	std::string type;



		//Read each row
	while (periodFile.is_open())
	{

		string runName;
		string runNo;
		stringstream runNameStream;
						

		if (periodFile >> currentLine){
			//If we are reading the runs from a month, the file names are read from the month.txt files

			if (currentLine == "%"){
				periodFile >> type;
			}
			else{
				if (type == "Muons"){
					energy=0;
				}
				else{
					istringstream (type) >> energy;
				}
					
					//Process the file
				if (!(type == "Muons" && !doMuons)){


					runNo=currentLine;
					runName = std::string("run")+currentLine;
					std::stringstream fileNStream;
					if (month!="MC"){
						fileNStream << "/home/inpactest/Yusen/dhcal_TB/ascii/" << month << "/" <<"ascii/"<< runName<<".beam.txt";
					}

					fileName = fileNStream.str();
					noProcessed++;
					cout<<endl;
					cout << noProcessed << ". " << month << " " << runName<<"\n";
					cout << fileName <<"\n";

					Run *theRun = new Run(energy, month, runName, runNo);
					theRun->process(fileName,discardDouble, discardTiming, discardActive, discardMultiple, calibrate, PID, currentDirectoryStr, ispion); 
					if ((theRun->getNoEvents())>10){
						listOfRuns.push_back(theRun);	
					}

				}

			}
		}
		else {

			periodFile.close();

			return 0;
		}
	}

	return 0;
}
