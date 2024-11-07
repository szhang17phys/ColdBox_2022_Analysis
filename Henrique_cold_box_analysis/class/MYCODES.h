#ifndef memorydepth
#define memorydepth 5000
#endif

#ifndef MYCODES
#define MYCODES

#include <fstream>
#include <iostream>
#include "Riostream.h"
#include <time.h>       // time_t, struct tm, difftime, time, mktime
#include <cmath>
#include <numeric>



#include "TROOT.h"
#include "TLatex.h"
#include <TMinuit.h>
#include <TStyle.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <THStack.h>
#include <TMath.h>
#include <TTree.h>
#include <TRandom3.h>
#include <vector>
#include <TGraph2D.h>
#include <TPolyLine3D.h>
#include <TLine.h>
#include <TTimeStamp.h>




using namespace std;

#include "variables.C"
#include "SignalProcessing.C"
#include "timeReader.C"
#include "readingCodes.C"
#include "wiener_filter.C"
#include "calibrationCodes.C"
#include "alpha_analysis.C"

#include "calibrationCodes_CT_validation.C"


#endif


