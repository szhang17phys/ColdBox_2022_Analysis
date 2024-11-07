// This is a class created to manage waveforms



#include "TH1.h"

// using namespace std;
class Signal
{

private:
    const Double_t step = 2.;// This step might change depending on the DIGITIZER you are using
    // time should always be 0, 2, 4, 6, 8, etc... 
    
    const Int_t nbits = 14; // DIGITIZER bits
    vector<Double_t> temp2;
    
    
    Bool_t debug = false;
    
public:
    
    vector<Double_t> signal;
    vector<Double_t> time;
    vector<Double_t> stddev;
    vector<Double_t> mean;
    
    
    
    //___________________________________________________________________________________________________
    void smoothWithMovingAvarage(Signal *smooth,  Int_t interactions = 35, Double_t cut = 0){
        if(debug)cout<<"\nMoving Avarage...\n"<<endl;
        if(checkSize(this)); // check if signal and time have the same size
        else return;
        
        Int_t n = this->signal.size();
        smooth->signal.clear();
        smooth->time.clear();
        
        Int_t midpoint = 0;
        Int_t width = 0;
        Double_t sum = 0;
        
        if(interactions%2==0){ // if it is even, we insert the middle point, e.g. 8 interactions takes 4 before, mid, 4 later
            midpoint = interactions/2+1;    //midpoint will be 5 here
            width = interactions+1;
        }
        else{
            midpoint = (interactions-1)/2 + 1; // e.g. 9 interactions the midpoint will be 5
            width = interactions;
        }
        
        for(Int_t i = 0; i < n; i++){
            
            if(i<midpoint || i>(n-midpoint)){ // make it to start at i = 5 and finish at i = (3000-5) = 2995
                smooth->signal.push_back(this->signal.at(i));
            }
            else if(i<(n-cut/step)){                
                for(Int_t j = (i-midpoint); j < (i+midpoint); j++) { //first one: from j = (5-5); j<(5+5)
                    sum = sum+this->signal.at(j);
                }
                smooth->signal.push_back(sum/width);
            }
            else{
                smooth->signal.push_back(0);
            }
                        
            smooth->time.push_back(i*step);
            
            sum=0;
        }        
    }
    
    void vec_statistics(const vector<Double_t> v, Double_t &mean, Double_t &dev){
        Double_t sum = std::accumulate(v.begin(), v.end(), 0.0);
        mean = sum / v.size();
    
        std::vector<Double_t> diff(v.size());
        std::transform(v.begin(), v.end(), diff.begin(), std::bind2nd(std::minus<Double_t>(), mean));
        Double_t sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
        dev = std::sqrt(sq_sum / (v.size()-1));
    }
    
    
    //___________________________________________________________________________________________________
    void peakSearch(Signal *peaks,  Int_t lag = 100, Double_t tolerance = 5, Double_t influence = 0){
        // lag in time domain = time to start looking data
        // tolerance of xtimes std deviation
        // influence from 0 to 1. 0 more robust
        
        //this method was implemented here
        // https://stackoverflow.com/questions/22583391/peak-signal-detection-in-realtime-timeseries-data
        if(debug)cout<<"\nSearching peak with method from overflow...\n"<<endl;
        if(checkSize(this)); // check if signal and time have the same size
        else return;
        
        peaks->signal.clear();
        peaks->time.clear();
        peaks->stddev.clear();
        peaks->mean.clear();
        
        
        TH1D *hbase = new TH1D("hbase","finding baseline",2*TMath::Power(2,nbits),-TMath::Power(2,nbits),TMath::Power(2,nbits));

        Int_t n = this->signal.size();
        
        vector<Double_t> filtered = this->signal;
        vector<Double_t> zeroVec(n);
        vector<Double_t> statistics(lag/step); // initiate with same points but zero

        peaks->signal = zeroVec;
        peaks->stddev = zeroVec;
        peaks->mean = zeroVec;
        
        peaks->time = this->time;
        
        Bool_t firstHit = true;
        
        Double_t mean_st = 0;
        Double_t dev_st = 0;
        
        for(Int_t i = 0; i<n; i++){
            
            
            if(i>=lag/step){ // "ignoring" the first 'lag' points
                if(firstHit){ // just the first time ... dont need this but saves time...
                    vec_statistics(statistics,mean_st,dev_st);
                    peaks->mean.at(i-1) = mean_st;
                    peaks->stddev.at(i-1) = tolerance*dev_st;
                    firstHit = false;
                }
                
                if(abs(this->signal.at(i)-peaks->mean.at(i-1))>peaks->stddev.at(i-1)){ // if its above tolerance
                    peaks->signal.at(i) = (this->signal.at(i)>peaks->mean.at(i-1))? 10 : -10; // set the signal marking an event
                    
                    filtered.at(i) = this->signal.at(i)*influence + (1-influence)*filtered.at(i-1); //Make influence lower
                }
                else{
                    peaks->signal.at(i) = 0; // no signal
                    
                    filtered.at(i) = this->signal.at(i);
                }
                // rotetate vector, so first item will be last. Then change last item
                std::rotate(statistics.begin(),statistics.begin()+1,statistics.end());
                statistics.back() = filtered.at(i);        
                vec_statistics(statistics,mean_st,dev_st);
                peaks->mean.at(i) = mean_st;
                peaks->stddev.at(i) = tolerance*dev_st;

            }
            else{ // before the lag is done
                statistics.at(i) = (this->signal.at(i));
            }
        }
        
        
        
        delete hbase;
    }
    
    
    //___________________________________________________________________________________________________
    Double_t fillSampleVector(Double_t from, Double_t to, Signal *sample, Bool_t kUnfill = false){
        if(debug)cout<<"\nFilling sample vector...\n"<<endl;
        if(checkSize(this)); // check if signal and time have the same size
        else return 0;
        
        // check conditions for from and to. return false if from > to or from<0
        // correct if from%step or to%step != 0
        // correct last point to match time limit
        if(checkFromTo(from,to));
        else return 0;
//         
        Double_t integral = 0;
        
        if(!kUnfill){
            for(Int_t i = from/step; i<=to/step; i++){
                integral += this->signal.at(i);
                sample->signal.push_back(this->signal.at(i));
                sample->time.push_back(this->time.at(i));
            }
        }
        else{
            for(Int_t i = from/step; i<to/step; i++){
                sample->signal.pop_back();
                sample->time.pop_back();
            }
        }
        
        return integral*step;
        
    }
    
    TGraph getGraph(){
        TGraph *g = new TGraph(this->time.size(), &this->time[0], &this->signal[0]);
        return *g;
    }
        
    
    //___________________________________________________________________________________________________
    Double_t getBaseline(Double_t from, Double_t to = -1, Bool_t kAdvanced = false){
        if(debug)cout<<"\nGetting baseline...\n"<<endl;
        if(to == -1){
            to = this->time.back();
        }
        
        if(checkSize(this)); // check if signal and time have the same size
        else return 0;
        
        // check conditions for from and to. return false if from > to or from<0
        // correct if from%step or to%step != 0
        // correct last point to match time limit
        if(checkFromTo(from,to));
        else return 0;
        
        Double_t baseline = 0;
        
        TH1D *hbase = new TH1D("hbase","finding baseline",2*TMath::Power(2,nbits),-TMath::Power(2,nbits),TMath::Power(2,nbits));
        
        if(kAdvanced==false){
            for(Int_t i = from/step; i<=to/step; i++){
                hbase->Fill(this->signal.at(i));
//                 cout << this->time.at(i) << " " << this->signal.at(i) << endl;
            }
            baseline = hbase->GetMean();
        }
        else{
            cout << "This option is not yet implemented >:) ... evaluationg normal baseline" << endl;
            delete hbase;
            baseline = this->getBaseline(from,to);
        }
        
        delete hbase;
        return baseline;
        
    }
    
    //___________________________________________________________________________________________________
    Double_t integrateSignal(Double_t from, Double_t to, Double_t offset = 0){  // offset is 0 in case t0 = 0
                                                                                // Use offset!=0 otherwise

        // from and to are in time domain !
        
        if(debug)cout<<"\nIntegrating baseline...\n"<<endl;
        from = from + offset;
        
        if(checkSize(this)); // check if signal and time have the same size
        else return 0;
        
        // check conditions for from and to. return false if from > to or from<0
        // correct if from%step or to%step != 0
        // correct last point to match time limit
        if(checkFromTo(from,to));
        else return 0;
        
        Double_t inte = 0;
        
        if(from!=this->time.at(static_cast<Int_t>(from/step))){
            cout << "Make sure the time is correct evaluated!!!" << endl;
        }
        
        for(Int_t i = from/step; i<=to/step; i++){
            inte += this->signal.at(i);
        }
        return inte*step;
        
    }
    
    
    //___________________________________________________________________________________________________
    Double_t integrateSignal_withCut(Double_t from, Double_t to, vector<Double_t> smoothed){ // offset is in case t0 != 0
        // from and to are in time domain !
        
        if(debug)cout<<"\nIntegrating baseline...\n"<<endl;
        
        if(checkSize(this)); // check if signal and time have the same size
        else return 0;
        
        // check conditions for from and to. return false if from > to or from<0
        // correct if from%step or to%step != 0
        // correct last point to match time limit
        if(checkFromTo(from,to));
        else return 0;
        
        Double_t inte = 0;
        
        if(from!=this->time.at(static_cast<Int_t>(from/step))){
            cout << "Make sure the time is correct evaluated!!!" << endl;
        }
        
        Int_t strike = 0;
        
        for(Int_t i = from/step; i<=to/step; i++){
            inte += this->signal.at(i);
            if(smoothed.at(i)<-10){
                strike++;
            }
        }
        if(strike>=15){
            return -11111;
        }
        else{
            return inte*step;
        }
    }
        
    //___________________________________________________________________________________________________
    void sumUpWaveForms(const Signal s1, const Signal s2){
        if(debug)cout<<"\nSuming up waveforms...\n"<<endl;
        this->signal.clear();
        
        if(checkSize(s1,s2));//all good
        else return;
    
        for(Int_t i = 0; i<s1.signal.size(); i++){
            this->signal.push_back(s1.signal.at(i)+s2.signal.at(i));
        }
        
        return;
    }
    
    //___________________________________________________________________________________________________
    void subtractWaveForms(const Signal s1, const Signal s2){
        if(debug)cout<<"\nSubtracting waveforms...\n"<<endl;
        this->signal.clear();
        
        if(checkSize(s1,s2));//all good
        else return;
        
        for(Int_t i = 0; i<s1.signal.size(); i++){
            this->signal.push_back(s1.signal.at(i)-s2.signal.at(i));
        }
        
        return;
    }
    
    void createTime(Int_t npoints = 12000){
        if(debug)cout<<"\nCreating time vector...\n"<<endl;
        npoints = (this->signal.size()==0) ? npoints : this->signal.size();
        this->time.clear();
        for(Int_t i = 0; i<npoints; i++){
            this->time.push_back(i*step);
        }
    }
    
    //___________________________________________________________________________________________________
    void printSignal(){
        Int_t nentries = this->signal.size();
        cout << "\nPrinting signal with " << nentries << " points " << endl;
        
        if(checkSize(this)){
            for(Int_t i = 0; i<nentries; i++)
            {
                cout << this->time.at(i) << "\t" << this->signal.at(i) << endl;
            }
        }
        else{
            for(Int_t i = 0; i<nentries; i++)
            {
                cout << this->signal.at(i) << endl;
            }
        }
        return;
    }
    
    //___________________________________________________________________________________________________
    void copyTime(const Signal s){
        if(debug)cout<<"\nCopying time...\n"<<endl;
        this->time =s.time;

        return;
        
    }
    
    //___________________________________________________________________________________________________
    void copySignal(const Signal s){
        if(debug)cout<<"\nCopying signal...\n"<<endl;
        this->signal = s.signal;
        return;
    }
        
    //___________________________________________________________________________________________________    
    Bool_t checkSize(const Signal *s){ // check if the vectors signal and time from s have the same size
        // return true if they have the same size, and false otherwise
        if(debug)cout<<"\nChecking size of time vector and signal...\n"<<endl;
        if(this->signal.size()!=this->time.size()){
            cout << "\n\nError !!!                        vectors have different size!!!\n\n" << endl;
            cout << "Signal = " << this->signal.size() << " and  Time = " << this->time.size() << endl;
            return false;
        }
        else{
            return true;
        }
        
    }
    Bool_t checkSize(const Signal s1, const Signal s2){ // check if signal vectors have the same size
        if(debug)cout<<"\nChecking size between two waveforms...\n"<<endl;
        if(s1.signal.size()!=s2.signal.size()){
            cout << "\n\nError !!!                        vectors have different size!!!\n\n" << endl;
            cout << "Signal = " << this->signal.size() << " and  Time = " << this->time.size() << endl;

            return false;
        }
        else{
            return true;
        }
        
    }
    
    //___________________________________________________________________________________________________
     Bool_t checkFromTo(Double_t &from, Double_t &to){
         if(debug)cout<<"\nCheck from and to parameters...\n"<<endl;
        if(from>=to || from<0){
            cout << "\n\nIntegration limits are not orthodox... " << endl;
            return false;
        }
        if(static_cast<Int_t>(from)%static_cast<Int_t>(step)!=0){
            from = from + static_cast<Int_t>(from)%static_cast<Int_t>(step);
        }
        if(static_cast<Int_t>(to)%static_cast<Int_t>(step)!=0){
            to = to + static_cast<Int_t>(to)%static_cast<Int_t>(step);
        }
        
        
        if(this->time.size()>0 && to > this->time.back()){
            cout << "\n\nIntegration higher limit is out of range... " << endl;
            to = this->time.back();
        }
        else if(this->time.size()<=0){
            cout << "\n\nTime vector was not defined!!!" << endl;
            return false;
        }
        
        return true;
     }



    
    
    
};
