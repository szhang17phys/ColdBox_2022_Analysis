// #include "MYCODES.h"


class TimeREADER{

private:  
  
  string year = "";
  string month = "";
  string day = "";
  
  string hour = "";
  string minute = "";
  string second = "";
  TTimeStamp da;


public:

  

  Double_t timeRead(string date, string time, string format_date="dd-mmm-yyyy", string format_time="hh:mm:ss"){

    Int_t n_dates = format_date.size();
    Int_t n_times = format_time.size();
    char *dateFragment = new char[n_dates];
    char *timeFragment = new char[n_times];

    
    year = "";
    month = "";
    day = "";
    
    hour = "";
    minute = "";
    second = "";

    strcpy(dateFragment,date.c_str());
    strcpy(timeFragment,time.c_str());

    // cout << "here is good " << dateFragment << endl;
    // cout << "here is good " << timeFragment << endl;

    for(Int_t i = 0; i<n_dates; i++){
      // cout << format_date[i] << " ";
      if(format_date[i]=='d'){
        // cout << dateFragment[i] << endl;
	day.append(1,dateFragment[i]);
      }
      else if(format_date[i]=='m'){
	month.append(1,dateFragment[i]);
      }
      else if(format_date[i]=='y'){
	year.append(1,dateFragment[i]);
      }
    }

    for(Int_t i = 0; i<n_times; i++){
      if(format_time[i]=='h'){
	hour.append(1,timeFragment[i]);
      }
      else if(format_time[i]=='m'){
	minute.append(1,timeFragment[i]);
      }
      else if(format_time[i]=='s'){
	second.append(1,timeFragment[i]);
      }
    }

    // cout << year << " " << month << " " << day << " " << hour << " " << minute << " " << second << endl;

    if(month.size()>2){
      for(int i = 0; i < month.size(); i++){
        // cout << "..." << endl;
	month[i] = (char)(tolower(month[i]));
      }
      rewrite_month(month);
    }

    // cout << year << "-" << month << "-" << day << " " << hour << ":" << minute << ":" << second << endl;
    da = TTimeStamp(stoul(year),stoul(month),stoul(day),stoul(hour),stoul(minute),stoul(second));
    return da.AsDouble();
 
  }


  void rewrite_month(string &month){
    if(month=="jan"){month="01"; return;}
    if(month=="feb"){month="02"; return;}
    if(month=="mar"){month="03"; return;}
    if(month=="apr"){month="04"; return;}
    if(month=="may"){month="05"; return;}
    if(month=="jun"){month="06"; return;}
    if(month=="jul"){month="07"; return;}
    if(month=="aug"){month="08"; return;}
    if(month=="set"){month="09"; return;}
    if(month=="out"){month="10"; return;}
    if(month=="nov"){month="11"; return;}
    if(month=="dec"){month="12"; return;}
  }



};
