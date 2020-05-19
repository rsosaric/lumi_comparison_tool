
#include <iostream>
#include <fstream>
#include <sstream>

#include <time.h>

using namespace std;

long calcTimestamp (int inYear, int inMon, int inDay, int inHour, int inMin, int inSec) {
  
  struct tm *timeinfo;
  time_t rawtime;
  
  rawtime = time(0);
  timeinfo = localtime (&rawtime);
  
  timeinfo->tm_year = inYear - 1900;
  timeinfo->tm_mon  = inMon - 1;
  timeinfo->tm_mday = inDay;
  timeinfo->tm_hour = inHour;
  timeinfo->tm_min  = inMin;
  timeinfo->tm_sec  = inSec;
  
  return timegm(timeinfo);
}

int main(int argc, char **argv) {
  
    int Nfiles = 1;
    int Nvars = 8;
    int NdataMax = 100000;
  
    if (argc != (Nfiles + 1)) {
        cout << "ERROR: Expecting " << Nfiles << " Timber output " << ((Nfiles == 1) ? "file." : "files.") << endl;
        exit(1);
    }
  
    string filename[Nfiles];
    for (int i=0; i<Nfiles; i++) filename[i] = string(argv[i+1]);
  
    string stVarName[Nvars];
    int iNdata[Nvars];
    for (int i=0; i<Nvars; i++) iNdata[i] = -1;
    int iTimestamp[Nvars][NdataMax];
    float dbVarval[Nvars][NdataMax];
  
    ifstream fin[Nfiles];
    string stfiline;
    int j, iNdataArr;
    int iYear, iMonth, iDay, iHour, iMin, iSec;
    for (int i=0; i<Nfiles; i++) {
        fin[i].open(filename[i].c_str());
    
        j = 0;
        iNdataArr = -1;
        while(getline(fin[i], stfiline)) {
            if (stfiline.substr(0,8) == string("VARIABLE")) {
                iNdataArr++;
                stVarName[iNdataArr] = stfiline.substr(10,26);
//                cout << "Variable name: " << stVarName[iNdataArr] << endl;
            } else if (stfiline.substr(0,2) == string("20")) {
                iYear = 0; iMonth = 0; iDay = 0; iHour = 0; iMin = 0; iSec = 0;
                iNdata[iNdataArr]++;
          
                iYear  = stoi(stfiline.substr(0,4));
                iMonth = stoi(stfiline.substr(5,2));
                iDay   = stoi(stfiline.substr(8,2));
                iHour  = stoi(stfiline.substr(11,2));
                iMin   = stoi(stfiline.substr(14,2));
                iSec   = stoi(stfiline.substr(17,2));
                if (stoi(stfiline.substr(20,1)) > 4) iSec += 1;
                iTimestamp[iNdataArr][iNdata[iNdataArr]] = calcTimestamp(iYear,iMonth,iDay,iHour,iMin,iSec);

                dbVarval[iNdataArr][iNdata[iNdataArr]] = stod(stfiline.substr(24));
          
            }
            j++;
        } //end loop over lines in the datafile
        iNdataArr++;
    } //end loop over datafiles
  
  
  
    int iDataCurr[Nvars];
    for (int i=0; i<Nvars; i++) iDataCurr[i] = 0;
    bool blTimeSame, blBreak;

//  cout << "test 05: iNdataArr = " << iNdataArr << endl;
    int iNdataOut = 0;
    ofstream flOutput;
    flOutput.open( filename[0].substr(0,filename[0].length() - 4) + string("_merged.csv") );
    flOutput << "sec,";
    for (int i=0; i<(iNdataArr); i++) {
        flOutput << stVarName[i];
        if (i != (iNdataArr - 1)) {
            flOutput << ",";
        } else {
            flOutput << endl;
        }
    }
    
    
    int iT_prev = 0;
    int iT_curr = iTimestamp[0][iDataCurr[0]];
    blBreak = false;
    
    while(!blBreak) {
        // search for the next minimum available timestamp
        for (int i=0; i<(iNdataArr); i++) {
            if ( iDataCurr[i] < iNdata[i] && iTimestamp[i][iDataCurr[i] + 1] < iT_curr && iTimestamp[i][iDataCurr[i] + 1] > iT_prev ) {
                iT_curr = iTimestamp[i][iDataCurr[i] + 1];
            }
        }
        
        // serch for the closest (less equal) timestamp in each data array
        for (int i=0; i<(iNdataArr); i++) {
            while( iDataCurr[i] < iNdata[i] && iTimestamp[i][iDataCurr[i] + 1] <= iT_curr) iDataCurr[i]++;
        }
        
        // write the data
        iNdataOut++;
        flOutput << iT_curr << ",";
        for (int i=0; i<(iNdataArr); i++) {
/*
            if ( iTimestamp[i][iDataCurr[i]] == iT_curr ) {
                flOutput << dbVarval[i][iDataCurr[i]];
            } else {
                flOutput << " === ";
            }
*/
            flOutput << dbVarval[i][iDataCurr[i]];
            if (i != (iNdataArr - 1)) {
                flOutput << ",";
            } else {
                flOutput << endl;
            }
        }
        iT_prev = iT_curr;
        
        // check if there is any timestamp above current
        for (int i=0; i<(iNdataArr); i++) {
            if ( iDataCurr[i] < iNdata[i] && iTimestamp[i][iDataCurr[i] + 1] > iT_curr ) {
                iT_curr = iTimestamp[i][iDataCurr[i] + 1];
                break;
            }
        }
        if (iT_curr == iT_prev) blBreak = true;
    }
    
    /*
    if (iNdataArr > 0) {
        blBreak = false;
        while(true) {
            blTimeSame = true;
            for (int i=1; i<(iNdataArr); i++) {
                if ( iTimestamp[i][iDataCurr[i]] > iTimestamp[0][iDataCurr[0]] ) {
                    iDataCurr[0]++;
                    if (iDataCurr[0] > iNdata[0]) {
                        blBreak = true;
                        break;
                    }
                    blTimeSame = false;
                } else if ( iTimestamp[i][iDataCurr[i]] < iTimestamp[0][iDataCurr[0]] ) {
                    iDataCurr[i]++;
                    if (iDataCurr[i] > iNdata[i]) {
                        blBreak = true;
                        break;
                    }
                    blTimeSame = false;
                }
            }
      
            if (blBreak) break;
      
            if (blTimeSame) {
                iNdataOut++;
                flOutput << iTimestamp[0][iDataCurr[0]] << "," << dbVarval[0][iDataCurr[0]] << ",";
                for (int i=1; i<(iNdataArr); i++) {
                    flOutput << dbVarval[i][iDataCurr[i]];
                    if (i != (iNdataArr - 1)) {
                        flOutput << ",";
                    } else {
                        flOutput << endl;
                    }
                }
                iDataCurr[0]++;
            }
        }
    }
    */
    cout << endl << "Output data points: " << iNdataOut << endl;
  
    flOutput.close();
  
    return 0;
}

