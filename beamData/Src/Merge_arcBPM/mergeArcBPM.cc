
#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>

#include "TMath.h"

using namespace std;

  const int NinCol     =     4; //B1_L, B1_R, B2_L, B2_R
  const int NinSubCol  =     3; //T, H, V
  const int NdataFmax  =    50;
  const int NdataMax   = 30000;
  const double scaleT  = 0.001;

double sciToDub(const string& str) {

   stringstream ss(str);
   double d = 0;
   ss >> d;

   if (ss.fail()) {
      string s = "Unable to format ";
      s += str;
      s += " as a number!";
      throw (s);
   }

   return (d);
}

long sciToLong(const string& str) {
   return (sciToDub(str));
}
  
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

void readDataSteer (string stFile, string stInColName[NinCol], string stInSubColName[NinSubCol], int iInColId[NinCol][NinSubCol], int iNdataF[NinCol], string stDataF[NinCol][NdataFmax], string& stOutFname) {
  
  const int Npar           = 3;
  const int NsubParsMax    = NinCol;
  const int NsubSubParsMax = NinSubCol;
  
  for (int j=0; j<NsubParsMax; j++) {
    for (int k=0; k<NsubSubParsMax; k++) {
      iInColId[j][k] = -1;
    }
  }
  
  for (int j=0; j<NsubParsMax; j++) {
    iNdataF[j] = 0;
  }
  
  string stPars[Npar]                = {"Input Columns:", "Data files:", "Output file:"};
  int NsubPars[Npar]                 = {       4,               4,              1};
  int NsubSubPars[Npar][NsubParsMax] = {{3, 3, 3, 3},     {0, 0, 0, 0},  {0, 0, 0, 0}};
  
  string stOFsubPars[1]              = {"Name:"};
  
  int iCurrPar, iCurrSubPar;
  bool blSubParRead[NsubParsMax], blSubSubParRead[NsubSubParsMax], blTmp;
  
  ifstream fin;
  string stfiline;
  int iSubstr;
  
  fin.open(stFile.c_str());
  
  iCurrPar = -1; iCurrSubPar = -1;
  while(getline(fin, stfiline)) {
    
    if (iCurrPar >= 0 && iCurrPar <= 1) {
      for (int j=0; j<NsubPars[iCurrPar]; j++) {
        iSubstr = stfiline.find(stInColName[j]);
        if (iSubstr < stfiline.size()) {
          if (blSubParRead[j]) {
            cout << "WARNING: Found subparameter '" << stInColName[j] << "' in parameter '" << stPars[iCurrPar] << "' while it is already marked as 'read'." << endl;  
          }
          if (iCurrSubPar != -1) {
            blTmp = blSubSubParRead[0];
            for (int k=1; k<NsubSubPars[iCurrPar][iCurrSubPar]; k++) blTmp = blTmp && blSubSubParRead[k];
            if (!blTmp) {
              cout << "WARNING: New subparameter '" << stInColName[j] << "' found while not all subsubpars for '" << stInColName[iCurrSubPar] << "' in parameter '" << stPars[iCurrPar] << "' have been read." << endl;
            }
          }
          iCurrSubPar = j;
          for (int k=0; k<NsubSubParsMax; k++) blSubSubParRead[k] = false;
        }
      } //end of loop over SubPars
    } else if (iCurrPar == 2) {
      iCurrSubPar = 0;
    }
    for (int i=0; i<Npar; i++) {
      iSubstr = stfiline.find(stPars[i]);
      if (iSubstr < stfiline.size()) {
        if (iCurrPar != -1) {
          blTmp = blSubParRead[0];
          for (int j=1; j<NsubPars[iCurrPar]; j++) blTmp = blTmp && blSubParRead[j];
          if (!blTmp) {
            cout << "WARNING: New parameter '" << stPars[i] << "' found while not all subpars for '" << stPars[iCurrPar] << "' have been read." << endl;
          }
        }
        iCurrPar = i;
        iCurrSubPar = -1;
        for (int j=0; j<NsubParsMax; j++) blSubParRead[j] = false;
        for (int k=0; k<NsubSubParsMax; k++) blSubSubParRead[k] = false;
      }
    } //end of loop over Pars
    
    if (iCurrPar == 0) {
      if (iCurrSubPar >= 0 && iCurrSubPar <= 3) {
        for (int k=0; k<NsubSubPars[iCurrPar][iCurrSubPar]; k++) {
          iSubstr = stfiline.find(stInSubColName[k]);
          if (iSubstr < stfiline.size()) {
            iInColId[iCurrSubPar][k] = stoi(stfiline.substr(iSubstr + stInSubColName[k].size() + 0));
            blSubSubParRead[k] = true;
          }
        }
      }
    } else if (iCurrPar == 1) {
      if (iCurrSubPar >= 0 && iCurrSubPar <= 3) {
        if (stfiline.find(stInColName[iCurrSubPar]) >= stfiline.size()) { //start reading only after subpar name line
          if (stfiline.find_first_not_of(' ') < stfiline.size()) { //don't read empty lines
            stDataF[iCurrSubPar][iNdataF[iCurrSubPar]] = stfiline.substr( stfiline.find_first_not_of(' ') );
            iNdataF[iCurrSubPar]++;
            blSubSubParRead[0] = true;
          }
        }
      }
    } else if (iCurrPar == 2) {
      for (int k=0; k<NsubPars[iCurrPar]; k++) {
        iSubstr = stfiline.find(stOFsubPars[k]);
        if (iSubstr < stfiline.size()) {
          stOutFname = stfiline.substr(iSubstr + stOFsubPars[k].size() + stfiline.substr(iSubstr + stOFsubPars[k].size()).find_first_not_of(' ') );
          blSubSubParRead[k] = true;
        }
      }
    }
    
    blSubParRead[iCurrSubPar] = blSubSubParRead[0];
    for (int k=1; k<NsubSubPars[iCurrPar][iCurrSubPar]; k++) blSubParRead[iCurrSubPar] = blSubParRead[iCurrSubPar] && blSubSubParRead[k];
    
  } //end of loop over steering file lines
  
  fin.close();
  
  return;
}

void readDataFile (string stFile, string stInSubColName[NinSubCol], int iInColId[NinSubCol], int& NreadPar, long lreadT[NdataMax], Double_t dbreadVar[NinSubCol-1][NdataMax]) {
  
    ifstream fin;
    string stfiline;
    int iLine;
  
    long lT_tmp;
    Double_t dbVar_tmp[NinSubCol-1];
    int id_tmp;
    bool bl_tmp;
  
    fin.open(stFile.c_str());
    int iSubstr, iSubstr_next, iCurrCol, iCurrId;
  
    int iDFcSign[NinSubCol];
    for (int i=0; i<NinSubCol; i++) {
        if (iInColId[i] < 0) {
            iDFcSign[i] = -1;
        } else {
            iDFcSign[i] = 1;
        }
    }
  
    iLine = -2;
    while(getline(fin, stfiline)) {
        if (iLine >= 0) {
            iCurrCol = 1;
            iSubstr = 0;
            iSubstr_next = 0;
            iCurrId = 0;
            while (iSubstr_next > -1) {
                iSubstr_next = stfiline.find(string(","), iSubstr);
                for (int i=0; i<NinSubCol; i++) {
                    if ((iInColId[i]*iDFcSign[i]) == iCurrCol) {
                        if (stInSubColName[i] == string("T:")) {
//                            lreadT[NreadPar] = sciToLong(stfiline.substr(iSubstr,iSubstr_next-iSubstr));
                            lT_tmp = sciToLong(stfiline.substr(iSubstr,iSubstr_next-iSubstr));
                        } else {
//                            dbreadVar[iCurrId][NreadPar] = stod(stfiline.substr(iSubstr,iSubstr_next-iSubstr)) * double(iDFcSign[i]);
                            dbVar_tmp[iCurrId] = stod(stfiline.substr(iSubstr,iSubstr_next-iSubstr)) * double(iDFcSign[i]);
                            iCurrId++;
                        }
                    }
                }
                iSubstr = iSubstr_next+1;
                iCurrCol++;
            } //end of one string reading
            
            id_tmp = -1;
            for (int j=0; j<NreadPar; j++) {
                if (lreadT[j] == lT_tmp) id_tmp = j;
            }
            if (id_tmp < 0) {
                lreadT[NreadPar] = lT_tmp;
                for (int i=0; i<(NinSubCol-1); i++) {
                    dbreadVar[i][NreadPar] = dbVar_tmp[i];
                }
                NreadPar++;
            } else {
                bl_tmp = true;
                for (int i=0; i<(NinSubCol-1); i++) bl_tmp = bl_tmp && ( (dbreadVar[i][id_tmp]/dbVar_tmp[i] - 1)*(dbreadVar[i][id_tmp]/dbVar_tmp[i] - 1) < 0.0000000001 );
                if (!bl_tmp) {
                    cout << "WARNING: In '" << stFile << "' second line with already registered timestamp and different values:" << endl;
                    cout << lreadT[id_tmp];  for (int i=0; i<(NinSubCol-1); i++) cout << "   " << dbreadVar[i][id_tmp];  cout << endl;
                    cout << lT_tmp;          for (int i=0; i<(NinSubCol-1); i++) cout << "   " << dbVar_tmp[i];          cout << endl << endl;
                }
            }
        }
        iLine++;
    } //end of file reading
  
  fin.close();
  
  return;
}

void sortData (int Npar, long lreadT[NdataMax], Double_t dbreadVar[NinSubCol-1][NdataMax]) {
  
  long lTmp;
  Double_t dbTmp;
  bool blIsSwap = true;
  
  while(blIsSwap) {
    blIsSwap = false;
    
    for (int i=0; i<(Npar-1); i++) {
      if (lreadT[i] > lreadT[i+1]) {
        lTmp        = lreadT[i];
        lreadT[i]   = lreadT[i+1];
        lreadT[i+1] = lTmp;
        for (int j=0; j<(NinSubCol-1); j++) {
          dbTmp             = dbreadVar[j][i];
          dbreadVar[j][i]   = dbreadVar[j][i+1];
          dbreadVar[j][i+1] = dbTmp;
        }
        blIsSwap = true;
      }
    }
  }
  
  return;  
}

void mergeData (int NparIn[NinCol], long lTIn[NinCol][NdataMax], Double_t dbVarIn[NinCol][NinSubCol-1][NdataMax], int iOutVarId[NinCol*(NinSubCol-1)][2], int& iOutData, long lOutDataVarT[NdataMax], Double_t dbOutDataVar[NinCol*(NinSubCol-1)][NdataMax]) {
  
  int iIter[NinCol];
  for (int i=0; i<NinCol; i++) iIter[i] = 0;
  bool blIsClose;
  
  iOutData = 0;
  
  for (iIter[0]=0; iIter[0]<NparIn[0]; iIter[0]++) {
    for (int j=1; j<NinCol; j++) {
      while ( long(lTIn[0][iIter[0]]*scaleT) > long(lTIn[j][iIter[j]]*scaleT) ) {
        if (iIter[j] < (NparIn[j]-1)) { 
          iIter[j]++;
        } else {
          break;
        }
      }
    }
    
    blIsClose = true;
    for (int j=1; j<NinCol; j++) blIsClose = blIsClose && (long(lTIn[0][iIter[0]]*scaleT) == long(lTIn[j][iIter[j]]*scaleT));
    
    if (blIsClose) {
      lOutDataVarT[iOutData] = long(double(lTIn[0][iIter[0]])*scaleT + 0.5);
      for (int k=0; k<(NinCol*(NinSubCol-1)); k++) {
        dbOutDataVar[k][iOutData] = dbVarIn[iOutVarId[k][0]][iOutVarId[k][1]][iIter[iOutVarId[k][0]]];
      }
      iOutData++;
    }
  }
  
  return;
}


int main(int argc, char **argv) {
  
  int Nfiles = 1;
  int Nvars = 8;
  
  int iDebug = 0;
  
  if (argc < Nfiles + 1) {
    printf("program usage:\n mergeArcBPM [arc BPM steering file]\n");
    return -1;
  }
    
  string filename[Nfiles];
  for (int i=0; i<Nfiles; i++) filename[i] = string(argv[i+1]);

  string stInSubParName[NinCol] = {"B1_L:", "B1_R:", "B2_L:", "B2_R:"};
  string stInSubSubParName[NinSubCol] = {"T:", "H:", "V:"};
  int iInParColId[NinCol][NinSubCol];
  
  string stOutParT = string("Timestamp");
  string stOutParName[NinCol*(NinSubCol-1)]  = {"B1_H_L", "B1_H_R", "B1_V_L", "B1_V_R", "B2_H_L", "B2_H_R", "B2_V_L", "B2_V_R"};
  int iOutParId[NinCol*(NinSubCol-1) + 1][2] = {{0, 0},   {1, 0},   {0, 1},   {1, 1},   {2, 0},   {3, 0},   {2, 1},   {3, 1}}; //{SubPar [0..NinCol], SubSubPar [0..NinSubCol]}
  int NOutVar;
  long lOutVarT[NdataMax];
  Double_t dbOutVar[NinCol*(NinSubCol-1)][NdataMax];
  
  int iNdataFiles[NinCol];
  string stDataFiles[NinCol][NdataFmax];
  string stOutFileName;
  
  readDataSteer(filename[0], stInSubParName, stInSubSubParName, iInParColId, iNdataFiles, stDataFiles, stOutFileName);

  if (iDebug == 1) {
    cout << "      T  H  V" <<endl;
    for (int j=0; j<NinCol; j++) {
      cout << stInSubParName[j] << " ";
      for (int k=0; k<NinSubCol; k++) {
        cout << iInParColId[j][k] << "  ";
      }
      cout << endl;
    }
    
    cout << endl;
    for (int j=0; j<NinCol; j++) {
      cout << stInSubParName[j] << endl;
      for (int k=0; k<iNdataFiles[j]; k++) {
        cout << stDataFiles[j][k] << endl;
      }
    }
  }
  
  int iNreadVar[NinCol];
  long lreadVarT[NinCol][NdataMax];
  Double_t dbreadVar[NinCol][NinSubCol-1][NdataMax];
  
  for (int i=0; i<NinCol; i++) {
    iNreadVar[i] = 0;
    for (int j=0; j<iNdataFiles[i]; j++) {
      readDataFile (stDataFiles[i][j], stInSubSubParName, iInParColId[i], iNreadVar[i], lreadVarT[i], dbreadVar[i]);
    }
    sortData(iNreadVar[i], lreadVarT[i], dbreadVar[i]);
  }
  
  if (iDebug == 2) {
    for (int i=0; i<NinCol; i++) {
      cout << stInSubParName[i] << endl;
      for (int k=0; k<iNreadVar[i]; k++) {
        cout << lreadVarT[i][k] << "   ";
        for (int j=0; j<(NinSubCol-1); j++) { 
          cout << dbreadVar[i][j][k] << "   ";
        }
        cout << endl;
      }
    }
  }
  
  mergeData (iNreadVar, lreadVarT, dbreadVar, iOutParId, NOutVar, lOutVarT, dbOutVar);
  
  if (iDebug == 3) {
    cout << stOutParT;
    for (int i=0; i<(NinCol*(NinSubCol-1)); i++) cout << "," << stOutParName[i];
    cout << endl;
    for (int j=0; j<NOutVar; j++) {
      cout << lOutVarT[j];
      for (int k=0; k<(NinCol*(NinSubCol-1)); k++) cout << "," << dbOutVar[k][j];
      cout << endl;
    }
  }
  
  ofstream flOutput;
  flOutput.open( stOutFileName );
  
  flOutput << stOutParT;
  for (int i=0; i<(NinCol*(NinSubCol-1)); i++) flOutput << "," << stOutParName[i];
  flOutput << endl;
  for (int j=0; j<NOutVar; j++) {
    flOutput << lOutVarT[j];
    for (int k=0; k<(NinCol*(NinSubCol-1)); k++) flOutput << "," << dbOutVar[k][j];
    flOutput << endl;
  }
  
  flOutput.close();
  
  return 0;
}

