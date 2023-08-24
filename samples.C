// ________________________________________ //
// Author: Henrique Souza
// Filename: samples.C
// Created: 2021
// ________________________________________ //
#include "MYCODES.h"



class SAMPLE: public ANALYZER{
  public:
    Int_t n_points = memorydepth;
    string file = "analyzed.root";
    string tree = "t1";






    SAMPLE(string m_myname = "s", string m_filename = "analyzed.root"){
      myname = m_myname;
      setAnalyzer(m_filename);

      cout << "@@@@@@@@@@@@@@@ ____________________________________ @@@@@@@@@@@@@@@" << endl;
      cout << "@@@@@@@@@@@@@@@ This class is not being used anymore @@@@@@@@@@@@@@@" << endl;
      cout << "@@@@@@@@@@@@@@@ ____________________________________ @@@@@@@@@@@@@@@" << endl;

    }
};
