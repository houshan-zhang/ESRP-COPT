#ifndef __VARIABLES_HEADER__
#define __VARIABLES_HEADER__

#include "../ThirdParty/copt/include/copt.h"

#include <stdlib.h>
#include <stdio.h>
#include <sstream>
#include <vector>
#include <math.h>
#include <string.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <random>
#include <algorithm>

using namespace std;

#define BENDERS_EPSILON             1e-3
#define EPSILON                     1e-5
#define PrintErr()                  printf("ERROR! [%s:%d]\n", __FILE__, __LINE__)

struct SampStruct;
struct SuppStruct;
struct CustStruct;
struct EdgeStruct;

/* Edge structure in graph */
struct EdgeStruct
{
   int id;
   int pos; /* Position in subproblem */
   CustStruct* cust; /* customer */
   SuppStruct* supp; /* supplier */
};

/* Customer structure in graph */
struct CustStruct
{
   int id;
   int pos; /* Position in subproblem */
   double x;
   double y;
   int degree;
   vector<EdgeStruct*> Edge;
   vector<double> demand;
   vector<double> cost;
};

/* Supplier structure in graph */
struct SuppStruct
{
   int id;
   int pos = -1; /* Position in master problem */
   double x;
   double y;
   int degree;
   double liveprob;
   vector<EdgeStruct*> Edge;
   int bandwidth;
   double working;
   double recover;
   vector<bool> alive;
   double status; /* status = alive[t] && (working || recover) */
};

/* Sample structure in graph */
struct SampStruct
{
   int id;
   int pos; /* Position in master problem */
   double* Cost;
   double* Demand;
   double SumDemand;
};

/* Problem InstStruct structure */
struct InstStruct
{
   /* Settings */
   int QoS; /* Quality of Service */
   int SolvingSetting; /* Solving settings (1 for (SNA), 2 for (SNA, SCNA) or 3 for (SNA,SCNA, and INA)) */
   int DisasterType;
   double TimeLimit; /* Time limit */
   //TODO
   double MemLimit; /* Memory limit */
   double MinGap;
   int Budget; /* Cardinality restriction */
   int ProduceSeed; /* Random seed of producing live-arc graphs (default: 1) */

   /* Graph structures */
   int NumSamp; /* Number of sampling */
   int NumCust; /* Nmber of custormer */
   int NumSupp; /* Number of supplier */
   int NumWork = 0; /* Number of working supplier */
   int NumReco; /* Number of supplier may be recovered (NumReco = NumSupp-NumWork) */
   int NumEdge = 0; /* Number of egde */
   vector<SampStruct*> Samp; /* Sample information */
   vector<CustStruct*> Cust; /* Custormer information */
   vector<SuppStruct*> Supp; /* Supplier information */
   vector<SuppStruct*> Reco; /* Damaged supplier information */
   vector<EdgeStruct*> Edge; /* Edge information */
   int* CustPos; /* COL Customers' position in MILP model of COPT */
   int* SuppPos; /* ROW Suppliers' position in MILP model of COPT */
   double Prob1 = 0.5;
   double Prob2 = -1.0;
   int NumIter = 0;

   /* Variables and parameters in solving */
   int NumBinCut = 0; /* Number of Benders cuts for binary relaxed solution */
   int NumFracCut = 0; /* Number of Benders cuts for fractional relaxed solution */
   int LastDepth = 0; /* Last branch-and-bound node depth */
   double LastRelaxObj = -1e+15; /* Last Lp value */
   double CutRhs = 0; /* Right hand side (constant term) of Benders cut */
   int NumThread = 1; /* Number of threads */
   int NumNode; /* Number of branch-and-bound nodes */
   double OptObj; /* Optimal objective value */
   double RelaxObj; /* Relaxed objective value */
   double* CutVal; /* Coefficients of Benders cut */
   int* CutInd; /* Variables indexes of Benders cut */
   vector<SuppStruct*> OneSupp;
   double* Pi; /* dual values */
   double VioVal; /* Violated value of Benders cut */
   double* RelaxSol; /* Relaxed solution in branch-and-bound tree */
   double BuildTime;

   /* COPT */
   copt_prob*  COPTLp; /* COPT Lp point for the Benders model */
   copt_env* COPTEnv; /* COPT environment point for the Benders model */
   copt_env* COPTEnvMaster; /* COPT environment point for the Benders model */
   copt_prob*  COPTLpMaster; /* COPT Lp point for the Benders model */
   copt_env* COPTEnvSub; /* COPT environment point for the Benders model */
   copt_prob* COPTLpSub; /* COPT Lp point for the Benders model */
   double COPTSolvingTime; /* COPT solving time */
   int COPTSeed = -1; /* COPT random seed */
   int COPTStatus; /* COPT status */
};


#endif
