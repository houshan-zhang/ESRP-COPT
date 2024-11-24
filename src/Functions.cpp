#include "Functions.h"
#include "Variables.h"

bool IsEq( double a, double b )
{
   return abs(a-b) <= EPSILON;
}

void GenerateGraph(InstStruct* Inst)
{
   srand(Inst->ProduceSeed);
   Inst->Cust = vector<CustStruct*>(Inst->NumCust);
   Inst->CustPos = new int[Inst->NumCust];
   Inst->SuppPos = new int[Inst->NumSupp];
   Inst->Supp = vector<SuppStruct*>(Inst->NumSupp);
   int counter = 0;
   /* Generate customer */
   for( int i = 0; i < Inst->NumCust; i++ )
   {
      Inst->Cust[i] = new CustStruct;
      CustStruct* cust = Inst->Cust[i];
      cust->id = i;
      cust->pos = i;
      Inst->CustPos[i] = i;
      cust->x = 100*(rand()/double(RAND_MAX));
      cust->y = 100*(rand()/double(RAND_MAX));
   }
   /* Generate supplier */
   for( int j = 0; j < Inst->NumSupp; j++ )
   {
      Inst->Supp[j] = new SuppStruct;
      SuppStruct* supp = Inst->Supp[j];
      supp->id = j;
      Inst->SuppPos[j] = Inst->NumCust + j;
      supp->x = 100*(rand()/double(RAND_MAX));
      supp->y = 100*(rand()/double(RAND_MAX));
      supp->bandwidth = rand()%40000+10000;
      double tmp = rand()/double(RAND_MAX);
      /* ---------- */
      supp->working = 1.0;
      if( Inst->DisasterType == 1 && sqrt(pow(supp->x-50, 2) + pow(supp->y-50, 2 )) <= sqrt(5000.0/3.1415926)) { supp->working = 0.0; }
      if( Inst->DisasterType == 2 )
      {
         if( sqrt(pow(supp->x-25, 2) + pow(supp->y-25, 2 )) <= sqrt(1250.0/3.1415926) )
            supp->working = 0.0;
         else if( sqrt(pow(supp->x-75, 2) + pow(supp->y-25, 2 )) <= sqrt(1250.0/3.1415926) )
            supp->working = 0.0;
         else if( sqrt(pow(supp->x-25, 2) + pow(supp->y-75, 2 )) <= sqrt(1250.0/3.1415926) )
            supp->working = 0.0;
         else if( sqrt(pow(supp->x-75, 2) + pow(supp->y-75, 2 )) <= sqrt(1250.0/3.1415926) )
            supp->working = 0.0;
      }
      if( Inst->DisasterType == 3 )
      {
         if( supp->x >= 5 && supp->x <= 10 ) { supp->working = 0.0; }
         if( supp->x >= 15 && supp->x <= 20 ) { supp->working = 0.0; }
         if( supp->x >= 25 && supp->x <= 30 ) { supp->working = 0.0; }
         if( supp->x >= 35 && supp->x <= 40 ) { supp->working = 0.0; }
         if( supp->x >= 45 && supp->x <= 50 ) { supp->working = 0.0; }
         if( supp->x >= 55 && supp->x <= 60 ) { supp->working = 0.0; }
         if( supp->x >= 65 && supp->x <= 70 ) { supp->working = 0.0; }
         if( supp->x >= 75 && supp->x <= 80 ) { supp->working = 0.0; }
         if( supp->x >= 85 && supp->x <= 90 ) { supp->working = 0.0; }
         if( supp->x >= 96.25 && supp->x <= 98.75 ) { supp->working = 0.0; }
      }
      if( Inst->DisasterType == 4 ) { supp->working = tmp < Inst->Prob1 ? 1.0 : 0.0; }
      /* ---------- */
      if( IsEq(supp->working, 1.0) )
         Inst->NumWork++;
      else
      {
         supp->pos = counter;
         Inst->Reco.push_back(supp);
         counter++;
      }
      supp->recover = 0.0;
   }
   Inst->NumReco = Inst->NumSupp - Inst->NumWork;
   /* Generate edge */
   for( int i = 0; i < Inst->NumCust; i++ )
   {
      CustStruct* cust = Inst->Cust[i];
      for( int j = 0; j < Inst->NumSupp; j++ )
      {
         SuppStruct* supp = Inst->Supp[j];
         if( sqrt(pow(cust->x-supp->x, 2) + pow(cust->y-supp->y, 2)) <= Inst->QoS )
         {
            EdgeStruct* edge = new EdgeStruct;
            Inst->Edge.push_back(edge);
            cust->Edge.push_back(edge);
            supp->Edge.push_back(edge);
            edge->id = Inst->NumEdge;
            edge->pos = Inst->NumCust + Inst->NumEdge;
            edge->cust = cust;
            edge->supp = supp;
            Inst->NumEdge++;
         }
      }
   }
   /* Write degree */
   for( int i = 0; i < Inst->NumCust; i++ )
      Inst->Cust[i]->degree = Inst->Cust[i]->Edge.size();
   for( int j = 0; j < Inst->NumSupp; j++ )
      Inst->Supp[j]->degree = Inst->Supp[j]->Edge.size();
}

void Sample(InstStruct* Inst)
{
   /* Generate sample */
   Inst->Samp = vector<SampStruct*>(Inst->NumSamp);
   for( int t = 0; t < Inst->NumSamp; t++ )
   {
      SampStruct* samp = new SampStruct;
      Inst->Samp[t] = samp;
      samp->id = t;
      samp->pos = Inst->NumReco + t;
      samp->Cost = new double[Inst->NumCust];
      samp->Demand = new double[Inst->NumCust];
      samp->SumDemand = 0;
   }
   for( int i = 0; i < Inst->NumCust; i++ )
   {
      CustStruct* cust = Inst->Cust[i];
      cust->demand = vector<double>(Inst->NumSamp);
      cust->cost = vector<double>(Inst->NumSamp);
      for( int t = 0; t < Inst->NumSamp; t++ )
      {
         cust->demand[t] = rand()%3000;
         cust->cost[t] = 1.0/Inst->NumSamp;
         Inst->Samp[t]->Cost[i] = cust->cost[t];
         Inst->Samp[t]->Demand[i] = cust->demand[t];
         Inst->Samp[t]->SumDemand += cust->demand[t];
      }
   }
   for( int j = 0; j < Inst->NumSupp; j++ )
   {
      SuppStruct* supp = Inst->Supp[j];
      supp->alive = vector<bool>(Inst->NumSamp);
      if( IsEq( Inst->Prob2, -1.0 ) )
         supp->liveprob = 0.5 + 0.5*(rand()/double(RAND_MAX));
      else
         supp->liveprob = Inst->Prob2;
      for( int t = 0; t < Inst->NumSamp; t++ )
      {
         double tmp = rand()/double(RAND_MAX);
         supp->alive[t] = tmp < supp->liveprob ? true : false;
      }
   }
}

void SetSupplier(InstStruct* Inst, vector<SuppStruct*> Supp, double sign)
{
   for( unsigned int j = 0; j < Supp.size(); j++ )
   {
      SuppStruct* supp = Supp[j];
      if( !IsEq(supp->working, 0.0) ) { PrintErr(); }
      supp->recover = sign;
   }
}

void SetSupplierStatus(InstStruct* Inst, vector<SuppStruct*> Supp, int t)
{
   int row_pos;
   for( unsigned int j = 0; j < Supp.size(); j++ )
   {
      SuppStruct* supp = Supp[j];
      if( supp->alive[t] )
         supp->status = IsEq(supp->working, 1.0) ? supp->working : supp->recover;
      else
         supp->status = 0.0;
      double rhs = supp->bandwidth * supp->status;
      row_pos = Inst->NumCust + supp->id;
      Inst->COPTStatus = COPT_SetRowUpper(Inst->COPTLpSub, 1, &row_pos, &rhs);
      if( Inst->COPTStatus ) { PrintErr(); }
   }
}

void ChgCOPTCoef(InstStruct* Inst, int t)
{
   Inst->COPTStatus = COPT_SetColObj(Inst->COPTLpSub, Inst->NumCust, Inst->CustPos, Inst->Samp[t]->Cost);
   if( Inst->COPTStatus ) { PrintErr(); }
   Inst->COPTStatus = COPT_SetRowLower(Inst->COPTLpSub, Inst->NumCust, Inst->CustPos, Inst->Samp[t]->Demand);
   Inst->COPTStatus = COPT_SetRowUpper(Inst->COPTLpSub, Inst->NumCust, Inst->CustPos, Inst->Samp[t]->Demand);
   if( Inst->COPTStatus ) { PrintErr(); }
}

/* Build the original problem formulation */
void BuildCOPTModel(InstStruct* Inst)
{
   /* Open Env and Lp */
   Inst->COPTStatus = COPT_CreateEnv(&Inst->COPTEnv);
   if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }

   Inst->COPTStatus = COPT_CreateProb(Inst->COPTEnv, &Inst->COPTLp);
   if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }

   /* Create and add variables to the  problem */
   int NumVar = Inst->NumReco + Inst->NumSamp * Inst->NumCust + Inst->NumSamp * Inst->NumEdge;
   double* Obj = new double[NumVar]; /* Objective Coefficients */
   double* Lb = new double[NumVar]; /* Variables low bounds */
   double* Ub = new double[NumVar]; /* Variables upper bounds */
   char* TypeVar = new char[NumVar]; /* Variables types */
   char** NameVar = new char*[NumVar]; /* Variables names */

   for( int i = 0; i < NumVar; i++ )
      NameVar[i] = new char[100];

   /* ( Col - 1 ) Create supplier variables */
   int counter = 0;
   for( int j = 0; j < Inst->NumSupp; j++ )
   {
      if( IsEq(Inst->Supp[j]->working, 0.0) )
      {
         Obj[counter] = 0.0;
         Lb[counter] = 0.0;
         Ub[counter] = 1.0;
         TypeVar[counter] = COPT_BINARY;
         sprintf(NameVar[counter], "supplier-%d", j);
         counter++;
      }
   }
   if( counter != Inst->NumReco ) { PrintErr(); }

   /* ( Col - 2 ) Create customer variables */
   for( int t = 0; t < Inst->NumSamp; t++ )
   {
      for( int j = 0; j < Inst->NumCust; j++ )
      {
         Obj[counter] = Inst->Samp[t]->Cost[j];
         Lb[counter] = 0.0;
         Ub[counter] = 1e+15;
         TypeVar[counter] = COPT_CONTINUOUS;
         sprintf(NameVar[counter], "loss-%d-%d", t, j);
         counter++;
      }
   }
   if( counter != Inst->NumReco + Inst->NumSamp * Inst->NumCust ) { PrintErr(); }

   /* ( Col - 3 ) Create flow variables */
   for( int t = 0; t < Inst->NumSamp; t++ )
   {
      for( int i = 0; i < Inst->NumCust; i++ )
      {
         CustStruct* cust = Inst->Cust[i];
         for( int j = 0; j < cust->degree; j++ )
         {
            Obj[counter] = 0.0;
            Lb[counter] = 0.0;
            Ub[counter] = 1e+15;
            TypeVar[counter] = COPT_CONTINUOUS;
            sprintf(NameVar[counter], "flow-%d-%d-%d", t, i, cust->Edge[j]->supp->id);
            counter++;
         }
      }
   }
   if( counter != NumVar ) { PrintErr(); }

   /* Add variables information to COPT */
   Inst->COPTStatus = COPT_AddCols(Inst->COPTLp, NumVar, Obj, NULL, NULL, NULL, NULL, TypeVar, Lb, Ub, NameVar);
   if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }

   delete []Obj;
   delete []Lb;
   delete []Ub;
   delete []TypeVar;

   for( int i = 0; i < NumVar; i++ )
   {
      delete []NameVar[i];
   }
   delete []NameVar;

   COPT_SetObjSense(Inst->COPTLp, COPT_MINIMIZE); /* Set optimization direction to be "MIN" */

   /* Add the constraints of the  problem */
   int* ConsInd = new int[max(Inst->NumEdge, Inst->NumReco) + 1];
   double* ConsVal = new double[max(Inst->NumEdge, Inst->NumReco) + 1];
   char sense;
   /* ( Row - 1 ) Creat the Cardinality (budget) constraint */
   sense = COPT_LESS_EQUAL; /* Set the constraint "<=" */
   counter = 0;
   for( int j = 0; j < Inst->NumSupp; j++ )
   {
      if( IsEq(Inst->Supp[j]->working, 0.0) )
      {
         ConsVal[counter] = 1;
         ConsInd[counter] = counter;
         counter++;
      }
   }
   double rhs = double(Inst->Budget);
   Inst->COPTStatus = COPT_AddRow(Inst->COPTLp, Inst->NumReco, ConsInd, ConsVal, sense, rhs, 0.0, NULL);
   if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }

   /* ( Row - 2 ) Creat the demand constraints */
   sense = COPT_EQUAL; /* Set the constraint "==" */
   for( int t = 0; t < Inst->NumSamp; t++ )
   {
      for( int i = 0; i < Inst->NumCust; i++ )
      {
         CustStruct* cust = Inst->Cust[i];
         for( int j = 0; j < cust->degree; j++ )
         {
            EdgeStruct* edge = cust->Edge[j];
            ConsVal[j] = 1;
            ConsInd[j] = Inst->NumReco + Inst->NumSamp*Inst->NumCust + t*Inst->NumEdge + edge->id;
         }
         ConsVal[cust->degree] = 1;
         ConsInd[cust->degree] = Inst->NumReco + t*Inst->NumCust + i;
         rhs = Inst->Samp[t]->Demand[i];
         Inst->COPTStatus = COPT_AddRow(Inst->COPTLp, cust->degree+1, ConsInd, ConsVal, sense, rhs, 0.0, NULL);
         if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }
      }
   }
   /* ( Row - 3 ) Creat the Capacity constraint */
   sense = COPT_LESS_EQUAL; /* Set the constraint "<=" */
   for( int t = 0; t < Inst->NumSamp; t++ )
   {
      for( int j = 0; j < Inst->NumSupp; j++ )
      {
         SuppStruct* supp = Inst->Supp[j];
         for( int i = 0; i < supp->degree; i++ )
         {
            EdgeStruct* edge = supp->Edge[i];
            ConsVal[i] = 1;
            ConsInd[i] = Inst->NumReco + Inst->NumSamp*Inst->NumCust + t*Inst->NumEdge + edge->id;
         }
         if( IsEq(supp->working, 1.0) && supp->alive[t] )
         {
            rhs = supp->bandwidth;
            Inst->COPTStatus = COPT_AddRow(Inst->COPTLp, supp->degree, ConsInd, ConsVal, sense, rhs, 0.0, NULL);
            if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }
         }
         else if( IsEq(supp->working, 0.0) && supp->alive[t] )
         {
            rhs = 0.0;
            ConsVal[supp->degree] = -supp->bandwidth;
            ConsInd[supp->degree] = supp->pos;
            Inst->COPTStatus = COPT_AddRow(Inst->COPTLp, supp->degree+1, ConsInd, ConsVal, sense, rhs, 0.0, NULL);
            if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }
         }
         else
         {
            rhs = 0.0;
            Inst->COPTStatus = COPT_AddRow(Inst->COPTLp, supp->degree, ConsInd, ConsVal, sense, rhs, 0.0, NULL);
            if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }
         }
      }
   }

   delete []ConsVal;
   delete []ConsInd;
}

void ChangeCOPTModel(InstStruct* Inst)
{
   vector<SuppStruct*> prior(Inst->NumReco);
   for( int i = 0; i < Inst->NumReco; i++ )
      prior[i] = Inst->Reco[i];
   if( Inst->SolvingSetting == -1 )
   {
      random_device rd;
      mt19937 g(rd());
      shuffle(prior.begin(), prior.end(), g);
   }
   else if( Inst->SolvingSetting == -2 )
      sort(prior.begin(), prior.end(), [](SuppStruct* a, SuppStruct* b) { return a->bandwidth > b->bandwidth; });
   else if( Inst->SolvingSetting == -3 )
      sort(prior.begin(), prior.end(), [](SuppStruct* a, SuppStruct* b) { return a->degree > b->degree; });
   else if( Inst->SolvingSetting == -4 )
      sort(prior.begin(), prior.end(), [](SuppStruct* a, SuppStruct* b) { return a->liveprob > b->liveprob; });
   double value = 1.0;
   for( int i = 0; i < min(Inst->Budget, Inst->NumReco); i++ )
   {
      Inst->COPTStatus = COPT_SetColLower(Inst->COPTLp, 1, &prior[i]->pos, &value);
      Inst->COPTStatus = COPT_SetColUpper(Inst->COPTLp, 1, &prior[i]->pos, &value);
   }
}

void SolveCOPTModel(InstStruct* Inst)
{
   /* Output information to the screen */
   Inst->COPTStatus = COPT_SetIntParam(Inst->COPTLp, COPT_INTPARAM_LOGTOCONSOLE, 1);

   /* Set the number of Threads in the calculation */
   Inst->COPTStatus = COPT_SetIntParam(Inst->COPTLp, COPT_INTPARAM_THREADS, Inst->NumThread);
   if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }

   /* Set the timelimit of the COPT (i.e., total timelimit minus the time cost in presolving) */
   cout << "COPT Solving Time Limit: " << Inst->TimeLimit << endl;
   Inst->COPTStatus = COPT_SetDblParam(Inst->COPTLp, COPT_DBLPARAM_TIMELIMIT, Inst->TimeLimit);
   if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }

   /* Set the absolute tolerance of COPT to 0.0 */
   Inst->COPTStatus = COPT_SetDblParam(Inst->COPTLp, COPT_DBLPARAM_ABSGAP, 0.0);
   if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }

   /* Set the relative tolerance of COPT to 0.0 */
   Inst->COPTStatus = COPT_SetDblParam(Inst->COPTLp, COPT_DBLPARAM_RELGAP, 0.0);

   /* COPT optimization process */
   clock_t TimeStart=clock();
   Inst->COPTStatus = COPT_Solve(Inst->COPTLp);
   if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }
   clock_t TimeEnd=clock();
   Inst->COPTSolvingTime = (double)(TimeEnd - TimeStart)/(double)CLOCKS_PER_SEC;

   /* Get the best feasible objective value within timelimit*/
   Inst->COPTStatus = COPT_GetDblAttr(Inst->COPTLp, COPT_DBLATTR_BESTOBJ, &Inst->OptObj);
   if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }

   /* Get the best relaxed objective value within timelimit*/
   Inst->COPTStatus = COPT_GetDblAttr(Inst->COPTLp, COPT_DBLATTR_BESTBND, &Inst->RelaxObj);
   if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }

   Inst->COPTStatus = COPT_GetIntAttr(Inst->COPTLp, COPT_INTATTR_NODECNT, &Inst->NumNode); /* Get the number of branch-and-bound nodes */

   /* Output the statistic information */
   cout<<endl<<endl;
   cout << "Results statistics:"<<endl;
   cout << "\t(a) Optimal objective value:\t\t\t\t" << Inst->OptObj <<endl;
   cout << "\t(b) Relaxed objective value:\t\t\t\t" << Inst->RelaxObj <<endl;
   cout << "\t(c) Number of branch-and-bound nodes:\t\t\t" << Inst->NumNode <<endl;
   cout << "\t(d) COPT solving time:\t\t\t\t\t" << Inst->COPTSolvingTime << endl;
   cout << "\t(e) Total time:\t\t\t\t\t\t" << Inst->COPTSolvingTime << endl;
   cout << "\t(f) Number of Benders cuts for binary solution:\t\t" << 0 << endl;
   cout << "\t(g) Number of Benders cuts for fractional solution:\t"<< 0 << endl;
   cout << "\t(h) Number of Benders iterations:\t\t\t"<< 0 << endl;
}


/* read parameters from a string */
void ReadParam(InstStruct* Inst, char* input)
{
   char* name;
   char* value;
   char* setfile;
   char* token;
   setfile = NULL;
   token = strtok(input, " ");
   while( token != NULL )
   {
      char *equalsign = strchr(token, '=');
      if( equalsign != NULL )
      {
         *equalsign = '\0';
         name = token;
         value = equalsign + 1;
         /* set parameter values */
         if( !strcmp(name, "MODE") ) /* Solving settings (1 for (COPT), 2 for (Benders) or 3 for (Benders + Initial Cuts)) */
            Inst->SolvingSetting = atoi(value);
         else if( !strcmp(name, "CUST") ) /* Number of customers */
            Inst->NumCust = atoi(value);
         else if( !strcmp(name, "SUPP") ) /* Number of suppliers */
            Inst->NumSupp = atoi(value);
         else if( !strcmp(name, "BUDGET") ) /* Budget (i.e. the Cardinality restriction) */
            Inst->Budget = atoi(value);
         else if( !strcmp(name, "SAMP") ) /* Number of the Monte-Carlo simulations */
            Inst->NumSamp = atoi(value);
         else if( !strcmp(name, "QOS") ) /* Limit of the quality of service */
            Inst->QoS = atoi(value);
         else if( !strcmp(name, "TIME") ) /* Time limitation (s) */
            Inst->TimeLimit = atof(value);
         else if( !strcmp(name, "MEM") ) /* Memory limitation (MB) */
            Inst->MemLimit = atoi(value);
         else if( !strcmp(name, "SEED") ) /* Random seed of sampling */
            Inst->ProduceSeed = atoi(value);
         else if( !strcmp(name, "MINGAP") ) /* Minimal gap */
            Inst->MinGap = atof(value);
         else if( !strcmp(name, "TYPE") )
            Inst->DisasterType = atoi(value);
         else if( !strcmp(name, "PROB1") )
            Inst->Prob1 = atof(value);
         else if( !strcmp(name, "PROB2") )
            Inst->Prob2 = atof(value);
         else if( !strcmp(name, "SETTING") )
         {
            if( CheckEnds(value) )
            {
               setfile = new char[200];
               strcpy(setfile, value);
            }
            else
               printf("unidentified settings file - %s\n", value);
         }
         else
            printf("undefined parameter - %s\n", name);
      }
      token = strtok(NULL, " ");
   }
   if( setfile != NULL )
   {
      ReadFile(Inst, setfile);
      delete []setfile;
   }
}

/* read parameters from a file */
void ReadFile(InstStruct* Inst, char* paramfile)
{
   char* name;
   char* value;
   FILE* file = fopen(paramfile, "r");
   if( file != NULL )
   {
      char line[200];
      while( fgets(line, sizeof(line), file) )
      {
         char* delimiter = strchr(line, '=');
         if( delimiter != NULL )
         {
            *delimiter = '\0';
            /* parameter name */
            name = line;
            DeleteSpace(name);
            /* parameter value */
            value = delimiter + 1;
            DeleteSpace(value);
            /* set parameter values */
            if( !strcmp(name, "MODE") ) /* Solving settings (1 for (COPT), 2 for (Benders) or 3 for (Benders + Initial Cuts)) */
               Inst->SolvingSetting = atoi(value);
            else if( !strcmp(name, "CUST") ) /* Number of customers */
               Inst->NumCust = atoi(value);
            else if( !strcmp(name, "SUPP") ) /* Number of suppliers */
               Inst->NumSupp = atoi(value);
            else if( !strcmp(name, "BUDGET") ) /* Budget (i.e. the Cardinality restriction) */
               Inst->Budget = atoi(value);
            else if( !strcmp(name, "SAMP") ) /* Number of the Monte-Carlo simulations */
               Inst->NumSamp = atoi(value);
            else if( !strcmp(name, "QOS") ) /* Limit of the quality of service */
               Inst->QoS = atoi(value);
            else if( !strcmp(name, "TIME") ) /* Time limitation (s) */
               Inst->TimeLimit = atof(value);
            else if( !strcmp(name, "MEM") ) /* Memory limitation (MB) */
               Inst->MemLimit = atoi(value);
            else if( !strcmp(name, "SEED") ) /* Random seed of sampling */
               Inst->ProduceSeed = atoi(value);
            else if( !strcmp(name, "MINGAP") ) /* Minimal gap */
               Inst->MinGap = atof(value);
            else if( !strcmp(name, "TYPE") )
               Inst->DisasterType = atoi(value);
            else if( !strcmp(name, "PROB1") )
               Inst->Prob1 = atof(value);
            else if( !strcmp(name, "PROB2") )
               Inst->Prob2 = atof(value);
            else
               printf("undefined parameter - %s\n", name);
         }
      }
      fclose(file);
   }
   else
   {
      printf("parameter file not found - using default parameters\n");
   }
}

/* check whether the input parameter is a .set file */
bool CheckEnds(char* input)
{
   size_t length = strlen(input);
   if( length < strlen(".set") )
      return false;
   char* ends = input + (length - strlen(".set"));
   if( !strcmp(ends, ".set") )
      return true;
   else
      return false;
}

/* remove space at the head and tail of the string */
void DeleteSpace(char* str)
{
   char *start = str;
   char *end = str + strlen(str) - 1;

   while( *start == ' ' )
      start++;

   while( end > start && *end == ' ' )
      end--;

   /* add the end flag */
   *(end + 1) = '\0';

   /* change the starting position */
   memmove(str, start, strlen(start) + 1);
}

/* Free and close COPT environment */
void CleanCOPTModel(InstStruct* Inst)
{
   /* Clena  problem */
   COPT_DeleteProb(&Inst->COPTLp);
   if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }
   COPT_DeleteEnv(&Inst->COPTEnv);
   if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }
}

void FreeMemory(InstStruct* Inst)
{
   delete []Inst->CustPos;
   delete []Inst->SuppPos;
   for( int i = 0; i < Inst->NumCust; i++ )
      delete Inst->Cust[i];
   for( int j = 0; j < Inst->NumSupp; j++ )
      delete Inst->Supp[j];
   for( int e = 0; e < Inst->NumEdge; e++ )
      delete Inst->Edge[e];
   for( int t = 0; t < Inst->NumSamp; t++ )
   {
      delete []Inst->Samp[t]->Cost;
      delete []Inst->Samp[t]->Demand;
      delete Inst->Samp[t];
   }
}
