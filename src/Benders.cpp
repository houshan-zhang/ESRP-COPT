#include "Benders.h"



/* COPT CALLBACK of adding Benders optimality cuts for binary relaxed solution */
int COPT_CALL BinCallback(copt_prob* LP, void* cb_data, int cbctx, void* cb_handle)
{
   InstStruct* Inst = (InstStruct*)cb_handle;

   Inst->CutRhs = 0.0;
   for( int j = 0; j < Inst->NumReco; j++ )
      Inst->CutVal[j] = 0.0; /* Initialize the cut Coefficients */
   Inst->CutVal[Inst->NumReco] = 1;
   double mincost;
   bool isviolate = false;

   /* Get relaxed solution */
   Inst->COPTStatus = COPT_GetCallbackInfo(cb_data, COPT_CBINFO_MIPCANDIDATE, Inst->RelaxSol);
   if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }

   /* Get the indexes where the solution takes the value of one */
   Inst->OneSupp.resize(0);
   for( int j = 0; j < Inst->NumReco; j++ )
   {
      SuppStruct* supp = Inst->Reco[j];
      if( Inst->RelaxSol[j] >= 0.999 )
         Inst->OneSupp.push_back(supp);
      else if( Inst->RelaxSol[j] >= 0.001 )
         PrintErr();
   }

   /* Open Suppliers */
   SetSupplier(Inst, Inst->OneSupp, 1.0);

   /* Calculate and add possible Benders cuts for each scenario */
   for( int t = 0; t < Inst->NumSamp; t++ )
   {
      for( int j = 0; j < Inst->NumReco; j++ )
         Inst->CutVal[j] = 0.0; /* Initialize the cut Coefficients */
      /* Change coeffcients */
      ChgCOPTCoef(Inst, t);

      /* Set suppliers status */
      SetSupplierStatus(Inst, Inst->Supp, t);

      /* Solve subproblem */
      Inst->COPTStatus = COPT_Solve(Inst->COPTLpSub);
      if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }
      Inst->COPTStatus = COPT_GetDblAttr(Inst->COPTLpSub, COPT_DBLATTR_LPOBJVAL, &mincost);
      if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }
      Inst->CutRhs = mincost;

      /* Get dual values */
      Inst->COPTStatus = COPT_GetRowInfo(Inst->COPTLpSub, COPT_DBLINFO_DUAL, Inst->NumSupp, Inst->SuppPos, Inst->Pi);
      if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }

      /* Get the right hand side of Benders cut */
      for( unsigned j = 0; j < Inst->OneSupp.size(); j++ )
      {
         SuppStruct* supp = Inst->OneSupp[j];
         if( supp->alive[t] )
            Inst->CutRhs -= supp->bandwidth * Inst->Pi[supp->id];
      }

      /* Get coefficients of Benders cut */
      for( int j = 0; j < Inst->NumReco; j++ )
      {
         SuppStruct* supp = Inst->Reco[j];
         if( supp->alive[t] )
         {
            Inst->CutVal[j] = -supp->bandwidth * Inst->Pi[supp->id];
         }
      }

      /* Calculate the violated value of cut */
      Inst->VioVal = mincost - Inst->RelaxSol[Inst->NumReco + t];

      /* If the violated value is greater than epsilon, then add the Benders cut by callback function */
      if( Inst->VioVal > BENDERS_EPSILON )
      {
         if( !isviolate )
         {
            isviolate = true;
            Inst->NumIter++;
         }
         Inst->CutInd[Inst->NumReco] = Inst->NumReco + t;
         Inst->COPTStatus = COPT_AddCallbackLazyConstr(cb_data, Inst->NumReco + 1, Inst->CutInd, Inst->CutVal, COPT_GREATER_EQUAL, Inst->CutRhs);
         if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }
         Inst->NumBinCut++; /* Number of Benders cuts for binary solution plus one */
      }
   }

   /* Close Suppliers */
   SetSupplier(Inst, Inst->OneSupp, 0.0);

   return 0;
}

//TODO
/* COPT CALLBACK of adding Benders optimality cuts for fractional relaxed solution */
int COPT_CALL FracCallback(copt_prob* LP, void* cb_data, int cbtx, void* cb_handle)
{
   InstStruct* Inst = (InstStruct*)cb_handle;

   Inst->CutRhs = 0.0;
   double mincost;
   double FracEPS;
   bool isviolate = false;
   for( int j = 0; j < Inst->NumReco; j++ )
      Inst->CutVal[j] = 0.0; /* Initialize the cut Coefficients */
   Inst->CutVal[Inst->NumReco] = 1;

   FracEPS = BENDERS_EPSILON;

   /* Get relaxed solution */
   Inst->COPTStatus = COPT_GetCallbackInfo(cb_data, COPT_CBINFO_RELAXSOLUTION, Inst->RelaxSol);
   if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }

   /* Open suppliers */
   for( int j = 0; j < Inst->NumReco; j++ )
   {
      SuppStruct* supp = Inst->Reco[j];
      SetSupplier(Inst, {supp}, Inst->RelaxSol[j]);
   }

   /* Calculate and add possible Benders cuts for each scenario */
   for( int t = 0; t < Inst->NumSamp; t++ )
   {
      for( int j = 0; j < Inst->NumReco; j++ )
         Inst->CutVal[j] = 0.0; /* Initialize the cut Coefficients */
      GetBendersCut(Inst, t, mincost);
      Inst->CutInd[Inst->NumReco] = Inst->NumReco + t;
      /* Calculate the violated value of cut */
      Inst->VioVal = mincost - Inst->RelaxSol[Inst->NumReco + t];

      /* If the violated value is greater than epsilon, then add the Benders cut by callback function */
      if( Inst->VioVal > FracEPS )
      {
         if( !isviolate )
         {
            isviolate = true;
            Inst->NumIter++;
         }
         Inst->COPTStatus = COPT_AddCallbackUserCut(cb_data, Inst->NumReco + 1, Inst->CutInd, Inst->CutVal, COPT_GREATER_EQUAL, Inst->CutRhs);
         if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }
         Inst->NumFracCut++; /* Number of Benders cuts for binary solution plus one */
      }
   }

   /* Close suppliers */
   SetSupplier(Inst, Inst->Reco, 0.0);

   return 0;
}

//TODO
void GetBendersCut(InstStruct* Inst, int t, double& mincost)
{
   /* Change coeffcients */
   ChgCOPTCoef(Inst, t);

   /* Set suppliers status */
   SetSupplierStatus(Inst, Inst->Supp, t);

   /* Solve subproblem */
   Inst->COPTStatus = COPT_Solve(Inst->COPTLpSub);
   if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }
   Inst->COPTStatus = COPT_GetDblAttr(Inst->COPTLpSub, COPT_DBLATTR_LPOBJVAL, &mincost);
   if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }
   Inst->CutRhs = mincost;

   /* Get dual values */
   Inst->COPTStatus = COPT_GetDblAttr(Inst->COPTLpSub, COPT_DBLINFO_DUAL, Inst->Pi);
   if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }

   /* Get coefficients of Benders cut */
   for( int j = 0; j < Inst->NumReco; j++ )
   {
      SuppStruct* supp = Inst->Reco[j];
      if( supp->alive[t] )
      {
         Inst->CutRhs -= supp->bandwidth * Inst->Pi[supp->id] * Inst->RelaxSol[supp->pos];
         Inst->CutVal[j] = -supp->bandwidth * Inst->Pi[supp->id];
      }
   }
}

/* Build the Benders subproblem formulation */
void BuildSubModel(InstStruct* Inst)
{
   /* Open COPT environment */
   Inst->COPTStatus = COPT_CreateEnv(&Inst->COPTEnvSub);
   if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }
   /* Create NRP InstStruct */
   Inst->COPTStatus = COPT_CreateProb(Inst->COPTEnvSub, &Inst->COPTLpSub);
   if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }
   /* COPT_OFF is default value of ScreenOutput in callback mode */
   Inst->COPTStatus = COPT_SetIntParam(Inst->COPTLpSub, COPT_INTPARAM_LOGTOCONSOLE, 0);
   if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }
   /* set threads is critical */
   Inst->COPTStatus = COPT_SetIntParam(Inst->COPTLpSub, COPT_INTPARAM_THREADS, 1);
   if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }
   /* Select dual simplex method optimizer */
   Inst->COPTStatus = COPT_SetIntParam(Inst->COPTLpSub, COPT_INTPARAM_LPMETHOD, 1);
   if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }

   /* Create and add variables to the Benders subproblem */
   int NumVar = Inst->NumCust + Inst->NumEdge;
   double* Obj = new double[NumVar]; /* Objective Coefficients */
   double* Lb = new double[NumVar]; /* Variables low bounds */
   double* Ub = new double[NumVar]; /* Variables upper bounds */
   char** NameVar = new char*[NumVar]; /* Variables names */

   for( int i = 0; i < NumVar; i++ )
      NameVar[i] = new char[100];

   /* ( Col - 1 ) Create customer variables */
   int counter = 0;
   for( int j = 0; j < Inst->NumCust; j++ )
   {
      Obj[counter] = 0.0;
      Lb[counter] = 0.0;
      Ub[counter] = 1e+15;
      sprintf(NameVar[counter], "loss-%d", j);
      counter++;
   }

   /* ( Col - 2 ) Create flow variables */
   for( int i = 0; i < Inst->NumCust; i++ )
   {
      CustStruct* cust = Inst->Cust[i];
      for( int j = 0; j < cust->degree; j++ )
      {
         Obj[counter] = 0.0;
         Lb[counter] = 0.0;
         Ub[counter] = cust->Edge[j]->supp->bandwidth;
         sprintf(NameVar[counter], "flow-%d-%d", i, cust->Edge[j]->supp->id);
         counter++;
      }
   }

   /* Add variables information to COPT */
   Inst->COPTStatus = COPT_AddCols(Inst->COPTLpSub, NumVar, Obj, NULL, NULL, NULL, NULL, NULL, Lb, Ub, NameVar);
   if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }

   delete []Obj;
   delete []Lb;
   delete []Ub;
   for( int i = 0; i < NumVar; i++ )
   {
      delete []NameVar[i];
   }
   delete []NameVar;

   COPT_SetObjSense(Inst->COPTLpSub, COPT_MINIMIZE); /* Set optimization direction to be "MIN" */

   /* Add the constraints of the subproblem */
   int* ConsInd = new int[Inst->NumEdge + 1];
   double* ConsVal = new double[Inst->NumEdge + 1];
   double rhs;
   char sense;

   /* ( Row - 1 ) Creat the demand constraints */
   sense = COPT_EQUAL; /* Set the constraint "==" */
   for( int i = 0; i < Inst->NumCust; i++ )
   {
      CustStruct* cust = Inst->Cust[i];
      for( int j = 0; j < cust->degree; j++ )
      {
         EdgeStruct* edge = cust->Edge[j];
         ConsVal[j] = 1;
         ConsInd[j] = edge->pos;
      }
      ConsVal[cust->degree] = 1;
      ConsInd[cust->degree] = cust->pos;
      rhs = 0.0;
      Inst->COPTStatus = COPT_AddRow(Inst->COPTLpSub, cust->degree+1, ConsInd, ConsVal, sense, rhs, 0.0, NULL);
      if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }
   }

   /* ( Row - 2 ) Creat the Capacity constraint */
   sense = COPT_LESS_EQUAL; /* Set the constraint "<=" */
   for( int j = 0; j < Inst->NumSupp; j++ )
   {
      SuppStruct* supp = Inst->Supp[j];
      for( int i = 0; i < supp->degree; i++ )
      {
         EdgeStruct* edge = supp->Edge[i];
         ConsVal[i] = 1;
         ConsInd[i] = edge->pos;
      }
      rhs = 0.0;
      Inst->COPTStatus = COPT_AddRow(Inst->COPTLpSub, supp->degree, ConsInd, ConsVal, sense, rhs, 0.0, NULL);
      if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }
   }

   delete []ConsVal;
   delete []ConsInd;
}

/* Build the Benders master formulation */
void BuildMasterModel(InstStruct* Inst)
{
   /* Initialization */
   Inst->CutInd = new int[Inst->NumReco+1];
   for( int j = 0; j < Inst->NumReco; j++ )
      Inst->CutInd[j] = j;
   Inst->CutVal = new double[Inst->NumReco+1];
   /* Relaxed solution during solving process */
   Inst->RelaxSol = new double[Inst->NumReco+Inst->NumSamp];
   Inst->Pi = new double[Inst->NumSupp];

   /* Open COPT environment */
   Inst->COPTStatus = COPT_CreateEnv(&Inst->COPTEnvMaster);
   if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }

   /* Create NRP InstStruct */
   Inst->COPTStatus = COPT_CreateProb(Inst->COPTEnvMaster, &Inst->COPTLpMaster);
   if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }

   /* Create and add variables to the master problem */
   int NumVar;
   NumVar = Inst->NumReco + Inst->NumSamp;
   double* Obj = new double[NumVar]; /* Objective Coefficients */
   double* Lb = new double[NumVar]; /* Variables low bounds */
   double* Ub = new double[NumVar]; /* Variables upper bounds */
   char* TypeVar = new char[NumVar]; /* Variables types */
   char** NameVar = new char*[NumVar]; /* Variables names */

   for( int i = 0; i < NumVar; i++ )
      NameVar[i] = new char[100];

   /* ( Col - 1 ) Create binary variables for each supplier*/
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

   /* ( Col - 2 ) Create continuous variables for each scenario */
   for( int t = 0; t < Inst->NumSamp; t++ )
   {
      Obj[counter] = 1;
      Lb[counter] = 0.0;
      Ub[counter] = Inst->Samp[t]->SumDemand;
      TypeVar[counter] = COPT_CONTINUOUS;
      sprintf(NameVar[counter], "scenario-%d", t);
      counter++;
   }

   /* Add variables information to COPT */
   Inst->COPTStatus = COPT_AddCols(Inst->COPTLpMaster, NumVar, Obj, NULL, NULL, NULL, NULL, TypeVar, Lb, Ub, NameVar);
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

   COPT_SetObjSense(Inst->COPTLpMaster, COPT_MINIMIZE); /* Set optimization direction to be "MIN" */

   /* Add the constraints of the master problem */
   int* ConsInd = new int[Inst->NumReco + 1];
   double* ConsVal = new double[Inst->NumReco + 1];
   char sense;
   double mincost;
   /* ( Row - 1 ) Creat the Cardinality (budget) constraint */
   sense = COPT_LESS_EQUAL; /* Set the constraint "<=" */
   counter = 0;
   for( int j = 0; j < Inst->NumSupp; j++ )
   {
      if( IsEq(Inst->Supp[j]->working, 0.0) )
      {
         ConsVal[counter] = 1;
         ConsInd[counter] = Inst->Supp[j]->pos;
         counter++;
      }
   }
   double rhs = double(Inst->Budget);
   Inst->COPTStatus = COPT_AddRow(Inst->COPTLpMaster, Inst->NumReco, ConsInd, ConsVal, sense, rhs, 0.0, NULL);
   if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }

   if( Inst->SolvingSetting == 3 )
   {
      /* ( Row - 2 ) Creat the submodular constraints */
      sense = COPT_GREATER_EQUAL; /* Set the constraint ">=" */
      for( int j = 0; j < Inst->NumReco; j++ )
         ConsVal[j] = 0.0; /* Initialize the cut Coefficients */
      for( int t = 0; t < Inst->NumSamp; t++ )
      {
         for( int j = 0; j < Inst->NumReco; j++ )
            ConsVal[j] = 0.0; /* Initialize the cut Coefficients */
         counter = 0;
         SampStruct* samp = Inst->Samp[t];
         SetSupplier(Inst, Inst->Reco, 0.0);
         SetSupplierStatus(Inst, Inst->Supp, t);
         ChgCOPTCoef(Inst, t);
         Inst->COPTStatus = COPT_Solve(Inst->COPTLpSub);
         if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }
         Inst->COPTStatus = COPT_GetDblAttr(Inst->COPTLpSub, COPT_DBLATTR_LPOBJVAL, &rhs);
         if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }
         for( int j = 0; j < Inst->NumSupp; j++ )
         {
            SuppStruct* supp = Inst->Supp[j];
            if( IsEq(supp->working, 0.0) )
            {
               SetSupplier(Inst, {supp}, 1.0);
               SetSupplierStatus(Inst, {supp}, t);
               if( supp->alive[t] )
               {
                  Inst->COPTStatus = COPT_Solve(Inst->COPTLpSub);
                  if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }
                  Inst->COPTStatus = COPT_GetDblAttr(Inst->COPTLpSub, COPT_DBLATTR_LPOBJVAL, &mincost);
                  if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }
               }
               else
                  mincost = rhs;
               ConsInd[counter] = supp->pos;
               ConsVal[counter] += (rhs - mincost);
               SetSupplier(Inst, {supp}, 0.0);
               SetSupplierStatus(Inst, {supp}, t);
               counter++;
            }
         }
         if( counter != Inst->NumReco ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }
         ConsInd[Inst->NumReco] = samp->pos;
         ConsVal[Inst->NumReco] = 1;
         Inst->COPTStatus = COPT_AddRow(Inst->COPTLpMaster, Inst->NumReco+1, ConsInd, ConsVal, sense, rhs, 0.0, NULL);
         if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }
      }
#if 1
      /* ( Row - 3 ) Creat the aggregated capacity constraints */
      sense = COPT_GREATER_EQUAL; /* Set the constraint ">=" */
      for( int j = 0; j < Inst->NumReco; j++ )
         ConsVal[j] = 0.0; /* Initialize the cut Coefficients */
      for( int t = 0; t < Inst->NumSamp; t++ )
      {
         for( int j = 0; j < Inst->NumReco; j++ )
            ConsVal[j] = 0.0; /* Initialize the cut Coefficients */
         counter = 0;
         SampStruct* samp = Inst->Samp[t];
         rhs = 0.0;
         for( int i = 0; i < Inst->NumCust; i++ )
         {
            CustStruct* cust = Inst->Cust[i];
            rhs += cust->demand[t];
         }
         for( int j = 0; j < Inst->NumSupp; j++ )
         {
            SuppStruct* supp = Inst->Supp[j];
            if( IsEq(supp->working, 1.0) && supp->alive[t] )
               rhs -= supp->bandwidth;
         }
         rhs = double(rhs)/Inst->NumSamp;
         for( int j = 0; j < Inst->NumSupp; j++ )
         {
            SuppStruct* supp = Inst->Supp[j];
            if( IsEq(supp->working, 0.0) )
            {
               ConsInd[counter] = supp->pos;
               if( supp->alive[t] )
                  ConsVal[counter] = supp->bandwidth;
               else
                  ConsVal[counter] = 0.0;
               counter++;
            }
         }
         if( counter != Inst->NumReco ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }
         ConsInd[Inst->NumReco] = samp->pos;
         ConsVal[Inst->NumReco] = 1;
         Inst->COPTStatus = COPT_AddRow(Inst->COPTLpSub, Inst->NumReco+1, ConsInd, ConsVal, sense, rhs, 0.0, NULL);
         if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }
      }
#endif
   }


   delete []ConsVal;
   delete []ConsInd;
}

/* Solve the NRP under BD framework */
void SolveMasterModel(InstStruct* Inst)
{
   /* Output information to the screen */
   Inst->COPTStatus = COPT_SetIntParam(Inst->COPTLpMaster, COPT_INTPARAM_LOGTOCONSOLE, 1);

   /* Set the number of Threads in the calculation */
   Inst->COPTStatus = COPT_SetIntParam(Inst->COPTLpMaster, COPT_INTPARAM_THREADS, Inst->NumThread);
   if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }

   /* Set the timelimit of the COPT (i.e., total timelimit minus the time cost in presolving) */
   double COPTTimeLimit = 0.0;
   if (Inst->TimeLimit <= Inst->BuildTime)
      COPTTimeLimit = 0.0;
   else
      COPTTimeLimit = Inst->TimeLimit - Inst->BuildTime;
   cout << "COPT Solving Time Limit: " << COPTTimeLimit << endl;
   Inst->COPTStatus = COPT_SetDblParam(Inst->COPTLpMaster, COPT_DBLPARAM_TIMELIMIT, COPTTimeLimit);
   if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }

   /* Set the absolute tolerance of COPT to 0.0 */
   Inst->COPTStatus = COPT_SetDblParam(Inst->COPTLpMaster, COPT_DBLPARAM_ABSGAP, 0.0);
   if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }

   /* Set the relative tolerance of COPT to 0.0 */
   Inst->COPTStatus = COPT_SetDblParam(Inst->COPTLpMaster, COPT_DBLPARAM_RELGAP, 0.0);
   if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }

   /* * solving the MIP model */
   clock_t TimeStart=clock();

   /* Set presolving keeping row structures and the callback functions working on the original model */
   Inst->COPTStatus = COPT_SetIntParam(Inst->COPTLpMaster, COPT_INTPARAM_PRESOLVE, 0);

   /* Set callback function of COPT for those binary relaxed solution */
   Inst->COPTStatus = COPT_SetCallback(Inst->COPTLpMaster, BinCallback, COPT_CBCONTEXT_MIPSOL, (void*)Inst);
   if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }

   // COPT 一次性只支持注册一个回调函数 !!!
   /* Set callback function of COPT for those fractional relaxed solution */
   //Inst->COPTStatus = COPT_SetCallback(Inst->COPTLpMaster, FracCallback, COPT_CBCONTEXT_MIPRELAX, (void*)Inst);
   //if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }

   /* COPT optimization process */
   Inst->COPTStatus = COPT_Solve(Inst->COPTLpMaster);
   if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }

   clock_t TimeEnd=clock();
   Inst->COPTSolvingTime = (double)(TimeEnd - TimeStart)/(double)CLOCKS_PER_SEC;

   /* Get the best feasible objective value within timelimit*/
   Inst->COPTStatus = COPT_GetDblAttr(Inst->COPTLpMaster, COPT_DBLATTR_BESTOBJ, &Inst->OptObj);
   if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }

   /* Get the best relaxed objective value within timelimit*/
   Inst->COPTStatus = COPT_GetDblAttr(Inst->COPTLpMaster, COPT_DBLATTR_BESTBND, &Inst->RelaxObj);
   if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }

   Inst->COPTStatus = COPT_GetIntAttr(Inst->COPTLpMaster, COPT_INTATTR_NODECNT, &Inst->NumNode); /* Get the number of branch-and-bound nodes */
   if( Inst->COPTStatus ) { PrintErr(); cout << "COPTStatus:" << Inst->COPTStatus << endl; }

   /* Output the statistic information */
   cout<<endl<<endl;
   printf("\t(a) Optimal objective value:\t\t\t\t%.2f\n", Inst->OptObj);
   printf("\t(b) Relaxed objective value:\t\t\t\t%.2f\n", Inst->RelaxObj);
   printf("\t(c) Number of branch-and-bound nodes:\t\t\t%d\n", Inst->NumNode);
   printf("\t(d) COPT solving time:\t\t\t\t\t%.2f\n", Inst->COPTSolvingTime);
   printf("\t(e) Total time:\t\t\t\t\t\t%.2f\n", Inst->BuildTime + Inst->COPTSolvingTime);
   printf("\t(f) Number of Benders cuts for binary solution:\t\t%d\n", Inst->NumBinCut);
   printf("\t(g) Number of Benders cuts for fractional solution:\t%d\n", Inst->NumFracCut);
   printf("\t(h) Number of Benders iterations:\t\t\t%d\n", Inst->NumIter);
}

/* Free and close COPT environment */
void CleanBendersModel(InstStruct* Inst)
{
   /* Clena master problem */
   COPT_DeleteProb(&Inst->COPTLpMaster);
   COPT_DeleteEnv(&Inst->COPTEnvMaster);

   /* Clean subproblem */
   COPT_DeleteProb(&Inst->COPTLpSub);
   COPT_DeleteEnv(&Inst->COPTEnvSub);

   delete []Inst->CutVal;
   delete []Inst->CutInd;
   delete []Inst->RelaxSol;
   delete []Inst->Pi;
}
