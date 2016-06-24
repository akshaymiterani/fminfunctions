// Copyright (C) 2016 - IIT Bombay - FOSSEE
//
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
// Author: Harpreet Singh, Pranav Deshpande and Akshay Miterani
// Organization: FOSSEE, IIT Bombay
// Email: toolbox@scilab.in

#define __USE_DEPRECATED_STACK_FUNCTIONS__

#include <iomanip>
#include <fstream>
#include <iostream>
#include "CoinPragma.hpp"
#include "CoinTime.hpp"
#include "CoinError.hpp"

#include "BonOsiTMINLPInterface.hpp"
#include "BonIpoptSolver.hpp"
#include "minuncTMINLP.hpp"
#include "BonCbc.hpp"
#include "BonBonminSetup.hpp"

#include "BonOACutGenerator2.hpp"
#include "BonEcpCuts.hpp"
#include "BonOaNlpOptim.hpp"

#include "sci_iofunc.hpp"
extern  "C"
{
#include "call_scilab.h"
#include <api_scilab.h>
#include <Scierror.h>
#include <BOOL.h>
#include <localization.h>
#include <sciprint.h>

int cpp_intfminunc(char *fname)
{
  using namespace Ipopt;
  using namespace Bonmin;

	CheckInputArgument(pvApiCtx, 11, 11); // We need total 12 input arguments.
	CheckOutputArgument(pvApiCtx, 3, 3);  // 3 output arguments
	
	//Function pointers, input matrix(Starting point) pointer, flag variable 
	int* funptr=NULL;
	int* gradhesptr=NULL;
	double* x0ptr=NULL;
	double flag1,flag2;
        
	// Input arguments
	double *cpu_time=NULL,*max_iter=NULL;
	static unsigned int nVars = 0,nCons = 0;
	unsigned int temp1 = 0,temp2 = 0, iret = 0;
	int x0_rows, x0_cols;
	double *intcon = NULL,*options=NULL, *ifval=NULL;
	int intconSize;
	
	// Output arguments
	double *fX = NULL, ObjVal=0,iteration=0,cpuTime=0,fobj_eval=0;
	double dual_inf, constr_viol, complementarity, kkt_error;
	double *fGrad=  NULL;
	double *fHess=  NULL;
	int rstatus = 0;
	int int_fobj_eval, int_constr_eval, int_fobj_grad_eval, int_constr_jac_eval, int_hess_eval;
	
	
  	//Objective Function
	if(getFunctionFromScilab(1,&funptr))
	{
		return 1;
	}

 	//Function for gradient and hessian
	if(getFunctionFromScilab(2,&gradhesptr))
	{
		return 1;
	}
	
	//Flag for Gradient from Scilab
	if(getDoubleFromScilab(3, &flag1))
	{
		return 1;
	}
	
	//Flag for Hessian from Scilab
	if(getDoubleFromScilab(5, &flag2))
	{
		return 1;
	}

	//x0(starting point) matrix from scilab
	if(getDoubleMatrixFromScilab(7, &x0_rows, &x0_cols, &x0ptr))
	{
		return 1;
	}
	
	if(getIntFromScilab(8, &intconSize))
	{
	  return 1;
	}
	
	temp1 = 1;
	temp2 = intconSize;
	// Getting intcon
	if (getDoubleMatrixFromScilab(9,&temp1,&temp2,&intcon))
	{
		return 1;
	}

    //Initialization of parameters
	nVars=x0_cols;
	nCons=0;
	
	temp1=1;
	temp2=5;
	if (getFixedSizeDoubleMatrixFromScilab(10,temp1,temp2,&options))
	{
		return 1;
	}
	
	temp1=1;
	temp2=5;
	if (getFixedSizeDoubleMatrixFromScilab(11,temp1,temp2,&ifval))
	{
		return 1;
	}
	


	
	SmartPtr<minuncTMINLP> tminlp = new minuncTMINLP(nVars, nCons, x0ptr, flag1, flag2, intconSize, intcon);

  BonminSetup bonmin;
  bonmin.initializeOptionsAndJournalist();
  //Register an additional option
   
  // Here we can change the default value of some Bonmin or Ipopt option
  bonmin.options()->SetStringValue("mu_oracle","loqo");
  
    //Register an additional option
	if((int)ifval[0])
            bonmin.options()->SetNumericValue("bonmin.integer_tolerance", (options[0]));
    if((int)ifval[1])
            bonmin.options()->SetIntegerValue("bonmin.node_limit", (int)(options[1]));
    if((int)ifval[2])
            bonmin.options()->SetNumericValue("bonmin.time_limit", (options[2]));
    if((int)ifval[3])
            bonmin.options()->SetNumericValue("bonmin.allowable_gap", (options[3]));
    if((int)ifval[4])
            bonmin.options()->SetIntegerValue("bonmin.iteration_limit", (int)(options[4]));
  
  //Now initialize from tminlp
  bonmin.initialize(GetRawPtr(tminlp));
  
  //Set up done, now let's branch and bound
  try {
    Bab bb;
    bb(bonmin);//process parameter file using Ipopt and do branch and bound using Cbc
  }
  catch(TNLPSolver::UnsolvedError *E) {
    //There has been a failure to solve a problem with Ipopt.
    std::cerr<<"Ipopt has failed to solve a problem!"<<std::endl;
    sciprint(999, "\nIpopt has failed to solve the problem!\n");
  }
  catch(OsiTMINLPInterface::SimpleError &E) {
    std::cerr<<E.className()<<"::"<<E.methodName()
	     <<std::endl
	     <<E.message()<<std::endl;
	  sciprint(999, "\nFailed to solve a problem!\n");
	}
  catch(CoinError &E) {
    std::cerr<<E.className()<<"::"<<E.methodName()
	     <<std::endl
	     <<E.message()<<std::endl;
	  sciprint(999, "\nFailed to solve a problem!\n");
	}
	rstatus=tminlp->returnStatus();
	if (rstatus >= 0 | rstatus <= 5){
		fX = tminlp->getX();
		ObjVal = tminlp->getObjVal();
		if (returnDoubleMatrixToScilab(1, 1, nVars, fX))
		{
			return 1;
		}

		if (returnDoubleMatrixToScilab(2, 1, 1, &ObjVal))
		{
			return 1;
		}

		if (returnIntegerMatrixToScilab(3, 1, 1, &rstatus))
		{
			return 1;
		}
		
	}
	else
	{
		if (returnDoubleMatrixToScilab(1, 0, 0, fX))
		{
			return 1;
		}

		if (returnDoubleMatrixToScilab(2, 1, 1, &ObjVal))
		{
			return 1;
		}

		if (returnIntegerMatrixToScilab(3, 1, 1, &rstatus))
		{
			return 1;
		}
		
		sciprint(999, "\nThe problem could not be solved!\n");
  }

	// As the SmartPtrs go out of scope, the reference count
	// will be decremented and the objects will automatically
	// be deleted(No memory leakage). 
	
  return 0;
}
}

