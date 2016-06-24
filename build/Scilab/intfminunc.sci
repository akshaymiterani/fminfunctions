// Copyright (C) 2015 - IIT Bombay - FOSSEE
//
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
// Author: Harpreet Singh, Pranav Deshpande and Akshay Miterani
// Organization: FOSSEE, IIT Bombay
// Email: toolbox@scilab.in

function [xopt,fopt,exitflag,gradient,hessian] = intfminunc (varargin)
  	// Solves a multi-variable unconstrainted optimization problem
  	//
  	//   Calling Sequence
  	//   xopt = fminunc(f,x0)
  	//   xopt = fminunc(f,x0,intcon)
  	//   xopt = fminunc(f,x0,intcon,options)
  	//   [xopt,fopt] = fminunc(.....)
  	//   [xopt,fopt,exitflag]= fminunc(.....)
  	//   [xopt,fopt,exitflag,gradient,hessian]= fminunc(.....)
  	//
  	//   Parameters
  	//   f : a function, representing the objective function of the problem 
  	//   x0 : a vector of doubles, containing the starting of variables.
  	//	 intcon : a vector of integers, represents which variables are constrained to be integers
  	//   options: a list, containing the option for user to specify. See below for details.
  	//   xopt : a vector of doubles, the computed solution of the optimization problem.
  	//   fopt : a scalar of double, the function value at x. 
  	//   exitflag : a scalar of integer, containing the flag which denotes the reason for termination of algorithm. See below for details.
  	//
  	//   Description
  	//   Search the minimum of an unconstrained optimization problem specified by :
  	//   Find the minimum of f(x) such that 
  	//
  	//   <latex>
  	//    \begin{eqnarray}
  	//    &\mbox{min}_{x}
  	//    & f(x)\\
  	//    \end{eqnarray}
  	//   </latex>
  	//
  	//   The routine calls Bonmin for solving the Un-constrained Optimization problem, Bonmin is a library written in C++.
  	//
  	// The options allows the user to set various parameters of the Optimization problem. 
  	// It should be defined as type "list" and contains the following fields.
	// <itemizedlist>
	//   <listitem>Syntax : options= list("IntegerTolerance", [---], "MaxNodes", [---], "MaxTime", ---, "AllowableGap", ---, "MaxIterations", ---);</listitem>
	//   <listitem>IntegerTolerance : a Scalar, containing the Integer tolerance value that the solver should take.</listitem>
	//   <listitem>MaxNodes : a Scalar, containing the maximum nodes that the solver should make.</listitem>
	//	 <listitem>MaxIterations : a Scalar, containing the Maximum Number of Iteration that the solver should take.</listitem>
	//   <listitem>AllowableGap : a Scalar, containing the allowable gap value that the solver should take.</listitem>
	//   <listitem>MaxTime : a Scalar, containing the Maximum amount of CPU Time that the solver should take.</listitem>
	// </itemizedlist>
	//
	// The exitflag allows to know the status of the optimization which is given back by Bonmin.
	// <itemizedlist>
	//   <listitem>exitflag=0 : Optimal Solution Found. </listitem>
	//   <listitem>exitflag=1 : InFeasible Solution.</listitem>
	//   <listitem>exitflag=2 : Output is Continuous Unbounded.</listitem>
	//   <listitem>exitflag=3 : Limit Exceeded.</listitem>
	//   <listitem>exitflag=4 : User Interrupt.</listitem>
	//   <listitem>exitflag=5 : MINLP Error.</listitem>
	// </itemizedlist>
	//
	// For more details on exitflag see the Bonmin page, go to http://www.coin-or.org/Bonmin
	//
	//
  	// Examples
  	//     //Find x in R^2 such that it minimizes the Rosenbrock function 
  	//     //f = 100*(x2 - x1^2)^2 + (1-x1)^2
  	//     //Objective function to be minimised
  	//     function y= f(x)
  	//   	    y= 100*(x(2) - x(1)^2)^2 + (1-x(1))^2;
  	//     endfunction
  	//	 //Starting point  
  	//     x0=[-1,2];
  	//	 //Gradient of objective function
  	//     function y= fGrad(x)
  	//   	     y= [-400*x(1)*x(2) + 400*x(1)^3 + 2*x(1)-2, 200*(x(2)-x(1)^2)];
  	//     endfunction
  	//	 //Hessian of Objective Function
  	//     function y= fHess(x)
  	//   	     y= [1200*x(1)^2- 400*x(2) + 2, -400*x(1);-400*x(1), 200 ];
  	//     endfunction
  	//     //Options
  	//     options=list("MaxIter", [1500], "CpuTime", [500], "Gradient", fGrad, "Hessian", fHess);
  	//     //Calling Ipopt
  	//     [xopt,fopt,exitflag,output,gradient,hessian]=fminunc(f,x0,options)
	// // Press ENTER to continue
  	//
  	// Examples
  	//      //Find x in R^2 such that the below function is minimum
  	//      //f = x1^2 + x2^2
  	//      //Objective function to be minimised
  	//      function y= f(x)
  	//   	     y= x(1)^2 + x(2)^2;
  	//      endfunction
  	//	  //Starting point  
  	//      x0=[2,1];
  	//      //Calling Ipopt  
  	//      [xopt,fopt]=fminunc(f,x0)
	// // Press ENTER to continue
  	//
  	// Examples
  	//     //The below problem is an unbounded problem:
  	//     //Find x in R^2 such that the below function is minimum
  	//     //f = - x1^2 - x2^2
  	//     //Objective function to be minimised
  	//     function y= f(x)
  	//        y= -x(1)^2 - x(2)^2;
  	//     endfunction
  	//	 //Starting point  
  	//     x0=[2,1];
  	//	 //Gradient of objective function
  	//     function y= fGrad(x)
  	//   	     y= [-2*x(1),-2*x(2)];
  	//     endfunction
  	//	 //Hessian of Objective Function
  	//     function y= fHess(x)
  	//   	     y= [-2,0;0,-2];
  	//     endfunction
    //     //Options
  	//     options=list();
  	//    //Calling Ipopt  
  	//    [xopt,fopt,exitflag,gradient,hessian]=fminunc(f,x0,options)


	//To check the number of input and output arguments
   	[lhs , rhs] = argn();
	
	//To check the number of arguments given by the user
   	if ( rhs<2 | rhs>4 ) then
    		errmsg = msprintf(gettext("%s: Unexpected number of input arguments : %d provided while should be int [2 3 4] "), "intfminunc", rhs);
    		error(errmsg);
   	end
 
	//Storing the 1st and 2nd Input Parameters  
   	fun = varargin(1);
   	x0 = varargin(2);
      
	//To check whether the 1st Input argument(fun) is a function or not
   	if (type(fun) ~= 13 & type(fun) ~= 11) then
   		errmsg = msprintf(gettext("%s: Expected function for Objective "), "intfminunc");
   		error(errmsg);
   	end
   
	//To check whether the 2nd Input argument(x0) is a vector/scalar
   	if (type(x0) ~= 1) then
   		errmsg = msprintf(gettext("%s: Expected Vector/Scalar for Starting Point"), "intfminunc");
   		error(errmsg);
  	end
   
	//To check and convert the 2nd Input argument(x0) to a row vector 
   	if((size(x0,1)~=1) & (size(x0,2)~=1)) then
   		errmsg = msprintf(gettext("%s: Expected Row Vector or Column Vector for x0 (Initial Value) "), "intfminunc");
   		error(errmsg);
   	else
   		if(size(x0,2)==1) then
   			x0=x0';		//Converting x0 to a row vector, if it is a column vector
   		else 
   	 		x0=x0;		//Retaining the same, if it is already a row vector
   		end   	 	
        	s=size(x0);	
   	end
   

  	//To check the match between f (1st Parameter) and x0 (2nd Parameter)
   	if(execstr('init=fun(x0)','errcatch')==21) then
		errmsg = msprintf(gettext("%s: Objective function and x0 did not match"), "intfminunc");
   		error(errmsg);
	end
	
	//Converting the User defined Objective function into Required form (Error Detectable)
   	function [y,check] = f(x)
   		if(execstr('y=fun(x)','errcatch')==32 | execstr('y=fun(x)','errcatch')==27)
			y=0;
			check=1;
		else
			y=fun(x);
			if (isreal(y)==%F) then
				y=0;
				check=1;
  			else
				check=0;
			end
		end
	endfunction
   
   //To add intcon
   intcon=[];
   if ( rhs >=3 ) then
			intcon = varargin(3);
	end
	
	//Error Checks for intcon
	 
	for i=1:size(intcon,2)
      if(intcon(i)>s) then
        errmsg = msprintf(gettext("%s: The values inside intcon should be less than the number of variables"), "intqpipopt");
        error(errmsg);
      end

      if (intcon(i)<0) then
        errmsg = msprintf(gettext("%s: The values inside intcon should be greater than 0 "), "intqpipopt");
        error(errmsg);
      end

      if(modulo(intcon(i),1)) then
        errmsg = msprintf(gettext("%s: The values inside intcon should be an integer "), "intqpipopt");
        error(errmsg);
      end
	end
   
	//To check whether options has been entered by user   
   	if ( rhs<4  ) then
      		param = list();
       else
      		param =varargin(4); //Storing the 3rd Input Parameter in an intermediate list named 'param'
    
   	end
   
	//If options has been entered, then check its type for 'list'   
   	if (type(param) ~= 15) then
   		errmsg = msprintf(gettext("%s: 3rd Input parameter should be a list (ie. Options) "), "intfminunc");
   		error(errmsg);
   	end
   
	//If options has been entered, then check whether an even number of entires has been entered   
   	if (modulo(size(param),2)) then
		errmsg = msprintf(gettext("%s: Size of parameters should be even"), "intfminunc");
		error(errmsg);
   	end
   	gradient = []
   	//Defining a function to calculate Gradient or Hessian if the respective user entry is OFF 
   	function [y,check]=gradhess(x,t)
		if t==1 then	//To return Gradient
			if(execstr('y=numderivative(fun,x)','errcatch')==10000)
				y=0;
				check=1;
			else
				y=numderivative(fun,x);
				if (isreal(y)==%F) then
					y=0;
					check=1;
  				else
					check=0;
				end
			end			
		else		//To return Hessian
			if(execstr('[grad,y]=numderivative(fun,x)','errcatch')==10000)
				y=0;
				check=1;
			else
				[grad,y]=numderivative(fun,x);
				if (isreal(y)==%F) then
					y=0;
					check=1;
  				else
					check=0;
				end	
			end
		end
   	endfunction
   	
   	options=param;

	//Flags to check whether Gradient is "ON"/"OFF" and Hessian is "ON"/"OFF" 
   	flag1=0;
   	flag2=0;
   	fGrad=[];
   	fGrad1=[];
   	fHess=[];
   	fHess1=[];
 
   
  //To check for correct input of Gradient and Hessian functions from the user	     	
   if (flag1==1) then
   		if (type(fGrad) ~= 13 & type(fGrad) ~= 11) then
  			errmsg = msprintf(gettext("%s: Expected function for Gradient of Objective"), "intfminunc");
   			error(errmsg);
   		end
   		
   		if(execstr('samplefGrad=fGrad(x0)','errcatch')==21)
			errmsg = msprintf(gettext("%s: Gradient function of Objective and x0 did not match"), "intfminunc");
   			error(errmsg);
		end
		
		samplefGrad=fGrad(x0);
		
		if (size(samplefGrad,1)==s(2) & size(samplefGrad,2)==1) then
		elseif (size(samplefGrad,1)==1 & size(samplefGrad,2)==s(2)) then
		else
   			errmsg = msprintf(gettext("%s: Wrong Input for Objective Gradient function(3rd Parameter)---->Row Vector function of size [1 X %d] is Expected"), "intfminunc",s(2));
   			error(errmsg);
   		end
   		
   		function [y,check] = fGrad1(x)
   			if(execstr('y=fGrad(x)','errcatch')==32 | execstr('y=fGrad(x)','errcatch')==27) then
				y = 0;
				check=1;
			else
				y=fGrad(x);
				if (isreal(y)==%F) then
					y = 0;
					check=1;
  				else
					check=0;
				end
			end
  		endfunction
   	end
   	if (flag2==1) then
   		if (type(fHess) ~= 13 & type(fHess) ~= 11) then
  			errmsg = msprintf(gettext("%s: Expected function for Hessian of Objective"), "intfminunc");
   			error(errmsg);
   		end
   		
   		if(execstr('samplefHess=fHess(x0)','errcatch')==21)
			errmsg = msprintf(gettext("%s: Hessian function of Objective and x0 did not match"), "intfminunc");
   			error(errmsg);
		end
		
		samplefHess=fHess(x0);
		
   		if(size(samplefHess,1)~=s(2) | size(samplefHess,2)~=s(2)) then
   			errmsg = msprintf(gettext("%s: Wrong Input for Objective Hessian function(3rd Parameter)---->Symmetric Matrix function of size [%d X %d] is Expected "), "intfminunc",s(2),s(2));
   			error(errmsg);
   		end
   		
   		function [y,check] = fHess1(x)
   			if(execstr('y=fHess(x)','errcatch')==32 | execstr('y=fHess(x)','errcatch')==27)
				y = 0;
				check=1;
			else
				y=fHess(x);
				if (isreal(y)==%F) then
					y = 0;
					check=1;
  				else
					check=0;
				end
			end
  		endfunction
   	end
	intconSize = length(intcon);
	
	//Pushing options as required to a double array
    optval = [];
    ifval =[];
    if length(options) == 0 then
        optval = [0 0 0 0 0];
        ifval = [0 0 0 0 0];
    else
        optval = [0 0 0 0 0];
        ifval = [0 0 0 0 0];
        for i=1:2:length(options)
            select options(i)
                case 'IntegerTolerance' then
					ifval(1) = 1;
                    optval(1) = options(i+1);
                case 'MaxNodes' then
					ifval(2) = 1;
                    optval(2) = options(i+1);
                case 'MaxTime' then
					ifval(3) = 1;
                    optval(3) = options(i+1);
                case 'AllowableGap' then
					ifval(4) = 1;
                    optval(4) = options(i+1);
                case 'MaxIterations' then
					ifval(5) = 1;
                    optval(5) = options(i+1);
                else
                    error(999, 'Unknown string argument passed.');
                end
        end
    end
    
    //Calling the Ipopt function for solving the above problem
	[xopt,fopt,status] = inter_fminunc(f,gradhess,flag1,fGrad1,flag2,fHess1,x0,intconSize,intcon,optval,ifval);
   
	//Calculating the values for output
   	xopt = xopt';
   	exitflag = status;

	//In the cases of the problem not being solved, return NULL to the output matrices
	if( status~=0 & status~=1 & status~=2 & status~=3 & status~=4 & status~=5 ) then
		xopt=[]
		fopt=[]
	end
	
	[ gradient, hessian] = numderivative(f, xopt, [], [], "blockmat");

	
    //To print output message
    select status
    
    case 0 then
        printf("\nOptimal Solution Found.\n");
    case 1 then
        printf("\nInFeasible Solution.\n");
    case 2 then
        printf("\nOutput is Continuous Unbounded.s\n");
    case 3 then
        printf("\nTime Limit Exceeded.\n");
    case 4 then
        printf("\nUser Interrupt.\n");
    case 5 then
        printf("\nMINLP Error.\n");
    else
        printf("\nInvalid status returned. Notify the Toolbox authors\n");
        break;
    end
    	

endfunction
