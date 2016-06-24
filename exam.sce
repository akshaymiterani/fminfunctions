function y = f(x) 
  y = 3*x(1)^2 + 2*x(1)*x(2) + x(2)^2 - 4*x(1) + 5*x(2) ;
endfunction

x0 = [0,0];
options = list("MaxNodes",1);

[xval, fval, status, gradient, hessian] = intfminunc(f, x0,[1],options)
