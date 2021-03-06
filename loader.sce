// This file is released under the 3-clause BSD license. See COPYING-BSD.
// Generated by builder.sce : Please, do not edit this file
// ----------------------------------------------------------------------------
//
libFOSSEE_Scilab_in_path = get_absolute_file_path('loader.sce');
//
// ulink previous function with same name
[bOK, ilib] = c_link('libFOSSEE_Scilab_intqpipopt');
if bOK then
  ulink(ilib);
end
//
list_functions = [ 'inter_fminunc';
];
addinter(libFOSSEE_Scilab_in_path + filesep() + 'libFOSSEE_Scilab_intqpipopt' + getdynlibext(), 'libFOSSEE_Scilab_intqpipopt', list_functions);
// remove temp. variables on stack
clear libFOSSEE_Scilab_in_path;
clear bOK;
clear ilib;
clear list_functions;
// ----------------------------------------------------------------------------
