# RKHS
Fortran90 implementation of the RKHS method



1) Compile the source code by opening a terminal and writing
> make

NOTE: You might change the Fortran compiler specified in the first line of the Makefile.

2) Run the code by typing 
> ./example.x

3) The code should produce the following output

> Slow evaluation at point    1.5000000000000000       0.50000000000000000       gives    1.6807426393791589     
 Slow evaluation at point    1.5000000000000000       0.50000000000000000       gives derivative of f with respect to x(1): -1.7620891821210760     
 Slow evaluation at point    1.5000000000000000       0.50000000000000000       gives derivative of f with respect to x(2):    2.9379071414204212     
 Slow evaluation at point    1.5000000000000000       0.50000000000000000       gives the Hessian: 
   2.0514523884667937        0.0000000000000000     
  -2.9275790292041037        5.5628495784683452     
 Fast evaluation at point    1.5000000000000000       0.50000000000000000       gives    1.6807426363836484     
 Fast evaluation at point    1.5000000000000000       0.50000000000000000       gives derivative of f with respect to x(1):   -1.7620891821013807     
 Fast evaluation at point    1.5000000000000000       0.50000000000000000       gives derivative of f with respect to x(2):    2.9379071414564648     
 Fast evaluation at point    1.5000000000000000       0.50000000000000000       gives the Hessian: 
   2.0514523888919172       -2.9275790292041037     
  -2.9275790292507025        5.5628495784669223     
 Fast evaluation at point    1.5000000000000000       0.50000000000000000       gives    1.6807426363836484     
 Fast evaluation at point    1.5000000000000000       0.50000000000000000       gives derivative of f with respect to x(1):   -1.7620891821013807     
 Fast evaluation at point    1.5000000000000000       0.50000000000000000       gives derivative of f with respect to x(2):    2.9379071414564648     
 Fast evaluation at point    1.5000000000000000       0.50000000000000000       gives the Hessian: 
   2.0514523888919172       -2.9275790292507025     
  -2.9275790292507025        5.5628495784669223     

> Evaluation with supplied mask:
 Fast evaluation at point    1.5000000000000000       0.50000000000000000       gives    1.6807426363836484     
 Fast evaluation at point    1.5000000000000000       0.50000000000000000       gives derivative of f with respect to x(1):   -1.7620891821013807     
 Fast evaluation at point    1.5000000000000000       0.50000000000000000       gives derivative of f with respect to x(2):    0.0000000000000000     
 Fast evaluation at point    1.5000000000000000       0.50000000000000000       gives the Hessian: 
   2.0514523888919172        0.0000000000000000     
   0.0000000000000000        5.5628495784669223  
 
 
 and should generate a binary file called "test.kernel" and a .csv file called "multidimensional-grid-RECOVERED.csv".
 
 TROUBLESHOOTING:
 If the code produces a different output, the most likely problem is that "multidimensional-grid.csv" could not be read properly.
 If available, check the contents of "multidimensional-grid-RECOVERED.csv" and see whether it matches to the contents of 
 "multidimensional-grid.csv". Most often the issue can be resolved by converting the line endings in "multidimensional-grid.csv" to
 the control character that is appropriate for your system.
 
4) For a tutorial on how to incorporate the RKHS module into your own code, go to src/example.f90. There you will find a step by step
   tutorial how kernels are initialized and evaluated. It also gives details on the format required for training data.
 
 
5) For a detailed step-by-step tutorial on how to use the RKHS module to construct PES, download PES_Tutorial.zip and execute the example codes.
 

 
 
