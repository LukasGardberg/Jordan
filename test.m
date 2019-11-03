A = [1 2 3;
     4 5 6;
     7 8 9];
 
A = testmatris(3);

% A = [-2 -1 1;
%       1 0 1;
%       0 1 1];
%   
%   A = [3 -2;
%        4 -1];
 
[ev, mult] = heltalsev(A);

A = [-1  0  1  0 -1  1;
      5  0 -2 -2  3 -2;
      1  0  0  0  1 -1;
     -7  0  2  2 -5  6;
      4  0 -3  0  4 -2;
      3  0 -2  0  3 -1];
 
  
A = [5 4 2 1;
     0 1 -1 -1;
     -1 -1 3 0;
     1 1 -1 2];
 
 jordanmatris(A)
 
 %%
  
A = compan(poly(1:20));
heltalsev(A)
 
jordanmatris(A)