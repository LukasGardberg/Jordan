function [ev,mult] = heltalsev(A, tol)
% Calculates the integer eigenvalues of A and their multiplicities.

% OUTPUT: 
% ev: Column matrix of found eigenvalues
% mult: Multiplicity of corresponding eigenvalues

% INPUT:
% A: Integer matrix
% tol: Rounding tolerance
% Anger hur stor avvikelse vi kan ha f?r att ?nd? r?kna ett egenv?rde som
% ett heltal. Heltal N +- tol r?knas allts? som heltal (?).

if nargin < 2
    tol = 0.01;
end

% Max tolerance, sqrt(2)/2, passes all eigenvalues

% Gets eigenvalues
eigs = eig(A);

% Rounded eigenvalues for comparison
eig_rounded = round(eigs);

% Checks if any egeinvalue is outisde the threshhold for rounding
distances = abs(eigs - eig_rounded);
eig_logical = distances > tol;

% Checks if any eigenvalue is not within our threshhold. If so, throw error
if any(eig_logical)
    error('The matrix has non-integer eigenvalues, exiting');
end

% eig_rounded now has all of our integer eigenvalues

ev = unique(eig_rounded);
mult = zeros(size(ev));

for i = 1:length(ev)
    mult(i) = sum(eig_rounded == eig_rounded(i));
end

end