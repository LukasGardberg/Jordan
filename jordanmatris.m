function J = jordanmatris(A,tol)

% Calculates the Jordan matrix for a square integer matrix A.

% Tolerance for rounding MATLAB precision
if nargin < 2
    tol = 0.01;
end

% Calculates the matrix eigenvalues and its multiplicities
[ev,mult] = heltalsev(A,tol);

size_A = length(A); % Assume square

eigs = cell(length(ev), 1);

nbr_blocks = 0;

% For each eigenvalue:

for j = 1:length(ev)
% Initialize / reset values
r = zeros(size_A, 1); % Rank
p = zeros(size_A, 1); % Kernel dim
b = zeros(size_A, 1);
n = zeros(size_A, 1); % Block multiplicites

    % Calculate the block multiplicities
    i = 1;
    while i <= mult(j)
    
        r(i) = rank((A - ev(j)*eye(size_A))^(i));
        p(i) = size_A - r(i);
        
        if(i == 1)
            b(i) = p(1);
        else
        
            b(i) = p(i) - p(i-1);
        
            % Break early since the rest will be zero
            if(r(i) == r(i-1))
                break
            end
        end
        
    i = i + 1;
    end
    
    size_length = i;
    
    for k = 1:(length(b) - 1)
        
        n(k) = b(k) - b(k + 1);
        
    end
    
    % Save block sizes and multipl. The array block_mult(a, b) denotes b blocks of size a.
    
    block_mult = zeros(size_length, 1);
    
    block_mult(:, 1) = n(1:size_length);
    nbr_blocks = nbr_blocks + sum(block_mult(:,1));
    
    % Store block multiplicities corresponding to eigenvalue ev(j).
    eigs{j} = block_mult;

end

% eigs now contain all info about our blocks. We can now create them.

% To store blocks
blocks = cell(nbr_blocks, 1);

% Keep track of nbr of blocks left to add
nbr_blocks_left = nbr_blocks;

% For each eigenvalue, add its jordanblocks to 'blocks'
for k = 1:length(ev)
    
    block_vector = eigs{k};
    
    for l = 1:length(block_vector)
        if(block_vector(l) ~= 0)

            nbr_blocks_left = nbr_blocks_left - 1;
            blocks{nbr_blocks_left + 1} = ev(k)*eye(l) + diag(ones(l-1,1),1);
            
        end
        
    end
    
    
end

J = blkdiag(blocks{:});

end

