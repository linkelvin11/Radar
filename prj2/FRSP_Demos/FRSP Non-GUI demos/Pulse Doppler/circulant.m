function circ = circulant( column_1, num_cols)
%CIRCULANT                                            [jMc 2/89]
%    circulant(C,P) makes a circulant matrix with P columns
%        and first column equal to C.
%    circulant(C) makes a square (M X M) circulant matrix,
%        where M = length(C).
%  see also CONVOLM, CONVMTX, HANKEL, TOEPLITZ

if nargin == 1
   num_cols = length(column_1);
end
column_1 = column_1(:);
circ = zeros(length(column_1),num_cols);           % allocate [ circ ]
for i=1:num_cols
   circ(:,i) = git_rotate(column_1,i-1);
end
