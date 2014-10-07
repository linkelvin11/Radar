function [ y ] = randgmix( a,m,s2,n )
%RANDGMIX Generates a block of iid random variables with a Gaussian mixture
%distribution.

% a - weight of each distribution
% mean - mean of each distribution
% s2 - variance of each distribution
% n - number of samples to be generated

if (nargin<4)
   if (nargin == 3)
       n = 1e5;
       disp('Only 3 inputs recieved; Setting n = 1e5');
   else
       disp('not enough inputs');
   end
end

if((size(a,1) ~= size(m,1))||(size(a,1)~=size(s2,1)))   %check input sizes  
    disp('you are stupid');
    y = 0;
    return
end

if (sum(a) ~= 1)      % check that probabilites add to 1
    disp('you are even more stupid');
    y = 0;
    return
end


if isreal(m)
%     disp('signal is real')
    y = arrayfun(@(a,m,s2) m + sqrt(s2)*randn(a*n,1),a,m,s2,'UniformOutput',false);
    y = cell2mat(reshape(y,[],1));
else
%     disp('signal is cplx')
    y = arrayfun(@(a,m,s2) m + sqrt(s2/2)*randn(a*n,1) + ...
                               1i*sqrt(s2/2)*randn(a*n,1),a,m,s2,'UniformOutput',false);
    y = cell2mat(reshape(y,[],1));
end

end

%  