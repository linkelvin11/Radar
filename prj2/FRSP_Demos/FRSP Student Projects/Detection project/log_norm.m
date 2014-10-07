function xln=log_norm2(M,m,SD)
%
% Function to create a log-normal random column vector with specified
% linear scale mean (m) and standard deviation (SD).  The vector length is M.
%
%

% Convert linear scale mean, std dev to log-scale equivalents

slog = sqrt(log(exp(log(SD^2) - 2*log(m)) + 1));
mlog = log(m) - (slog^2)/2;

% Generate log-normal variates to use for cell means; first, generate
% log-space sequence
%
	xlog=slog.*randn(M,1)+mlog;
	mxlog=mean(xlog);
	sxlog=std(xlog);

%
% convert to log-normal; compute sample mean, std. dev. and
% compare to expected
%
 xln=exp(xlog);
	mxln=mean(xln);
	stdxln=std(xln);


% Show 100-bin histogram and also a theoretical log-normal
% curve with the desired parameters.

%[counts bins]=hist(xln,100);
%bar(bins,counts/sum(counts))
%xlabel('x')
%ylabel('density')
end

