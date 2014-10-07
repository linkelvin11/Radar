function q=marcum(a,t,prec)
% For (x,y) ~ N(0,a,1,1,0) marcum(a,t,prec) computes the 
%    probability that x^2+y^2 > t^2 using a method from 
%    Brennan & Reed, IEEE Trans. Info Theory, 1965, 
%    pp. 312-313
%
% Original MATLAB code obtained from P.F. Swaszek
% at URL http://www.ele.uri.edu/Courses/ele510/Marcum.m.
%
% Modified by M. Richards July 2002
% Modified again by M. Richards in March 2004 to avoid infinite loops
%   for very small probabiliities ( q < 1e-6), which occur for large values
%   of 'a', which occur when the SNR in the calling function
%   is large
%
% prec=precision of computation. prec=0.001 => compute Q to 0.1% accuracy.
%

ct = .5*t^2;
ca = .5*a^2;
et = exp(-ct);
g  = 1 - et;
k  = exp(-ca);
gnew = et;
q = 1 - g*k;

n=0;
while ( n < a*t/2) | ( g*k/(1-(.5*a*t/n)^2) > prec*q )
    if (q < 1e-6)
        return
    end
	n = n + 1;
	gnew = gnew * ct/n;
	g = g - gnew;
	
	k = k * ca/n;
	q = q - g*k;
end
