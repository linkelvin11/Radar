function [amp,del_index] = peakinterp(z)

%PEAKINTERP
%  peakinterp performs a quadratic interpolation of the
%  peak defined by three consecutive data values in the
%  three-element vector z. The middle value, z(2), must be
%  the largest element.  No checking is done of this
%  requirement.  Real data is assumed.  There is no
%  guarding against edge effects, either.
%  peakinterp returns the interpolated peak amplitude and
%  the interpolated peak location relative to the center
%  element in range bins, e.g. index = -0.4 means the
%  interpolated peak lies 0.4 samples to the "left" of the
%  center element, i.e. at sample #1.6, so to speak.
%
%  Mark Richards, February 1997

del_index = -0.5*(z(3) - z(1))/(z(1) - 2*z(2) + z(3));

k=del_index;
amp = 0.5*((k-1)*k*z(1) - 2*(k-1)*(k+1)*z(2) + (k+1)*k*z(3));
