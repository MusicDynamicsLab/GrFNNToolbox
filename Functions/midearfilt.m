%% midearfilt
%  z = midearfilt(s,varargin)
%
%  Creates a transfer function for the middle ear and filters the data in
%  vector s with it to output filtered data z. Optional second argument is 
%  sampling frequency, else sampling frequency is assumed to be 100,000 Hz.

function z = midearfilt(s,varargin)

% Sampling Frequency (this code should work/stable for Fs > 70e3)
Fs = 100e3;
if nargin==2, Fs=varargin{1};end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cat_or_human = 2; % Choose 1 for cat and 2 for human


tdres = 1/Fs;
fp    = 1e3;  % prewarping frequency 1 kHz 
TWOPI = 2*pi;
C     = TWOPI*fp/tan((TWOPI/2)*fp*tdres);

if cat_or_human == 1  % cat
    
    %%% cat middle ear function %%%%%%%%%%%%
    
    m11 = C/(C + 693.48); 
    m12 = (693.48 - C)/C; 
    m21 = 1/(power(C,2) + 11053*C + 1.163e8);
    m22 = -2*power(C,2) + 2.326e8;
    m23 = power(C,2) - 11053*C + 1.163e8;
    m24 = power(C,2) + 1356.3*C + 7.4417e8;
    m25 = -2*power(C,2) + 14.8834e8;
    m26 = power(C,2) - 1356.3*C + 7.4417e8; 
    m31 = 1/(power(C,2) + 4620*C + 909059944);
    m32 = -2*power(C,2) + 2*909059944; 
    m33 = power(C,2) - 4620*C + 909059944; 
    m34 = 5.7585e5*C + 7.1665e7; 
    m35 = 14.333e7; 
    m36 = 7.1665e7 - 5.7585e5*C;

    sos = [1 -1 0 1 m11*m12 0;...
           m21*m24 m21*m25 m21*m26 1 m21*m22 m21*m23;...
           m31*m34 m31*m35 m31*m36 1 m31*m32 m31*m33];
    %%% end of cat middle ear function

elseif cat_or_human == 2 % human

    %%% Human middle ear function %%%%%%%%%%%%%%%%%%%%

    m11 = 1/(power(C,2) + 5.9761e+003*C + 2.5255e+007);
    m12 = (-2*power(C,2) + 2*2.5255e+007);
    m13 = (power(C,2) - 5.9761e+003*C + 2.5255e+007);
    m14 = (power(C,2) + 5.6665e+003*C);
    m15 = -2*power(C,2);
    m16 = (power(C,2) - 5.6665e+003*C);
    
    m21 = 1/(power(C,2) + 6.4255e+003*C + 1.3975e+008);
    m22 = (-2*power(C,2) + 2*1.3975e+008);
    m23 = (power(C,2) - 6.4255e+003*C + 1.3975e+008);
    m24 = (power(C,2) + 5.8934e+003*C + 1.7926e+008);
    m25 = (-2*power(C,2) + 2*1.7926e+008);
    m26 = (power(C,2)-5.8934e+003*C+1.7926e+008);
    
    m31 = 1/(power(C,2) + 2.4891e+004*C+1.2700e+009);
    m32 = (-2*power(C,2) + 2*1.2700e+009);
    m33 = (power(C,2) - 2.4891e+004*C + 1.2700e+009);
    m34 = (3.1137e+003*C + 6.9768e+008);
    m35 = 2*6.9768e+008;
    m36 = (-3.1137e+003*C + 6.9768e+008);    
    
    sos = [m11*m14 m11*m15 m11*m16 1 m11*m12 m11*m13;...
           m21*m24 m21*m25 m21*m26 1 m21*m22 m21*m23;...
           m31*m34 m31*m35 m31*m36 1 m31*m32 m31*m33];

    %%% end of Human middle ear function %%%%%%%%%%%%%%
end


[num, den] = sos2tf(sos, 1); % gain of 1 used
z = filter(num, den, s);
