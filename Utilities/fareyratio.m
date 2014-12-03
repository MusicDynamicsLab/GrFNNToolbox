function [num, denom, level, error] = fareyratio(fractions, pertol)

if nargin < 2
   pertol = .01;
end

num = [];
denom = [];
level = [];
error = [];

for f = fractions

   if f > 1
      frac = 1/f;
   else
      frac = f;
   end

   frac*pertol;

   ln = 0;
   ld = 1;

   rn = 1;
   rd = 1;

   l = 1;

   if (abs(frac - ln/ld) <= frac*pertol)
      n = ln;
      d = ld;
      e = abs(frac - ln/ld);

   elseif (abs(frac - rn/rd) <= frac*pertol)
      n = rn;
      d = rd;
      e = abs(frac - rn/rd);

   else

     cn = ln+rn;
     cd = ld+rd;
     l  = l + 1;

     while (abs(frac - cn/cd) > frac*pertol)

	if frac > cn/cd
	   ln=cn;ld=cd;
	else
	   rn=cn;rd=cd;
	end

	cn = ln+rn;
	cd = ld+rd;
        l  = l + 1;

     end

     n = cn;
     d = cd;
     e = abs(frac - cn/cd);
   end

   if f > 1
      tmp = n;
      n = d;
      d = tmp;
   end
%   e = abs(f - n/d);

   num   = [num   n];
   denom = [denom d];
   level = [level l];
   error = [error e];
end
