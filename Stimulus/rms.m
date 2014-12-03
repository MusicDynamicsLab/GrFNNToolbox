function RMS = rms(x)

if length(x) == sum(size(x))-1 && size(x,2) > 1
    x = x.';
end

RMS = zeros(1,size(x,2));

for i = 1:length(RMS)
    RMS(i) = sqrt(sum(real(x(:,i)).^2)/size(x,1));
end