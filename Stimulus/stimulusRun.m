%% stimulusRun
%  y = stimulusRun(s, tx)
%
%  Interpolates stimulus s for given time index tx

%%
function y = stimulusRun(s, tx)

rm = tx - floor(tx); 
if rm == 0  % integer, so valid index
    y = s.x(:,tx);
else        % We're between indices, so linear interpolation
    y = s.x(:,tx-rm) * (1-rm) + s.x(:,tx+1-rm) * rm;
end