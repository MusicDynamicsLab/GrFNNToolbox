%% Interpolate stimulus value for given t
function y = stimulusRun(s, tx)

% tx is an index into stimulus data
% First check if index is integer or not
rm = tx - floor(tx); % NOTE don't name this 'rem'! That's a func and will slow things!
if rm == 0
    % integer, so valid index
    y = s.x(:,tx);
else
    % We're between indices, so linear interpolation
    y = s.x(:,tx-rm) * (1-rm) + s.x(:,tx+1-rm) * rm;
    % (t-rm) & (t+1-rm) are much faster than floor(t) & ceil(t)
end