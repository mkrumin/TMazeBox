function changePatchColor(h)

% this is a recursive function, which will change the properties of all the
% patches in the figure/axis h

if isequal(h.Type, 'patch')
    h.EdgeColor = 'k';
    % add other properties to change here
else
    ch = h.Children;
    for iCh = 1:length(ch)
        changePatchColor(ch(iCh));
    end
end