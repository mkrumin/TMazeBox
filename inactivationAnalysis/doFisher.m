function out = doFisher(data)

nTests = length(data.cc{1});
nr = data.nr;
nl = cellfun(@minus, data.nn, data.nr, 'UniformOutput', false);
for iTest = 1:nTests
    x = [nr{1}(iTest), nr{2}(iTest);...
        nl{1}(iTest), nl{2}(iTest)];
    [out.h(iTest), out.p(iTest), out.stats(iTest)] = fishertest(x);
end
