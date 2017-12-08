function dataOut = fitDebiasedPC(data)

alpha = 0.1;
idx = data.finished & data.random;
optiStim = data.optiStim;
contrast = data.contrast;
behavior = data.behavior;
session = data.iSession;

idxNone = idx & ~optiStim(:, 1) & ~optiStim(:,2);
idxLeft = idx & optiStim(:, 1) & ~optiStim(:,2);
idxRight = idx & ~optiStim(:, 1) & optiStim(:,2);
idxBoth = idx & optiStim(:, 1) & optiStim(:,2);

idx = {idxNone; idxLeft; idxRight; idxBoth};
groupName = {'none'; 'left'; 'right'; 'both'};
LineStyle = {'.k'; '.r'; '.b'; '.m'};
pcLineStyle = {'k'; 'r'; 'b'; 'm'};
nGroups = length(idx);
nSessions = max(session);

%%
for iGroup = 1:nGroups
    for iSession = 1:nSessions
        idxSession = idx{iGroup} & session == iSession;
        cc{iGroup, iSession} = unique(contrast(idxSession));
        nn{iGroup, iSession} = nan(size(cc{iGroup, iSession}));
        nr{iGroup, iSession} = nan(size(cc{iGroup, iSession}));
        for iC = 1:length(cc{iGroup, iSession})
            idxC = idxSession & contrast == cc{iGroup, iSession}(iC);
            nn{iGroup, iSession}(iC) = sum(idxC);
            nr{iGroup, iSession}(iC) = sum(behavior(idxC) == 'R');
        end
        %     [pp{iGroup}, ci{iGroup}]= binofit(nr{iGroup}, nn{iGroup}, alpha);
    end
    ccAll{iGroup} = unique(contrast(idx{iGroup}));
    nnAll{iGroup} = nan(size(ccAll{iGroup}));
    nrAll{iGroup} = nan(size(ccAll{iGroup}));
    for iC = 1:length(ccAll{iGroup})
        idxC = idx{iGroup} & contrast == ccAll{iGroup}(iC);
        nnAll{iGroup}(iC) = sum(idxC);
        nrAll{iGroup}(iC) = sum(behavior(idxC) == 'R');
    end
    
end

%% Fit logit psychometric curve to pooled data

x0 = [0, 0.2, 0.1, 0.1];
figure
for iGroup = 1:nGroups
    ccFit = ccAll{iGroup};
    nnFit = nnAll{iGroup};
    nrFit = nrAll{iGroup};
    pOut{iGroup} = fminsearch(@logLikFun, x0);
    xx = -50:50;
    yy = pcCurve(xx, pOut{iGroup});
    plot(ccFit, nrFit./nnFit, LineStyle{iGroup}, 'MarkerSize', 20);
    hold on;
    plot(xx, yy, pcLineStyle{iGroup}, 'LineWidth', 3);
end

%% Fit logit psychometric curve session-by-session correcting for bias

% figure
for iGroup = 1:nGroups
    x0 = [pOut{iGroup}(1)*ones(1, nSessions), pOut{iGroup}(2:4)];
    ccFit = cc(iGroup, :);
    nnFit = nn(iGroup, :);
    nrFit = nr(iGroup, :);
    pOutBiased{iGroup} = fminsearch(@logLikBiasedFun, x0);
%     xx = -50:50;
%     yy = pcCurve(xx, pOut{iGroup});
%     plot(ccFit, nrFit./nnFit, LineStyle{iGroup}, 'MarkerSize', 20);
%     hold on;
%     plot(xx, yy, pcLineStyle{iGroup}, 'LineWidth', 3);
end

%% plot debiased results

figure;

for iGroup = 1:nGroups
    xx = -50:50;
    yy = pcCurve(xx, pOut{iGroup});
    plot(xx, yy, pcLineStyle{iGroup}, 'LineWidth', 3);
    hold on;
    yy = pcCurve(xx, [0, pOutBiased{iGroup}(end-2:end)]);
    plot(xx, yy, pcLineStyle{iGroup}, 'LineWidth', 1);
    xlim([-50 50]);
    ylim([0 1]);
    plot(xlim, [0.5 0.5], 'k:', [0 0], ylim, 'k:');
end
    

plotDebiased();

%%
dataOut = pOut;

    function logL = logLikFun(p)
        
        f = pcCurve(ccFit, p);
        logL = sum(log(f).*nrFit + log(1-f).*(nnFit-nrFit));
        logL = -logL; % find minimum, not maximum
        
    end

    function logL = logLikBiasedFun(pAll)
        
        logL = 0;
        for iS = 1:length(ccFit)
            p = [pAll(iS), pAll(end-2:end)];
            f = pcCurve(ccFit{iS}, p);
            logL = logL + sum(log(f).*nrFit{iS} + log(1-f).*(nnFit{iS}-nrFit{iS}));
        end
        logL = -logL; % find minimum, not maximum
        
    end

end

function out = pcCurve(x, p)

b = p(1);
s = p(2);
upL = p(3);
lowL = p(4);

out = (1-upL-lowL)./(1+exp(-s*(x-b)))+lowL;

end