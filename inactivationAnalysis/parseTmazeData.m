function out = parseTmazeData(data)

nSessions = length(data);
for iSession = 1:nSessions
    SESSION = data(iSession).SESSION;
    EXP = data(iSession).EXP;
    if (iSession == 1)
        out = getSessionRes(EXP, SESSION);
    else
        out(iSession) = getSessionRes(EXP, SESSION);
    end
end
out = out(:);
