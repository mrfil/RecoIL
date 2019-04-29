function coilRank = chooseCoilRank(sig,energyLevel)
sizeSig = size(sig);
nc = sizeSig(1);
for ii = 1:nc
    if ((sum(sig(1:ii) - sig(end))) > energyLevel*sum(sig - sig(end)))
        break;
    end
end
coilRank = ii;
end