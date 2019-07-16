function index = samplePD(PD)

cPD = cumsum(PD);
r = rand();
index = find(r<=cPD,1,'first');

if isempty(index)
    PD
    cPD
    r
end

end