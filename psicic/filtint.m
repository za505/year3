function outdata = filtint(indata)

channel = 1;

outdata = indata;

if indata.CellCount > 0
source = indata.SourceImage(:,:,channel);

for cellcount = 1:length(indata.cells)
    temp = roipoly(indata.SourceImage,indata.cells(cellcount).border(1,:),indata.cells(cellcount).border(2,:));
    temp2 = source(temp);
    mins(cellcount) = min(temp2(:));
    maxs(cellcount) = max(temp2(:));
    means(cellcount) = mean(temp2(:));
    clear temp temp2
end


%bad = find(mins > 1.4*min(indata.SourceImage(:)));
bad = find(means > 0.98*indata.parameters.contourlevel);

outdata.cells(bad) = [];
outdata.CellCount = length(outdata.cells);
end
