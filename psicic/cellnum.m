function cellnum(indata,colour)
%Label all cells in the current image.
%
%   CELLNUM(DATA,COLOUR)

if ~exist('colour') || isempty(colour), colour = 'y';end

for cellcount = 1:length(indata.cells)
    text(indata.cells(cellcount).mid(1,floor(indata.parameters.npoints/2)),indata.cells(cellcount).mid(2,floor(indata.parameters.npoints/2)),num2str(cellcount),'Color',colour)
end