
dirListing = dir('*.fig'); 
for d = 1:length(dirListing) 
fileName = dirListing(d).name; 
openfig(fileName); 
plot2svg('fileName.svg'); 
end