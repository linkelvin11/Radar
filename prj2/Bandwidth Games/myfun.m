
dirListing = dir('*.fig'); 
for d = 1:length(dirListing) 
fileName = dirListing(d).name; 
openfig(fileName); 
plot2svg(strcat(fileName,'.svg')); 
end