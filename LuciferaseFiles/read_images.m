%% 
% Analysis code from "Coordination of robust single cell rhythms in  the Arabidopsis
% circadian clock via spatial waves of gene expression" (BioRxiv)
% doi https://doi.org/10.1101/208900s 
%
% Created by Mark Greenwood (U. of Cambridge). 09/02/2018.
%
% Reads all .tiff images in a directory into a 3d array
%%

function imageStack=read_images(fol)

contents = dir(sprintf('%s/*.tif', fol));

for i = 1:numel(contents)
    filename = contents(i).name;
    [~, name] = fileparts(filename);
    imageStack(:, :, i) = imread(sprintf('%s/%s.tif', fol, name)); 
end

end


