% Reads all .tiff images in a directory into a 3d array

function imageStack=read_images(fol)

contents = dir(sprintf('%s/*.tif', fol));

for i = 1:numel(contents)
    filename = contents(i).name;
    [~, name] = fileparts(filename);
    imageStack(:, :, i) = imread(sprintf('%s/%s.tif', fol, name)); 
end

end


