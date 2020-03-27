%convert 3D wavelet transform tensor to vector
function output = coef_to_vec(coef, start)
    output = zeros(256*256*(10-start),1);
    for j = 1:10-start %for each image/layer
        layer = coef(:,:,j+start-1); %get jth frequency sub-band
        layer = layer(:); %convert to a vector
        startIndex = 256*256*(j-1) + 1;
        endIndex = 256*256*j; 
        output(startIndex:endIndex) = layer; %add to position in output
    end
end