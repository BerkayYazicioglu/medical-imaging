function dataset = zubal_brain(disp)
    % Dataset of brain
    load brain.txt brain
    brain = normalize(brain,'range',[0,1]);
    if disp
        imshow(brain);
    end
    % Converting the dataset into 3D array
    pixels = 256;
    dim = size(brain);
    dataset = zeros(pixels,pixels,dim(1)*dim(2)/pixels^2);
    for i = 1:dim(2)/pixels-1
        for j = 1:dim(1)/pixels-1
            dataset(:,:,i*j) = brain((j-1)*pixels+1:j*pixels,(i-1)*pixels+1:i*pixels);
        end
    end
end