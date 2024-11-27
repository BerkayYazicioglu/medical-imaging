function dataset = mcgill_brain(normal)
    if normal
        temp = double(loadminc('phantom_1.0mm_normal_crisp.mnc'));
    else
        temp = double(loadminc('phantom_1.0mm_msles2_crisp.mnc'));
    end
    dim = size(temp);
    % zero-padding to have equal dimensions
    offset = dim(1)-dim(2);
    dataset = zeros(dim(1),dim(1),dim(1));
    dataset(:,offset/2:dim(1)-offset/2-1,offset/2:dim(1)-offset/2-1) = temp;
    dim = size(dataset);
    if mod(dim(1),2) == 0
        temp = dataset;
        dataset = zeros(size(temp)+1);
        dataset(1:dim(1),1:dim(2),1:dim(3)) = temp;
    end
end


