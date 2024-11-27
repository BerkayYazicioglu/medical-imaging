% Each element on the discrete dataset is replaced with 
% corresponding [T1,T2,T2*,PD] values of the tissue in this order
function tissues = mcgill_decoder(slice)
    MRvals = [[   0   0   0  0   ];... % background
              [2569 329  58  1   ];... % CSF
              [ 833  83  69  0.86];... % gray matter
              [ 500  70  61  0.77];... % white matter
              [ 350  70  58  1   ];... % fat
              [ 900  47  30  1   ];... % muscle/skin
              [2569 329  58  1   ];... % skin
              [   0   0   0  0   ];... % skull
              [ 833  83  69  0.86];... % glial matter
              [ 500  70  61  0.77];... % meat
              [ 752 237 204  0.76]];   % ms lesion
    dims = size(slice);
    tissues = zeros(dims(1),dims(2),4);
    for i = 1:dims(1)
        for j = 1:dims(2)
            tissues(i,j,:) = MRvals(slice(i,j)+1,:);
        end
    end
end