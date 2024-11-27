%%
% returns 3D magnetic field in Tesla
function B0 = magField(dataset,val,mode)
    dimensions = size(dataset);
    B0_shift = zeros(dimensions); 
    std = val/10^6;
    % chemical shift
    if mode == 3
        disp('Generating magnetic field with chemical shift');
        fat_locations = find(dataset == 4);
        glial_locations = find(dataset == 8);
        ms_locations = find(dataset == 10);

        B0_shift(fat_locations) = -3.35*10^-6;
        B0_shift(glial_locations) = -2*10^-6;
        B0_shift(ms_locations) = -10*10^-6;
    % random noise
    elseif mode == 2
        fprintf('Generating magnetic field with std = %.2f ppm random noise\n',std*10^6/val);
        std = 100*val/10^6;
    else
        fprintf('Generating magnetic field with std = %.2f ppm random noise\n',std*10^6/val);  
    end
    noise = double(-std+2*std.*rand(dimensions(1),dimensions(2)));
    B0 = val*(ones(dimensions)+repmat(noise,[1,1,dimensions(3)])+B0_shift); 
end