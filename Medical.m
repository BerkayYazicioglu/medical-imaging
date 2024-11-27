clear, clc, close all;

%------------ IMAGING PARAMETERS -------------
% dataset selection (0,1) 0->ms lesion 1->normal
data_select = 0;
% inhomogeneity mode (1,2,3) 1->homogeneous 2->noisy 3->chemical shift
mode = 1;
% echo mode (0,1) 0->spin 1->gradient
echo = 0;
% slice indices
coronal_ind = 0; % for coronal selection only
sagittal_ind = 0; % for sagittal selection only
axial_ind = 108; % for axial selection only
% for arbitrary selection only
vector = [0 1 2];
normal = vector./sqrt(sum(vector.^2));
point = [0 0 108]; % offset parameters of slice
% sequence parameters
TS = 50;   % ms
TE = 200;  % ms
TR = 6000; % ms
% coronal , sagittal , axial , arbitrary
select = [0 0 0 1];   
%---------------------------------------------

% Getting the dataset 
dataset = mcgill_brain(data_select);
pixels = size(dataset);
N = pixels(1);
%% coronal
if select(1)
    disp('---getting coronal image---');
    [readings,mags,images] = MRI(3,dataset,90,[0 1 0],[0 coronal_ind 0],TS,TE,TR,mode,echo);
    % getting different readings 
    coronal_signal = rescale(readings(:,1:N));
    coronal_mags = rescale(mags(:,1:N));
    coronal_image = rescale(images(:,1:N));
    figure; ax=axes;
    montage({coronal_mags,coronal_signal,coronal_image},...
    	  'BorderSize',[0 10],'BackgroundColor',[1 1 1],'Size',[1 3]);
    title(ax,'Coronal MRI with M_{xy}(TE) values, k-space readings and reconstructed image');
end
%% sagittal
if select(2)
    disp('---getting sagittal image---');
    [readings,mags,images] = MRI(3,dataset,90,[1 0 0],[sagittal_ind 0 0],TS,TE,TR,mode,echo);
    % getting different readings 
    sagittal_signal = rescale(readings(:,N+1:2*N));
    sagittal_mags = rescale(mags(:,N+1:2*N));
    sagittal_image = rescale(images(:,N+1:2*N));
    figure; ax=axes;
    montage({sagittal_mags,sagittal_signal,sagittal_image},...
    	  'BorderSize',[0 10],'BackgroundColor',[1 1 1],'Size',[1 3]);
    title(ax,'Sagittal MRI with M_{xy}(TE) values, k-space readings and reconstructed image');
end
%% axial
if select(3)
    disp('---getting axial image---');
    [readings,mags,images] = MRI(3,dataset,90,[0 0 1],[0 0 axial_ind],TS,TE,TR,mode,echo);
    % getting different readings 
    axial_signal = rescale(readings(:,2*N+1:3*N));
    axial_mags = rescale(mags(:,2*N+1:3*N));
    axial_image = rescale(images(:,2*N+1:3*N));
    figure; ax=axes;
    montage({axial_mags,axial_signal,axial_image},...
    	  'BorderSize',[0 10],'BackgroundColor',[1 1 1],'Size',[1 3]);
    title(ax,'Axial MRI with M_{xy}(TE) values, k-space readings and reconstructed image');
end
%% arbitrary
if select(4)
    disp('---getting arbitrary image---');
    [readings,mags,images] = MRI(3,dataset,90,normal,point,TS,TE,TR,mode,echo);
    % getting different readings 
    coronal_signal = rescale(readings(:,1:N));
    coronal_mags = rescale(mags(:,1:N));
    coronal_image = rescale(images(:,1:N));
    sagittal_signal = rescale(readings(:,N+1:2*N));
    sagittal_mags = rescale(mags(:,N+1:2*N));
    sagittal_image = rescale(images(:,N+1:2*N));
    axial_signal = rescale(readings(:,2*N+1:3*N));
    axial_mags = rescale(mags(:,2*N+1:3*N));
    axial_image = rescale(images(:,2*N+1:3*N));
    % projection scaling and interpolation for correct arbitrary imaging
    xzscale = 1/cosd(atand(normal(3)/normal(2)));
    xyscale = 1/cosd(atand(normal(2)/normal(3)));
    yzscale = 1/cosd(atand(normal(3)/normal(1)));
    yxscale = 1/cosd(atand(normal(1)/normal(3)));
    zyscale = 1/cosd(atand(normal(1)/normal(2)));
    zxscale = 1/cosd(atand(normal(2)/normal(1)));
    scales = [xzscale xyscale yzscale yxscale zyscale zxscale];
    for i = 1:length(scales)
        if scales(i) == Inf
            scales(i) = 1;
        end
    end
    % scaling to project from x y z readings
    sag_scale = affine2d([scales(6) 0 0; 0 scales(3) 0; 0 0 1]);
    cor_scale = affine2d([scales(5) 0 0; 0 scales(1) 0; 0 0 1]);
    axi_scale = affine2d([scales(4) 0 0; 0 scales(2) 0; 0 0 1]);
    sagittal_proj = imwarp(sagittal_image,sag_scale,'cubic'); 
    coronal_proj = imwarp(coronal_image,cor_scale,'cubic');  
    axial_proj = imwarp(axial_image,axi_scale,'cubic'); 
    % plotting
    figure; ax=axes;
    montage({coronal_mags,coronal_signal,coronal_image,coronal_proj,...
             sagittal_mags,sagittal_signal,sagittal_image,sagittal_proj,... 
             axial_mags,axial_signal,axial_image,axial_proj},...
    	  'BorderSize',[5 10],'BackgroundColor',[1 1 1],'Size',[3 4]);
    title(ax,...
    'Arbitrary slice MRI with M_{xy}(TE) values, k-space readings,reconstructed and projected images');
end