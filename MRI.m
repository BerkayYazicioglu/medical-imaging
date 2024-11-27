%%
% dataset: NxNxN discrete values
function [readings,mags,images] = MRI(T,dataset,alpha,n,P,TS,TE,TR,mode,echo)
    gama_hat = 4.258; % kHz/G
    G_max = 0.5;       % G/mm
    % producing magnetic field
    B0 = magField(dataset,T,mode);
    % slice selection
    dr = [1 1 1]; % mm
    [slices,v0s,dv,v_hat] = slice_selection(dataset,T,B0,dr,n,P,G_max);
    
    % decoding tissue parameters on selected slice  
    samples = mcgill_decoder(slices);
    % M0 magnetization of selected slice
    B0s = v0s/gama_hat;
    M0 = bulkMag(samples,B0s);
    
    % forming pulse sequences and acquiring data
    TS_min = 1/(gama_hat*G_max);
    fprintf('Minimum data acquisition window: %f ms\n',TS_min);
    phi = 0;   % degrees
    params = [TE,TR,phi,alpha,TS,G_max,dv,v_hat,T*10^4];
    [readings,mags,images] = pulse_sequence(samples,M0,v0s,params,echo);
    readings = abs(readings);
    mags = abs(mags);
    images = abs(images);
    disp('Imaging finished')
end
