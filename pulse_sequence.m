%%
% inputs
% samples = 2D slice with T1,T2,T2*,PD parameters at each voxel
% M0 = initial longitudinal magnetization of selected 2D slice
% v0 = initial larmor frequency array of selected 2D slice
% params:
% TE = echo time in ms
% TR = repetition time in ms
% phi = initial phase in degrees
% alpha = tip angle in degrees
% Ts: data acquisition window in terms of ms
% G_max: maximum gradient strength in terms of G/mm
% bandwidth: bandwidth of RF pulse in kHz
% v_hat: carrier frequency of RF pulse in kHz
% base_mag: field strength of the device in G
function [readings,mags,images]= pulse_sequence(samples,M0,v0,params,echo)
    % initializations
    gama_hat = 4.258; % kHz/G
    T1 = samples(:,:,1);
    T2 = samples(:,:,2);
    T2star = samples(:,:,3);
    dims = size(samples);
    % parameters
    TE = params(1);
    TR = params(2);
    phi = params(3);
    alpha = params(4);
    Ts = params(5);
    G_max = params(6);
    bandwidth = params(7);
    v_hat = params(8);
    base_mag = params(9);
    % calculating reconstruction parameters 
    % since we have a cubic database,following calculations are simplified 
    N = dims(1); % since resolution is 1 mm
    T = Ts/N;  % ms
    du = 1/N;  % 1/mm
    dv = 1/N;  % 1/mm
    Gu = du/(gama_hat*T); % G/mm
    Gv = linspace(G_max,-G_max,N); % G/mm
    Tp = dv*0.5*N/(gama_hat*G_max); % ms
    % k-space traversal parameters
    % RF is a sinc with two side lobs
    tau_p = 6/bandwidth;
    maxTs = TE-tau_p/2;
    fprintf('Maximum data acquisition window: %f ms\n',maxTs);
    if Ts > maxTs
        disp('Given Ts exceeds allowable value, setting it to maximum available')
        Ts = floor(maxTs);
    end
    % relaxation equations after RF pulse with refocusing
    syms t Mz0p
    % Magnetization values for different echo sequences 
    Mxy = @(t,Mz0p) Mz0p.*sind(alpha).*exp(-(1i*2*pi*v0*t-deg2rad(phi))).*exp(-t./T2); 
    Mz = @(t,Mz0p) M0.*(1-exp(-t./T1))+ Mz0p.*cosd(alpha).*exp(-t./T1);
    fxy = @(t,Mz0p) Mz0p.*sind(alpha).*exp(-t./T2)*exp(deg2rad(phi));
    % readout time samples
    t0_adc = TE-Ts/2-tau_p;
    t1_adc = TE+Ts/2-tau_p;
    t_adc = linspace(t0_adc,t1_adc,N);
    if echo
        Mxy = @(t,Mz0p) Mz0p.*sind(alpha).*exp(-(1i*2*pi*v0*t-deg2rad(phi))).*exp(-t./T2star); 
        Mz = @(t,Mz0p) M0.*(1-exp(-t./T1))+ Mz0p.*cosd(alpha).*exp(-t./T1);
        fxy = @(t,Mz0p) Mz0p.*sind(alpha).*exp(-t./T2star)*exp(deg2rad(phi));
        % readout time samples
        t0_adc = TE-Ts/2-3*tau_p/2;
        t1_adc = TE+Ts/2-3*tau_p/2;
        t_adc = linspace(t0_adc,t1_adc,N);
        disp('Using gradient echoes');
    else
        disp('Using spin echoes');
    end
    % coil reading arrays 
    % k space matrices
    coilx = zeros(N);
    coily = zeros(N);
    coilz = zeros(N);
    offset = floor(N/2);
    % iterating over v axis on frequency domain
    initMz = M0;    
    fprintf('\t---iterating on k-space---\n');
    u_pos = (-offset:offset)';
    v1_pos = 0:N-1;
    v2_pos = -offset:offset;
    exp_u = zeros(N,length(t_adc));
    exp_v1 = exp_u;
    exp_v2 = exp_v1;
    % actual modulation of the signal
    % if demod and received match they should cancel out 
    received = zeros(dims(1),dims(2),N);
    % demodulation at receiver size
    demod = exp(1i*2*pi*gama_hat*base_mag*t_adc);
    % preprocessing variables to increase speed
    for i = 1:N
        received(:,:,i) = exp(-1i*2*pi*v0*t_adc(i));
        exp_u(:,i) = exp(-1i*2*pi*gama_hat*Gu*u_pos*t_adc(i));
        exp_v1(i,:) = exp(-1i*2*pi*gama_hat*Gv(i)*Tp*v1_pos);
        exp_v2(i,:) = exp(-1i*2*pi*gama_hat*Gv(i)*Tp*v2_pos);
    end
    for i = 1:N
        % iterating over u axis on frequency domain
        if mod(i,40)==1
            fprintf('\t---k-space line:%d---\n',i);
        end
        for j = 1:N
            % single RF pulse reading
            received_x = received(:,1:N,j);
            received_y = received(:,N+1:2*N,j);
            received_z = received(:,2*N+1:3*N,j);
            f = fxy(t_adc(i),initMz);
            coilx(i,j) = sum(demod(j)*received_x.*sum(f(:,1:N).*(exp_u(:,j)*exp_v1(i,:)),'all'),'all');
            coily(i,j) = sum(demod(j)*received_y.*sum(f(:,N+1:2*N).*(exp_u(:,j)*exp_v1(i,:)),'all'),'all');
            coilz(i,j) = sum(demod(j)*received_z.*sum(f(:,2*N+1:3*N).*(exp_u(:,j)*exp_v2(i,:)),'all'),'all');
            if (isnan(coilx(i,j))||isnan(coily(i,j))||isnan(coilz(i,j)))
                disp('nan error');
                return
            end
        end
        initMz = Mz(TR,initMz); % for next iteration
    end
    % outputs
    readings = [coilx,coily,coilz]; % actual signals placed wrt pulse sequences
    mags = Mxy(TE,initMz);  % magnetization matrices that we want to reconstruct
    imx = ifftshift(ifft2(coilx));
    imx_temp = imx(1:floor(N/2),:);
    imx(1:ceil(N/2),:) = imx(ceil(N/2):end,:);
    imx(ceil(N/2)+1:end,:) = imx_temp;
    imx = rot90(imx,3);
    imy = ifftshift(ifft2(coily));
    imy_temp = imy(1:floor(N/2),:);
    imy(1:ceil(N/2),:) = imy(ceil(N/2):end,:);
    imy(ceil(N/2)+1:end,:) = imy_temp;
    imy = rot90(imy,3);
    imz = rot90(ifftshift(ifft2(coilz)),3);
    images = [imx,imy,imz]; % reconstructed images
end
