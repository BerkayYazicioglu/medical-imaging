%%
% inputs
% slice: the selected 2D slice consisting of MR tissue parameters
% B_slice: 2D magnetic field corresponding to the slice parameter selected
% output
% M0: bulk magnetization of selected parameters
function M0 = bulkMag(slice,B_slice)
    h = 6.626068e-34*10^3;             % mm^2*kg/ms
    k = physconst('Boltzmann');        % J/K
    gama_hat = 4.258;                  % kHz/G
    T = 298;                           % K, room temperature
    M0 = B_slice .* slice(:,:,4);      % slice(:,:,4) = PD
    M0 = M0 * gama_hat^2 * h^2/(4*k*T);% actual value
end