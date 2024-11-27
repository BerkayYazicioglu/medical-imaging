%%
% inputs
% dr: 1x3 slice thickness in terms of mm
% n: 1x3 vector in terms of mm, normal to the selected plane
% p: 1x3 scalar in terms of mm, offset of selected plane
% T: homogenous scalar magnetic field value in terms of Tesla
% B0: 3D magnetic field in terms of Tesla 
% dataset: 3D discrete tissue array 
% maxG: maximum value of available gradients
% outputs
% out: selected 2D slice on x,y,z coils
% v0: initial 2D larmor frequencies of selected slices
% dv: bandwidth of RF pulse in kHz
% v_hat: center frequency of RF pulse in kHz
function [out,v0,dv,v_hat] = slice_selection(dataset,T,B0,dr,n,p,maxG)
    T = 10^4*T;           % G
    B0 = 10^4*B0;         % G
    dims = size(dataset); % X,Y,Z
    gama_hat = 4.258;     % kHz/G
    % offsets
    x = floor(dims(1)/2);  
    y = floor(dims(2)/2);  
    % finding gradients for selecting given slice
    G = n.*maxG; % G/mm
    fprintf('Slice selection gradients: Gx=%.2f G/cm | Gy=%.2f G/cm | Gz=%.2f G/cm\n',...
                              G(1)*10,G(2)*10,G(3)*10);
    % finding bandwidth of selected slice thickness
    dv = gama_hat*(G*dr');  % kHz
    fprintf('RF pulse bandwidth: %.2f kHz\n',dv);
    % finding center frequency
    v_hat = gama_hat*(T+G*p'); % kHz
    fprintf('RF pulse frequency: %.3f MHz\n',v_hat/10^3);
    % 3D excited larmor frequency array
    v = zeros(dims);
    coilx = zeros(dims(2),dims(3));
    coily = zeros(dims(3),dims(1));
    coilz = zeros(dims(1),dims(2));
    %------------
    xsurf = zeros(dims(2),dims(3));
    ysurf = zeros(dims(3),dims(1));
    zsurf = zeros(dims(1),dims(2));
    %------------
    % forming V0 matrices seen by coils
    v0x = gama_hat*T*ones(size(coilx));
    v0y = gama_hat*T*ones(size(coily));
    v0z = gama_hat*T*ones(size(coilz));
    ind = 1;
    for i = 1:dims(1)
        for j = 1:dims(2)
            for k = 1:dims(3)
                pos = [(j-1)-x ; y-i+1 ; k-1];
                v(i,j,k) = gama_hat*(B0(i,j,k)+G*pos);
                % getting data that falls into selected bandwidth
                if ((v(i,j,k) < v_hat+dv/2) && (v(i,j,k) > v_hat-dv/2))
                    coilx(j,k) = dataset(i,j,k);
                    coily(i,k) = dataset(i,j,k);
                    coilz(i,j) = dataset(i,j,k);
                    v0x(j,k) =  gama_hat*(B0(i,j,k));
                    v0y(i,k) =  gama_hat*(B0(i,j,k));
                    v0z(i,j) =  gama_hat*(B0(i,j,k));
                    if ind <= dims(1)*dims(2)
                        xsurf(ind) = pos(1);
                        ysurf(ind) = pos(2);
                        zsurf(ind) = pos(3);
                    end
                    ind = ind+1;
                end
            end
        end
    end
    % correcting orientations
    coily = rot90(flip(coily));
    v0y = rot90(flip(v0y));
    coilx = rot90(coilx);
    v0x = rot90(v0x);
    % Plotting selected slice over larmor frequency field on 3D 
    [X,Y,Z] = meshgrid([-x:x],[y:-1:-y],[0:dims(3)-1]);
    if ind > dims(1)*dims(2)+1
        [xsurf,ysurf] = meshgrid([-x:x],[y:-1:-y]);
        if n(3) == 0
            zsurf = [];
            if n(2) == 0
                ysurf = [];
                xsurf = p(1);
            else
                ysurf = p(2)-(n(1)*(xsurf-p(1)))/n(2);
                zsurf = meshgrid([0:dims(3)-1])';
            end
        else
            zsurf = p(3)-(n(1)*(xsurf-p(1))+n(2)*(ysurf-p(2)))/n(3);
        end
    end
    % plotting selected slice,larmor frequency and 3D dataset together
    figure;
    Ds = smooth3(dataset);
    hiso = patch('XData',X,'YData',Y,'ZData',Z,...
           isosurface(X,Y,Z,Ds,5),'FaceColor',[1,.75,.65],'EdgeColor','none');
    hold on;
    isonormals(Ds,hiso)
    hcap = patch('XData',X,'YData',Y,'ZData',Z,...
           isocaps(X,Y,Z,dataset,5),'FaceColor','interp','EdgeColor','none');
    hold on;
    daspect([1,1,1])
    view(45,45);
    light('Position',[-1 1 1],'Style','local')
    lightangle(45,45);
    lighting flat
    hcap.AmbientStrength = 0.6;
    hiso.SpecularColorReflectance = 0;
    hiso.SpecularExponent = 50;hold on;
    s = slice(X,Y,Z,v,-x,y,0);hold on;
    set(s,'EdgeColor','none','SpecularStrength',0);hold on;
    slice(X,Y,Z,v,xsurf,ysurf,zsurf);colorbar;
    xlabel("x");
    ylabel("y");
    zlabel("z");
    xlim([-x,x]);
    ylim([-y,y]);
    zlim([0,dims(3)-1]);
    title({'Larmor frequencies of 3D dataset',''});
    % concatenating coil outputs
    out = [coilx,coily,coilz];
    v0 = [v0x,v0y,v0z];
end 
    