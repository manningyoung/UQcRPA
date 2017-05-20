function [ Fij , gindex , ekin] = calc_aux(unk,kpoints,kweights,qpoints,Nfft,gvecs,iband,jband)
% Calculates the auxiliary function F_{i,j} using the FFT method.
%
% Also returns gindex for mapping G-vectors in gvecs to the corresponding Fij.
% And ekin = |q+G|^2 (useful for plotting/debugging).
%
% FIXME: assumes spatially symmetric FFT (Nfft(1)=Nfft(2)=Nfft(3)=N).

fprintf('Calculating F_{%d,%d}...\n',iband,jband);

if (Nfft(1) == Nfft(2)) && (Nfft(2) == Nfft(3))
    N = Nfft(1);
else
    error('FFT not spatially symmetric - check Nfft.')
end
scale = 1/N^3;

% Setup the real and reciprocal space intervals along one dimension
% so we can determine the coordinates (G-vectors) of the FFT output.
% The real space interval is normalised so that dg=1.
x_total = 1;
dx = x_total/N;
g_max = 1/(2*dx);
dg = 1/x_total;
if mod(N,2)==0 % even FFT
    x = -x_total/2:dx:x_total/2-dx;
    g = -g_max:dg:g_max-dg;
else % odd FFT
    x = -x_total/2:dx:x_total/2;
    g = -g_max:dg:g_max;
end

for iq = 1:length(qpoints)
    
    fprintf('Doing FFTs at q-point %d of %d...',iq,length(qpoints));
    
    % Determine which G-vectors the linear indices map to in reciprocal space
    for ig = 1:N^3
        [ix,iy,iz] = ind2sub([N N N],ig);
        gvec = [g(ix) g(iy) g(iz)];
        [found,idx] = ismember(gvec,gvecs{iq},'rows');
        if found
            gindex{iq}(idx) = ig; % maps an index in gvecs to an Fij
        end
        ekin{iq}(ig) = norm(qpoints(iq,:)+gvec)^2; % kinetic energy |q+G|^2
    end
    
    for ik = 1:length(kpoints) % FFT each k-point
        
        % Determine the index of the k+q point
        kqpoint = mod(kpoints(ik,:)+qpoints(iq,:),1); % k+q (mod BZ)
        kqpoint = round(kqpoint,1); % avoids round-off error - CAUTION
        [found,ikq] = ismember(round(kqpoint,1),kpoints,'rows');
        if ~found
            error('k+q point is not an existing k-point: [%f %f %f].',kqpoint(1),kqpoint(2),kqpoint(3));
        end
        
        % Create the data array
        data1d = conj(unk{ik,iband}.*unk{ikq,jband});
        
        % Transfer data to 3D array for FFT
        for ir = 1:length(data1d)
            [ix,iy,iz] = ind2sub([N N N],ir);
            data3d(ix,iy,iz) = data1d(ir);
        end

        % Do the FFT
        % Note: usage of fftshift is tricky - see e.g. MathWorks forums
        fftdata3d = scale*fftshift(fftn(ifftshift(data3d),[N N N]));
        
        % Since we know which G-vector the linear index maps to,
        % we can transfer the fftdata to a 1D array.
        % Note: ind2sub <-> sub2ind works with ix varying fastest.
        for iz = 1:N
            for iy = 1:N
                for ix = 1:N
                    ig = sub2ind([N N N],ix,iy,iz);
                    fftdata1d{iq,ik}(ig) = fftdata3d(ix,iy,iz);
                end
            end
        end
        
    end
    
    fprintf('Done.\n');
    
end

% Do the k-space integration (sum) for each G-vector and q-point

fprintf('Summing over k-points...');

for iq = 1:length(qpoints)
    for ig = 1:N^3
        ksum = 0;
        for ik = 1:length(kpoints)
            ksum = ksum + fftdata1d{iq,ik}(ig)*kweights(ik);
        end
        Fij{iq}(ig) = ksum/length(kpoints); % normalise by 1/N_k
    end
end

fprintf('Done.\n');

end