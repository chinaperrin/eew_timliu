function [compFeatures,iE] = get_compFeatures_growingSnippet(sz,sy,sx,sr,snpLength,nsnp)
%Compute various waveform statistics on increasingly long signal snippets
%of vector sums of the three components X, Y, Z in three unit dimensions
%(acc, vel & dsp)

% nsnp                                  No. of snippets over which statistics are computed
% snpLength                             length of each snippet in sec
nspersnippet = snpLength*sr;            % No. of samples per snippet
iE           = (1:nsnp)'*nspersnippet;  % Indices of End of Waveform snippet
iE           = iE(iE<=numel(sz));
nsnp         = numel(iE);

dt   = 1/sr;
fnyq = sr/2;

pcc_zy = zeros(nsnp,1);
pcc_zx = zeros(nsnp,1);
pcc_yx = zeros(nsnp,1);

for isnp = 1:nsnp

    z = sz(1:iE(isnp));
    y = sy(1:iE(isnp));
    x = sx(1:iE(isnp));
    n = numel(z);
    
    % Pearsons linear correlation coefficient
    pcc_zy(isnp) = corr(z,y);
    pcc_zx(isnp) = corr(z,x);
    pcc_yx(isnp) = corr(y,x);

    %     % Fourier transforms
    %     Z = dt*fft(z,n);
    %     nf = n/2 + 1;
    %     f  = linspace(0,fnyq,nf);
    %     %plot(f,abs(Z(1:nf)),'k','lineWidth',2);

    %         fft_X = fft(X) # fft computing
        %         fft_Y = fft(Y)
        %         fft_Z = fft(Z)
        %         fft_vector_sum = fft(acc_vector_sum)
        %     
        %         fft_amp_X = np.abs(fft_X[arange(1, n/2)])
        %         fft_amp_Y = np.abs(fft_Y[arange(1, n/2)])
        %         fft_amp_Z = np.abs(fft_Z[arange(1, n/2)])
        %         fft_amp_vector_sum = np.abs(fft_vector_sum[range(1, n/2)])
        %     
        %         energy_x = np.sum(fft_amp_X**2) / len(fft_amp_X)
        %         energy_y = np.sum(fft_amp_Y**2) / len(fft_amp_Y)
        %         energy_z = np.sum(fft_amp_Z**2) / len(fft_amp_Z)
        %         energy_all = energy_x + energy_y + energy_z / 3 
        %         energy_acc_vector_sum = np.sum(fft_amp_vector_sum**2) / len(fft_amp_vector_sum)
    

end
compFeatures.pcc_zy = pcc_zy;
compFeatures.pcc_zx = pcc_zx;
compFeatures.pcc_yx = pcc_yx;



end