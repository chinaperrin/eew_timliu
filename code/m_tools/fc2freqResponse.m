function fc2freqResponse(fc,fOrder,figNum)

% Computes and plots combined and individual frequency responses for all 
% passbands of a filterbank with corner frequency matrix <fc>.
%
% mam, 121221


fprintf(1,'\n-----------------------------------------------------------\n')
addpath(genpath('../../../matlab/sac/'))

ftSize = 15;
sr     = 100;
fNyq   = sr/2;
wNyq   = fNyq*2*pi;

if (nargin < 2)
    figNum = 555;
end

nb = numel(fc(:,1));

Y            = cell(nb,1);
legendString = cell(nb+1,1);


% % Load example waveform
% velFullName = '/scratch/memeier/data/socal/scsn_100101_111231/cms/14718764/14718764.CI.BAR.HHZ.sac';
% 
% [out] = read_sac_trace3(velFullName);
%  72     sr    = out.sr;
%  73     t     = out.t;
%  74     S.raw = out.sraw/100;  
%  
% [vel,meta]  = read_any_trace(velFullName);
% tvel        = meta.t;
% sr          = meta.sr;


% B. Filter   -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
figure(figNum)
clf; box on; hold on

%y_sum    = zeros(length(vel),1);
C        = linspace(0,1,nb);
resp_env = zeros(4096,1);

for ib = 1:nb
    fLow = fc(ib,1);
    fUp  = fc(ib,2);
    
    wLo         = fc(ib,1)*2*pi;                        % Circular corner frequencies
    wUp         = fc(ib,2)*2*pi;
    [B,A]       = butter(fOrder,[wLo/wNyq wUp/wNyq],'bandpass');              % Compute coefficients
    %[y ff]      = freqz(B,A,4096,sr);
    [h,resp,ff] = plot_freqResp(A,B,sr,fUp,figNum,0);

    
    resp_env = resp_env + resp;
    
    colr = [1-C(ib),1-C(ib),C(ib)];
    plot(ff,abs(resp),'color',colr,'linewidth',2)
    
    hl = line([fLow fLow],[0.4 1],'Color',colr,'lineWidth',2); hasbehavior(hl,'legend',false);
    hl = line([fUp  fUp], [0.4 1],'Color',colr,'lineWidth',2); hasbehavior(hl,'legend',false);
    
    legendString{ib} = [num2str(fLow),' - ',num2str(fUp),' Hz'];
    
end

legendString{nb+1} = 'Envelope';


% Design plot 
plot(ff,abs(resp_env),'r','lineWidth',2)
grid on
set(gca,'xScale','log','fontSize',ftSize)
set(gca,'yScale','log')
set(gca,'xLim',[0.07,170])
set(gca,'YLim',[1e-1,2.1])
set(gca,'ytick',[1e-1,0.707,1], 'YAxisLocation', 'right')
hleg = legend(legendString,'Location','SouthWest','fontSize',ftSize);
set(hleg,'fontSize',ftSize)
title('Frequency Response', 'fontSize',ftSize)
xlabel('Frequency', 'fontSize',ftSize)
ylabel('Amplitude', 'fontSize',ftSize)

[xlims] = get(gca,'xlim');
line([xlims(1),xlims(2)],[1 1],'color','k')
%line([xlims(1),xlims(2)],[1.5 1.5],'color','r')
line([xlims(1),xlims(2)],[0.707 0.707],'color','r')




%% Example waveform
fprintf(1,'\nREDO EXAMPLE WAVEFORM PART\n')
% fLow = fc(1,1);
% fUp  = fc(end,2);
% 
% fprintf(1,'\nPassing example waveform through filter ...\n')
% [yExample,~,~] = butter_pass_tdomain_f(vel,fLow,fUp,sr,realisationType,0);
% 
% subplot(3,1,3)
% hold on
% plot(tvel,yExample,'k','lineWidth',2);
% plot(tvel,y_sum,'r','lineWidth',2);
% title('Example waveform','fontSize',ftSize)
% hleg = legend('original','summed up passband-outputs','Location','SouthEast');
% set(hleg,'fontSize',ftSize)
% set(gca,'fontSize',ftSize)
% xlabel('Time [sec]', 'fontSize',ftSize)
% ylabel('Amplitude [m/s]', 'fontSize',ftSize)
% set(gca,'xlim',[5 20])