function Fitting_Gpufit()
%% unfold data
 data = importdata('data.mat');

        xdata  = data(:,1);
        ydata = data(:,2);
 
    
%% Multiplies data by #fits, needed for performance testing
    number_fits = 1;
    ydata_mat = repmat(single(ydata), [1, number_fits]);

%% calculate fft (need to feed frequencies as initial parameters)
    samplerate = length(xdata) / max(xdata);                      % sample rate (Hz)
    nsamples = length(xdata);
    yfft = abs(fft(ydata));
    yfft = yfft(1:floor(nsamples / 2));                           % discard half of points due to nyquist
    frange = samplerate / nsamples * (0:nsamples / 2 - 1);        % generate frequency range in MHz

    yfft = yfft(frange > 0.00);                                   % discard frequencies below 1 MHz
    frange = frange(frange > 0.00);

%% estimate initial values    
    
    % get doublet hyperfine frequencies from fft
    minpeakheigth = (max(yfft) - min(yfft)) / 5;
    [peaks locs] = findpeaks(yfft,'SORTSTR','descend','MINPEAKHEIGHT',minpeakheigth);   
    
    freq = frange(locs(1:2))
    freqamp = peaks(1:2);   
     
    % Estimates for amplitudes, offset and phases; t2* currently hardcoded
    amp = (max(ydata) - min(ydata)) / 5;
    offset = (max(ydata) + min(ydata)) / 2; 

    phase = 2.8e-08 ./ (2.*pi) .*10.^6;  % Chinese Paper relation, crudest approximation: phi = detuning * pulse length
    
    t2star1 =0.75e-06;
%% fitting

    % fit 1 with p = 1
    % model = 'exp(-(x./t2star)^1)*(A1*cos(2*pi*f1*(x - x1)) + A2*cos(2*pi*f2*(x-x2))) + c';
    
    % parameters: [A1 A2 c f1 f2 t2star x1 x2]
    startpoints = [amp amp offset freq(1) freq(2) t2star1 phase phase];
    startpoints_mat = repmat(single(startpoints'), [1, number_fits]) %  .*(1 + 0.1.*sin(rand(8, number_fits)));
    
    [x1, states, chi, ite, time] = gpufit(single(ydata_mat),[], ModelID.RAMSEY_FIXED_P, single(startpoints_mat), 1e-19, 10000,[1 1 1 1 1 1 1 1]', EstimatorID.LSE, single(xdata) );
    fitf1 = struct('A1', x1(1, :), 'A2', x1(2, :), 'c', x1(3, :), 'f1', x1(4, :), 'f2', x1(5, :), 't2star', x1(6, :), 'x1', x1(7, :), 'x2', x1(8, :));

    % fit 2, with p being optimised as well. Use values of fit 1 as startparameters.
    % model = 'exp(-(x./t2star)^p)*(A1*cos(2*pi*f1*(x - x1)) + A2*cos(2*pi*f2*(x-x2))) + c';
    
    % parameters: [A1 A2 c f1 f2 p t2star x1 x2]
    p=ones(1,number_fits);
    startpoints2 = [fitf1.A1; fitf1.A2; fitf1.c; fitf1.f1; fitf1.f2; p; fitf1.t2star; fitf1.x1; fitf1.x2];
  
    [x2, states, chi, ite, time] = gpufit(single(ydata_mat),[], ModelID.RAMSEY_VAR_P, single(startpoints2), 1e-19, 10000,[1 1 1 1 1 1 1 1 1]', EstimatorID.LSE, single(xdata) );
    fitf = struct('A1', x2(1,:), 'A2', x2(2,:), 'c', x2(3,:), 'f1', x2(4,:), 'f2', x2(5,:), 'p', x2(6,:), 't2star', x2(7,:), 'x1', x2(8,:), 'x2', x2(9,:));
    
number_converged = sum(states == 0);
fprintf('\nratio converged         %6.2f %%\n', number_converged / number_fits * 100);

%% Results and plotting
%%% plotting

epilog = char(['T2* = ' num2str(round(fitf.t2star*10^6,2)) ' us p=' num2str(round(fitf.p,2))]);

% set figure
if ~isempty(findobj('name','exp fit'))
    
    fh =  findobj('name','exp fit');
    set(0, 'currentfigure', fh);  % for figures
    clf;
    
    else
    
    figure('name','exp fit','Position',[350 200 800 800/(2*1.618)]);    %define figure handle
    
end

    % plot data + fit
    subplot(1,2,1);          
    plot(xdata *10^6, ydata,'-', 'MarkerSize', 10);
    text(1, min(ydata) + 0.2*abs(min(ydata)), epilog);
    
    hold on;
    pfit1 = @(x) exp(-(x./fitf1.t2star).^1).*((fitf1.A1.*cos(2.*pi.*fitf1.f1.*(x - fitf1.x1)) + fitf1.A2.*cos(2.*pi.*fitf1.f2.*(x-fitf1.x2))))+fitf1.c;
    pfit2 = @(x) exp(-(x./fitf.t2star).^fitf.p).*((fitf.A1.*cos(2.*pi.*fitf.f1.*(x - fitf.x1)) + fitf.A2.*cos(2.*pi.*fitf.f2.*(x-fitf.x2))))+fitf.c;
    % plot var p
    plot(xdata .*10^6, pfit2(xdata), 'r-');
    %plot p=1
    plot(xdata .*10^6, pfit1(xdata), 'k--');
      
    xlabel('\tau (us)');
    ylabel('Ramsey Contrast');

    axis([0 15 -inf inf]);    
    axis tight;
    
    legend off;
    
    % plot fft
    subplot(1,2,2);
    plot(frange/10^6, yfft); % in MHz
    xlabel('Frequency (Mhz)');
    ylabel('FFT (arb. units)');
%     
    hold on;
    plot(freq*10^-6, freqamp, 'k^')
 
    
 
 
end
