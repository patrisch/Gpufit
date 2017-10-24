function Fitting_CPU()
%% unfold data
 data = importdata('data.mat');

        xdata  = data(:,1);
        ydata = data(:,2);

%% calculate fft 

%     if showfft

    samplerate = length(xdata)/max(xdata);                     % sample rate (Hz)
    nsamples = length(xdata);
    yfft = abs(fft(ydata));
    yfft = yfft(1:floor(nsamples/2));                                 % discard half of points due to nyquist
    frange = samplerate/nsamples*(0:nsamples/2-1);        % generate frequency range in MHz

    yfft = yfft(frange>0.00);                                  % discard frequencies below 1 MHz
    frange = frange(frange>0.00);

%% estimate initial values    
    
    % get doublet hyperfine frequencies from fft
    minpeakheigth = (max(yfft)-min(yfft))/5;
    [peaks locs] = findpeaks(yfft,'SORTSTR','descend','MINPEAKHEIGHT',minpeakheigth);   
    
    freq = frange(locs(1:2))
    freqamp = peaks(1:2);   
    
    % Estimates for amplitudes, offset and phases; t2* currently hardcoded
    amp = (max(ydata)-min(ydata))/5;
    offset = (max(ydata)+min(ydata))/2; 
    
    phase = 2.8e-08 ./ (2.*pi) .*10.^6;  % Chinese Paper relation, crudest approximation: phi = detuning * pulse length
    
    t2star1 =0.75e-06;
    
%% fitting

    % fit 1 with p = 1
    % model = 'exp(-(x./t2star)^1)*(A1*cos(2*pi*f1*(x - x1)) + A2*cos(2*pi*f2*(x-x2))) + c';
    
    % parameters: [A1 A2 c f1 f2 t2star x1 x2]
    startpoints = [amp amp offset freq(1) freq(2) t2star1 phase phase];

    options1 = optimoptions('lsqcurvefit', ...
        'Jacobian', 'on',...
        'Algorithm', 'levenberg-marquardt',...
        'MaxIter', 10000,...
        'MaxFunEvals', 10000,...
        'TolFun', 1e-19,...
        'TolX', 1e-19,...
        'UseParallel', true);
    
    [x1, resnorm_f1] = lsqcurvefit(@(x,xdata) model1(x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8), xdata), startpoints, xdata, ydata,[],[], options1);

    fitf1 = struct('A1', x1(1), 'A2', x1(2), 'c', x1(3), 'f1', x1(4), 'f2', x1(5), 't2star', x1(6), 'x1', x1(7), 'x2', x1(8))


    % fit 2, with p being optimised as well. Use values of fit 1 as startparameters.
    % model = 'exp(-(x./t2star)^p)*(A1*cos(2*pi*f1*(x - x1)) + A2*cos(2*pi*f2*(x-x2))) + c';
    
    % parameters: [A1 A2 c f1 f2 p t2star x1 x2]
    p=1.0;
    startpoints2 = [fitf1.A1 fitf1.A2 fitf1.c fitf1.f1 fitf1.f2 p fitf1.t2star fitf1.x1 fitf1.x2];
    
    options2 = optimoptions('lsqcurvefit',...
        'Jacobian', 'on',...
        'Algorithm', 'levenberg-marquardt',...
        'MaxIter', 10000,...
        'MaxFunEvals', 10000,...
        'TolFun', 1e-19,...
        'TolX', 1e-19,...
        'UseParallel', true);
    
    [x2, resnorm_f] = lsqcurvefit(@(x,xdata) model2(x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8), x(9), xdata), startpoints2, xdata, ydata,[],[], options2);
    
    fitf = struct('A1', x2(1), 'A2', x2(2), 'c', x2(3), 'f1', x2(4), 'f2', x2(5), 'p', x2(6), 't2star', x2(7), 'x1', x2(8), 'x2', x2(9))

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
    % plot fit 2, p varied
    plot(xdata *10^6, pfit2(xdata),'r-');
    % plot fit 1, p=1
    plot(xdata *10^6, pfit1(xdata),'k--');

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
    


function [F, J] = model1(A1, A2, c, f1, f2, t2star, x1, x2, x)
    F = (c+exp(-x./t2star).*((A1.*cos(f1.*pi.*(x-x1).*2.0)+A2.*cos(f2.*pi.*(x-x2).*2.0))));
    J = [(cos(f1.*pi.*(x-x1).*2.0).*exp(-x./t2star)),...
        (cos(f2.*pi.*(x-x2).*2.0).*exp(-x./t2star)),...
        (ones(1,size(x,1))'),...
        (A1.*pi.*sin(f1.*pi.*(x-x1).*2.0).*exp(-x./t2star).*(x-x1).*-2.0),...
        (A2.*pi.*sin(f2.*pi.*(x-x2).*2.0).*exp(-x./t2star).*(x-x2).*-2.0),...
        (1.0./t2star.^2.*x.*exp(-x./t2star).*(A1.*cos(f1.*pi.*(x-x1).*2.0)+A2.*cos(f2.*pi.*(x-x2).*2.0))),...
        (A1.*f1.*pi.*sin(f1.*pi.*(x-x1).*2.0).*exp(-x./t2star).*2.0),...
        (A2.*f2.*pi.*sin(f2.*pi.*(x-x2).*2.0).*exp(-x./t2star).*2.0)];
end

function [F2, J2] = model2(A1, A2, c, f1, f2, p, t2star, x1, x2, x)
    F2 = (exp(-(x./t2star).^p).*((A1.*cos(2.*pi.*f1.*(x - x1)) + A2.*cos(2.*pi.*f2.*(x-x2)))) + c);
    J2 = [(cos(f1.*pi.*(x-x1).*2.0).*exp(-(x./t2star).^p)),...
        (cos(f2.*pi.*(x-x2).*2.0).*exp(-(x./t2star).^p)),...
        (ones(1,size(x,1))'),...
        (A1.*pi.*sin(f1.*pi.*(x-x1).*2.0).*exp(-(x./t2star).^p).*(x-x1).*-2.0),...
        (A2.*pi.*sin(f2.*pi.*(x-x2).*2.0).*exp(-(x./t2star).^p).*(x-x2).*-2.0),...
        (-log(x./t2star).*exp(-(x./t2star).^p).*(x./t2star).^p.*(A1.*cos(f1.*pi.*(x-x1).*2.0)+A2.*cos(f2.*pi.*(x-x2).*2.0))),...
        (p.*1.0./t2star.^2.*x.*exp(-(x./t2star).^p).*(x./t2star).^(p-1.0).*(A1.*cos(f1.*pi.*(x-x1).*2.0)+A2.*cos(f2.*pi.*(x-x2).*2.0))),...
        (A1.*f1.*pi.*sin(f1.*pi.*(x-x1).*2.0).*exp(-(x./t2star).^p).*2.0),...
        (A2.*f2.*pi.*sin(f2.*pi.*(x-x2).*2.0).*exp(-(x./t2star).^p).*2.0)];
end
 
    
 
 
end
