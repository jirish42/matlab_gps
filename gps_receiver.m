%% Load Data

fs = 2046000;
df = loadFile('gps_1520.dat')'; % replace with SDR sample filename
%% Acquisition Stage
% Determine visible satellites and estimate code chip phase and doppler shift
% Performs a 2d search for maximum correlation between signal and code replica
% over code phase and doppler frequency

d = df(1:20000);
t = [0:length(d)-1]*(1/fs);

generate_plots = [6,7]; % choose which prn numbers to visualize acquisition for
plots = [];
acq_data = [];

fstep = 100;
frange = -10000:fstep:10000; % doppler shift search range

for i=1:37
    % prepare c/a code replica
    prn = cacode(i,2);
    c = repmat(-1*prn + ~prn,1,round(length(d)/length(prn)));
    C = fft(c,length(d));
    
    peaks = [];
    xcorrs = [];
    for f=frange
        % shift data in frequency then correlate with code replica
        d_shifted = real(d).*cos(2*pi*f*t) + imag(d).*sin(2*pi*f*t);
        D = fft(d_shifted,length(d));
        r = ifft(conj(C).*D);
        r = r(1:2046); % code is periodic over 2046 samples, so no need to consider repeated correlations
        [peak, s_offset] = max(abs(r)); 
        peaks = [peaks; peak, s_offset]; 
        if ismember(i,generate_plots)
            xcorrs = [xcorrs; abs(r(1:2046))];
        end

    end
    if ismember(i,generate_plots)
        plots = [plots xcorrs];
    end
    [fpeak, f_index] = max(peaks(:,1));
    acq_data = [acq_data; fpeak, frange(f_index), peaks(f_index,2)]; % [peak, freq offset, sample offset]
    
end

visible_satellites = [];
noise = mean(acq_data(33:37,1)); % PRN's 33-37 are unused, so will give a good estimate of the noise floor
for i=1:32
    if acq_data(i,1)/noise > 2
        visible_satellites = [visible_satellites; [i acq_data(i,:)]];
    end
end

stem(acq_data(:,1)/noise);
title("Acquisition Results");
axis([1,37,0,4.5]);
xlabel("PRN Number");
ylabel("Peak Correlation/Noise");
%% Tracking Stage/Navigation Data Recovery

i = 6; % choose any locked satellite
prn = cacode(i,2);

% initialize frequency and phase shifts with data from acquisition stage
f_offset_init = -1200;
s_offset_init = 625;
pseudosymbols = [];
f_offset = f_offset_init;
s_offset = s_offset_init;

num_symbols=15; % number of navigation bits to record
for s=0:num_symbols-1
    % recalibrate freq and phase shifts every 20 ms
    d = df(s_offset+s*20*2046:s_offset+(s+1)*20*2046 - 1);
    t = [s*length(d):(s+1)*length(d) - 1]*(1/fs);
    c = repmat(-1*prn + ~prn,1,round(length(d)/length(prn)));
    max_corr = abs(dot(c,d));
    f_offset_refined = f_offset;
    s_offset_refined = s_offset;
    
    for delta_f=-50:50
        for delta_s=-1:1
            d = df(s_offset+delta_s+s*20*2046:s_offset+delta_s+(s+1)*20*2046 - 1);
            d_shifted = real(d).*cos(2*pi*(f_offset+delta_f)*t) + imag(d).*sin(2*pi*(f_offset+delta_f)*t);
            corr = abs(dot(c,d_shifted));
            if corr>max_corr
                max_corr = corr;
                f_offset_refined = f_offset+delta_f;
                s_offset_refined = s_offset+delta_s;
            end
        end
    end
    f_offset = f_offset_refined;
    s_offset = s_offset_refined;
    for ps=0:19
        d = df(s_offset+(20*s+ps)*2046:s_offset+(20*s+ps+1)*2046 - 1);
        t = [(20*s+ps)*length(d):(20*s+ps+1)*length(d) - 1]*(1/fs);
        d_shifted = real(d).*cos(2*pi*f_offset*t) + imag(d).*sin(2*pi*f_offset*t);

        c = repmat(-1*prn + ~prn,1,round(length(d)/length(prn)));
        r = dot(c,d_shifted);
        pseudosymbols = [pseudosymbols r<0];
    end
   
end

stem(pseudosymbols);
title("Navigation Message Pseudosymbols (PRN 6)");
xlabel("Code Periods (1ms)");
ylabel("Pseudosymbol");
%% Average Pseudosymbols into Navigation Message

symbol_start = 41; % pseudosymbol offset of first navigation bit boundary
symbols = [];
for n=symbol_start:20:length(pseudosymbols)
    symbols = [symbols mean(pseudosymbols(n:n+19))>0.5];
end
waveform = conv(upsample(symbols,20),ones(1,20));
figure(2);
plot(waveform);
axis([0, length(waveform)-20, -0.5, 1.5]);
title("Navigation Message Waveform (PRN 6)");
xlabel("Time (ms)");
ylabel("Navigation Message");

