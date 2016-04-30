MATLAB Code:

Main Code :

videoFReader = vision.VideoFileReader('C:\Users\User\Desktop\Spring\DSP\Project 1\iphone.MOV');
videoPlayer = vision.VideoPlayer;
 while ~isDone(videoFReader)
   videoFrame = step(videoFReader);
   step(videoPlayer, videoFrame);
   pause(0.04);
 end
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Brightness Calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 brightness=zeros(1,390);
 videoPlayer=VideoReader('C:\Users\User\Desktop\Spring\DSP\Project 1\iphone.MOV');
 for i = 1:390	
frame = read(videoPlayer, i);  % frame-by-frame reading
redPlane = frame(:,:,1);       % redplane gives intensity information
brightness(i) = sum(sum(redPlane)) / (size(frame, 1) * size(frame, 2)); 
 end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Butterworth Filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BPM_L=40;
BPM_H=230;
FrameRate=20;
[b, a] = butter(2, [((BPM_L / 60) / FrameRate * 2), ((BPM_H/ 60) / FrameRate * 2)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Decision of end parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  fprintf('Values of b are \n');
%  disp(b);
%  fprintf('Values of a are \n');
%  disp(a);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Brightness Filtering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 filtBrightness = filter(b, a, brightness);
 figure(1);
 stem(filtBrightness);
 xlabel('Frame Numbers');
 ylabel('Brightness');
 title('Framerate=60');

 fftbright=abs(fft2(filtBrightness)); 
 figure(2);
 plot(fftbright);
 xlabel('Frequency');
 ylabel('Amplitude');
 title('Framerate=60');
 %%%%%%%%%%%%%%%%%%%%%%%%% Sliding Window%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 slidingWindow=6;
 hannWindow=hann(size(slidingWindow,2));
 fftMagnitude=(fft(slidingWindow.*hannWindow));
 display(fftMagnitude);
 %%%%%%%%%%%%%%%%%%%%%%%%% Peak Detection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 rangeOfInterest=((BPM_L:BPM_H)/60)*(size(fftMagnitude,2)/7.667)+1; 
 fftRange=fft(rangeOfInterest);
  [peaksValues, peakIndices] = findpeaks(fftMagnitude(1));
% [peaksValues, peakIndices] = findpeaks(fftMagnitude);
 display(peaksValues);
 % Find the highest peak
[maxPeakValue, maxPeakIndex] = max(peaksValues);
% Translate the peak index to an FFT vector index
bpmFreqIndex = rangeOfInterest(peakIndices(maxPeakIndex));
display(bpmFreqIndex);
% Get the frequency in bpm that corresponds to the highest  peak
bpmPeak = (bpmFreqIndex - 1) * (7.667 /size(fftMagnitude, 2)) * 60;
display(bpmPeak);
figure(3);
 stem(rangeOfInterest);
 xlabel('Frame Numbers');
 ylabel('Amplitude');
 title('RangeOfInterest');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Smoothing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

WINDOW_SIZE=0.01;
fftResolution = 1 / WINDOW_SIZE;
lowFreq = bpmPeak / 60 - 0.5 * fftResolution;
smoothingResolution = 120 / 60;
testFreqs = round(fftResolution / smoothingResolution);
display(testFreqs);
power = zeros(1, testFreqs);
display(size(smoothingResolution+lowFreq));
display(size(0:testFreqs - 1)');
% freqs = (0:testFreqs - 1) .* smoothingResolution+lowFreq;
for freqs=1:1:50;
for k = 1:testFreqs,
re = 0; im = 0;
for j = 0:(size(b, 2) - 1),
phi = 2 * pi * freqs * (j / 7.667);
re = re + b(j+1) * cos(phi);

im = im + b(j+1) * sin(phi);

end
end
end

[maxPeakValue, maxPeakIndex] = max(power);
smoothedBpm = 60 * freqs(maxPeakIndex);
display(smoothedBpm);

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Quantization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  uniquan(filtBrightness,100);
%   alaw(filtBrightness);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Uniform Quantization Function: 

function out_seq=uniquan(fftbright,bits)
 in_max=max(abs(fftbright));           % fftbright is final fft signal of Project 1
in_min=0;
step=(in_max-in_min)/(2^bits);
fftbright=floor((fftbright-in_min)/step);
delta=in_min+(step/2);
L=length(fftbright);
out_seq = zeros(size(fftbright), 'double');
for y=1:L
out_seq(y) = ((double(fftbright(y)) * step) + delta)*100;
figure(4);
hold on;
plot(y,out_seq(y));
xlabel('Frequency');
 ylabel('Amplitude');
 title('Quantized Output, Number of bits=100');
display(out_seq(y));
display(y);
end
 end
