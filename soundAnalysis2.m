function [EOD,Spike,EODR]=soundAnalysis2(sound)


% [value,sample]=findpeaks(diff([sound(1); sound]),'MINPEAKHEIGHT',0.5,'MINPEAKDISTANCE',150);
[value,sample]=findpeaks( sound,'MINPEAKHEIGHT',.4,'MINPEAKDISTANCE',150);
EOD=zeros([length(sound),1]); % new timevector
for j=1:length(sample);
    EOD(sample(j))=1;
end

for i=1:length(sample);
    if i==1;
        Spike(i)=sample(i);
    elseif sample(i)>=sample(i-1)+10;
        Spike=[Spike;sample(i)];
    else
    end
end
EODR=zeros(size(EOD));
for i=2:length(Spike);
    EODR(Spike(i-1):Spike(i))=1/((Spike(i)-Spike(i-1))/10000);
end
% times=10000*1:10000*120;
% plot(times(sample),10000./diff([Spike(1,:) ;Spike]),'r')
% hold on
% plot(times,EODR)
% plot(times,sound,'r')
% hold on
% plot(times(sample),1,'ok')
end