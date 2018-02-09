% function this function detrends the time-series and its amplitude.
% bandwith is defunct.

function z = mirela_amp_detrend(data_col, bandwidth, gind)

%bandwidth is always set to 10 in local_regression
%[times data_col] = textread('data.txt', '%f\t%f');
times = 1:length(data_col);
times = times';


fv1 = local_regression(times, data_col, bandwidth, 1);
x=data_col-fv1;
y=abs(x);
% i=x./y;
i=sign(x);
R=local_regression(times, y, bandwidth, 1);
mR=mean(R);
z=y./R;
z=mR*(z.*i)+mean(fv1);
%if gind > 0
   %plot(times, z);
%end
%if gind < 0
 %subplot(211), plot(times,data_col,'--mo'); hold on; plot(times,fv1); hold off;
 %subplot(212), plot(times,z,'--ro');
%end

end

