
function fv = local_regression(times, data_col)

bandwidth=10; % empirical bandwidth

n=length(data_col);     % length of each data series
bs=zeros(n,2);
te=zeros(2,1);
ov=ones(n,1);
ones_times=[ov times]; % column of ones added to times column
wcoeff=(1/(2*pi*bandwidth^2))^0.5; % coeff in front of kernel
hv=2*bandwidth^2;
for ia=1:n
    w(:,ia)=wcoeff*exp(-(times-times(ia)).^2/hv); % w(j,:) is the jth set of weights
end
%w a square vector
for ia=1:n
    W=diag(w(ia,:));
    te=ones_times'*W;        % te*ones_times gives column [S,S_x S_x,S_xx] which is 2 x n : see below
    % te2 = ones_times.*w(i,:);
    bs(ia,:)=(inv(te*ones_times)*(te*data_col))';
end

fv=bs(:,1)+bs(:,2).*times;

end
