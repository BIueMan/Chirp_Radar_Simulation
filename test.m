% create a random matrix of size m x n
m = 4;
n = 5;
X = rand(m, n);

% perform fft on every row using matlab with a specified output length
n = 10; % output length
Y = fft(X,n,2);

Y2 = zeros(m, n);

% perform fft on every row using a loop with a specified output length
for i = 1:m
    Y2(i,:) = fft(X(i,:),n);
end