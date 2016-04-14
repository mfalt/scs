n = 10;
a = randn(n,1);
b = randn(n,1);
c = randn(n,1);
d = randn(n,1);
u = randn(n,1);
l = min(randn(n,1),u);

f = @(t) c+t*d+abs(a+t*d);
ts = -1:0.01:1;
fs = zeros(length(ts),1);
for i = 1:length(ts)
    fs(i) = norm(f(ts(i)));
end
plot(ts,fs)