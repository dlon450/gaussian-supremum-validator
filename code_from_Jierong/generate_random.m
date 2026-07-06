function [c,b,ksi] = generate_random(n,d)
c = randn(d,1);
b = rand(1,1,'double');
ksi = randn(n,d);
