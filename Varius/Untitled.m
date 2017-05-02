n=10;
a=rand(n);
bara=mean(a,2);
sum=zeros(n);
for i=1:n
    sum=sum+(a(:,i)-bara)*(a(:,i)-bara)';
end
sums=1/(n)*sum
b=cov(a)
cond(b)
cond(sums)