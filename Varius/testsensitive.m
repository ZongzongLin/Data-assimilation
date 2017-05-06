n=4;
m=10;
anum=1;
bnum=1;
Iteration=100;
x0=rand(n,1);
L=gegraph(n);
Lv=mtov(L,n);
alpha=rand(anum,1);
beta=rand(bnum,1);
arange=linspace(0,1,m+1);
brange=linspace(0,1,m+1);
ya=zeros(n,Iteration,m);
yb=zeros(n,Iteration,m);
for k=1:m+1
    xa=x0;
    xb=x0;
    ya(:,1,k)=x0;
    yb(:,1,k)=x0;
    for i=1:Iteration-1
        tempa=virusdynamic(xa,Lv,arange(k),beta,n);
        ya(:,i+1,k)=tempa(1:n);
        xa=ya(:,i+1,k);
        tempb=virusdynamic(xb,Lv,alpha,brange(k),n);
        yb(:,i+1,k)=tempb(1:n);
        xb=yb(:,i+1,k);
    end
end
figure(1)
for k=1:m
    plot(ya(1,:,k))
    hold on
end
title('change of alpha')
hold off
figure(2)
for k=1:m
    plot(yb(1,:,k))
    hold on
end
title('change of beta')
hold off