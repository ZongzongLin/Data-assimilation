function y=virusdynamic(x,matv,alpha,beta,n)
y=zeros(n,1);
L=zeros(n);
p=1;
for i=1:n
    for j=1:n
        if j~=1
            L(i,j)=matv(p);
            p=p+1;
        end
    end
end
for i=1:n
    multiply=1;
    for k=1:n
        if L(i,k)~=0
            multiply=multiply*(1-L(i,k)*x(k));
        end
    end
    if length(alpha)==1
        y(i)=(1-(1-alpha)*multiply)*(1-x(i))+(1-beta)*x(i);
    else
        y(i)=(1-(1-alpha(i))*multiply)*(1-x(i))+(1-beta(i))*x(i);
    end
end
y=[y;matv;alpha;beta];
