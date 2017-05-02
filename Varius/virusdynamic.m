function y=virusdynamic(x,matv,alpha,beta,n,pnum,SG,sigma)
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
    y(i)=(1-(1-alpha(i))*multiply)*(1-x(i))+(1-beta(i))*x(i);
end
y=[y;matv;alpha;beta];
if strcmp(SG,'Gaussian')
        y=y+sigma*randn(pnum,1);
elseif strcmp(SG,'Normal')
        y=y+sigma*(2*rand(pnum,1)-1);
end 
y=cutoff(y);
