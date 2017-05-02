function matv=mtov(mat,n)
matv=zeros(n*(n-1),1);
p=1;
for i=1:n
    for j=1:n
        if j~=i
            matv(p)=mat(i,j);
            p=p+1;
        end
    end
end
            
        