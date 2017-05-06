function L=gegraph(n)
flag=0;
while(flag==0)
    L=ones(n)-eye(n);
    for i=1:n
        for j=1:n
            if j~=i
                if rand<0.2
                    L(i,j)=0;
                end
            end
        end
    end
    P=L;
    for k=1:n
        P=P+L^(n);
    end
    flag=1;
    for i=1:n
        for j=1:n
            if P(i,j)==0
                flag=0;
            end
        end
    end
end
