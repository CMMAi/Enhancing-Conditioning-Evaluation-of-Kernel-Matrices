function x=gauss_jordan(A,force)

n=size(A,1);
norm2=sqrt(sum((A).^2));% New change
A=A./repmat(norm2,n,1);% New change
for i=1:n-1
    if A(i,i)==0
        [max_in_column,k]=max(abs(A(i+1:end,i)));
        A(i,:)=A(i,:)+A(i+k,:);
        force(i)=force(i)+force(i+k);
    end
    
    for j=i+1:n
        CC=A(j,i)/A(i,i);
        A(j,:)=A(j,:)-CC*A(i,:);
        force(j)=force(j)-CC*force(i);
    end
end

x=zeros(n,1);

for i=n:-1:1
    S=0;
    for j=i+1:n
        S=S+A(i,j)*x(j);
    end
    S=force(i)-S;
    x(i)=S/A(i,i);
end
x = x./norm2'; % New change

end
