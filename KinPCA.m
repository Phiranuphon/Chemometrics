function [PCA]=KinPCA(data)
[m n b]=size(data);

for x=1:b
[test.V test.D] = eig(cov(data(:,:,x)));
PCA.V(:,:,x)=test.V;
PCA.D(:,:,x)=test.D;

%plot eigenvector 1-5
figure
hold on
i=1;
for j=n:-1:n-4
subplot(5,1,i)
plot(PCA.V(:,j,x))
i=i+1;
end
xlabel(x)
end

end