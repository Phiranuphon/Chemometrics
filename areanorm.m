function [normdata]=areanorm(data,start,stop);

%data=data';
[n m]=size(data);
data1=data(:,start:stop);
[n1 m1]=size(data1);
normdata=zeros(n,m);
    for x=1:n
        testmin(1,1)=min(data1(x,[1:15]));
        testmin(1,2)=min(data1(x,[m1-15:m1])) ;
        testmin(1,3)=mean(testmin(1,1:2));

        testmin1=find(data1(x,:)==testmin(1,1));
        testmin2=find(data1(x,:)==testmin(1,2));

        test_area(1,:)=sum(data1(x,testmin1(1,1):testmin2(1,1)))-(length(data1(x,testmin1(1,1):testmin2(1,1)))*testmin(1,3));
        normdata(x,:)=data(x,:)/test_area(1,1);
    end


end