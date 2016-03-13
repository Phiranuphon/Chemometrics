function [IWRMSEallspec,OLS_percent]= IWOLS_check_spec(sample,lipidlibrary,compnum,range,step);
WNstr=80;
WNstp=WNstr+step;
WNstp1=WNstr+step;

count=1;

[n m]=size(sample)
j=compnum;
Cmix=eye(j);

while(WNstp<range)

 if compnum==4
        [C,A,R,ssq]=ALS(sample,lipidlibrary,j,[1:4],WNstr:WNstp);
    elseif compnum==3
        [C,A,R,ssq]=ALS(sample,lipidlibrary,j,[2:4],WNstr:WNstp);%2:4 only for POL, C10 must use 1:4 and j must be 4  
    end

OLS(1:n,:)=C(j+1:n+j,:);%%create OLS conc. matrix

    for y=1:j
        for x=1:n
          OLS_percent(x,y)=OLS(x,y)/sum(OLS(x,1:j));%re-calculate conc. into %
        end
    end
 
RMSEspec =zeros(n,m);
      
    if j==3
        for x=1:n
            syn_spec(x,:) = OLS_percent(x,1)*lipidlibrary(2,:)+OLS_percent(x,2)*lipidlibrary(3,:)+OLS_percent(x,3)*lipidlibrary(4,:);%for 3components
            syn_spec_norm(x,:)=syn_spec(x,:)/max(syn_spec(x,WNstr:WNstp));
            sample_norm(x,:)=sample(x,:)/max(sample(x,WNstr:WNstp));
        end
    elseif j==4
        for x=1:n
            syn_spec(x,:) = OLS_percent(x,1)*lipidlibrary(1,:)+OLS_percent(x,2)*lipidlibrary(2,:)+OLS_percent(x,3)*lipidlibrary(3,:)+OLS_percent(x,4)*lipidlibrary(4,:);%for 4components
            syn_spec_norm(x,:)=syn_spec(x,:)/max(syn_spec(x,WNstr:WNstp));
            sample_norm(x,:)=sample(x,:)/max(sample(x,WNstr:WNstp));
        end
    end
    
    for x=1:n
        for y=1:m
            subt_spec(x,y)=sample_norm(x,y)-syn_spec_norm(x,y);
            subt_spec_sqr(x,y)=sqrt(subt_spec(x,y)^2);
        end
    end
    for x=1:n
        RMSEspec(x,1)=sum(subt_spec_sqr(x,:))/m;
    end
    
RMSEall(1:n,count)=RMSEspec(:,1);
RMSEall(n+1,count)=WNstp;

count=count+1;
WNstp=WNstp1+(step*count);


end

    for x=1:n
        testmin(1,1)=min(RMSEall(x,:));
        [testa1,testb1]=find(RMSEall(x,:)==testmin(1,1));
        RMSEall(x+n+1,1)=RMSEall(x,testb1);
        RMSEall(x+n+1,2)=RMSEall(n+1,testb1);

    end


IWRMSEallspec=RMSEall;

end