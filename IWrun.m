function [IWRMSEall,OLS_percent,RMSEGC,reglist,allC]=IWrun(GCdata,sample,lipidlibrary,compnum,range,step,type);

if type==1
[IWRMSEall,OLS_percent,RMSEGC]= IWOLS_check(GCdata,sample,lipidlibrary,compnum,range,step);
elseif type==2
[IWRMSEall,OLS_percent,RMSEGC]=IWOLS_check_meanC(GCdata,sample,lipidlibrary,compnum,range,step);
end

[n m]=size(sample);
allC=zeros(n,compnum);
    if type==1
        reglist(:,1)=IWRMSEall(n+3:2*n+2,2);
    
         if compnum==4
             for x=1:n
                 [C,A,R,ssq]=ALS(sample(x,:),lipidlibrary,compnum,1:4,[81:reglist(x,1)]);
                    allC(x,:)=C(7,1:4);
             end   
         elseif compnum==3
                for x=1:n
                    [C,A,R,ssq]=ALS(sample(x,:),lipidlibrary,compnum,2:4,[81:reglist(x,1)]);    
                allC(x,:)=C(8,1:3);
                end
         end
 elseif type==2
         reglist(1,1)=IWRMSEall(3,2);
          if compnum==4
             for x=1:n
                 [C,A,R,ssq]=ALS(sample(x,:),lipidlibrary,compnum,1:4,[81:reglist(1,1)]);
                    allC(x,:)=C(7,1:4);
             end   
         elseif compnum==3
                for x=1:n
                    [C,A,R,ssq]=ALS(sample(x,:),lipidlibrary,compnum,2:4,[81:reglist(1,1)]);    
                allC(x,:)=C(8,1:3);
                end
          end     
         
    end
       
    if compnum==4        
        figure
        ylim([0 0.6])
        p=plot(GCdata(1,:)','b');
        hold on
        [Cn Cm]=size(allC);
        allC(Cn+1,:)=mean(allC);
        plot(allC(Cn+1,:)','r')
        set(get(p,'Parent'), 'XTick', [1:4] );
        set(get(p,'Parent'), 'XTickLabel', {'DA','OA', 'PA', 'LA'} );
        set(get(p,'Parent'),'FontSize',30,'fontWeight','bold')
  figure
plot(IWRMSEall(2,:),IWRMSEall(1,:),'ob')
hold on
plot(IWRMSEall(3,2),IWRMSEall(3,1),'or')
        
elseif compnum==3
        figure
        ylim([0 0.6])
        p=plot(GCdata(1,2:4)','b');
        hold on
        [Cn Cm]=size(allC);
        allC(Cn+1,:)=mean(allC);
        plot(allC(Cn+1,:)','r')% plot(allC(Cn+1,:)','r')
        set(get(p,'Parent'), 'XTick', [1:3] );
        set(get(p,'Parent'), 'XTickLabel', {'OA', 'PA', 'LA'} );
        set(get(p,'Parent'),'FontSize',30,'fontWeight','bold')
figure
plot(IWRMSEall(2,:),IWRMSEall(1,:),'ob')
hold on
plot(IWRMSEall(3,2),IWRMSEall(3,1),'or')
    end

end


%check 0 entry
%for x=1:15
%[n m]=find(allC(x,:)==0)
%if isempty(n)
%test(x,:)=allC(x,:)
%end
%end





  