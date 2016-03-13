function [bslcrt,bslcrtminus,WN,a] =baseline_FTIR(file,header,curfol);


fid = fopen(file);
b=textscan(fid,'%n%n','headerlines',18,'delimiter','\t');
fclose(fid);
a(:,1)=b{1,1};
a(:,2)=b{1,2};



WN(:,1)=a(:,1);
bsl(:,1)=a(:,2);
[m n]=size(a);

%%first loop
for i=1;   
y(:,i) = polyfit(WN,bsl(:,i),5);

for x=1:m
    bslcrt(x,i)=bsl(x,i)-((WN(x,1)^5*y(1,i))+(WN(x,1)^4*y(2,i))+(WN(x,1)^3*y(3,i))+(WN(x,1)^2*y(4,i))+(WN(x,1)*y(5,i))+(y(6,i)));

     if(bslcrt(x,i)>0)
       bslcrtminus(x,i)=0;
     else
       bslcrtminus(x,i)=bslcrt(x,i);
     end

end
end

%%2nd loop and so on
for i=2:10;   
y(:,i) = polyfit(WN,bslcrtminus(:,i-1),5);

for x=1:m
    bslcrt(x,i)=bslcrt(x,i-1)-((WN(x,1)^5*y(1,i))+(WN(x,1)^4*y(2,i))+(WN(x,1)^3*y(3,i))+(WN(x,1)^2*y(4,i))+(WN(x,1)*y(5,i))+(y(6,i)));

     if(bslcrt(x,i)>0)
       bslcrtminus(x,i)=0;
     else
       bslcrtminus(x,i)=bslcrt(x,i);
     end

end

end

new(:,1)=header(1:18,1);
new(1:11,2)=header(1:11,2);
new(19:2019,1)=num2cell(WN(:,1));
new(19:2019,2)=num2cell(bslcrt(:,10));

new(12,2)=num2cell(a(1,1));%first x
new(13,2)=num2cell(a(2001,1));%last x
%new(12,2)=num2cell(m);%size
new(14,2)=num2cell(2001);
new(15,2)=num2cell(a(1,2));%first y
%new(14,2)=num2cell(a(1024,2));%last y
new(16,2)=num2cell(max(a(1:2001,2)));%max y 
new(17,2)=num2cell(min(a(1:2001,2)));%min y



dlmcell('C:\Users\Kin\Desktop\MATLAB\test.txt',new,'\t','w');%% change back to text
 

newname=[curfol,'\baseline\',file,'.txt'];%%new destination
oldname= ('C:\Users\Kin\Desktop\MATLAB\test.txt');%%old destination
copyfile(oldname,newname);

end

 