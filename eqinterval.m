function eqinterval(str);

curfol=cd;
for j={str}%change the number of folder
%foldernum=num2str(j);
dishfol= [curfol,'\',j{:}];%,foldernum
cd(dishfol);    

mkdir ('eqinterval')
imported=dir('*.txt');   
[m n]=size(imported);
for i=1:m
raw=imported(i,1).name;
fid = fopen(imported(i,1).name);
tscan= textscan(fid,'%n%n','headerlines',3,'delimiter','\t');

c(1,:)=tscan{1,1}';
c(2,:)=tscan{1,2}';
fclose(fid);

%%wavenumber
newWN(1,1)=c(1,1);
flnum=floor(c(1,1024));
decinum=c(1,1)-floor(c(1,1));
flnum=flnum+decinum;
for x=1:1073
if newWN<flnum
newWN(1,x+1)=newWN(1,x)+1;
end
end


%interpol by Lagrange

yi(1,1) = c(2,1);
yi(1,2) = LagrangeInter(c(1,1:3),c(2,1:3),newWN(1,2));
yi(1,3) = LagrangeInter(c(1,1:3),c(2,1:3),newWN(1,3));
for x=4:1074
for y=1:1024
flsel(1,y)=abs(newWN(1,x)-c(1,y));
end
[idx idx] = min(flsel(1,:));
flselidx(1,x)=idx;
idx2=idx-2;
yi(1,x) = LagrangeInter(c(1,idx2:idx),c(2,idx2:idx),newWN(1,x));
end
newWN(2,:)=yi;


prtfile(:,1)=num2cell(newWN(1,:));
prtfile(:,2)=num2cell(newWN(2,:));
dlmcell('C:\Users\Kin\Desktop\test interval\test.txt',prtfile,'\t','w');%% change back to text
newname=[dishfol,'\eqinterval\',raw,'.txt'];%%new destination
oldname= ('C:\Users\Kin\Desktop\test interval\test.txt');%%old destination
copyfile(oldname,newname);


end
end
end

