clc,close all,clear all
%% read excel file
 [num1,txt1,input] = xlsread('KRAKOW - FULL DATA.xlsx','AllConnectionsInCurrSolution');
%[num1,txt1,input] = xlsread('KRAKOW - FULL DATA.xlsx','TRAFFIC-DEMAND');

%% get sort matrix
%%
x1=num1(:,1);x2=num1(:,2);x3=num1(:,3);
[B,I] = sort(x1);
x1=B;x2=x2(I);x3=x3(I);
num1=[x1,x2,x3];
%%
y1=diff(x1);
[ind,y2]=find(y1>0);

y3=num1(ind);
y3(length(y3)+1)=num1(ind(end)+1);
j=length(y3);
%% make output 

z=zeros(j+1);
z(1,2:end)=y3;
z(2:end,1)=y3';
output=cell(j+1);
output(1,2:end)=num2cell(y3);
output(2:end,1)=num2cell(y3');
  
 ind(end+1) = length(x1);
 inds = 1;
 for i = 1:length(y3)
    inde = ind(i);

    for j = inds:inde
        if(x2(j)==452)
            aa =10;
        end
        [indy,yy]=find(z(:,1)==x2(j));
        z(i+1,indy)=x3(j);
        if(x3(j)==0)
            z(i+1,indy)=inf;
        end
        
    end
    inds = inde;
 end
  xlswrite('output2',z)