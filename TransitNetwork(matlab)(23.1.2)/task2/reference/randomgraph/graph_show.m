%%{
clear
close all 
clc


%% Load matricies 
load('KRAKOW_data_final.mat')
KRAKOWTravelTimes = load('cracow.mat').cracow;
geration=0;Parent=[];
min_TF=inf;
MandlTravelTimes=KRAKOWTravelTimes;
MandlDemand=KRAKOWDemand;
total_passnger=sum(sum(MandlDemand));
weight=lt(KRAKOWTravelTimes,1000);    %make weight matrix to detrmine the links between nodes

%% Graph Construction 
G=graph(weight,'omitselfloops');
hh= plot((G),'Markersize',4);
%% Travel Time and Demand extraction for all existente links
% for a = 1 : height(G.Edges)
%     garbage.s = G.Edges.EndNodes(a,1);
%     garbage.t = G.Edges.EndNodes(a,2);
%     garbage.Z(a,:) = KRAKOWTravelTimes(garbage.s,garbage.t); % get
%     garbage.D(a,:) = KRAKOWTravelTimes(garbage.s,garbage.t);
%  end
% G.Edges.Time = garbage.Z(:,1);
% G.Edges.Demand = garbage.D(:,1);
%,'XData',(No_name(:,1)),'YData',(No_name(:,2)));  %plot the graph

