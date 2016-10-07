function [theGrid,gridNoNoise] = areaGenerationNewVarLin( timeseriesFile )

%%grid size and time series size.
dimX = 50;
dimY = 70;

sMin = 2.0;
%%Step 1. generate initial time series.
timeseriesStruct = load(timeseriesFile);
ts1 = timeseriesStruct.ts1; ts1 = ts1(:);
ts2 = timeseriesStruct.ts2; ts2 = ts2(:);
ts3 = timeseriesStruct.ts3; ts3 = ts3(:);
ts4 = timeseriesStruct.ts4; ts4 = ts4(:);
ts5 = timeseriesStruct.ts5; ts5 = ts5(:);


dimT = length(ts1);
%%Step 2. generate connections between areas.
m1 = mean(ts1); s1 = std(ts1); ts1 = (ts1-m1)./s1;
m2 = mean(ts2); s2 = std(ts2); ts2 = (ts2-m2)./s2;
m3 = mean(ts3); s3 = std(ts3); ts3 = (ts3-m3)./s3;
m4 = mean(ts4); s4 = std(ts4); ts4 = (ts4-m4)./s4;
m5 = mean(ts5); s5 = std(ts5); ts5 = (ts5-m5)./s5;



%%connect areas 4 and 5 at zero lag
a45 = 0.25; 
%a45 = 0.0;
ts4 = (1-a45).*ts4+a45.*ts5;

%%connect areas 1 and 3, with area 1 connecting at a lag of 10.
a13 = 0.3; lag13 = 0; %a = connection strength
%a13 = 0.0; lag13 = 0;
temp1 = circshift(ts1,[lag13,0]);

%%connect areas 5 and 3 with area 5 connecting at a lag of 5 to area 3
a53 = 0.4; lag53 = 0; %a = connection strength
%a53 = 0.0; lag53 = 0;
temp5 = circshift(ts5,[lag53,0]);
ts3 = (1-a13-a53)*ts3+a53*temp5+a13*temp1 ;


%%normalize the time series to zero mean and unit var (so we can control
%%the variance when we add them to an area).
m1 = mean(ts1); s1 = std(ts1); ts1 = (ts1-m1)./s1;
m2 = mean(ts2); s2 = std(ts2); ts2 = (ts2-m2)./s2;
m3 = mean(ts3); s3 = std(ts3); ts3 = (ts3-m3)./s3;
m4 = mean(ts4); s4 = std(ts4); ts4 = (ts4-m4)./s4;
m5 = mean(ts5); s5 = std(ts5); ts5 = (ts5-m5)./s5;

%%plot them
figure(); 
subplot(5,1,1); plot(ts1); title('Ts1');
subplot(5,1,2); plot(ts2); title('Ts2');
subplot(5,1,3); plot(ts3); title('Ts3');
subplot(5,1,4); plot(ts4); title('Ts4');
subplot(5,1,5); plot(ts5); title('Ts5');

%%check the connections here
 figure(); crosscorr(ts1,ts2,50); title('ts1,ts2');
 figure(); crosscorr(ts1,ts3,50); title('ts1,ts3');
 figure(); crosscorr(ts1,ts4,50); title('ts1,ts4');
 figure(); crosscorr(ts1,ts5,50); title('ts1,ts5');
 figure(); crosscorr(ts2,ts3,50); title('ts2,ts3');
 figure(); crosscorr(ts2,ts4,50); title('ts2,ts4');
 figure(); crosscorr(ts2,ts5,50); title('ts2,ts5');
 figure(); crosscorr(ts3,ts4,50); title('ts3,ts4');
 figure(); crosscorr(ts3,ts5,50); title('ts3,ts5');
 figure(); crosscorr(ts4,ts5,50); title('ts4,ts5');

%%generate areas (make also areas 1-3 to overlap)
sFactor1 = 4;
sFactor2 = 9;
sFactor3 = 16;
sFactor4 = 9;
sFactor5 = 4;



[fooGrid1,area1] = generateArea(30,10,2,8,ts1,sFactor1,dimX,dimY,dimT,sMin);
[fooGrid2,area2] = generateArea(30,26,4,12,ts2,sFactor2,dimX,dimY,dimT,sMin);
[fooGrid3,area3] = generateArea(30,42,2,8,ts3,sFactor3,dimX,dimY,dimT,sMin);
[fooGrid4,area4] = generateArea(10,20,0.5,5,ts4,sFactor4,dimX,dimY,dimT,sMin);
[fooGrid5,area5] = generateArea(10,29,1,5,ts5,sFactor5,dimX,dimY,dimT,sMin);
linePlot(area1,area2,area3);
%%display the area sizes
areaSize = find(fooGrid1~=0); disp(['Area 1 size: ',num2str(length(areaSize))]);
areaSize = find(fooGrid2~=0); disp(['Area 2 size: ',num2str(length(areaSize))]);
areaSize = find(fooGrid3~=0); disp(['Area 3 size: ',num2str(length(areaSize))]);
areaSize = find(fooGrid4~=0); disp(['Area 4 size: ',num2str(length(areaSize))]);
areaSize = find(fooGrid5~=0); disp(['Area 5 size: ',num2str(length(areaSize))]);
%%display the area core sizes
areaSize = find(fooGrid1==1); disp(['Area 1 core size: ',num2str(length(areaSize))]);
areaSize = find(fooGrid2==1); disp(['Area 2 core size: ',num2str(length(areaSize))]);
areaSize = find(fooGrid3==1); disp(['Area 3 core size: ',num2str(length(areaSize))]);
areaSize = find(fooGrid4==1); disp(['Area 4 core size: ',num2str(length(areaSize))]);
areaSize = find(fooGrid5==1); disp(['Area 5 core size: ',num2str(length(areaSize))]);
%%plot the overlap map.
plotOverlapMap(fooGrid1,fooGrid2,fooGrid3,fooGrid4,fooGrid5);
%%superimpose the area time series to the grid.
theGrid = superimposeAreas2(area1,area2,area3,area4,area5);
linePlotGrid(theGrid);
gridNoNoise = theGrid;
% 
% %%plot the link maps (ground truth)
 mLag = 0; alpha = 0.001;
 lmap1 = getNetworkLink(fooGrid1, fooGrid2, fooGrid3, fooGrid4, fooGrid5, theGrid, mLag, alpha);
 figure(); pcolor(flipdim(lmap1,1)); shading flat; colorbar();
 lmap2 = getNetworkLink(fooGrid2, fooGrid1, fooGrid3, fooGrid4, fooGrid5, theGrid, mLag, alpha);
 figure(); pcolor(flipdim(lmap2,1)); shading flat; colorbar();
 lmap3 = getNetworkLink(fooGrid3, fooGrid1, fooGrid2, fooGrid4, fooGrid5, theGrid, mLag, alpha);
 figure(); pcolor(flipdim(lmap3,1)); shading flat; colorbar();
 lmap4 = getNetworkLink(fooGrid4, fooGrid1, fooGrid2, fooGrid3, fooGrid5, theGrid, mLag, alpha);
 figure(); pcolor(flipdim(lmap4,1)); shading flat; colorbar();
 figure(); lmap5 = getNetworkLink(fooGrid5, fooGrid1, fooGrid2, fooGrid3, fooGrid4, theGrid, mLag, alpha);
 pcolor(flipdim(lmap5,1)); shading flat; colorbar();

%%plot the variance of the time series
plotVarMap(theGrid); title('Variance without noise');

%final step is to add noise
for i = 1:dimX
    for j = 1:dimY
        wgn = randn(1,dimT); wgn = wgn(:);
        ts = theGrid(i,j,:); ts = ts(:);
        if(isnan(ts(1)))%%no area time series
            ts = wgn;
        else
            ts = ts+wgn;
        end;
        theGrid(i,j,:) = ts;
    end;
end;


%%variane plot again
plotVarMap(theGrid); title('Variance with noise');
linePlotGrid(theGrid);
plotCorrelationMap(30,10,theGrid);title('Area 1');
plotCorrelationMap(30,26,theGrid);title('Area 2');
plotCorrelationMap(30,42,theGrid);title('Area 3');
plotCorrelationMap(10,20,theGrid);title('Area 4');
plotCorrelationMap(10,29,theGrid);title('Area 5');

testFinalCorr(fooGrid1,fooGrid2,fooGrid3,fooGrid4,fooGrid5,theGrid);

end


function [fooGrid,area] = generateArea(x_c,y_c,rCore,rPeriphery,ts,sFactor,dimX,dimY,dimT,sMin)
%x_c, y_c: the area's central point
%rCore: radius of the area's periphery
%rPeriphery: area radius (rPeriphery > rCore)
%tsVar: variance of the area's time series.
%dimX,dimY: grid dimensions
if(rCore >= rPeriphery)
    disp('Warning: area radius is less or equal than area core');
end;
%step 1. Setup a foo grid and mark the area's points in there
fooGrid = zeros(dimX,dimY);
%%mark the area's core and periphery
for i = 1:dimX
    for j = 1:dimY
        if(euclidianDistance(i,j,x_c,y_c) <= rCore)
            fooGrid(i,j) = 1;
        elseif(euclidianDistance(i,j,x_c,y_c) <= rPeriphery)
            fooGrid(i,j) = 2;
        end;
    end;
end;

%%set the variance of the central time series
tsCenter = sqrt(sFactor).*ts;
%%construct an empty grid [dimX,dimY,dimT] that will hold the area.
area = zeros(dimX,dimY,dimT);
%%add the area's time series
minDecay = 10000;
decayMap = zeros(dimX,dimY);
for i = 1:dimX
    for j = 1:dimY
        if(fooGrid(i,j) == 1)%%grid cell in core, no decay
            area(i,j,:) = tsCenter;
        end;
        if(fooGrid(i,j) == 2)%%grid cell in periphery, need to add decay
            distance = euclidianDistance(x_c,y_c,i,j);
            decay = varianceDecay(distance,sFactor,sMin,rPeriphery,rCore);            
            if(minDecay > decay)
                minDecay = decay;
            end;
            decayMap(i,j) = decay;
            decayTs = (decay).*ts;
            area(i,j,:) = decayTs;
        end;
    end;
end;
figure();imagesc(decayMap); title('Decay factor, s'); colorbar();
disp(['Min Decay: ',num2str(minDecay)]);

varMap = zeros(dimX,dimY);
for i = 1:dimX
    for j = 1:dimY
        varMap(i,j) = var(area(i,j,:));
    end;
end;
figure(); imagesc(fooGrid);
figure(); imagesc(varMap(:,:,1));colorbar();


end


function d = varianceDecay(x,s2,smin2,rp,rc)

d = sqrt( (s2-smin2)/(rc-rp)*x + (smin2*rc -s2*rp)/(rc-rp) );

end


function dist = euclidianDistance(x1,y1,x2,y2)
dist = sqrt( (x1-x2)^2+(y1-y2)^2);
end


function theGrid = superimposeAreas2(area1,area2,area3,area4,area5)

dimX = size(area1,1);
dimY = size(area1,2);
dimT = size(area1,3);

theGrid = zeros(dimX,dimY,dimT);

for i = 1:dimX
    for j = 1:dimY
        ts1 = area1(i,j,:); ts1 = ts1(:);
        ts2 = area2(i,j,:); ts2 = ts2(:);
        ts3 = area3(i,j,:); ts3 = ts3(:);
        ts4 = area4(i,j,:); ts4 = ts4(:);
        ts5 = area5(i,j,:); ts5 = ts5(:);
        theGrid(i,j,:) = ts1+ts2+ts3+ts4+ts5;
    end;
end;

end

function theGrid = superimposeAreas(area1,area2,area3,area4,area5)

dimX = size(area1,1);
dimY = size(area1,2);
dimT = size(area1,3);

theGrid = zeros(dimX,dimY,dimT);

for i = 1:dimX
    for j = 1:dimY
        ts1 = area1(i,j,:); ts1 = ts1(:);
        ts2 = area2(i,j,:); ts2 = ts2(:);
        ts3 = area3(i,j,:); ts3 = ts3(:);
        ts4 = area4(i,j,:); ts4 = ts4(:);
        ts5 = area5(i,j,:); ts5 = ts5(:);
        
        c1 = getCFactor(ts1,ts2,ts3,ts4,ts5);
        c2 = getCFactor(ts2,ts1,ts3,ts4,ts5);
        c3 = getCFactor(ts3,ts1,ts2,ts4,ts5);
        c4 = getCFactor(ts4,ts1,ts2,ts3,ts5);
        c5 = getCFactor(ts5,ts1,ts2,ts3,ts4);
        
        newTs = c1.*ts1+ c2.*ts2+ c3.*ts3+ c4.*ts4+ c5.*ts5;
        if(isnan(newTs(1)))
            newTs = zeros(1,dimT); newTs = newTs(:);
        end;
        theGrid(i,j,:)=newTs;
    end;
end;

end

function cx = getCFactor(tsx,ts1,ts2,ts3,ts4)

cx= sqrt( var(tsx)/( var(tsx) + var(ts1) + var(ts2) + var(ts3) + var(ts4) ) );

end




function plotVarMap(theGrid)
dimX = size(theGrid,1);
dimY = size(theGrid,2);
varMap = zeros(dimX,dimY);
for i = 1:dimX
    for j = 1:dimY
        ts = theGrid(i,j,:); ts = ts(:);
        varMap(i,j) = var(ts);
    end;
end;
figure();
imagesc(varMap); colorbar();
end

function omap =  plotOverlapMap(area1,area2,area3,area4,area5)
dimX = size(area1,1);
dimY = size(area2,2);
omap = zeros(dimX,dimY);
for i = 1:dimX
    for j = 1:dimY
        if(area1(i,j)~=0)
            if(omap(i,j)==0)
                omap(i,j) = area1(i,j);
            else
                omap(i,j) = -1;
            end;
        end;
    
        if(area2(i,j)~=0)
            if(omap(i,j)==0)
                omap(i,j) = area2(i,j);
            else
                omap(i,j) = -1;
            end;
        end;
        
        if(area3(i,j)~=0)
            if(omap(i,j)==0)
                omap(i,j) = area3(i,j);
            else
                omap(i,j) = -1;
            end;
        end;
        
        if(area4(i,j)~=0)
            if(omap(i,j)==0)
                omap(i,j) = area4(i,j);
            else
                omap(i,j) = -1;
            end;
        end;
        
        if(area5(i,j)~=0)
            if(omap(i,j)==0)
                omap(i,j) = area5(i,j);
            else
                omap(i,j) = -1;
            end;
        end;
    end;
end;
figure();imagesc(omap);

end

function plotCorrelationMap(x,y,theGrid)
dimX = size(theGrid,1);
dimY = size(theGrid,2);
corrmap = zeros(dimX,dimY);
tsFrom = theGrid(x,y,:);
tsFrom = tsFrom(:);
for i = 1:dimX
    for j = 1:dimY
        if(i~=x || j~=y)
            tsTo = theGrid(i,j,:); tsTo = tsTo(:);
            corrmap(i,j) = corr(tsFrom,tsTo);
        end;
    end;
end;
figure();
imagesc(corrmap); colorbar(); colormap(jet(10));
end

function testFinalCorr(foo1,foo2,foo3,foo4,foo5,theGrid)
dimX = size(theGrid,1);
dimY = size(theGrid,2);
dimT = size(theGrid,3);
ts1 = zeros(1,dimT);ts1 = ts1(:);
ts2 = zeros(1,dimT);ts2 = ts2(:);
ts3 = zeros(1,dimT);ts3 = ts3(:);
ts4 = zeros(1,dimT);ts4 = ts4(:);
ts5 = zeros(1,dimT);ts5 = ts5(:);
cc1 = 0;
        cc2 = 0; 
        cc3 = 0;
        cc4 = 0;
        cc5 = 0;
for i = 1:dimX
    for j = 1:dimY
        
        if(foo1(i,j)~=0)
            ts = theGrid(i,j,:);ts = ts(:);
            ts1 = ts1+ts;
            cc1 = cc1+1;
        end;
        if(foo2(i,j)~=0)
            ts = theGrid(i,j,:); ts = ts(:);
            ts2 = ts2+ts;
            cc2 = cc2+1;
        end;
        if(foo3(i,j)~=0)
            ts = theGrid(i,j,:); ts = ts(:);
            ts3 = ts3+ts;
            cc3 = cc3+1;
        end;
        if(foo4(i,j)~=0)
            ts = theGrid(i,j,:); ts = ts(:);
            ts4 = ts4+ts;
            cc4 = cc4+1;
        end;
        if(foo5(i,j)~=0)
            ts = theGrid(i,j,:); ts = ts(:);
            ts5 = ts5+ts;
            cc5=cc5+1;
        end;
    end;
end;

figure();crosscorr(ts1,ts3,50); title('ts1,ts3');
figure();crosscorr(ts3,ts5,50); title('ts3,ts5');
figure();crosscorr(ts1,ts2,50); title('ts1,ts2');
figure();crosscorr(ts2,ts3,50); title('ts2,ts3');
figure();crosscorr(ts4,ts5,50); title('ts4,ts5');

end

function linkMap = getNetworkLink(areaFrom,area2,area3,area4,area5,theGrid,maxLag,alpha)

dimX = size(theGrid,1);
dimY = size(theGrid,2);
dimT = size(theGrid,3);

%%construct area cumulative anomalies
areaCaFrom = zeros(1,dimT); areaCaFrom = areaCaFrom(:);
areaCa2 = zeros(1,dimT); areaCa2 = areaCa2(:);
areaCa3 = zeros(1,dimT); areaCa3 = areaCa3(:);
areaCa4 = zeros(1,dimT); areaCa4 = areaCa4(:);
areaCa5 = zeros(1,dimT); areaCa5 = areaCa5(:);
for i = 1:dimX
    for j = 1:dimY
        if(areaFrom(i,j)~=0)
            ts = theGrid(i,j,:); ts = ts(:);
            areaCaFrom = areaCaFrom+ts;
        end; 
        if(area2(i,j)~=0)
            ts = theGrid(i,j,:); ts = ts(:);
            areaCa2 = areaCa2+ts;
        end;
        
        if(area3(i,j)~=0)
            ts = theGrid(i,j,:); ts = ts(:);
            areaCa3 = areaCa3+ts;
        end;
        
        if(area4(i,j)~=0)
            ts = theGrid(i,j,:); ts = ts(:);
            areaCa4 = areaCa4+ts;
        end;        
        if(area5(i,j)~=0)
            ts = theGrid(i,j,:); ts = ts(:);
            areaCa5 = areaCa5+ts;
        end;
    end;
end;
%%done now compute the links from areaFrom to all other areas.
%%areaFrom to area2
sigCorrs = getCorrelogram(areaCaFrom,areaCa2,maxLag,alpha,1);
ind = find(~isnan(sigCorrs));
if(~isempty(ind))
    corrs = sigCorrs(ind);
    bestCor = max(corrs);
    lw2 = std(areaCaFrom)*std(areaCa2)*bestCor;
else
    lw2 = 0;
end;
%to area 3...
sigCorrs = getCorrelogram(areaCaFrom,areaCa3,maxLag,alpha,1);
ind = find(~isnan(sigCorrs));
if(~isempty(ind))
    corrs = sigCorrs(ind);
    bestCor = max(corrs);
    lw3 = std(areaCaFrom)*std(areaCa3)*bestCor;
else
    lw3 = 0;
end;

%%to area 4
sigCorrs = getCorrelogram(areaCaFrom,areaCa4,maxLag,alpha,1);
ind = find(~isnan(sigCorrs));
if(~isempty(ind))
    corrs = sigCorrs(ind);
    bestCor = max(corrs);
    lw4 = std(areaCaFrom)*std(areaCa4)*bestCor;
else
    lw4 = 0;
end;

%%and to area 5
sigCorrs = getCorrelogram(areaCaFrom,areaCa5,maxLag,alpha,1);
ind = find(~isnan(sigCorrs));
if(~isempty(ind))
    corrs = sigCorrs(ind);
    bestCor = max(corrs);
    lw5 = std(areaCaFrom)*std(areaCa5)*bestCor;
else
    lw5 = 0;
end;

%%finally construct the link map
linkMap = zeros(dimX,dimY);
if(lw2 > 0)
    for i = 1:dimX
        for j = 1:dimY
            if(area2(i,j)~=0)
                linkMap(i,j) = lw2;
            end;
        end;
    end;
end;

if(lw3 > 0)
    for i = 1:dimX
        for j = 1:dimY
            if(area3(i,j)~=0)
                linkMap(i,j) = lw3;
            end;
        end;
    end;
end;

if(lw4 > 0)
    for i = 1:dimX
        for j = 1:dimY
            if(area4(i,j)~=0)
                linkMap(i,j) = lw4;
            end;
        end;
    end;
end;

if(lw5 > 0)
    for i = 1:dimX
        for j = 1:dimY
            if(area5(i,j)~=0)
                linkMap(i,j) = lw5;
            end;
        end;
    end;
end;

for i = 1:dimX
    for j = 1:dimY
        if(areaFrom(i,j)~=0)
            linkMap(i,j) = NaN;
        end;
    end;
end;

end

function linePlotGrid(grid)
dimY = size(grid,2);
myX = 30;
mvar = zeros(1,dimY);
for i = 1:dimY
    ts = grid(myX,i,:);
    mvar(i) = var(ts);
end;
figure(); plot(mvar);
end

function linePlot(area1,area2,area3)

dimY = size(area1,2);
myX = 30;
var1 = zeros(1,dimY);
var2 = zeros(1,dimY);
var3 = zeros(1,dimY);

for i = 1:dimY
    ts1 = area1(myX,i,:); 
    myVar = var(ts1);
    var1(i) = myVar;
    ts2 = area2(myX,i,:);
    myVar = var(ts2);
    var2(i) = myVar;
    ts3 = area3(myX,i,:);
    myVar = var(ts3);
    var3(i) = myVar;
end;
figure(); 
hold on;
% ind = var1 == 0;
% var1(ind) = NaN;
% ind = var2 == 0;
% var2(ind) = NaN;
% ind = var3 == 0;
% var3(ind) = NaN;

plot(var1,'+-b', 'LineWidth',2);
plot(var2,'+-r', 'LineWidth',2);
plot(var3,'+-k', 'LineWidth',2);
ylabel('Variance');
legend({'Area 1','Area 2', 'Area 3'});




end

