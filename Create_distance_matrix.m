
%% Creates a grid of points that form the backbone of the distance calculations
DA = 0;
if DA == 1
   clear all
   load crcb_domain_Solomons reef_centroid Landraw LandOutline
   load RegionOutline G*
   
   ReefsInStudyArea = inpolygon(reef_centroid(:,1),reef_centroid(:,2),G2(:,1),G2(:,2));
   
   figure(1), clf, hold on, box on
   plot(reef_centroid(:,1),reef_centroid(:,2),'g.')
   plot(reef_centroid(ReefsInStudyArea,1),reef_centroid(ReefsInStudyArea,2),'r.')
   plot(Landraw(:,1),Landraw(:,2),'b')
   
   % Create a coarse-scale grid for the whole region
   gr = 40;
   Lon = linspace(min(reef_centroid(:,1)),max(reef_centroid(:,1)),gr);
   Lat = linspace(min(reef_centroid(:,2)),max(reef_centroid(:,2)),gr+1);
   [LN,LT] = meshgrid(Lon,Lat); LN = LN(:); LT = LT(:);
   INSIDE = inpolygon(LN,LT,G(:,1),G(:,2));
   LT(INSIDE==0) = []; LN(INSIDE==0) = [];
   INSIDE = inpolygon(LN,LT,G2(:,1),G2(:,2));
   LT(INSIDE) = []; LN(INSIDE) = [];
   Coarse = ones(length(LT),1);
   
   % Distance along diagonal between grid points
   LongH = sqrt((Lon(2)-Lon(1)).^2 + (Lat(2)-Lat(1)).^2);
   LongO = sqrt((Lat(2)-Lat(1)).^2);
   LongA = sqrt((Lon(2)-Lon(1)).^2);
   
   % Create a fine-scale grid for the sample region
   gr2 = 55;
   Lon2 = linspace(157.7478,159.2126,gr2);
   Lat2 = linspace(-7.8857,-7.1933,gr2+1);
   [LN2,LT2] = meshgrid(Lon2,Lat2); LN2 = LN2(:); LT2 = LT2(:);
   INSIDE2 = inpolygon(LN2,LT2,G2(:,1),G2(:,2));
   LT2(INSIDE2==0) = []; LN2(INSIDE2==0) = [];
   Coarse = [Coarse; zeros(length(LT2),1)];
   
   % Distance along diagonal between grid points
   ShortH = sqrt((Lon2(2)-Lon2(1)).^2 + (Lat2(2)-Lat2(1)).^2);
   ShortO = sqrt((Lat2(2)-Lat2(1)).^2);
   ShortA = sqrt((Lon2(2)-Lon2(1)).^2);
   
   LT = [LT; LT2]; LN = [LN; LN2];
   plot(LN,LT,'ro','markersize',3)
   axis tight
   
   %% Exclude any points that fall within an island outline
   for a = 1:length(LandOutline)
      %    if mod(a,100) == 0; disp(a/length(LandOutline)); end
      IN = inpolygon(LN,LT,LandOutline{a}(:,1),LandOutline{a}(:,2));
      if sum(IN) > 0
         LN(IN) = []; LT(IN) = []; Coarse(IN) = [];
      end
   end
   plot(LN,LT,'ko','markersize',3); % LN and LT are not on land
   
   % Go through the points and identify the closest and their distance in a network matrix
   % We need to do this cheaply, because there are so many grid points. We do this by only
   % using a set of closest points, and leaving the rest of the connections (either those
   % that go across land, or those that are too far away, as NaNs.
   GRD = inf.*ones(length(LN));
   for i = 1:length(LT)
      % How far from this grid point to all other grid points?
      Distances = pdist2([LN(i) LT(i)],[LN LT]);
      
      % Go through all nearby distances to find the shortest path
      
      % If we're in coarse space
      if Coarse(i) == 1
         Near = [find(abs(Distances - LongH)  < 5e-3) ...
            find(abs(Distances - LongO)  < 5e-3) ...
            find(abs(Distances - LongA)  < 5e-3)];
      else % If we're in fine space
         Near = [find(abs(Distances - ShortH) < 5e-3) ...
                 find(abs(Distances - ShortO) < 5e-3) ...
                 find(abs(Distances - ShortA) < 5e-3)];
         % If we're on the edge of fine space
         if LN(i) <= Lon2(2) | LN(i) >= Lon2(end-2) | LT(i) <= Lat2(2) | LT(i) >= Lat2(end-2)
            Near = find(Distances < LongH.*1.01);
         end
      end
      
      
      if isempty(Near) == 1
         keyboard
      end
      for j = 1:length(Near)
         GRD(i,Near(j)) = Distances(Near(j));
         GRD(Near(j),i) = Distances(Near(j));
      end
   end
   A = GRD < inf;
   
   % Use Dijkstra's algorithm to track shortest paths through this network
   tic; [costs,paths] = dijkstra(A,GRD); toc
   
   %% The previous results track distances among the grid points shown. We therefore need to find the
   %%     closest grid point to each reef, and add that extra distance to the total.
   
   NumReefs = length(reef_centroid);
   for i = 1:NumReefs
      DistanceToGrid = pdist2(reef_centroid(i,1:2),[LN LT]);
      [DistanceToClosestGridPoint(i),ClosestGridPoint(i)] = min(DistanceToGrid);
   end
   
   DistanceMatrix = zeros(NumReefs);
   for i = 1:NumReefs
      for j = i+1:NumReefs
         % I don't really care about distances that don't include at least one reef in the study area
         if ReefsInStudyArea(i) == 1 | ReefsInStudyArea(j) == 1
            
            DistanceMatrix(i,j) = DistanceToClosestGridPoint(i) ...
               + DistanceToClosestGridPoint(j) ...
               + costs(ClosestGridPoint(i),ClosestGridPoint(j));
            
            % If the direct distance between the two points is less than the fine grid size
            if sqrt(sum((reef_centroid(i,:)-reef_centroid(j,:)).^2)) < ShortH
               % Just head there directly
               DistanceMatrix(i,j) = sqrt(sum((reef_centroid(i,:)-reef_centroid(j,:)).^2));
            end
            
            DistanceMatrix(j,i) = DistanceMatrix(i,j);
         end
      end
   end
   
   save TEMP_distances
end

clear all
load TEMP_distances reef* Close* Num* L* paths G2

%% This section plots routes

for a = 1:4000
   Reef1 = ceil(rand*NumReefs); ST = ClosestGridPoint(Reef1);
   Reef2 = ceil(rand*NumReefs); ND = ClosestGridPoint(Reef2);
   if inpolygon(reef_centroid(Reef1,1),reef_centroid(Reef1,2),G2(:,1),G2(:,2)) == 1
      
      figure(1), clf, hold on, box on
      plot(Landraw(:,1),Landraw(:,2),'b')
      plot(LN,LT,'k.','markersize',5)
      plot(reef_centroid(:,1),reef_centroid(:,2),'g.','markersize',4)
      
      plot(LN(paths{ST,ND}),LT(paths{ST,ND}),'k','linewidth',3)
      plot(LN(paths{ST,ND}),LT(paths{ST,ND}),'m.','markersize',14)
      plot(reef_centroid(Reef1,1),reef_centroid(Reef1,2),'b.','markersize',20)
      plot(reef_centroid(Reef2,1),reef_centroid(Reef2,2),'r.','markersize',20)
      pause
   end
end
































