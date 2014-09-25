% Andy Gallagher

function [GraphIds,im_,GR, Res, Hits, GI, Blocks,time1,time2] = GreedyAssembly5R(greedyOrder,ap,SCO, nr,nc,mni)
% SHOULD ENFORCE THE SIZE OF THE PUZZLE! ENABLE SEEDING OF THE PIECES? 
% ENFORCE THE no-overlap policy.
% enforce single piece additions?
% Quads first, then singles?
% 
% locked means that I know how big the puzzle is... 
%   Remember, I need a "buffer row" on the top and bottom and left and
%   right. 
%
% givenPieces = []; Indexes of pieces that are in the correct postion 
% to start with. Assumes "locked" aspect ratio... 
%

%    addLots = 0; % instead of iterating, just add singles one after another
                % after the initial QuadGrid is formed. 
%    v1v2Flag = 1; %choose v1 for giving preference to 4surrounded neighbors, etc. 
                % v2 is just the most confidenct. (when adding next single)


Hits =[0 0]; % how many proposed block merges are rejected because of overlap?


% adjust the score matrix to be normalized. Each a,b score is divided by
% the second best match (to either A or B) on that side of the piece. This 
% will need to be changed for other rotations. 

normSCO = SCO; 

mirror = [ 3 4 1 2  16 13 14 15  9 10 11 12 6 7 8 5  ]; %what is the symmetry? %this will cut processing time in half. 


t = 0.000000000000001; 
% NEED TO SPEED THIS UP A LOT! 
% Normalize the score matrix so that each element is the ``confidence''
% of a match. (i.e. P(match| features a b). 
% for ii = 1:1:size(SCO,3) %over each possible arrangement.
%     ii
% for jj = 1:1:size(SCO,1) %over each row. 
% for kk = 1:1:size(SCO,2) %over each possible arrangement. 
%     value = SCO(jj,kk,ii); 
%     n1 = min(SCO( jj, list~=kk,ii)); %best nonsame 
%     n2 = min(SCO( list~=jj, kk,ii)); 
%     nval = (value+t)/(min(n1,n2)+t);    % closest match 
%     normSCO(jj,kk,ii) = nval; 
% end
% end
% end %the easy one to program, but slow
start_time = cputime;
for ii = 1:1:size(SCO,3) %over each possible arrangement.
    fprintf('Processing Scores Matrix %d\n',ii);
    %    ii
    [aaa,bbb] = sort(SCO(:,:,ii),2); % sorted over each row.
    rowmins = aaa(:,1:2); %2smallest over each row.
    rowminloc = bbb(:,1); %location of minimum in each row.
    
    [aaa,bbb] = sort(SCO(:,:,ii)); % sorted over each column.
    colmins = aaa(1:2,:); %2smallest over each row.
    colminloc = bbb(1,:); %location of minimum in each row.
    
    for jj = 1:1:size(SCO,1) %over each row.
        values = SCO(jj,:,ii);% the values in the row.
        n1 = values.*0 + rowmins(jj,1); %the minimum for that row...
        n1(rowminloc) = rowmins(jj,2); % the second lowest
        %each position can also be replaced by the smallest nonsame value in  the column...
        n2 = values.*0 + colmins(1,:); %the smallest value in each column.
        % but whereever the row is the same, we use the second lowest value instead
        n2(jj==colminloc)= colmins(2,jj==colminloc);
        nval = (values+t)./(min([n1;n2])+t);
        normSCO(jj,:,ii) = nval;
    end
end



Blocks = {}; %empty list of blocks for initialization
Rots =[]; %keep track of the rotation of pieces within the blocks. 
iters = 1; 
gogo = 1; 

%VECTORIZED = normSCO(:);

ST = 1.25; %the stopping threshold... 
pairlist = [];
%while gogo

%load('/home/hevc/Desktop/Lajan_Research/lajan_puzzle/mextest/greedyOrder');
ord = 1;
numPairs = size(greedyOrder,1);
%load('../mextest/greedyOrd');
%greedyOrder = greedyOrd;
%numPairs = 491;

tf = [4 2 1 3];
pathsDone = false;
%while gogo
scount = 0;
while ord <= numPairs

    if ord > numPairs
        pathsDone = true;
    end

    if pathsDone 
        [aa,BB] = min(normSCO(:)); % the most confident remaining match
        
        if(isnan(aa) || aa>ST )
            gogo = 0;
            break;
        end

        [R,C,How] = ind2sub(size(normSCO),BB);
    else    
        %[R,C,How] = ind2sub(size(normSCO),BB);
        pair = greedyOrder(ord,:);
        R = pair(1);
        C = pair(3);
        How = tf(pair(2));
        ord = ord + 1;
    end
    
    
    P1 = R; P2 = C;
    %R and C are the next suggested pieces to join.
    Rb = R;
    Cb = C;
    
    Rr = 1; Cr = 1;
    %find what blocks they belong to.
    rmems = zeros(numel(Blocks),1);
    cmems = zeros(numel(Blocks),1);
    findr = [];
    findc = [];
    for ii = 1:1:numel(Blocks);
        rmems(ii)= ismember(R, Blocks{ii});
        cmems(ii)= ismember(C, Blocks{ii});
    end
    if(sum(rmems))
        findr = find(rmems);
        Rb = Blocks{findr};
        Rr = Rots{findr};
    end
    if(sum(cmems))
        findc = find(cmems);
        Cb = Blocks{findc};
        Cr = Rots{findc};
    end
    

    % okay, now join the pieces to get a new piece:
    
    if(numel(findr==1) && numel(findc==1))
        if(findr~=findc)
            [b r  s]=joinPiecesR(Rb,Cb,Rr,Cr,R,C,How);
            BSize = sum(b(:)>0);
        else
            s=0; BSize=0;
        end
    else
        [b r s]=joinPiecesR(Rb,Cb,Rr,Cr,R,C,How);
        
        BSize = sum(b(:)>0);
    end
    Hits = Hits+[s==1 s==0];
    %b
    if(s==1)

        scount = scount + 1;
        pairlist(scount,:) = [R C];
        % b is the new puzzle piece.
        if(numel(findr))
            Blocks{findr} =b;
            Rots{findr} =r;
            if(numel(findc))
                Blocks{findc} = [];
                Rots{findc} =[];
            end
        elseif(numel(findc)) %the findr is empty (meaning first piece is a single)
            Blocks{findc} =b;
            Rots{findc} =r;
        else % both pieces were singles.
            Blocks{end+1} = b;
            Rots{end+1} =r;
        end
        
        %%% make a frame: 
        %Blocks
        %frame_image = makeFrameImage(Blocks,Rots,ap); 
        %imshow(frame_image);
        %pause(0.2);
        %if(WRITE)
        %    imwrite(uint8(frame_image),strcat('ThankYouSeq/V2_', num2str(scount), '.jpg'),'Quality',95); 
        %end
        %Blocks{end+1} = b; %this is the new piece.
        %now, update the SCO matrix for the "taken" pieces.
        normSCO(P1,:,How) = NaN;
        normSCO(:,P2,How) = NaN;
        
        %HowN = mod(How-1+2,4) +1;
        HowN  = mirror(How); %get the dual situation...
        normSCO(P2,:,HowN) = NaN;
        normSCO(:,P1,HowN) = NaN;
        
        
        %%% Nan out all piecewise combinations of pieces
        if((numel(Rb)>1)||(numel(Cb)>1))
            %%%
            Rb1 = Rb(Rb>0);
            Cb1 = Cb(Cb>0);
            
            for i1 = 1:1:numel(Rb1)
                iii1 = Rb1(i1);
                for i2 = 1:1:numel(Cb1)
                    iii2 = Cb1(i2);
                    normSCO(iii1,iii2,:) = NaN;
                    normSCO(iii2,iii1,:) = NaN;
                end
            end
        end
        % Nan out blocks that
        
        %delets the puzzle chunks that were used to form the new one.
        %    if(numel(findr))
        %       Blocks{findr} =[];
        %    end
        %    if(numel(findc))
        %        Blocks{findc} = [];
        %    end
    else  %unsuccessful match (e.g. pieces actually overlap)
        %HowN = mod(How-1+2,4) +1;
        HowN  = mirror(How); %get the dual situation...
        normSCO(P1,P2,How) = NaN;
        normSCO(P2,P1,HowN) = NaN;
        
    end
    %VECTORIZED(BB)=NaN;
    if(BSize> size(normSCO,1)-1)
        gogo=0;
    end
    
    iters=iters+1;
    
    if(floor(iters/100)*100 == iters)
        fprintf('%d\t%d %d\t%d %d %d %d %d %d %d %d \n',iters,  R ,C ,findr, findc, How, sum(Rb(:)>0), sum(Cb(:)>0), BSize, numel(Blocks), sum(normSCO(:)>0))
        %fprintf('%.2f \n',aa);
        %[iters  R C findr findc How sum(Rb(:)>0) sum(Cb(:)>0) BSize numel(Blocks) sum(normSCO(:)>0)]
        
        %[aa]
    end
    
%BlockSize = zeros(numel(Blocks),1);
%for iii = 1:1:numel(Blocks)
%    BlockSize(iii) = sum(Blocks{iii}(:)>0);     
%end
%[~,bb] = max(BlockSize); 
%GI = Blocks{bb};
%renderBlocks(GI,28,ap);
%pause
end
time1 = cputime - start_time;

% find the biggest block. This is the completed puzzle. 
BlockSize = zeros(numel(Blocks),1);
for iii = 1:1:numel(Blocks)
    BlockSize(iii) = sum(Blocks{iii}(:)>0);     
end
[~,bb] = max(BlockSize); 
GI = Blocks{bb};
GR = Rots{bb}; %the corresponding rotations.


start_time = cputime;
numOfParts = nr*nc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ct = 1;
[dests,paths] = findPathsLen(mni,numOfParts,3);
locns = dests;
numDests = size(dests,1);
%Parts = [Prts];
Parts = GI(find(GI))';
%occupied_locns = Ws;
[r c] = find(GI);
occupied_locns = [r c];
remaining_parts = setdiff(1:numOfParts,Parts);
%Locn = squeeze(paths(Parts,remaining_parts,:,:));
Locn = [];

locns = [];
for t = 1:size(occupied_locns,1)
  loc = occupied_locns(t,:);
  P = Parts(t);
  for i = 1:numDests
    %L = dests(i,:) + loc;
    L = loc - dests(i,:);
    if ismember(L,occupied_locns,'rows'), continue; end
    if ~isempty(locns) && ismember(L,locns,'rows')
      ind = find(ismember(locns,L,'rows'));
      Locn(:,ind,:) = Locn(:,ind,:) + shiftdim(paths(P,remaining_parts,i,:));
    else
      locns = [locns;L];
      Locn = cat(2,Locn,shiftdim(paths(P,remaining_parts,i,:)));
    end
  end
end

done = false(size(mni));
while true
  typePths = sum(Locn > 0,3);
  TypePths = typePths;
  Lc = sum(Locn,3);
  probs1 = Lc ./ repmat(sum(Lc,2),1,size(Lc,2));
  probs2 = Lc ./ repmat(sum(Lc,1),size(Lc,1),1);

  if max(typePths(:)) <= 1, break; end

  if(size(typePths,1) == 0) break; end
  msk = typePths == max(typePths(:));
  typePths = msk .* (probs1 .* probs2);
  %typePths = (probs1 .* probs2);

  [u,v] = max(typePths(:));

  [part_ind,loc_ind] = ind2sub(size(typePths),v);

  if TypePths(part_ind,loc_ind) <= 1, break; end
  %if ismember(P,Parts), Locn(P,loc_ind,:) = 0; continue; end
  P = remaining_parts(part_ind);
  loc = locns(loc_ind,:);
  Parts = [Parts P];
  occupied_locns = [occupied_locns;loc];

  rem_inds = remaining_parts ~= P;
  remaining_parts = remaining_parts(rem_inds);


  [dests,paths] = findPathsLen_omit(mni,numOfParts,3,Parts);
  n_dests = size(dests,1);
  paths = paths(Parts,:,:,:);
  occ_loc_enc = occupied_locns(:,1) + 2*numOfParts*occupied_locns(:,2);
  dests_enc = dests(:,1) + 2*numOfParts*dests(:,2);
  %map = repmat(occ_loc_enc,1,n_dests) + repmat(dests_enc',length(occ_loc_enc),1);
  map = repmat(occ_loc_enc,1,n_dests) - repmat(dests_enc',length(occ_loc_enc),1);
  rem_locs = unique(setdiff(map,occ_loc_enc));
  paths_re = permute(paths,[1 3 2 4]);
  paths_lin = reshape(paths_re,numel(map),numOfParts,4);
  map = map(:);
  locns = [rem_locs - 2*numOfParts*round(rem_locs/(2*numOfParts)) round(rem_locs/(2*numOfParts))];
  size(locns)
  Locn = zeros(length(rem_locs),numOfParts,4);
  for i = 1:length(rem_locs)
    Locn(i,:,:) = sum(paths_lin(find(map == rem_locs(i)),:,:),1);
  end
  Locn = Locn(:,remaining_parts,:);
  Locn = permute(Locn,[2 1 3]);

  %Locn = Locn(rem_inds,:,:);
  %id = ~ismember(locns,loc,'rows');
  %locns = locns(id,:);
  %Locn = Locn(:,id,:);
  %for i = 1:numDests
  %  L = loc - dests(i,:);% + loc;
  %  if ismember(L,occupied_locns,'rows'), continue; end
  %  if ismember(L,locns,'rows')
  %    ind = find(ismember(locns,L,'rows'));
  %    if (numel(remaining_parts) == 1), break; end
  %    Locn(:,ind,:) = Locn(:,ind,:) + shiftdim(paths(P,remaining_parts,i,:));
  %  else
  %    locns = [locns;L];
  %    Locn = cat(2,Locn,shiftdim(paths(P,remaining_parts,i,:)));
  %  end
  %end

  %if length(Parts) == 431, break; end
  %if rem(length(Parts),30) == 0,
  %  [mni,done,cost_ranks] = updatemni(mni,done,Parts,occupied_locns,cost_ranks);
  %  [~,paths] = findPathsLen(mni,numOfParts,4);
  %end
end

%Prts = Parts;
%Ws = occupied_locns;
%%W = repmat(max(Ws,[],1) + 1,size(Ws,1),1) - Ws;
%W = Ws - repmat(min(Ws,[],1),size(Ws,1),1) + repmat([1 1],size(Ws,1),1);
%GI = zeros(max(W,[],1));
%GI(sub2ind(size(GI),W(:,1),W(:,2))) = Prts;

time2 = cputime - start_time;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%for it = 1:60
%[dests,paths] = findPaths2(mni,1728,1);
%while true
%
%    assembled = GI(GI>0);
%    [r,c] = ind2sub(size(GI),find(GI));
%    loc = [r c];
%    remaining = setdiff(1:1728,assembled);
%    if isempty(remaining), break; end
%    Pths = paths(GI(find(GI)),remaining,:);
%    
%    locns = [];
%    Locn = [];
%    for i = 1:size(Pths,1)
%        for j = 1:size(Pths,3)
%            D = loc(i,:) - dests(j,:);
%            if ismember(D,loc,'rows'), continue; end
%            id = [];
%            if ~isempty(locns), id = find(ismember(locns,D,'rows')); end
%            if id, Locn(id,:) = Locn(id,:) + (Pths(i,:,j)>0);
%            else 
%                locns = [locns;D];
%                Locn(end+1,:) = Pths(i,:,j);
%            end
%        end
%    end
%    
%    [u,v] = max(Locn(:));
%    if u == 0, break; end
%    [p,q] = ind2sub(size(Locn),v);
%    %size(Locn)
%    %[p q]
%    %remaining(q)
%    l = locns(p,:)
%    %pause
%    %if l(1) > 0 && l(2) > 0 && l(1) <= size(GI,1) && l(2) <= size(GI,2)
%    %    GI(l(1),l(2)) = remaining(q);
%    %else
%        if l(1) <= 0,           GI = [zeros(1,size(GI,2));GI]; end
%        if l(2) <= 0,           GI = [zeros(size(GI,1),1) GI]; end
%        if l(1) > size(GI,1),   GI = [GI;zeros(1,size(GI,2))]; end
%        if l(2) > size(GI,2),   GI = [GI zeros(size(GI,1),1)]; end
%    %end
%    if l(1) == 0, l(1) = 1; end
%    if l(2) == 0, l(2) = 1; end
%    GI(l(1),l(2)) = remaining(q);
%
%end

GR = zeros(size(GI));
GR(GI>0) = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% can "derotate" the pieces here... 
%%% 
% derotate the pieces: 
            
for jj = 1:1:numel(GR)
    if(GR(jj)>0)
        rv = GR(jj); 
        rotateCCW = rv-1; 
        bid = GI(jj); 
        ap{bid} = imrotate(ap{bid},90*rotateCCW);
    end
end

[im_,GraphIds] = renderPiecesFromGraphIDs(ap,GI,0);
%evaluate the puzzle piece fitting together: 
[Res] = EvalPuzzleAssembly(GraphIds, nc, nr);

function renderBlocks(GI,PP,ap)
    im = uint16(zeros([size(GI)*PP 3]));
    for i = 1:size(GI,1)
        for j = 1:size(GI,2)
            if GI(i,j) == 0, continue; end
            im(PP*(i-1)+(1:PP),PP*(j-1)+(1:PP),:) = ap{GI(i,j)};
        end
    end
    %figure,imshow(im);
    imshow(im);
