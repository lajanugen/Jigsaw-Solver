function [cor,timeStats] = puzzle();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% This demo code creates and solves type 1 and type 2 jigsaw puzzles from
%%% an image. 
%%%
%%%
%%% Copyright (C) 2012, Andrew Gallagher
%%% This code is distributed with a non-commercial research license.
%%% Please see the license file license.txt included in the source directory.
%%%
%%% Please cite this paper if you use this code: 
%%% 
%%% 
%%% Jigsaw Pieces with Pieces of Unknown Orientation, Andrew Gallagher,
%%% CVPR 2012. 
%%% 
%%% Andrew Gallagher
%%% Aug. 15, 2012. 
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
imNum = 5;
cor = [];
timeStats = [];
for imNo = 15:20

%Demo Really Big Puzzle!
addpath('./PuzCode/');
%imlist= '../images/1.png';  % people watching waterfall, river
imlist = imread(['/home/hevc/Desktop/Lajan_Research/lajan_puzzle/datasets/mcgill/' num2str(imNo) '.jpg']);

%PP = 28; nc = 24; nr = 18;
PP = 14; nc = 54; nr = 40;

%nc = 24; % number of puzzle pieces across the puzzle
%nr = 18;  % number of puzzle pieces up and down

placementKey = reshape(1:nr*nc, nr,nc);
ScramblePositions =1;% 0 to keep pieces in original locations, 1 to scrable
ScrambleRotations =0;% 0 to keep upright orientation, 1 to scramble

%make the puzzle pieces:
[pieceMat,ap, pi] = PuzzleFun(imlist,PP, nc,nr);

%% SCRAMBLE THE PIECES:
%if(ScramblePositions)
%    ap2 = ap;
%    ppp = randperm(numel(ap));
%    for iii = 1:1:numel(ppp)
%        ap2{iii} = ap{ppp(iii)};
%    end
%    ap=ap2;
%else
%    ppp = 1:numel(ap);
%end
%if(ScrambleRotations)
%    randscram = floor(rand(numel(ap),1)*4)+1;
%    for jj = 1:1:numel(ap)
%        ap2{jj} = imrotate(ap{jj},90*(randscram(jj)-1));
%    end
%    ap=ap2;
%end
%
%
%% compute PuzzlePiece Atribute
%fprintf('Compute Attributes for Puzzle\n');
%pieceAttributes = ComputePieceAttributes(ap);  % Need to write this
%fprintf('Done with Compute Attributes for Puzzle\n');
%
%fprintf('Compute Pairwise Compatibility Scores for Puzzle, Method %d\n',7);
%
%SCO = ComparePiecePairA_ROT(ap, 7,pieceAttributes,ScrambleRotations);  


% SCRAMBLE THE PIECES:
if(ScramblePositions)
    ap2 = ap;
    ppp = 1:numel(ap);
    %ppp = randperm(numel(ap));
    for iii = 1:1:numel(ppp)
        ap2{iii} = ap{ppp(iii)};
    end
    ap=ap2;
else
    ppp = 1:numel(ap);
end
if(ScrambleRotations)
    randscram = floor(rand(numel(ap),1)*4)+1;
    for jj = 1:1:numel(ap)
        ap2{jj} = imrotate(ap{jj},90*(randscram(jj)-1));
    end
    ap=ap2;
end

pcs3 = [];
for i = 1:numel(ap)
	pcs3 = cat(3,pcs3,ap{i});
end
pcs3 = double(pcs3);

%% compute PuzzlePiece Atribute
%%fprintf('Compute Attributes for Puzzle\n');
%%pieceAttributes = mComputePieceAttributes(pcs3);  % Need to write this
%%fprintf('Done with Compute Attributes for Puzzle\n');
%
%%fprintf('Compute Pairwise Compatibility Scores for Puzzle, Method %d\n',7);
%
%%SCO = ComparePiecePairA_ROT(ap, 7,pcs3,ScrambleRotations);  
%
%%SCO = mComparePiecePairA_ROT(7,pcs3,ScrambleRotations);  
%%load('../mextest/sco');


start_time = cputime;

addpath('../mextest');
[greedyOrder,sco,origSco,mni,time0] = greedySolve2(pcs3,nr*nc,nc,nr);
SCO = single(origSco);
%SCO = mscore_mgc_orig(double(pcs3));
SCO = cat(3,SCO(:,:,3),SCO(:,:,2),SCO(:,:,4),SCO(:,:,1));
%SCO = max(max(max(SCO))) - SCO;
%tr = true(nr*nc);
%SCO(cat(3,tr,tr,tr,tr)) = inf;

fprintf('Done with Score Computation for Puzzle Method %d\n',7);
[GI, GR, im, Results, i_,time1,time2] = DoAllAssemblyOfPuzzle(greedyOrder,ap,SCO,nr,nc,1,ppp,mni);
%[GI, i_, Blocks] = DoAllAssemblyOfPuzzle(greedyOrder,ap,SCO,nr,nc,1,ppp,mni);
%save('GI','GI');
%save('i_','i_');
%save('mni','mni');
%save('Blocks','Blocks');


%load('Blocks');
%B = {};
%b = 1;
%for iii = 1:1:numel(Blocks)
%    BlockSize = sum(Blocks{iii}(:)>0);    
%    if BlockSize, 
%        B{b} = Blocks{iii};
%        b = b + 1; 
%    end
%end
%load('mni');
%load('GI');
%load('i_');
%%renderBlocks(B{1},PP,ap);
%%renderBlocks(B{2},PP,ap);
%for it = 1:60
%
%[dests,paths] = findPaths2(mni,432,1);
%assembled = GI(GI>0);
%[r,c] = ind2sub(size(GI),find(GI));
%loc = [r c];
%remaining = setdiff(1:432,assembled);
%Pths = paths(GI(find(GI)),remaining,:);
%
%locns = [];
%Locn = [];
%for i = 1:size(Pths,1)
%    for j = 1:size(Pths,3)
%        D = loc(i,:) - dests(j,:);
%        if ismember(D,loc,'rows'), continue; end
%        id = [];
%        if ~isempty(locns), id = find(ismember(locns,D,'rows')); end
%        if id, Locn(id,:) = Locn(id,:) + (Pths(i,:,j)>0);
%        else 
%            locns = [locns;D];
%            Locn(end+1,:) = Pths(i,:,j);
%        end
%    end
%end
%
%[u,v] = max(Locn(:));
%u
%if u == 0, break; end
%[p,q] = ind2sub(size(Locn),v);
%remaining(q)
%l = locns(p,:)
%%pause
%%if l(1) > 0 && l(2) > 0 && l(1) <= size(GI,1) && l(2) <= size(GI,2)
%%    GI(l(1),l(2)) = remaining(q);
%%else
%    if l(1) <= 0,           GI = [zeros(1,size(GI,2));GI]; end
%    if l(2) <= 0,           GI = [zeros(size(GI,1),1) GI]; end
%    if l(1) > size(GI,1),   GI = [GI;zeros(1,size(GI,2))]; end
%    if l(2) > size(GI,2),   GI = [GI zeros(size(GI,1),1)]; end
%%end
%if l(1) == 0, l(1) = 1; end
%if l(2) == 0, l(2) = 1; end
%GI(l(1),l(2)) = remaining(q);
%
%end
%
%renderBlocks(GI,PP,ap);

imwrite(i_,['a' num2str(imNo) '.jpg'],'Quality',99);
imwrite(im,['b' num2str(imNo) '.jpg'],'Quality',99);
%figure; 
%imagesc(im);axis image;axis off
%title('solved puzzle'); 
%imwrite(im,'SolvedDemo.jpg','Quality',99);

%imS  = renderPiecesFromGraphIDs(ap,placementKey,0);
%figure; 
%imagesc(imS);axis image;axis off
%title('original puzzle'); 
%imwrite(imS,'ScrambledDemo.jpg','Quality',99);



% evaluate the results:
G_safe = GI;
G_undo = (G_safe==0);
G_safe(G_safe==0) = 1;
G_temptemp = ppp(G_safe);
G_temptemp(G_undo) = 0;
[Res] = EvalPuzzleAssembly(G_temptemp, nc, nr);
% count the number with the right rotation:
%G_temprot = GR(ppp);
LUT = [1 4 3 2];
%sum(LUT(G_temprot)==randscram');
cloc = Res(3); K = nr*nc; 
fprintf('\n\nAccuracy: %.0f pieces of %.0f are exactly in the correct position.\n', cloc,K ); 
%fprintf('Orientation: %d pieces have correct orientation.\n\n\n', sum(LUT(G_temprot)==randscram')); 

cor = [cor;Res];

time = cputime - start_time;

timeStats = [timeStats;time0 time1 time2 time];

[cor timeStats]


end


end


function renderBlocks(GI,PP,ap)
    im = uint16(zeros([size(GI)*PP 3]));
    for i = 1:size(GI,1)
        for j = 1:size(GI,2)
            if GI(i,j) == 0, continue; end
            im(PP*(i-1)+(1:PP),PP*(j-1)+(1:PP),:) = ap{GI(i,j)};
        end
    end
    figure,imshow(im);
end
