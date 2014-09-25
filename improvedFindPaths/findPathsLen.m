function [dests,paths] = findPaths2(mni,numOfParts,len)

[dests,paths] = findPaths(mni,numOfParts,1);
Dests = dests; Paths = paths; % For case len = 1
for i = 2:len
  [Dests,Paths] = findPaths(mni,numOfParts,i);
  [exst,inds] = ismember(Dests,dests,'rows');
  newDests = ~exst;
  exstngDests = inds(inds > 0);
  %size(Paths(:,:,inds,:))
  %size(paths)
  paths(:,:,exstngDests,:) = paths(:,:,exstngDests,:) + Paths(:,:,exst,:);
  dests = [dests;Dests(~exst,:)];
  paths = cat(3,paths,Paths(:,:,newDests,:));
end

