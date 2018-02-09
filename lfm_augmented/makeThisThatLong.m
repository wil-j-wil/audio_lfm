function newSig = makeThisThatLong(oldSig,length)

% a function which makes 'oldSig' 'length' long

% Will Wilkinson - Jan 2018

warning('off','MATLAB:colon:nonIntegerIndex')
[oldDur,ind] = max(size(oldSig));
if ind == 1
    oldSig = oldSig';
end
if nargin < 2
    length = oldDur;
end


if oldDur == length % do nothing
    newSig = oldSig;
elseif oldDur > length % downsample
    DS = oldDur / length;
    newSig = oldSig(:,1:DS:end);
else % interpolate
    DS = length / oldDur;
    newSig = [];
    for i=1:size(oldSig,1)
        newSig = [newSig; interp(oldSig(i,:),ceil(DS))]; %#ok<AGROW>
    end
    newSig = newSig(:,1:ceil(DS)/DS:end);
end

if ind == 1
    newSig = newSig';
end
warning('on','MATLAB:colon:nonIntegerIndex')
end