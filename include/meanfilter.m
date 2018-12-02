function [ pfilt, pstd ] = meanfilter(p,num,dir)

% meanfilter
%
%   function to preform a n point mean filter on a vector
%   [ pfilt, pstd ] = meanfilter(p,n,dir)  where n is odd
%   dir can be 
%       'C' - centered (this is the default)
%       'F' - returned value is for n points FORWARD in time
%       'B' - returned value is for n points BACKWARD in time

if nargin < 3,
    dir='C';
end

nn = length(p);
num = ceil(num/2)*2 -1;
if num<3; num=3; end;
pfilt = zeros(length(p),1);
pstd = zeros(length(p),1);
nnn=floor(num/2);

if strcmp(dir,'C'),
  % calculate the edges
  for n = 1:nnn,
	  pfilt(n)=mean(p(1:n+nnn));
          pfilt(nn-(n-1))=mean(p(nn-(n-1)-nnn:nn));
          	  pstd(n)=std(p(1:n+nnn));
          pstd(nn-(n-1))=std(p(nn-(n-1)-nnn:nn));
  end
  % now fill in the rest
  for n = nnn+1:nn-nnn,
    pfilt(n) = mean(p(n-nnn:n+nnn));
       pstd(n) = std(p(n-nnn:n+nnn));
  end
end
if strcmp(dir,'F'),
  % calculate the edges
  for n = 1:num,
	 % pfilt(n)=mean(p(1:n+nnn));
  %        pfilt(nn-(n-1))=mean(p(nn-(n-1):nn));
         pfilt(nn-(n-1))=mean(p(nn-num:nn));
                  pstd(nn-(n-1))=std(p(nn-num:nn));
  end
  % now fill in the rest
  for n = 1:nn-num,
    pfilt(n) = mean(p(n:n+num));
       pstd(n) = std(p(n:n+num));
  end
end
if strcmp(dir,'B'),
  % calculate the edges
  for n = 1:num,
	%  pfilt(n)=mean(p(1:n));
    pfilt(n)=mean(p(1:num));
        pstd(n)=std(p(1:num));
      %    pfilt(nn-(n-1))=mean(p(nn-(n-1)-nnn:nn));
  end
  % now fill in the rest
  for n = num+1:nn,
    pfilt(n) = mean(p(n-num:n));
       pstd(n) = std(p(n-num:n));
  end
end
return;
