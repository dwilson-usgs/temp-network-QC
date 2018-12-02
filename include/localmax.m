function [ Imax, maxval ] = localmax(a)

%  [ Imax, maxval ] = localmax(a)
%
% find the local max, returns indicies and values

mincount=0;
l=length(a);
      Imax=[];
      maxval=[];
  if a(2) < a(1),
    mincount = mincount+1;
    Imax(mincount) = 1;
    maxval(mincount) = a(1);
  end

for n = 2:l-1,
  if a(n-1) < a(n) & a(n+1) < a(n),
    mincount = mincount+1;
    Imax(mincount) = n;
    maxval(mincount) = a(n);
  end
end

  if a(l-1) < a(l),
    mincount = mincount+1;
    Imax(mincount) = l;
    maxval(mincount) = a(l);
  end



 return;
