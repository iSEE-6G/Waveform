function z=cyconv2D(x,y)
%Non-Fourier domain cyclic convolution
%
%  z=cyconv(x,y)
 siz=num2cell(size(x));
 subs=cellfun(@(n)[2:n,1:n],siz,'uni',0);
 x=x(subs{:});
 z=convn(x,y,'valid');
end